#!/usr/bin/env bash
# Deep-sky galaxy pipeline: calibrate, stack, deconvolve, denoise, stretch.
#
#   ./process.sh sessions/<session>.conf [first step]
#
# Steps (start later with the step number, e.g. "5" to only redo the stretch):
#   1 masters   master dark, dark-flat, flat                      (Siril)
#   2 stack     calibrate, register, select 90 %, stack the lights  (Siril)
#   3 post      crop, gradient removal, plate solve               (Siril)
#   4 deconv    Richardson-Lucy on the galaxy                     (Python)
#   5 denoise   NL-Bayes, partial strength                        (Siril)
#   6 stretch   Gaia white balance, sky-only denoise, stretch     (Python)
# The Siril scripts are written to <WORK>/siril/ so you can read or rerun each one; logs are next to them.
set -euo pipefail

REPO=$(cd "$(dirname "$0")" && pwd)
CONF=${1:?usage: ./process.sh sessions/<session>.conf [first step]}
FIRST=${2:-1}
# shellcheck source=/dev/null
source "$CONF"

# --- settings: the session file sets what it needs, the rest falls back to these ---
: "${SIRIL:=$HOME/Downloads/Siril-1.4.3-x86_64.AppImage}"
: "${PYTHON:=$HOME/iss/.venv/bin/python}"
: "${NAME:?NAME missing}" "${LIGHTS:?LIGHTS missing}" "${WORK:?WORK missing}"
: "${RA:?RA missing}" "${DEC:?DEC missing}" "${MAJOR:?MAJOR missing}" "${MINOR:?MINOR missing}" "${PA:?PA missing}"
: "${DARKS:=}" "${MASTER_DARK:=}"            # darks: a folder of frames, or a ready master
: "${FLATS:=}" "${MASTER_FLAT:=}"            # flats: a folder of frames, or a ready master
: "${DARKFLATS:=}" "${FLAT_BIAS:=}"          # flats' offset: dark-flat frames, or a constant (e.g. 80)
: "${LIGHT_BIAS:=}"                          # no darks for the lights: subtract this constant (ADU, 16-bit) instead
: "${HOT_DARK:=}"                            # no darks: hot/cold pixel map from this master dark (any offset)
: "${FOCAL:=750}" "${PIXEL:=2.9}"
: "${CROP:=40 40 1856 1016}"                 # x y w h after registration (1936 x 1096 frames)
: "${KEEP_BEST:=90%}"                        # frames kept by weighted FWHM
: "${SOLVE:=online}"                         # online = Siril/VizieR first, then offline; offline = astrometry.net only
: "${WCS_FROM:=}"                            # last resort: a solved image of the same field (copies its WCS)
: "${DECONV_ITERS:=20}" "${DECONV_TV:=0.002}" "${DECONV_MAJOR:=$MAJOR}" "${DECONV_MINOR:=$MINOR}"
: "${DENOISE_MOD:=0.7}"
: "${STRETCH_OPTS:=--lumsmooth=0 --black=4 --sky=0.13 --bg=flat --skydn=1.0 --skychroma=3 --skyflat=0.15}"
: "${TAG:=final}"
: "${KEEP_FRAMES:=0}"                        # 1 = keep the per-frame calibrated/registered files

mkdir -p "$WORK/siril" "$WORK/masters"
S="$WORK/siril"

run_siril() {   # run_siril <name>: runs $S/<name>.ssf, log in $S/<name>.log
    echo "--- Siril: $1"
    if ! "$SIRIL" siril-cli -s "$S/$1.ssf" > "$S/$1.log" 2>&1; then
        grep -v '^progress' "$S/$1.log" | tail -25
        echo "Siril script $1 failed (full log: $S/$1.log)"; exit 1
    fi
    grep -E 'Total: |Integration of|Background noise value \(channel: #1\)|Resolution:|Image center|NL-Bayes exec' "$S/$1.log" || true
}
py() { "$PYTHON" -I "$@"; }

# ---------------------------------------------------------------- 1. masters
if [ "$FIRST" -le 1 ]; then
    {
        echo "requires 1.4.0"; echo "set16bits"; echo "setext fits"
        if [ -z "$MASTER_DARK" ] && [ -n "$DARKS" ]; then
            echo "# master dark: mean, winsorized 3/3 sigma rejection, no normalisation (keep the absolute level)"
            echo "cd \"$DARKS\""; echo "link dark -out=$WORK/masters/dark"
            echo "cd \"$WORK/masters/dark\""; echo "stack dark rej 3 3 -nonorm -out=../master_dark"
        fi
        if [ -z "$MASTER_FLAT" ] && [ -n "$FLATS" ]; then
            if [ -n "$DARKFLATS" ]; then
                echo "# master dark-flat"
                echo "cd \"$DARKFLATS\""; echo "link df -out=$WORK/masters/darkflat"
                echo "cd \"$WORK/masters/darkflat\""; echo "stack df rej 3 3 -nonorm -out=../master_darkflat"
                FLAT_OFFSET="-dark=$WORK/masters/master_darkflat"
            elif [ -n "$FLAT_BIAS" ]; then
                FLAT_OFFSET="-bias=\"=$FLAT_BIAS\""
            else
                echo "need DARKFLATS or FLAT_BIAS for the flats" >&2; exit 1
            fi
            echo "# master flat: offset removed, then multiplicative normalisation"
            echo "cd \"$FLATS\""; echo "link flat -out=$WORK/masters/flat"
            echo "cd \"$WORK/masters/flat\""; echo "calibrate flat $FLAT_OFFSET"
            echo "stack pp_flat rej 3 3 -norm=mul -out=../master_flat"
        fi
    } > "$S/1_masters.ssf"
    run_siril 1_masters
fi
[ -z "$MASTER_DARK" ] && [ -n "$DARKS" ] && MASTER_DARK="$WORK/masters/master_dark"
[ -z "$MASTER_FLAT" ] && [ -n "$FLATS" ] && MASTER_FLAT="$WORK/masters/master_flat"

# ---------------------------------------------------------------- 2. stack
if [ "$FIRST" -le 2 ]; then
    CAL=""
    PRE=""
    if [ -n "$MASTER_DARK" ]; then
        CAL="$CAL -dark=$MASTER_DARK -cc=dark"
    else
        [ -n "$LIGHT_BIAS" ] && CAL="$CAL -bias=\"=$LIGHT_BIAS\""
        if [ -n "$HOT_DARK" ]; then
            # hot pixels sit at the same place whatever the offset: map them from an older master dark
            PRE="load $HOT_DARK
find_hot $WORK/masters/hotpixels 3 3"
            CAL="$CAL -cc=bpm $WORK/masters/hotpixels.lst"
        fi
    fi
    [ -n "$MASTER_FLAT" ] && CAL="$CAL -flat=$MASTER_FLAT"
    cat > "$S/2_stack.ssf" <<EOF
requires 1.4.0
set16bits
setext fits
$PRE
cd "$LIGHTS"
link light -out=$WORK/lights
cd "$WORK/lights"
# dark, flat, hot/cold pixels from the dark, debayer (-equalize_cfa keeps the flat from tinting)
calibrate light$CAL -cfa -equalize_cfa -debayer
# few stars per frame (small field, bright sky): more sensitive star detection
setfindstar -sigma=0.5 -relax=on
register pp_light -2pass -minpairs=6
seqapplyreg pp_light -filter-wfwhm=$KEEP_BEST
set32bits
# winsorized sigma 3/3, additive+scale normalisation, weighted by FWHM
stack r_pp_light rej 3 3 -norm=addscale -output_norm -rgb_equal -weight=wfwhm -out=$WORK/${NAME}_stack
EOF
    run_siril 2_stack
    if [ "$KEEP_FRAMES" != 1 ]; then
        find "$WORK/lights" -maxdepth 1 \( -name 'pp_*' -o -name 'r_pp_*' -o -type l \) -delete
    fi
fi

# ---------------------------------------------------------------- 3. post
if [ "$FIRST" -le 3 ]; then
    cat > "$S/3_post.ssf" <<EOF
requires 1.4.0
cd "$WORK"
load ${NAME}_stack
# cut the ragged registration edges, remove the gradient (2nd-degree polynomial)
crop $CROP
subsky 2 -samples=20 -tolerance=1.0
save ${NAME}_linear
EOF
    run_siril 3_post
    # plate solve (online catalogue via VizieR); its own script so a server outage does not stop the run
    cat > "$S/3_solve.ssf" <<EOF
requires 1.4.0
cd "$WORK"
load ${NAME}_linear
platesolve $RA $DEC -focal=$FOCAL -pixelsize=$PIXEL -force
save ${NAME}_linear
EOF
    echo "--- plate solve ($SOLVE)"
    [ "$SOLVE" = offline ] && echo "online solve skipped (SOLVE=offline)" > "$S/3_solve.log"
    if [ "$SOLVE" != offline ] && "$SIRIL" siril-cli -s "$S/3_solve.ssf" > "$S/3_solve.log" 2>&1; then
        grep -E 'Resolution:|Image center' "$S/3_solve.log" || true
    elif echo "--- astrometry.net offline ($(grep -o -m1 -E 'skipped|libcurl error[^.]*|HTTP code [0-9]+|could not be aligned' "$S/3_solve.log"))" &&
         py "$REPO/scripts/solve_local.py" "$WORK/${NAME}_linear.fits" "$RA" "$DEC" \
            "$(awk "BEGIN {printf \"%.3f\", 206.265 * $PIXEL / $FOCAL}")" 2>&1 | grep -v -i warn | tail -3; then
        :   # solved offline (pipefail: the pipeline fails when solve_local.py does)
    elif [ -n "$WCS_FROM" ]; then
        echo "plate solve failed ($(grep -o -m1 -E 'libcurl error[^.]*|HTTP code [0-9]+|could not be aligned' "$S/3_solve.log")) -> WCS copied from $WCS_FROM"
        py "$REPO/scripts/align_wcs.py" "$WCS_FROM" "$WORK/${NAME}_linear.fits" "$WORK/${NAME}_linear.fits" 2>&1 | grep -v -i warn
    else
        grep -v '^progress' "$S/3_solve.log" | grep -i -E 'error|failed' | tail -5
        echo "Plate solve failed. If the catalogue server is down, set WCS_FROM=<a solved image of this field> in the"
        echo "session file and rerun from step 3:  ./process.sh $CONF 3"; exit 1
    fi
fi

# ---------------------------------------------------------------- 4. deconv
if [ "$FIRST" -le 4 ]; then
    echo "--- deconvolution"
    py "$REPO/scripts/deconv_galaxy.py" "$WORK/${NAME}_linear.fits" "$WORK/${NAME}_deconv.fits" \
        "$DECONV_ITERS" "$DECONV_TV" "--gal=$RA,$DEC,$DECONV_MAJOR,$DECONV_MINOR,$PA" 2>&1 | grep -v -i -E "warn|RADECSYS"
fi

# ---------------------------------------------------------------- 5. denoise
if [ "$FIRST" -le 5 ]; then
    cat > "$S/5_denoise.ssf" <<EOF
requires 1.4.0
cd "$WORK"
load ${NAME}_deconv
denoise -nocosmetic -mod=$DENOISE_MOD
save ${NAME}_deconv_dn
EOF
    run_siril 5_denoise
fi

# ---------------------------------------------------------------- 6. stretch
if [ "$FIRST" -le 6 ]; then
    echo "--- stretch"
    # shellcheck disable=SC2086
    py "$REPO/scripts/stretch_galaxy.py" "$WORK/${NAME}_deconv_dn.fits" "$NAME" "$MAJOR" "$MINOR" "$PA" "$RA" "$DEC" \
        $STRETCH_OPTS "--tag=$TAG" 2>&1 | grep -v -i -E "warn|RADECSYS|_ureduce|true_divide|out = np.where"
fi
echo "done: $WORK/final/"

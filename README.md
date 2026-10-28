# astro-processing

Deep-sky galaxy processing for the ASI290MC on the 150/750 Newtonian: from raw Ekos frames to a finished,
north-up image with one command. Classical math only (calibration, Richardson-Lucy, NL-Bayes, wavelets,
curves), no AI tools.

First result: NGC 6503, 25 Oct 2026 (`sessions/2026-10-25_ngc6503.conf`, notes in `docs/ngc6503-2026-10-25.md`).


## Quick start

```bash
cp sessions/2026-10-25_ngc6503.conf sessions/<date>_<target>.conf   # edit paths and the galaxy's numbers
./process.sh sessions/<date>_<target>.conf          # all steps, ~5-8 minutes
./process.sh sessions/<date>_<target>.conf 6        # only redo the stretch (e.g. after changing STRETCH_OPTS)
```

Result: `<WORK>/final/<NAME>_final_field.png` (north up, east left) and `<NAME>_final_crop_2x.png`.
The linear images after each step stay in `<WORK>` (`_stack`, `_linear`, `_deconv`, `_deconv_dn`), and the
Siril scripts the run wrote, with their logs, in `<WORK>/siril/`.


## What you need

- **Siril 1.4.3 AppImage** at `~/Downloads/Siril-1.4.3-x86_64.AppImage` (or set `SIRIL=` in the session file).
- **Python 3** with `requirements.txt` (numpy, scipy, astropy, pillow). Default: the ISS project's
  `~/iss/.venv/bin/python`; or `python3 -m venv .venv && .venv/bin/pip install -r requirements.txt` and `PYTHON=.venv/bin/python`.
- **Internet:** the plate solve uses VizieR (Siril), the white balance uses the Gaia archive at ESA (cached per
  target in `<WORK>/<NAME>_gaia.csv`). When VizieR is down, see `WCS_FROM` below.
- **Disk:** ~7 GB free while step 2 runs (250 frames). The per-frame files are deleted after stacking
  (`KEEP_FRAMES=1` keeps them).


## Capture checklist (so the processing has what it needs)

- **Lights:** 30 s, gain 100, **offset 30**. Dithering on: Ekos Guide → Options → Dither 3 px, every 3 frames,
  threshold 1.5 px, settle 5 s, timeout 60 s. Rename the sequence when you change target (Ekos does not).
- **Darks:** 30–50 × the lights' exposure, same gain and offset, at the end of the night outdoors
  (uncooled camera: temperature must match), nose capped.
- **Flats:** 30, Ekos *Flat* with *Duration: ADU* 30000 ± 1000, before touching the camera, focus or rotation.
  Light: `tools/flat_panel.html` full screen behind two layers of white T-shirt over the tube. Device backlight at
  maximum, dim with the page's ↑/↓ (a dimmed backlight flickers and bands short exposures).
- **Dark-flats:** 30 × the flats' exposure, same gain and offset as the flats, nose capped (copy the Flat job,
  change only the frame type to Dark).
- **Check what you shot** before processing: Ekos headers carry no pointing. Plate solve one frame per block of the
  night (gaps in DATE-OBS = slews). On 25 Oct, 196 NGC 6503 frames were saved as NGC 5866.

Put the lights of one target in one folder. Links keep the originals untouched:

```bash
mkdir -p ~/Downloads/NGC_6503/raw_lights
for i in $(seq 177 372); do ln -sf ../../NGC_5866/Light/NGC_5866_Light_$i.fits ~/Downloads/NGC_6503/raw_lights/a_$i.fits; done
for i in $(seq 1 58); do n=$(printf %03d $i); ln -sf ../Light/NGC_6503_Light_$n.fits ~/Downloads/NGC_6503/raw_lights/b_$n.fits; done
```


## Session file

A shell file with `NAME=value` lines (see `sessions/`). Required:

- `NAME`: short name for output files, e.g. `ngc6503`
- `LIGHTS`: folder with the lights (or links to them)
- `WORK`: output folder (created)
- `RA`, `DEC`: the galaxy's centre, e.g. `17:49:26.4`, `+70:08:40`
- `MAJOR`, `MINOR`, `PA`: its size in arcmin and position angle in degrees (NED, SIMBAD or Wikipedia)

Calibration (all optional, but use them):

- `DARKS` (folder) or `MASTER_DARK` (a master from an earlier run, path without `.fits`)
- `FLATS` (folder) or `MASTER_FLAT`
- `DARKFLATS` (folder), or `FLAT_BIAS=<ADU>` when you have none: a constant subtracted from the flats
- No darks for these lights (e.g. the offset changed): `LIGHT_BIAS=<ADU>` subtracts a constant instead
  (an error of tens of ADU is harmless: step 3 removes any leftover constant), and `HOT_DARK=<older master dark>`
  maps hot/cold pixels for cosmetic correction (they sit at the same places whatever the offset)

Tuning (defaults are the 25 Oct NGC 6503 result):

- `KEEP_BEST=90%`: frames kept, by weighted FWHM
- `CROP="40 40 1856 1016"`: x y w h cut after registration (ragged edges)
- `SOLVE=online`: plate solve with Siril (VizieR) first, astrometry.net offline when that fails;
  `SOLVE=offline` skips the online try (it costs ~3 min when VizieR is down)
- `WCS_FROM`: last resort, an already solved image of the same field (copies its WCS; handles flips, shifts and
  rotations up to 5°)
- `LIGHT_BIAS`, `HOT_DARK`: see above (lights without matching darks)
- `DECONV_ITERS=20`, `DECONV_TV=0.002`: Richardson-Lucy iterations and smoothing
- `DECONV_MINOR`: minor axis (arcmin) of the area left open for deconvolution; make it smaller when stars next to
  the galaxy get dark rings (NGC 5866 needed 1.0 instead of 2.8)
- `DENOISE_MOD=0.7`: NL-Bayes strength (1 = full)
- `STRETCH_OPTS`: options for `scripts/stretch_galaxy.py` (below)
- `TAG=final`: output name suffix; change it when trying variants so results are not overwritten


## What each step does

**1. Masters (Siril).** Master dark: mean of the darks with winsorized sigma rejection 3/3, no normalisation (a
dark must keep its absolute level). Master dark-flat the same way. Master flat: each flat minus the dark-flat, then
a mean with multiplicative normalisation (only the flat's shape matters).

**2. Stack (Siril).**
- `calibrate light -dark -cc=dark -flat -cfa -equalize_cfa -debayer`: subtract the dark, divide by the flat
  (dust donuts, vignetting), replace hot/cold pixels found in the dark, keep the flat from tinting, debayer RGGB.
- `setfindstar -sigma=0.5 -relax=on` and `register -2pass -minpairs=6`: these frames hold only 10–30 stars, so
  star detection is made more sensitive; two passes pick the best reference frame first.
- `seqapplyreg -filter-wfwhm=90%`: keep the sharpest 90 %.
- `stack rej 3 3 -norm=addscale -output_norm -rgb_equal -weight=wfwhm`: winsorized 3/3 rejection (satellites,
  planes), each frame's sky level and contrast matched, sharper frames weighted more, equal R/G/B background.

**3. Post (Siril).** Crop, then `subsky 2`: a 2nd-degree polynomial background (moonlight, town glow) removed per
channel. Then the plate solve writes sky coordinates (WCS) into the header; steps 4 and 6 need them. Siril's solver
(online catalogue) is tried first, then `scripts/solve_local.py`: astrometry.net's `solve-field` offline, with the
4205-4207 index series (all-sky, 1.3 GB in `~/.local/share/astrometry`, copied from KStars' folder on the Pi).
If both fail and `WCS_FROM` is set, `scripts/align_wcs.py` finds the flip and sub-pixel shift between the two images by
phase correlation, moves the new image onto the solved one's grid and copies its header. That shift leaves a
smeared strip (~16 px) on one edge of the final image; crop it off.
Siril's photometric colour calibration is not used: on this camera it left blue too strong.

**4. Deconvolution (`scripts/deconv_galaxy.py`).**
- PSF measured from the image: bright (>60σ), unsaturated stars with no neighbour within 12 px, centred to
  sub-pixel accuracy and median-combined (13 stars, FWHM ~9 px ≈ 7″ for NGC 6503).
- Richardson-Lucy with total-variation regularisation (Dey et al. 2006), on the luminance only.
- Applied only where the signal is strong (weight 0 at ≤3σ, 1 at ≥15σ), never on saturated stars, and every field
  star is masked with a disc sized to its brightness. Unmasked stars get dark rings, because RL pulls light from a
  star's halo into its core. Inside the galaxy's ellipse (0.45 × major by 0.6 × `DECONV_MINOR`) nothing is masked.
- Colour kept: every channel is multiplied by the same factor (deconvolved / original luminance).

**5. Denoise (Siril).** NL-Bayes (non-local Bayesian, Lebrun–Buades–Morel 2013) at `DENOISE_MOD` strength: the
result is 70 % denoised + 30 % original, which keeps a natural grain.

**6. Stretch (`scripts/stretch_galaxy.py`).**
1. White balance from Gaia DR3: stars G < 16 in the field, R/G/B flux measured in an aperture, a fit of colour
   ratio against Gaia BP−RP, scaled so a Sun-like star (BP−RP = 0.82) is white.
2. R and B shifted onto G by the measured offset (atmospheric dispersion, ~0.3 px at 35–42° altitude).
3. `--bg=flat`: one constant per channel subtracted (the gradient was removed in step 3). `--bg=grid` models the
   background on a 64 px grid: only for data without flats, it adds blotches otherwise.
4. Colour noise: the colour (not the brightness) smoothed with a 2.5 px Gaussian.
5. Sky only (where the signal is below ~2–8σ): `--skydn=K` soft-thresholds 4 starlet (à trous) wavelet scales at
   K × each scale's noise; `--skychroma=S` adds S px colour smoothing. 1.0 / 3 is light; 1.5 / 4 mild;
   2.5 / 6 strong, too smooth (faint stars disappear).
6. Stretch: black point `--black=4` sky-σ below the sky (2 clips noise dips to pure black grains), arcsinh curve
   with the white point above the galaxy core, sky at `--sky=0.13`, `--skyflat=0.15` lowers contrast just around
   the sky level. `--lumsmooth` (luminance blur in faint areas) is off; the wavelet step replaces it.
7. Saturation ×1.35 keeping hue, rotated to north up / east left from the WCS, field and 2× crop saved.

Also available: `--bin=2` (2×2 software binning). On this camera it reduced noise only 1.2–1.3×, because debayering
and registration correlate neighbouring pixels, so it is not used.


## Troubleshooting

- **Many frames fail to register** → keep `setfindstar -sigma=0.5 -relax=on` (step 2 does). On a finished stack
  use `setfindstar reset`, or noise is detected as thousands of stars.
- **Plate solve fails at the right coordinates** → are the frames really that target? Solve one at other candidates.
- **"libcurl error" / "HTTP code 503"** in `3_solve.log` → the VizieR server is down; the offline solve takes
  over (set `SOLVE=offline` to skip the wait).
- **Offline solve fails** → check `~/.local/share/astrometry` holds the index files; the system's astrometry.net
  helper scripts break with NumPy 2, so `solve_local.py` runs solve-field with `--no-remove-lines --uniformize 0`.
- **Donuts left or bright rings after flat-fielding** → the flats don't match the lights (camera moved, wrong set).
- **Dark rings around stars after deconvolution** → a smaller `DECONV_MINOR`.
- **`find_hot name`** (if you use it in your own scripts) writes `name.lst`: give the name without `.lst`.


## Layout

- `process.sh`: the pipeline
- `scripts/deconv_galaxy.py`, `scripts/stretch_galaxy.py`, `scripts/align_wcs.py`: the Python steps (each has
  its usage in its first lines and runs on its own too)
- `sessions/*.conf`: one file per session/target
- `tools/flat_panel.html`: white screen for flats
- `docs/`: notes per session

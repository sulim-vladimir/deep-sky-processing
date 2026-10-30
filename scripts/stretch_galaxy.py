"""Stretch a Siril linear, colour-calibrated galaxy stack (no flats) into display images.

Run:  ~/iss/.venv/bin/python stretch_galaxy.py <linear_cc.fits> <name> <major'> <minor'> <PA deg> [ra dec] [--bin=2]
      e.g. stretch_galaxy.py result_ngc5866_linear_cc.fits ngc5866 6.4 2.8 128 --bin=2
Output: final/<name>[_bin2]_field.png and final/<name>[_bin2]_crop_*x.png next to the input file.

Other options: --bg=grid|flat (grid = masked 64-px background model, flat = only a constant per channel,
relying on Siril's subsky), --black=2.2 (black point in sky sigmas below the sky), --sky=0.10 (sky level
in the output), --chroma=2.5 (colour-noise smoothing everywhere, px; 0 = off),
--lumsmooth=1.3 (luminance smoothing in faint areas, px; 0 = off), --tag=<suffix>,
--galmask=1.9,3.0 (galaxy mask = these multiples of the major, minor axis; kept off the sky and the white-balance
stars. Use ~1.0,1.3 for a galaxy that fills much of the field, e.g. M81).
Sky-only options (masked off the galaxy and stars): --skydn=K wavelet (starlet, 4 scales) soft threshold at
K x noise per scale; --skychroma=S extra colour smoothing (px); --skyflat=A (0..0.9) lowers the stretch
contrast right at the sky level so sky noise is shown smaller.

--bin=N averages NxN pixels after white balance / channel alignment (done at full resolution).
At 0.8"/px with 6-7" stars the data are heavily oversampled, so 2x2 costs no detail and cuts noise.
All pixel-based settings below are scaled by 1/N so they cover the same patch of sky.

Steps:
 1. background model on a fine grid with stars and the galaxy masked -> removes vignetting,
    dust blotches and gradients that flats would normally take out; background made neutral
 2. colour noise: chrominance smoothed; luminance smoothed only where the signal is faint
 3. arcsinh stretch, black point just under the sky, white point on the galaxy core
 4. hue-preserving saturation, rotate to north up / east left, galaxy crop
"""
import os
import sys

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import proj_plane_pixel_scales
from scipy import ndimage as nd
import astropy.units as au
from PIL import Image

BIN = 1
OPT = {"bg": "grid", "black": "2.2", "sky": "0.10", "lumsmooth": "1.3", "tag": "",
       "skydn": "0", "skychroma": "0", "skyflat": "0", "chroma": "2.5", "galmask": "1.9,3.0"}
args = []
for a in sys.argv[1:]:
    if a.startswith("--bin="):
        BIN = int(a.split("=")[1])
    elif a.startswith("--") and a[2:].split("=")[0] in OPT:
        k, v = a[2:].split("=", 1)
        OPT[k] = v
    else:
        args.append(a)
sys.argv[1:] = args
fn, name = sys.argv[1], sys.argv[2]
major, minor, pa = (float(v) for v in sys.argv[3:6])
GM_A, GM_B = (float(v) for v in OPT["galmask"].split(","))
tag = (name if BIN == 1 else f"{name}_bin{BIN}") + (f"_{OPT['tag']}" if OPT["tag"] else "")


def px(v):
    """a length given in full-resolution pixels -> pixels of the (binned) working image"""
    return v / BIN

HERE = os.path.dirname(os.path.abspath(fn))
OUT = os.path.join(HERE, "final")
os.makedirs(OUT, exist_ok=True)

hdr = fits.getheader(fn)
wcs = WCS(hdr, naxis=2)
raw = fits.getdata(fn)
rgb = np.moveaxis(raw.astype(np.float64), 0, -1)                # (H, W, 3), FITS row 0 = bottom
if np.issubdtype(raw.dtype, np.integer):                        # 16-bit stack (e.g. an external one): to 0..1
    rgb /= 65535.0
del raw
H, W, _ = rgb.shape
scale = proj_plane_pixel_scales(wcs).mean() * 3600               # arcsec/px

# galaxy position: from argv, or the built-in table
if len(sys.argv) > 7:
    gal = SkyCoord(sys.argv[6], sys.argv[7], unit=("hourangle", "deg"))
else:
    gal = SkyCoord(*{"ngc5866": ("15:06:29.5", "+55:45:48"),
                     "ngc6503": ("17:49:26.4", "+70:08:40")}[name], unit=("hourangle", "deg"))
gx, gy = wcs.world_to_pixel(gal)
# galaxy major-axis direction in pixels (PA is east of north)
ahead = gal.directional_offset_by(pa * au.deg, 0.01 * au.deg)
ax, ay = wcs.world_to_pixel(ahead)
theta = np.arctan2(ay - gy, ax - gx)
yy, xx = np.mgrid[0:H, 0:W]
u = (xx - gx) * np.cos(theta) + (yy - gy) * np.sin(theta)
v = -(xx - gx) * np.sin(theta) + (yy - gy) * np.cos(theta)
a_px, b_px = major * 60 / scale / 2, minor * 60 / scale / 2
gal_mask = (u / (GM_A * a_px)) ** 2 + (v / (GM_B * b_px)) ** 2 < 1
print(f"{name}: galaxy at x={gx:.0f} y={gy:.0f}, {scale:.3f}\"/px")

# --- 0a. white balance from Gaia star colours (Siril's PCC leaves blue too strong) ---
# Sun-like stars (BP-RP = 0.82) are made white. Gaia stars are cached next to the input.
gaia_cache = os.path.join(HERE, f"{name}_gaia.csv")
if not os.path.exists(gaia_cache):
    import urllib.parse
    import urllib.request
    centre = wcs.pixel_to_world(W / 2, H / 2)
    q = ("SELECT ra,dec,phot_g_mean_mag,bp_rp FROM gaiadr3.gaia_source WHERE "
         f"1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{centre.ra.deg:.5f},{centre.dec.deg:.5f},0.3)) "
         "AND phot_g_mean_mag<16 AND bp_rp IS NOT NULL")
    url = "https://gea.esac.esa.int/tap-server/tap/sync?" + urllib.parse.urlencode(
        {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": q})
    with open(gaia_cache, "wb") as fh:
        fh.write(urllib.request.urlopen(url, timeout=120).read())
cat = np.genfromtxt(gaia_cache, delimiter=",", names=True)
spx, spy = wcs.world_to_pixel(SkyCoord(cat["ra"], cat["dec"], unit="deg"))
R, R1, R2 = 10, 15, 22
yy2, xx2 = np.mgrid[-R2:R2 + 1, -R2:R2 + 1]
rr = np.hypot(yy2, xx2)
meas, cents = [], []
for x, y, c in zip(spx, spy, cat["bp_rp"]):
    xi, yi = int(round(x)), int(round(y))
    if not (R2 < xi < W - R2 and R2 < yi < H - R2) or gal_mask[yi, xi]:
        continue
    cut = rgb[yi - R2:yi + R2 + 1, xi - R2:xi + R2 + 1]
    if cut[R2 - 2:R2 + 3, R2 - 2:R2 + 3].max() > 0.9:
        continue                                                   # saturated
    sub = cut - np.median(cut[(rr >= R1) & (rr <= R2)], axis=0)
    fl = sub[rr <= R].sum(axis=0)
    if fl.min() > 0:
        meas.append((c, *fl))
        w8 = np.clip(sub[R2 - 8:R2 + 9, R2 - 8:R2 + 9], 0, None)
        gy8, gx8 = np.mgrid[-8:9, -8:9]
        cents.append([((w8[..., k] * gy8).sum() / w8[..., k].sum(), (w8[..., k] * gx8).sum() / w8[..., k].sum())
                      for k in range(3)])
m = np.array(meas)
p1 = np.polyfit(m[:, 0], np.log10(m[:, 1] / m[:, 2]), 1)
p2 = np.polyfit(m[:, 0], np.log10(m[:, 3] / m[:, 2]), 1)
wb = np.array([10 ** -np.polyval(p1, 0.82), 1.0, 10 ** -np.polyval(p2, 0.82)])
rgb *= wb
print(f"white balance from {len(m)} Gaia stars: R x{wb[0]:.3f}  B x{wb[2]:.3f}")

# --- 0b. align R and B onto G (atmospheric dispersion at low altitude) ---
cents = np.array(cents)                                            # (n, 3, 2) dy, dx
offs = np.median(cents - cents[:, 1:2, :], axis=0)
for k in (0, 2):
    rgb[..., k] = nd.shift(rgb[..., k], -offs[k], order=3, mode="nearest")
print("channel shift (dy, dx px): R", offs[0].round(2), " B", offs[2].round(2))


def sky_noise(img):
    """robust per-pixel noise of the luminance, measured on a smooth-subtracted copy"""
    l = img.mean(axis=2)
    r = l - nd.median_filter(l, size=9)
    return 1.4826 * np.median(np.abs(r - np.median(r)))


noise_full = sky_noise(rgb)

# --- 0c. software binning ---
if BIN > 1:
    H, W = H // BIN * BIN, W // BIN * BIN
    rgb = rgb[:H, :W].reshape(H // BIN, BIN, W // BIN, BIN, 3).mean(axis=(1, 3))
    wcs = wcs[::BIN, ::BIN]                                        # numpy order: (y, x)
    H, W = rgb.shape[:2]
    scale *= BIN
    gx, gy = wcs.world_to_pixel(gal)
    ax, ay = wcs.world_to_pixel(ahead)
    theta = np.arctan2(ay - gy, ax - gx)
    yy, xx = np.mgrid[0:H, 0:W]
    u = (xx - gx) * np.cos(theta) + (yy - gy) * np.sin(theta)
    v = -(xx - gx) * np.sin(theta) + (yy - gy) * np.cos(theta)
    a_px, b_px = a_px / BIN, b_px / BIN
    gal_mask = (u / (GM_A * a_px)) ** 2 + (v / (GM_B * b_px)) ** 2 < 1
    noise_bin = sky_noise(rgb)
    print(f"binned {BIN}x{BIN}: {W}x{H} px, {scale:.2f}\"/px, galaxy at x={gx:.0f} y={gy:.0f}")
    print(f"sky noise per pixel: full {noise_full:.3e} -> binned {noise_bin:.3e}  "
          f"(improvement x{noise_full / noise_bin:.2f}, ideal x{BIN:.0f})")
else:
    print(f"sky noise per pixel: {noise_full:.3e}")

# --- star mask ---
lum = rgb.mean(axis=2)
med = np.median(lum)
mad = 1.4826 * np.median(np.abs(lum - med))
sm = nd.gaussian_filter(lum, px(1.5))
local_bg = nd.median_filter(sm[::4, ::4], size=max(3, int(px(15))))
local_bg = nd.zoom(local_bg, (H / local_bg.shape[0], W / local_bg.shape[1]), order=1)[:H, :W]
stars = sm - local_bg > 4 * mad / 1.5 * (1 if BIN == 1 else 0.7)
lab, n = nd.label(stars)
# grow each star in proportion to its brightness so halos are masked too
peak = nd.maximum(sm - local_bg, lab, range(1, n + 1))
star_mask = np.zeros_like(stars)
for r in (3, 6, 10, 16):
    big = np.isin(lab, np.nonzero(peak > mad * {3: 0, 6: 30, 10: 150, 16: 800}[r])[0] + 1)
    star_mask |= nd.binary_dilation(big, iterations=max(1, round(px(r))))
mask = star_mask | gal_mask
print(f"stars masked: {n}, masked area {mask.mean() * 100:.0f} %")

# --- 1. background model (per channel) on a grid of boxes ---
box = int(px(64))
gh, gw = int(np.ceil(H / box)), int(np.ceil(W / box))
bg = np.empty((H, W, 3))
for c in range(3):
    grid = np.full((gh, gw), np.nan)
    ch = rgb[..., c]
    for i in range(gh):
        for j in range(gw):
            sl = (slice(i * box, (i + 1) * box), slice(j * box, (j + 1) * box))
            vals = ch[sl][~mask[sl]]
            if vals.size > 0.25 * box * box:
                m, s = np.median(vals), 1.4826 * np.median(np.abs(vals - np.median(vals)))
                vals = vals[np.abs(vals - m) < 2.5 * s]
                grid[i, j] = np.median(vals)
    # fill holes (galaxy, bright stars) from neighbours, then smooth lightly
    while np.isnan(grid).any():
        filled = nd.generic_filter(grid, np.nanmedian, size=3, mode="nearest")
        grid = np.where(np.isnan(grid), filled, grid)
    grid = nd.median_filter(grid, size=3, mode="nearest")
    grid = nd.gaussian_filter(grid, 0.8, mode="nearest")
    full = nd.zoom(grid, (box, box), order=3, mode="nearest")[:H, :W]
    bg[..., c] = full
if OPT["bg"] == "flat":                           # constant per channel only
    bg = np.broadcast_to(np.array([np.median(rgb[..., c][~mask]) for c in range(3)]), rgb.shape)
rgb = rgb - bg                                    # sky now ~0 and neutral in every channel

# --- noise level ---
lum = rgb.mean(axis=2)
sky = ~mask
sigma = 1.4826 * np.median(np.abs(lum[sky] - np.median(lum[sky])))

# --- 2. noise reduction ---
chroma = rgb - lum[..., None]
if float(OPT["chroma"]) > 0:                                     # --chroma=0: no colour-noise smoothing
    chroma = np.dstack([nd.gaussian_filter(chroma[..., c], px(float(OPT["chroma"]))) for c in range(3)])
snr = nd.gaussian_filter(lum, px(2)) / sigma
w_smooth = np.clip(1 - (snr - 1) / 6, 0, 1)                       # 1 in the sky, 0 on real signal
ls = float(OPT['lumsmooth'])
lum_s = (1 - w_smooth) * lum + w_smooth * nd.gaussian_filter(lum, px(ls)) if ls > 0 else lum
lum = lum_s

# sky-only noise reduction (weight 1 in empty sky, 0 on the galaxy and stars)
w_sky = nd.gaussian_filter(np.clip(1 - (snr - 2) / 6, 0, 1), px(2))
K = float(OPT["skydn"])
if K > 0:
    c = lum.copy()
    rec = np.zeros_like(lum)
    for j in range(4):
        step = 2 ** j
        kern = np.zeros(4 * step + 1)
        kern[::step] = np.array([1, 4, 6, 4, 1]) / 16
        cn = nd.convolve1d(nd.convolve1d(c, kern, axis=0, mode="reflect"), kern, axis=1, mode="reflect")
        wj = c - cn
        sj = 1.4826 * np.median(np.abs(wj[sky]))
        rec += np.sign(wj) * np.maximum(np.abs(wj) - K * sj, 0)
        c = cn
    lum_dn = c + rec
    lum = w_sky * lum_dn + (1 - w_sky) * lum
SC = float(OPT["skychroma"])
if SC > 0:
    chroma_s = np.dstack([nd.gaussian_filter(chroma[..., c_], px(SC)) for c_ in range(3)])
    chroma = w_sky[..., None] * chroma_s + (1 - w_sky[..., None]) * chroma
rgb = lum[..., None] + chroma

# --- 3. stretch ---
black = -float(OPT['black']) * sigma
x = np.clip(lum - black, 1e-12, None)
core = nd.gaussian_filter(lum, px(2))[gal_mask].max()                 # brightest part of the galaxy
top = max(core * 1.15, np.percentile(x, 99.9))


def asinh(vv, a, t):
    return np.arcsinh(vv / a) / np.arcsinh(t / a)


L = np.clip(asinh(x, 4 * sigma, top), 0, 1)
L_sky = np.median(L[sky])
# put the sky at ~0.10 with a gentle toe, keep the top fixed
L = np.clip((L - L_sky * 0.6) / (1 - L_sky * 0.6), 0, 1)
L_sky = np.median(L[sky])
gamma = np.log(float(OPT['sky'])) / np.log(max(L_sky, 1e-3))
L = L ** np.clip(gamma, 0.6, 1.6)
A = float(OPT["skyflat"])
if A > 0:
    L_sky = np.median(L[sky])
    sL = 1.4826 * np.median(np.abs(L[sky] - L_sky))
    uu = (L - L_sky) / (3 * sL)
    L = np.clip(L - min(A, 0.9) * (L - L_sky) * np.exp(-uu ** 2), 0, 1)

out = np.clip(rgb - black, 0, None) * (L / x)[..., None]
l2 = out.mean(axis=2, keepdims=True)
out = np.clip(l2 + 1.35 * (out - l2), 0, None)                   # saturation, hue kept
m = out.max(axis=2, keepdims=True)
out = np.where(m > 1, out / m, out)

# --- 4. orientation: north up, east left ---
img = out[::-1]                                                    # display rows top-down
gy_d = H - 1 - gy
pn = wcs.world_to_pixel(gal.directional_offset_by(0 * au.deg, 0.02 * au.deg))
pe = wcs.world_to_pixel(gal.directional_offset_by(90 * au.deg, 0.02 * au.deg))
north = np.array([pn[0] - gx, -(pn[1] - gy)])                      # (dx, dy) in display coords
east = np.array([pe[0] - gx, -(pe[1] - gy)])
cx, cy = gx, gy_d
best = None
for k in range(4):
    for flip in (False, True):
        def tr(vec, k=k, flip=flip):
            dx, dy = vec
            for _ in range(k):                                     # np.rot90 = counter-clockwise
                dx, dy = dy, -dx
            return np.array([-dx if flip else dx, dy])
        nv, ev = tr(north), tr(east)
        score = -nv[1] / np.linalg.norm(nv) - ev[0] / np.linalg.norm(ev)   # up = -y, left = -x
        if best is None or score > best[0]:
            best = (score, k, flip)
_, k, flip = best
Hd, Wd = img.shape[:2]
img = np.rot90(img, k)
for _ in range(k):
    cx, cy = cy, Wd - 1 - cx
    Hd, Wd = Wd, Hd
if flip:
    img = img[:, ::-1]
    cx = Wd - 1 - cx
print(f"orientation: rot90 x{k}, mirror {flip}")


def save(im, fname, upscale=1):
    p = Image.fromarray((np.clip(im, 0, 1) * 255 + 0.5).astype(np.uint8))
    if upscale != 1:
        p = p.resize((p.width * upscale, p.height * upscale), Image.LANCZOS)
    p.save(os.path.join(OUT, fname))
    print("saved", os.path.join(OUT, fname))


save(img, f"{tag}_field.png")
half = int(a_px * 1.6)
y0, x0 = int(cy) - half, int(cx) - half
y0, x0 = max(0, min(y0, img.shape[0] - 2 * half)), max(0, min(x0, img.shape[1] - 2 * half))
save(img[y0:y0 + 2 * half, x0:x0 + 2 * half], f"{tag}_crop_{2 * BIN}x.png", 2 * BIN)

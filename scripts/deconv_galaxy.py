"""Richardson-Lucy deconvolution of a linear colour stack, applied to luminance where the signal is strong.

Run:  ~/iss/.venv/bin/python -I deconv_galaxy.py <linear.fits> <out.fits> [iterations=20] [tv=0.002] [--gal=RA,DEC,major',minor',PA]
      e.g. deconv_galaxy.py result_ngc5866_flat_linear_cc.fits out.fits 20 0.002 --gal=15:06:29.5,+55:45:48,6.4,2.8,128
Input: linear RGB FITS (3, H, W); header copied unchanged. --gal needs a WCS in the header: peaks inside that
ellipse (0.45 x major by 0.6 x minor semi-axes) are NOT protected (galaxy knots and the core); every other star is.

Steps:
 1. PSF measured from isolated, unsaturated stars (sub-pixel centred, median-combined)
 2. RL with total-variation regularisation (Dey et al. 2006) on the luminance
 3. result blended in only where the smoothed signal is well above the sky noise, and never on
    saturated stars; field stars masked (deconvolving them leaves dark rings around them)
 4. colour kept: every channel is scaled by the same factor L_deconv / L
"""
import sys

import numpy as np
from astropy.io import fits
from scipy import ndimage as nd
from scipy.signal import fftconvolve

GAL = None
argv = []
for a in sys.argv[1:]:
    if a.startswith("--gal="):
        GAL = a[6:].split(",")
    else:
        argv.append(a)
fn, out_fn = argv[0], argv[1]
ITERS = int(argv[2]) if len(argv) > 2 else 20
TV = float(argv[3]) if len(argv) > 3 else 0.002

hdr = fits.getheader(fn)
rgb = fits.getdata(fn).astype(np.float64)           # (3, H, W)
L = rgb.mean(axis=0)
H, W = L.shape
bg = np.median(L)
sig = 1.4826 * np.median(np.abs(L - bg))

# --- 1. PSF ---
s = nd.gaussian_filter(L, 1.5)
sat_level = 0.85 * rgb.max()
saturated = nd.binary_dilation(rgb.max(axis=0) > sat_level, iterations=12)


def select_stars(R, iso, cen, thr):
    """Stamps (2R+1 square) of bright (> thr sky sigma), unsaturated stars with nothing brighter than 5 % of
    the peak beyond `iso` px from the centre, centred on their centroid (within `cen` px)."""
    peaks = (s == nd.maximum_filter(s, 2 * R + 1)) & (s - bg > thr * sig)
    yy, xx = np.mgrid[-R:R + 1, -R:R + 1]
    rr = np.hypot(yy, xx)
    out, where = [], []
    for y, x in zip(*np.nonzero(peaks)):
        if not (2 * R < y < H - 2 * R and 2 * R < x < W - 2 * R) or saturated[y, x]:
            continue
        cut = L[y - R:y + R + 1, x - R:x + R + 1] - bg
        if (cut[rr > iso] > 0.05 * cut[R, R]).any():
            continue
        w = np.clip(cut, 0, None) * (rr < cen)
        cy, cx = (w * yy).sum() / w.sum(), (w * xx).sum() / w.sum()
        cut = nd.shift(cut, (-cy, -cx), order=3, mode="constant")
        out.append(cut / cut.sum())
        where.append((y, x))
    return out, where


# standard selection (what the 25 Oct images used); if it finds too few stars - wider stars on a night of
# poorer seeing, or a shallower stack - size the stamps to the stars and lower the brightness threshold
R = 15
stamps, where = select_stars(R, 12, 10, 60)
if len(stamps) < 5:
    pk = (s == nd.maximum_filter(s, 31)) & (s - bg > 20 * sig) & ~saturated
    widths = []
    for y, x in sorted(zip(*np.nonzero(pk)), key=lambda p: -s[p])[:20]:
        if 20 < y < H - 20 and 20 < x < W - 20:
            c = np.clip(L[y - 15:y + 16, x - 15:x + 16] - bg, 0, None)
            g = np.mgrid[-15:16, -15:16]
            widths.append(2.3548 * np.sqrt((c * (g[0] ** 2 + g[1] ** 2)).sum() / c.sum() / 2))
    fw = float(np.median(widths))
    R = int(np.ceil(max(15, 2 * fw)))
    for thr in (60, 40, 25, 15):
        stamps, where = select_stars(R, max(12, 1.4 * fw), max(10, 1.2 * fw), thr)
        if len(stamps) >= 5:
            break
    print(f"few stars at the standard settings: star width ~{fw:.1f} px -> stamp {2 * R + 1} px, threshold {thr} sigma")
ys, xs = [p[0] for p in where], [p[1] for p in where]
if len(stamps) < 3:
    print(f"WARNING: PSF from only {len(stamps)} stars - check the result for artefacts")
psf = np.median(np.array(stamps), axis=0)
yy, xx = np.mgrid[-R:R + 1, -R:R + 1]
rr = np.hypot(yy, xx)
psf -= np.median(psf[rr > R - 2])
psf = np.clip(psf, 0, None) * (rr <= R - 1)
psf /= psf.sum()


def fwhm(img):
    """FWHM from the second moment of a centred stamp (Gaussian equivalent)."""
    m = np.clip(img, 0, None)
    r2 = (m * rr ** 2).sum() / m.sum()
    return 2.3548 * np.sqrt(r2 / 2)


print(f"PSF from {len(stamps)} stars, FWHM ~{fwhm(psf):.2f} px")

# --- 2. RL + TV on luminance ---
ped = max(0.0, -L.min()) + 3 * sig                   # keep everything positive
obs = L + ped
est = obs.copy()
psf_f = psf[::-1, ::-1]
eps = 1e-12
for i in range(ITERS):
    conv = fftconvolve(est, psf, mode="same")
    est = est * fftconvolve(obs / np.maximum(conv, eps), psf_f, mode="same")
    if TV > 0:
        gy, gx = np.gradient(est)
        norm = np.sqrt(gx ** 2 + gy ** 2) + 1e-3 * sig
        div = np.gradient(gx / norm, axis=1) + np.gradient(gy / norm, axis=0)
        est = est / np.maximum(1 - TV * div, 0.1)
dec = est - ped

# --- 3. where to use it ---
snr = nd.gaussian_filter(L - bg, 3) / sig
w = np.clip((snr - 3) / 12, 0, 1)
w = np.where(saturated, 0, w)

# star protection: a disc per star, out to where its profile falls to ~1.5 sky sigma
gal_xy = None
if GAL:
    from astropy.wcs import WCS
    from astropy.coordinates import SkyCoord
    from astropy.wcs.utils import proj_plane_pixel_scales
    wcs = WCS(hdr, naxis=2)
    import astropy.units as au
    c = SkyCoord(GAL[0], GAL[1], unit=("hourangle", "deg"))
    gx, gy = wcs.world_to_pixel(c)
    ax_, ay_ = wcs.world_to_pixel(c.directional_offset_by(float(GAL[4]) * au.deg, 0.01 * au.deg))
    th = np.arctan2(ay_ - gy, ax_ - gx)
    scale = proj_plane_pixel_scales(wcs).mean() * 3600
    ga, gb = 0.45 * float(GAL[2]) * 60 / scale, 0.6 * float(GAL[3]) * 60 / scale
    gal_xy = (gx, gy, ga, gb, th)
psf_sigma = fwhm(psf) / 2.3548
sp = (s == nd.maximum_filter(s, 9)) & (s - bg > 6 * sig)
star_mask = np.zeros_like(L)
yy_full, xx_full = np.ogrid[0:H, 0:W]
n_prot = 0
for y, x in zip(*np.nonzero(sp)):
    if gal_xy:
        dx, dy = x - gal_xy[0], y - gal_xy[1]
        u = dx * np.cos(gal_xy[4]) + dy * np.sin(gal_xy[4])
        v = -dx * np.sin(gal_xy[4]) + dy * np.cos(gal_xy[4])
        if (u / gal_xy[2]) ** 2 + (v / gal_xy[3]) ** 2 < 1:
            continue
    amp = s[y, x] - bg
    r = psf_sigma * np.sqrt(2 * np.log(max(amp / (1.5 * sig), 1.01))) + 3
    if saturated[y, x]:
        r *= 1.5
    r = min(r, 60)
    y0, y1, x0, x1 = max(0, int(y - r - 1)), min(H, int(y + r + 2)), max(0, int(x - r - 1)), min(W, int(x + r + 2))
    d = np.hypot(yy_full[y0:y1] - y, xx_full[:, x0:x1] - x)
    star_mask[y0:y1, x0:x1] = np.maximum(star_mask[y0:y1, x0:x1], (d <= r).astype(float))
    n_prot += 1
star_mask = np.clip(nd.gaussian_filter(star_mask, 2) * 1.5, 0, 1)
w = w * (1 - star_mask)
print(f"protected {n_prot} stars" + (f" (galaxy ellipse {gal_xy[2]:.0f} x {gal_xy[3]:.0f} px left open)" if gal_xy else ""))
w = nd.gaussian_filter(w, 2)
L_out = w * dec + (1 - w) * L

# --- 4. apply to colour ---
ratio = np.where(w > 0.01, (L_out - bg + 3 * sig) / np.maximum(L - bg + 3 * sig, eps), 1.0)
ratio = np.clip(ratio, 0.2, 8)
out = rgb * ratio[None]

# star sharpness before / after (same isolated stars)
before, after = [], []
for y, x in zip(ys, xs):
    if not (2 * R < y < H - 2 * R and 2 * R < x < W - 2 * R) or saturated[y, x] or w[y, x] < 0.9:
        continue
    for img, lst in ((L, before), (out.mean(axis=0), after)):
        c = img[y - R:y + R + 1, x - R:x + R + 1] - bg
        lst.append(fwhm(c * (rr < 10)))
if before:
    print(f"star FWHM in the deconvolved zone: {np.median(before):.2f} -> {np.median(after):.2f} px "
          f"({len(before)} stars)")
fits.writeto(out_fn, out.astype(np.float32), hdr, overwrite=True)
print("saved", out_fn)

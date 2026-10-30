"""Remove vignetting and gradient from a linear RGB stack taken without flats, keeping a large galaxy's faint halo.

Run:  ~/iss/.venv/bin/python -I vignette_bg.py <in.fits> <out.fits> --gal=RA,DEC,major',minor',PA [--galmask=1.3]
      [--cmax=0.4]  (how far the vignetting centre may sit from the frame centre, in frame widths; 0 = fixed)
Input: linear RGB FITS (3, H, W) with a WCS (plate solved). Output: float32 0..1, header kept.

Siril's polynomial subsky cannot follow strong vignetting and, when the galaxy fills much of the field, some of its
samples sit on the halo: it then digs a dark moat around the galaxy. Here the background is a physical model with
few free parameters, fitted on sky only:
    B(x, y) = c0 + c1*r^2 + c2*r^4 + c3*r^6 + gx*x + gy*y,   r = distance from the vignetting centre
The centre is fitted on the luminance, then each channel gets its own coefficients. Sky samples: 16 x 16 px block
medians outside the galaxy ellipse (galmask x the catalogue axes) with stars and residual halo removed by an
iterative 2.5-sigma clip. The model minus its median is subtracted, so the sky level stays where it was.
"""
import sys

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import proj_plane_pixel_scales
import astropy.units as au
from scipy.optimize import least_squares

opts = {"gal": None, "galmask": "1.3", "cmax": "0.4"}
args = []
for a in sys.argv[1:]:
    if a.startswith("--") and a[2:].split("=")[0] in opts:
        k, v = a[2:].split("=", 1)
        opts[k] = v
    else:
        args.append(a)
fn_in, fn_out = args
if not opts["gal"]:
    sys.exit("--gal=RA,DEC,major,minor,PA is required")

hdr = fits.getheader(fn_in)
raw = fits.getdata(fn_in)
img = raw.astype(np.float64)
if np.issubdtype(raw.dtype, np.integer):
    img /= 65535.0
del raw
_, H, W = img.shape

# --- galaxy ellipse in pixels
ra, dec, major, minor, pa = opts["gal"].split(",")
wcs = WCS(hdr, naxis=2)
scale = proj_plane_pixel_scales(wcs).mean() * 3600
gal = SkyCoord(ra, dec, unit=("hourangle", "deg"))
gx, gy = wcs.world_to_pixel(gal)
ax, ay = wcs.world_to_pixel(gal.directional_offset_by(float(pa) * au.deg, 0.01 * au.deg))
th = np.arctan2(ay - gy, ax - gx)
k = float(opts["galmask"])
a_px, b_px = k * float(major) * 60 / scale / 2, k * float(minor) * 60 / scale / 2

# --- block medians (16 px) as sky samples
B = 16
h, w = H // B, W // B
blocks = img[:, :h * B, :w * B].reshape(3, h, B, w, B)
med = np.median(blocks, axis=(2, 4))                      # (3, h, w)
yc, xc = (np.mgrid[0:h, 0:w] + 0.5) * B - 0.5
u = (xc - gx) * np.cos(th) + (yc - gy) * np.sin(th)
v = -(xc - gx) * np.sin(th) + (yc - gy) * np.cos(th)
edge = 2                                                  # stacking edges
ok = ((u / a_px) ** 2 + (v / b_px) ** 2 > 1)
ok[:edge] = ok[-edge:] = False
ok[:, :edge] = ok[:, -edge:] = False
lum = med.sum(0)
xs, ys = xc / W - 0.5, yc / H - 0.5                       # normalised coordinates


def design(cx, cy, x, y):
    r2 = ((x - cx) ** 2 + ((y - cy) * H / W) ** 2)
    return np.stack([np.ones_like(x), r2, r2 ** 2, r2 ** 3, x, y], axis=-1)


def fit_lin(cx, cy, data, m):
    A = design(cx, cy, xs[m], ys[m])
    c, *_ = np.linalg.lstsq(A, data[m], rcond=None)
    return c


# iterative clip on the luminance with the centre fixed at the frame centre, then fit the centre
m = ok.copy()
for _ in range(5):
    c = fit_lin(0.0, 0.0, lum, m)
    res = lum - design(0.0, 0.0, xs, ys) @ c
    s = 1.4826 * np.median(np.abs(res[m] - np.median(res[m])))
    m = ok & (np.abs(res) < 2.5 * s)
cmax = float(opts["cmax"])
if cmax > 0:
    f = least_squares(lambda p: lum[m] - design(p[0], p[1], xs[m], ys[m]) @ fit_lin(p[0], p[1], lum, m), [0.0, 0.0],
                      bounds=([-cmax, -cmax], [cmax, cmax]))
    cx, cy = f.x
else:
    cx, cy = 0.0, 0.0
for _ in range(3):                                        # re-clip around the final centre
    c = fit_lin(cx, cy, lum, m)
    res = lum - design(cx, cy, xs, ys) @ c
    s = 1.4826 * np.median(np.abs(res[m] - np.median(res[m])))
    m = ok & (np.abs(res) < 2.5 * s)
print(f"vignetting centre ({(cx + 0.5) * W:.0f}, {(cy + 0.5) * H:.0f}) px, {m.sum()} sky samples "
      f"({100 * m.sum() / m.size:.0f} % of blocks), galaxy ellipse {2 * a_px:.0f} x {2 * b_px:.0f} px")

# --- per channel model at full resolution, subtracted keeping the median level
yy, xx = np.mgrid[0:H, 0:W]
XF, YF = xx / W - 0.5, yy / H - 0.5
out = np.empty_like(img, dtype=np.float32)
for ch in range(3):
    c = fit_lin(cx, cy, med[ch], m)
    model = design(cx, cy, XF, YF) @ c
    res = med[ch][m] - design(cx, cy, xs[m], ys[m]) @ c
    print(f"  channel {ch}: model range {model.min():.5f}..{model.max():.5f}, residual rms {res.std():.2e}")
    out[ch] = img[ch] - (model - np.median(model))
hdr["HISTORY"] = "vignette_bg.py: radial vignetting + gradient model subtracted"
fits.writeto(fn_out, out, hdr, overwrite=True)
print(f"saved {fn_out}")

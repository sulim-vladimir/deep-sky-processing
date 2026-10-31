"""Second background stage: removes what the first one (Siril's subsky or vignette_bg.py) leaves, i.e. a steeper
fall-off at the edges and low-level patches, while keeping the galaxies' faint halos.

Run:  ~/iss/.venv/bin/python -I bg_residual.py <in.fits> <out.fits> --gal=RA,DEC,major',minor',PA [--gal=...]
      [--smooth=160] [--galmask=2.0]
Input: linear RGB FITS (3, H, W) with a WCS (plate solved). Output: float32, header kept.

Sky samples: 32 x 32 px block medians, stars masked (3 sigma, dilated 10 px), blocks inside any galaxy ellipse
(galmask x the catalogue axes, each semi-axis at least 0.75 x smooth) left out, blocks that stand out from their
neighbourhood clipped at 2.5 sigma (star halos, small galaxies). The gaps are filled by normalised Gaussian
convolution at --smooth px, so the model has no structure smaller than that: a halo is never flattened outside the
masks either. The model minus its median is subtracted per channel, so the sky level stays where it was.
On the Draco Trio (31 Oct 2026), after BG=radial: 3x3 sky zones from +-0.05 to +-0.01 pixel-sigma.
"""
import sys

import numpy as np
import astropy.units as au
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from scipy import ndimage as nd

opts = {"smooth": "160", "galmask": "2.0"}
gals, args = [], []
for a in sys.argv[1:]:
    if a.startswith("--gal="):
        gals.append(a[6:])
    elif a.startswith("--") and a[2:].split("=")[0] in opts:
        k, v = a[2:].split("=", 1)
        opts[k] = v
    else:
        args.append(a)
fn_in, fn_out = args
smooth, gm = float(opts["smooth"]), float(opts["galmask"])

hdr = fits.getheader(fn_in)
raw = fits.getdata(fn_in)
img = raw.astype(np.float32)
if np.issubdtype(raw.dtype, np.integer):
    img /= 65535.0
del raw
_, H, W = img.shape
wcs = WCS(hdr, naxis=2)
scale = proj_plane_pixel_scales(wcs).mean() * 3600

B = 32
gh, gw = H // B, W // B
yy, xx = np.mgrid[0:gh, 0:gw] * B + B / 2
gal = np.zeros((gh, gw), bool)
for g in gals:
    ra, dec, major, minor, pa = g.split(",")
    c = SkyCoord(ra, dec, unit=("hourangle", "deg"))
    cx, cy = wcs.world_to_pixel(c)
    ax, ay = wcs.world_to_pixel(c.directional_offset_by(float(pa) * au.deg, 0.01 * au.deg))
    t = np.arctan2(ay - cy, ax - cx)
    u = (xx - cx) * np.cos(t) + (yy - cy) * np.sin(t)
    v = -(xx - cx) * np.sin(t) + (yy - cy) * np.cos(t)
    sa = max(gm * float(major) * 30 / scale, 0.75 * smooth)   # semi-axes in px (a thin edge-on still gets a wide mask)
    sb = max(gm * float(minor) * 30 / scale, 0.75 * smooth)
    gal |= (u / sa) ** 2 + (v / sb) ** 2 < 1

L = img.mean(axis=0)
med = np.median(L)
sig = 1.4826 * np.median(np.abs(L - med))
star = nd.binary_dilation(nd.gaussian_filter(L, 2) > med + 3 * sig, iterations=10)
out = img.copy()
for k in range(3):
    ch = img[k]
    blocks = np.full((gh, gw), np.nan, np.float32)
    for j in range(gh):
        for i in range(gw):
            if gal[j, i]:
                continue
            s = ch[j * B:(j + 1) * B, i * B:(i + 1) * B][~star[j * B:(j + 1) * B, i * B:(i + 1) * B]]
            if s.size > 0.4 * B * B:
                blocks[j, i] = np.median(s)
    for _ in range(3):                                         # clip blocks that stand out from their surroundings
        ok = ~np.isnan(blocks)
        f = nd.gaussian_filter(np.where(ok, blocks, 0), 3, mode="nearest")
        n = nd.gaussian_filter(ok.astype(np.float32), 3, mode="nearest")
        dev = blocks - f / np.maximum(n, 1e-6)
        s = 1.4826 * np.nanmedian(np.abs(dev))
        blocks = np.where(np.abs(dev) > 2.5 * s, np.nan, blocks)
    ok = ~np.isnan(blocks)
    f = nd.gaussian_filter(np.where(ok, blocks, 0), smooth / B, mode="nearest")
    n = nd.gaussian_filter(ok.astype(np.float32), smooth / B, mode="nearest")
    model = f / np.maximum(n, 1e-6)
    full = nd.zoom(model, (H / gh, W / gw), order=3, mode="nearest")[:H, :W]
    full = np.pad(full, ((0, H - full.shape[0]), (0, W - full.shape[1])), mode="edge")
    out[k] = ch - (full - np.median(full))
    print(f"  channel {k}: {ok.sum()} sky blocks of {gh * gw}, model range "
          f"{(model.min() - np.median(model)) / sig:+.2f} .. {(model.max() - np.median(model)) / sig:+.2f} pixel-sigma")
hdr["HISTORY"] = f"bg_residual.py: smooth background residual ({smooth:.0f} px) subtracted"
fits.writeto(fn_out, out, hdr, overwrite=True)
print(f"saved {fn_out}")

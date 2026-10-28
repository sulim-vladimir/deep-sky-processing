"""Give a new (unsolved) image of the same field the WCS of an earlier solved one.

Finds the flip (none / vertical) and sub-pixel shift by phase correlation of the luminance,
resamples the new image onto the old pixel grid and saves it with the old header.
Run: ~/iss/.venv/bin/python -I align_wcs.py <solved_ref.fits> <new.fits> <out.fits>
"""
import sys
import numpy as np
from astropy.io import fits
from scipy import ndimage as nd

ref_fn, new_fn, out_fn = sys.argv[1:4]
ref_h = fits.getheader(ref_fn)
ref = fits.getdata(ref_fn).astype(np.float64)
new = fits.getdata(new_fn).astype(np.float64)


def prep(a):
    l = a.mean(axis=0)
    l = l - nd.median_filter(l, 31)
    return np.clip(l, 0, np.percentile(l, 99.9))


def xcorr(a, b):
    F = np.fft.fft2(a) * np.conj(np.fft.fft2(b))
    c = np.fft.ifft2(F / (np.abs(F) + 1e-12)).real
    y, x = np.unravel_index(c.argmax(), c.shape)
    peak = c.max() / c.std()
    # sub-pixel by parabola
    def sub(m, p, q):
        d = (m - 2 * p + q)
        return 0.5 * (m - q) / d if d != 0 else 0.0
    dy = y + sub(c[y - 1, x], c[y, x], c[(y + 1) % c.shape[0], x])
    dx = x + sub(c[y, x - 1], c[y, x], c[y, (x + 1) % c.shape[1]])
    if dy > c.shape[0] / 2: dy -= c.shape[0]
    if dx > c.shape[1] / 2: dx -= c.shape[1]
    return peak, dy, dx


R = prep(ref)
best = None
for flip in (False, True):
    N = prep(new[:, ::-1] if flip else new)
    pk, dy, dx = xcorr(R, N)
    print(f"flip={flip}: correlation peak {pk:.1f}, shift dy={dy:.2f} dx={dx:.2f}")
    if best is None or pk > best[0]:
        best = (pk, flip, dy, dx)
pk, flip, dy, dx = best
src = new[:, ::-1] if flip else new
out = np.stack([nd.shift(src[c], (dy, dx), order=3, mode="nearest") for c in range(src.shape[0])]).astype(np.float32)
# check: residual shift should be ~0
print("residual:", np.round(xcorr(R, prep(out))[1:], 2))
fits.writeto(out_fn, out, ref_h, overwrite=True)
print("saved", out_fn)

"""Give a new (unsolved) image of the same field the WCS of an earlier solved one.

Finds the flip (none / vertical), a small rotation (up to 5 deg, searched only when needed) and the sub-pixel
shift by phase correlation of the luminance,
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
        best = (pk, flip, 0.0, dy, dx)
# a weak peak means the field is also rotated (camera turned between nights): search small rotations,
# coarse then fine
if best[0] < 15:
    flip = best[1]
    N0 = prep(new[:, ::-1] if flip else new)
    for step, centre, span in ((0.25, 0.0, 5.0), (0.05, None, 0.25)):
        c0 = best[2] if centre is None else centre
        for ang in np.arange(c0 - span, c0 + span + 1e-9, step):
            pk, dy, dx = xcorr(R, nd.rotate(N0, ang, reshape=False, order=1))
            if pk > best[0]:
                best = (pk, flip, float(ang), dy, dx)
    print(f"rotation search: {best[2]:+.2f} deg, correlation peak {best[0]:.1f}, shift dy={best[3]:.2f} dx={best[4]:.2f}")
pk, flip, ang, dy, dx = best
src = new[:, ::-1] if flip else new
if ang:
    src = np.stack([nd.rotate(src[c], ang, reshape=False, order=3, mode="nearest") for c in range(src.shape[0])])
out = np.stack([nd.shift(src[c], (dy, dx), order=3, mode="nearest") for c in range(src.shape[0])]).astype(np.float32)
# check: residual shift should be ~0
print("residual:", np.round(xcorr(R, prep(out))[1:], 2))
fits.writeto(out_fn, out, ref_h, overwrite=True)
print("saved", out_fn)

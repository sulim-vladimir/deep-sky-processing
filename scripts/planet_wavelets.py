"""RegiStax-style finishing of a stacked planetary image (e.g. AutoStakkert's TIF/PNG output), from the command line.

Run:  ~/iss/.venv/bin/python -I planet_wavelets.py <stack.tif|png|fits> [options]
      e.g. planet_wavelets.py jupiter_AS_P25.tif --gains=5,5,3.2,2,1.3,1 --denoise=2.5,1.8,0.8
      Defaults tuned on Saturn, 10 Sep 2026 (ASI290MC, 4050 mm, 0.15"/px): they match the per-layer detail of
      the RegiStax result of that night. At this sampling the detail sits in layers 3-5 (4-16 px).
Output next to the input: <name>_wav.png (16-bit) and <name>_wav_noden.png (same, without wavelet denoise).
--grid writes <name>_grid.png: the same image with weaker / chosen / stronger sharpening side by side.

Steps (all classical, no AI):
 1. RGB align (--rgbalign, default on): R and B shifted onto G by sub-pixel phase correlation on the planet
    (atmospheric dispersion; what RegiStax's "RGB align" does).
 2. Wavelets: a trous B3-spline decomposition into 6 layers (scales 1, 2, 4 ... 32 px), like RegiStax's 6 sliders
    in "dyadic, default" mode. Each layer is multiplied by its gain (--gains, 1 = unchanged); before that, each
    layer is soft-thresholded at K x its own noise (--denoise, per layer from the fine end; RegiStax's "denoise"
    boxes). Applied to the luminance only by default (--lum=1) so colour noise is not amplified; --lum=0 does each
    channel separately, as RegiStax does.
 3. Colour balance (--balance=none|disk): none keeps the stack's colour (default); disk makes the planet's mean
    R, G, B equal (grey world), which turns Saturn and Jupiter too white.
 4. Levels: black point = sky background, white point = 99.95th percentile of the planet, then gamma (--gamma).
"""
import os
import sys

import numpy as np
from scipy import ndimage as nd
from PIL import Image

OPT = {"gains": "5,5,3.2,2.0,1.3,1.0", "denoise": "2.5,1.8,0.8,0,0,0", "lum": "1", "rgbalign": "1",
       "balance": "none", "gamma": "1.0", "grid": "0"}
args = []
for a in sys.argv[1:]:
    if a.startswith("--") and a[2:].split("=")[0] in OPT:
        k, v = (a[2:].split("=", 1) + ["1"])[:2]
        OPT[k] = v
    else:
        args.append(a)
if len(args) != 1:
    sys.exit(__doc__)
fn = args[0]
base = os.path.splitext(fn)[0]


def load(path):
    if path.lower().endswith((".fit", ".fits", ".fts")):
        from astropy.io import fits
        a = fits.getdata(path).astype(np.float64)
        if a.ndim == 3 and a.shape[0] in (3, 4):
            a = np.moveaxis(a[:3], 0, -1)
        a = a[::-1]                                       # FITS row 0 = bottom
    else:
        im = Image.open(path)
        a = np.asarray(im).astype(np.float64)
        if a.ndim == 3:
            a = a[..., :3]
    a -= a.min()
    return a / a.max()


img = load(fn)
mono = img.ndim == 2
if mono:
    img = img[..., None]
lum0 = img.mean(axis=2)

# planet mask: above a threshold between sky and disk, the largest blob, slightly grown
sky = np.percentile(lum0, 5)
thr = sky + 0.15 * (np.percentile(lum0, 99.5) - sky)
lab, n = nd.label(nd.gaussian_filter(lum0, 2) > thr)
if n == 0:
    sys.exit("no planet found")
disk = lab == (np.argmax(nd.sum(np.ones_like(lum0), lab, range(1, n + 1))) + 1)
disk = nd.binary_fill_holes(disk)
print(f"{os.path.basename(fn)}: {img.shape[1]} x {img.shape[0]}, {'mono' if mono else 'RGB'}, "
      f"planet {int(disk.sum() ** 0.5 * 1.128)} px across")


# ---- 1. RGB align
MAXSHIFT = 8                                            # px: dispersion is a few px; larger peaks are noise


def subpixel_shift(a, b):
    """shift that moves b onto a (dy, dx): cross-correlation of the edges, peak searched within MAXSHIFT px,
    refined with a parabola (plain, not phase, correlation: steadier on a noisy blue channel)"""
    def edges(z):
        z = nd.gaussian_filter(z, 1.5)
        return np.hypot(nd.sobel(z, 0), nd.sobel(z, 1))
    w = np.outer(np.hanning(a.shape[0]), np.hanning(a.shape[1]))
    ea, eb = edges(a), edges(b)
    A, B = np.fft.fft2((ea - ea.mean()) * w), np.fft.fft2((eb - eb.mean()) * w)
    r = np.real(np.fft.ifft2(A * np.conj(B)))
    H0, W0 = r.shape
    yy, xx = np.ogrid[:H0, :W0]
    far = (np.minimum(yy, H0 - yy) > MAXSHIFT) | (np.minimum(xx, W0 - xx) > MAXSHIFT)
    r = np.where(far, r.min(), r)
    y, x = np.unravel_index(np.argmax(r), r.shape)
    def par(m, c, p):
        d = m - 2 * c + p
        return 0.0 if d == 0 else 0.5 * (m - p) / d
    H, W = r.shape
    dy = y + par(r[(y - 1) % H, x], r[y, x], r[(y + 1) % H, x])
    dx = x + par(r[y, (x - 1) % W], r[y, x], r[y, (x + 1) % W])
    return (dy + H / 2) % H - H / 2, (dx + W / 2) % W - W / 2


if not mono and OPT["rgbalign"] == "1":
    ys, xs = np.where(disk)
    pad = 20
    sl = (slice(max(ys.min() - pad, 0), ys.max() + pad), slice(max(xs.min() - pad, 0), xs.max() + pad))
    for c, name in ((0, "R"), (2, "B")):
        dy, dx = subpixel_shift(img[sl][..., 1], img[sl][..., c])
        img[..., c] = nd.shift(img[..., c], (dy, dx), order=3, mode="nearest")
        print(f"  {name} shifted onto G by ({dy:+.2f}, {dx:+.2f}) px")


# ---- 2. wavelets
def atrous(a, nlev=6):
    k = np.array([1, 4, 6, 4, 1], float) / 16
    layers, c = [], a
    for j in range(nlev):
        kj = np.zeros(4 * 2 ** j + 1); kj[::2 ** j] = k
        s = nd.convolve1d(nd.convolve1d(c, kj, axis=0, mode="reflect"), kj, axis=1, mode="reflect")
        layers.append(c - s)
        c = s
    return layers, c


def sharpen(a, gains, den):
    layers, resid = atrous(a, len(gains))
    out = resid.copy()
    for j, (w, g) in enumerate(zip(layers, gains)):
        k = den[j] if j < len(den) else 0.0
        if k > 0:
            s = 1.4826 * np.median(np.abs(w[~disk_grown] - np.median(w[~disk_grown]))) if (~disk_grown).sum() > 500 \
                else 1.4826 * np.median(np.abs(w - np.median(w)))
            w = np.sign(w) * np.maximum(np.abs(w) - k * s, 0)
        out += g * w
    return out


disk_grown = nd.binary_dilation(disk, iterations=10)
gains = [float(v) for v in OPT["gains"].split(",")]
den = [float(v) for v in OPT["denoise"].split(",")]


def process(gains, den):
    x = img.copy()
    if OPT["lum"] == "1" and not mono:
        L = x.mean(axis=2)
        Ls = sharpen(L, gains, den)
        x = x * (np.clip(Ls, 1e-6, None) / np.clip(L, 1e-6, None))[..., None]
        x = np.where(L[..., None] > 1e-4, x, Ls[..., None])
    else:
        x = np.dstack([sharpen(x[..., c], gains, den) for c in range(x.shape[2])])
    # ---- 3. colour balance on the disk
    if not mono and OPT["balance"] == "disk":
        m = np.array([x[..., c][disk].mean() for c in range(3)])
        x = x * (m.mean() / m)
    # ---- 4. levels + gamma
    bg = np.median(x[~disk_grown]) if (~disk_grown).sum() > 500 else 0.0
    top = np.percentile(x[disk], 99.95)
    x = np.clip((x - bg) / (top - bg), 0, 1) ** (1 / float(OPT["gamma"]))
    return x[..., 0] if mono else x


def save(x, path):
    a = (np.clip(x, 0, 1) * 65535 + 0.5).astype(np.uint16)
    if a.ndim == 2:
        Image.fromarray(a).save(path)
    else:                                               # PIL cannot write 16-bit RGB PNG: use OpenCV if present
        try:
            import cv2
            cv2.imwrite(path, a[..., ::-1])
        except ImportError:
            Image.fromarray((a >> 8).astype(np.uint8)).save(path)
    print(f"  saved {path}")


save(process(gains, den), base + "_wav.png")
save(process(gains, [0.0] * len(gains)), base + "_wav_noden.png")

if OPT["grid"] == "1":
    from PIL import ImageDraw, ImageFont
    ys, xs = np.where(disk)
    r = int(max(ys.ptp(), xs.ptp()) * 0.65)
    cy, cx = int(ys.mean()), int(xs.mean())
    tiles = []
    for f, label in ((0.5, "weaker"), (1.0, "chosen"), (1.6, "stronger")):
        g = [1 + f * (v - 1) for v in gains]
        x = process(g, den)[max(cy - r, 0):cy + r, max(cx - r, 0):cx + r]
        t = Image.fromarray((np.clip(x, 0, 1) * 255).astype(np.uint8)).convert("RGB")
        t = t.resize((max(t.width, 400), max(t.height, 400) * t.height // t.height), Image.LANCZOS) if t.width < 400 else t
        d = ImageDraw.Draw(t)
        try:
            font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", max(16, t.width // 18))
        except OSError:
            font = ImageFont.load_default()
        d.text((8, 6), f"{label}: {','.join(f'{v:.2f}' for v in g)}", font=font, fill=(255, 230, 80))
        tiles.append(t)
    G = Image.new("RGB", (sum(t.width for t in tiles) + 20, max(t.height for t in tiles)), (255, 255, 255))
    xo = 0
    for t in tiles:
        G.paste(t, (xo, 0)); xo += t.width + 10
    G.save(base + "_grid.png"); print(f"  saved {base}_grid.png")

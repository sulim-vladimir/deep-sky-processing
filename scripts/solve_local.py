"""Plate solve a linear RGB FITS offline with astrometry.net and write the solution (WCS) into its header.

Run:  python -I solve_local.py <image.fits> <ra> <dec> [arcsec_per_px=0.80] [index_dir=~/.local/share/astrometry]
      e.g. solve_local.py ngc6503_linear.fits 17:49:26.4 +70:08:40

Needs `solve-field` (apt: astrometry.net) and index files for the field size: for the ASI290MC on 750 mm
(25.7' x 14.6') the 4205-4207 series (all-sky, ~1 GB; copied from KStars' folder on the Pi).
The image is solved from its luminance exactly as astropy reads the array, so the WCS fits every later step
that reads the file with astropy.
"""
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

fn, ra, dec = sys.argv[1:4]
scale = float(sys.argv[4]) if len(sys.argv) > 4 else 0.80
index_dir = os.path.expanduser(sys.argv[5] if len(sys.argv) > 5 else "~/.local/share/astrometry")
if not shutil.which("solve-field"):
    sys.exit("solve-field not installed (sudo apt install astrometry.net)")
if not os.path.isdir(index_dir) or not any(f.startswith("index-") for f in os.listdir(index_dir)):
    sys.exit(f"no astrometry.net index files in {index_dir}")

with fits.open(fn) as hd:
    hdr = hd[0].header.copy()
    data = hd[0].data.copy()
lum = data.mean(axis=0).astype(np.float32) if data.ndim == 3 else data.astype(np.float32)

with tempfile.TemporaryDirectory() as tmp:
    fits.writeto(os.path.join(tmp, "lum.fits"), lum, overwrite=True)
    cfg = os.path.join(tmp, "astrometry.cfg")
    with open(cfg, "w") as f:
        f.write(f"add_path {index_dir}\nautoindex\ninparallel\n")
    cmd = ["solve-field", "--config", cfg, "--overwrite", "--no-plots", "--no-verify", "--crpix-center",
           "--scale-units", "arcsecperpix", "--scale-low", f"{scale * 0.9:.3f}", "--scale-high", f"{scale * 1.1:.3f}",
           "--ra", ra, "--dec", dec, "--radius", "2", "--downsample", "2", "--cpulimit", "120",
           "--no-remove-lines", "--uniformize", "0",
           "--new-fits", "none", "--index-xyls", "none", "--rdls", "none", "--corr", "none", "--match", "none",
           os.path.join(tmp, "lum.fits")]
    r = subprocess.run(cmd, capture_output=True, text=True)
    wcs_file = os.path.join(tmp, "lum.wcs")
    if not os.path.exists(wcs_file):
        print(r.stdout[-1500:], r.stderr[-800:])
        sys.exit("astrometry.net could not solve the image")
    w = WCS(fits.getheader(wcs_file))

# replace any old WCS keys with the new solution
for key in list(hdr.keys()):
    if key.startswith(("CTYPE", "CRPIX", "CRVAL", "CDELT", "CD1_", "CD2_", "PC1_", "PC2_", "CUNIT", "CROTA",
                       "A_", "B_", "AP_", "BP_", "LONPOLE", "LATPOLE", "EQUINOX", "RADESYS", "WCSAXES",
                       "IMAGEW", "IMAGEH")):
        del hdr[key]
hdr.update(w.to_header())
fits.writeto(fn, data, hdr, overwrite=True)
c = w.pixel_to_world(lum.shape[1] / 2, lum.shape[0] / 2)
sc = np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix))) * 3600
print(f"solved offline: centre {c.to_string('hmsdms', precision=1)}, {sc:.3f}\"/px")

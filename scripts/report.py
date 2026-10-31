"""Short report of a processing run: Markdown text to post on a forum or use in an article, and the same as a web
page with the images (<NAME>_<TAG>_report.html, next to them in final/: open it in a browser).

Run by process.sh (step 7) with the session settings in the environment:
    report.py <WORK> <NAME> <TAG>     -> <WORK>/final/<NAME>_<TAG>_report.md and .html
Everything comes from what the run left behind: the lights' FITS headers, Siril's logs in <WORK>/siril/ (frames
stacked, per-frame star FWHM from registration), <WORK>/<NAME>_stars.lst (Siril findstar on the stack), the WCS of
<NAME>_linear.fits and the stretch log. Missing pieces (e.g. an external stack without lights) are left out.
REPORT_LANG=uk (Ukrainian, default) or en. Optional session settings for the text: TITLE (e.g. "NGC 2340 galaxy group"), ABOUT (a sentence or two about the
object), TELESCOPE (default "Ньютон 150/750" / "150/750 Newtonian"), CAMERA (when the headers have none, e.g. a DSLR).
"""
import glob
import os
import re
import sys
from datetime import datetime

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

work, name, tag = sys.argv[1:4]
env = os.environ.get
siril_dir = os.path.join(work, "siril")


def read(path):
    try:
        with open(path, errors="replace") as f:
            return f.read()
    except OSError:
        return ""


def fits_files(folder):
    if not folder or not os.path.isdir(folder):
        return []
    return sorted(f for f in glob.glob(os.path.join(folder, "*")) if re.search(r"\.(fits?|fts)$", f, re.I))


def raw_files(folder):
    if not folder or not os.path.isdir(folder):
        return []
    return sorted(f for f in glob.glob(os.path.join(folder, "*")) if os.path.isfile(f))


def fmt_range(a, b, unit="", nd=0):
    a, b = round(a, nd), round(b, nd)
    return f"{a:.{nd}f}{unit}" if a == b else f"{a:.{nd}f}–{b:.{nd}f}{unit}"


# ---- lights: headers
lights = fits_files(env("LIGHTS"))
hdrs = []
for f in lights:
    try:
        hdrs.append(fits.getheader(f))
    except OSError:
        pass
n_lights = len(lights) or len(raw_files(env("LIGHTS")))
skip = set()                                                     # LIGHTS_SKIP: positions in file-name order, from 1
for r in (env("LIGHTS_SKIP") or "").split():
    a, _, b = r.partition("-")
    skip.update(range(int(a), int(b or a) + 1))
exps = [h.get("EXPTIME") for h in hdrs]
used = [e for i, e in enumerate(exps, 1) if i not in skip and e]
exp = used[0] if used and len(set(used)) == 1 else None          # one exposure for the stacked frames
gain = offset = None
times, temps, alts = [], [], []
camera = env("CAMERA") or ""
for h in hdrs:
    gain = h.get("GAIN", gain)
    offset = h.get("OFFSET", offset)
    camera = camera or (h.get("INSTRUME") or "").replace("ZWO CCD ", "ZWO ")
    if h.get("DATE-OBS"):
        try:
            times.append(datetime.fromisoformat(h["DATE-OBS"][:19]))
        except ValueError:
            pass
    if h.get("CCD-TEMP") is not None:
        temps.append(float(h["CCD-TEMP"]))
    if h.get("OBJCTALT") is not None:
        alts.append(float(h["OBJCTALT"]))
focal = hdrs[0].get("FOCALLEN") if hdrs else None
focal = float(focal) if focal else float(env("FOCAL") or 0) or None

# ---- Siril logs
stack_log = read(os.path.join(siril_dir, "2_stack.log"))
siril_ver = re.search(r"siril (\d+\.\d+\.\d+)", stack_log + read(os.path.join(siril_dir, "3_post.log")))
m = re.findall(r"Integration of (\d+) images", stack_log)
n_stacked = int(m[-1]) if m else None
frame_fwhm = np.array([float(x) for x in re.findall(r"Gaussian profile stars in image, channel #1 \(FWHM ([\d.]+)\)",
                                                    stack_log)])
stretch_log = read(os.path.join(siril_dir, "6_stretch.log"))
m = re.search(r"white balance from (\d+) Gaia stars", stretch_log)
n_gaia = int(m.group(1)) if m else None

# ---- stack: WCS and Siril findstar
scale = field = None
lin = os.path.join(work, f"{name}_linear.fits")
if os.path.exists(lin):
    h = fits.getheader(lin)
    try:
        scale = proj_plane_pixel_scales(WCS(h, naxis=2)).mean() * 3600
        field = (h["NAXIS1"] * scale / 60, h["NAXIS2"] * scale / 60)
    except Exception:
        pass
stars = []
for line in read(os.path.join(work, f"{name}_stars.lst")).splitlines():
    p = line.split()
    if line.startswith("#") or len(p) < 15:
        continue
    fx, fy, sat = float(p[7]), float(p[8]), int(p[14])
    if not sat:
        stars.append((fx, fy))
stars = np.array(stars)

# ---- text (REPORT_LANG: uk = Ukrainian, default; en = English)
LANG = env("REPORT_LANG") or "uk"
T = {
    "en": dict(
        acq="Acquisition", date="Date", alt="altitude", tel="Telescope", cam="camera", sensor="sensor",
        newt="150/750 Newtonian", lights="Lights", s="s", minute="min", left_out="{n} left out",
        no_dark_exp=" (no darks at their exposure)", gain="gain", offset="offset", stacked="{n} stacked",
        cal="Calibration", darks="{n} {w}", w_dark=("dark", "darks", "darks"), mdark="master dark from an earlier night",
        nodark="no darks (constant offset subtracted{hot})", hot=", hot pixels mapped from an older dark",
        flats="{n} {w}", w_flat=("flat", "flats", "flats"), dflats="{n} {w}", w_dflat=("dark-flat", "dark-flats", "dark-flats"),
        mflat="master flat from another night", noflat="no flats",
        quality="Image quality", px="px", scale="Scale {s}″/px, field {w}′ × {h}′",
        fwhm1="Star FWHM, single frames (Siril): median {p} px = {a}″ (range {lo}–{hi}″)",
        fwhm2="Star FWHM, stack (Siril findstar, {n} {w}): {p} px = {a}″, roundness {r}", w_star=("star", "stars", "stars"),
        proc="Processing", tools=" (Siril {v} + Python; classical methods, no AI)",
        stack="Calibration, registration and stacking in Siril: winsorized sigma rejection, frames weighted by star "
              "sharpness", keep=", sharpest {k} kept",
        bg="Background", bg_poly="2nd-degree polynomial gradient removal (Siril)",
        bg_radial="vignetting model (radial) + linear gradient", bg_none="no gradient removal",
        bg_res=", then a smooth residual model ({p} px)", bg_mask="; the target{x} masked",
        bg_mask_x=" and the other galaxies in the field",
        solve="Plate solving", deconv="Deconvolution: Richardson-Lucy, {n} {w}, on the galaxy only", w_iter=("iteration", "iterations", "iterations"),
        colour="Colour: white balance from {n} Gaia DR3 {w} (Sun-like stars made white)", w_star_by=("star", "stars", "stars"),
        dn="Noise reduction: NL-Bayes at {p} % + wavelet smoothing of the sky only",
        stretch="Stretch: arcsinh curve{x}; north up, east left",
        stars_white="; bright star cores fade to white as in a photograph", clipped=" (centres clipped in the camera made neutral)",
        images="Images", nodn="no noise reduction",
        h_orient="North up, east left. Click an image for full size.", h_dn="Noise reduction", h_nodn="None",
        h_dn_group="Noise reduction", h_crop="Centre, enlarged 2×.", h_light="Light", h_dark="Dark",
        h_theme="Page theme", h_field="full field", h_crop_alt="centre at 2×",
        months=["January", "February", "March", "April", "May", "June", "July", "August", "September", "October",
                "November", "December"]),
    "uk": dict(
        acq="Зйомка", date="Дата", alt="висота", tel="Телескоп", cam="камера", sensor="сенсор",
        newt="Ньютон 150/750", lights="Кадри", s="с", minute="хв", left_out="{n} не використано",
        no_dark_exp=" (немає дарків для цієї витримки)", gain="гейн", offset="офсет", stacked="у стек {n}",
        cal="Калібрування", darks="{n} {w}", w_dark=("дарк", "дарки", "дарків"), mdark="майстер-дарк з іншої ночі",
        nodark="без дарків (віднято сталий рівень{hot})", hot=", гарячі пікселі за старшим дарком",
        flats="{n} {w}", w_flat=("флет", "флети", "флетів"), dflats="{n} {w}", w_dflat=("дарк-флет", "дарк-флети", "дарк-флетів"),
        mflat="майстер-флет з іншої ночі", noflat="без флетів",
        quality="Якість зображення", px="пкс", scale="Масштаб {s}″/пкс, поле {w}′ × {h}′",
        fwhm1="FWHM зірок, окремі кадри (Siril): медіана {p} пкс = {a}″ (від {lo} до {hi}″)",
        fwhm2="FWHM зірок, стек (Siril findstar, {n} {w}): {p} пкс = {a}″, округлість {r}", w_star=("зірка", "зірки", "зірок"),
        proc="Обробка", tools=" (Siril {v} + Python; класичні методи, без ШІ)",
        stack="Калібрування, вирівнювання та складання в Siril: вінзоризоване сигма-відкидання, кадри зважені за "
              "різкістю зірок", keep=", залишено {k} найрізкіших",
        bg="Фон", bg_poly="видалення градієнта поліномом 2-го ступеня (Siril)",
        bg_radial="модель віньєтування (радіальна) + лінійний градієнт", bg_none="без видалення градієнта",
        bg_res=", потім гладка модель залишку ({p} пкс)", bg_mask="; ціль{x} замасковано",
        bg_mask_x=" та інші галактики в полі",
        solve="Астрометрія", deconv="Деконволюція: Річардсон–Люсі, {n} {w}, лише на галактиці", w_iter=("ітерація", "ітерації", "ітерацій"),
        colour="Колір: баланс білого за {n} {w} Gaia DR3 (зірки сонячного типу — білі)", w_star_by=("зіркою", "зірками", "зірками"),
        dn="Шумозаглушення: NL-Bayes на {p} % + вейвлет-згладжування лише неба",
        stretch="Розтяг: крива arcsinh{x}; північ угорі, схід ліворуч",
        stars_white="; ядра яскравих зірок переходять у білий, як на фотографії",
        clipped=" (центри, обрізані в камері, зроблено нейтральними)",
        images="Зображення", nodn="без шумозаглушення",
        h_orient="Північ угорі, схід ліворуч. Натисніть на зображення, щоб відкрити повний розмір.",
        h_dn="Шумозаглушення", h_nodn="Без нього", h_dn_group="Шумозаглушення", h_crop="Центр, збільшено у 2 рази.",
        h_light="Світла", h_dark="Темна", h_theme="Тема сторінки", h_field="повне поле", h_crop_alt="центр, збільшено у 2 рази",
        months=["січня", "лютого", "березня", "квітня", "травня", "червня", "липня", "серпня", "вересня", "жовтня",
                "листопада", "грудня"]),
}[LANG]


def pl(n, forms):
    """the noun form for a count: Ukrainian has one / few (2-4) / many; English one / many"""
    n = int(n)
    if LANG == "uk":
        if n % 10 == 1 and n % 100 != 11:
            return forms[0]
        return forms[1] if 2 <= n % 10 <= 4 and not 12 <= n % 100 <= 14 else forms[2]
    return forms[0] if n == 1 else forms[2]


def num(v, nd=1):
    """a number with the language's decimal separator"""
    t = f"{v:.{nd}f}"
    return t.replace(".", ",") if LANG == "uk" else t


def day(t):
    return f"{t.day} {T['months'][t.month - 1]} {t.year}"


def rng(a, b, unit, nd=0):
    a, b = round(a, nd), round(b, nd)
    return f"{num(a, nd)}{unit}" if a == b else f"{num(a, nd)}–{num(b, nd)}{unit}"


title = env("TITLE") or name.upper().replace("NGC", "NGC ").replace("IC", "IC ").replace("  ", " ")
date_txt = day(times[0]) if times else ""
lines = [f"# {title}" + (f" — {date_txt}" if date_txt else ""), ""]
if env("ABOUT"):
    lines += [env("ABOUT"), ""]

lines.append(f"**{T['acq']}**")
if times:
    t = f"- {T['date']}: {day(times[0])}, {times[0]:%H:%M}–{max(times):%H:%M} UTC"
    if alts:
        t += f", {T['alt']} {alts[0]:.0f}° → {alts[-1]:.0f}°"
    lines.append(t)
tel = env("TELESCOPE") or T["newt"]
optics = f"- {T['tel']}: {tel}" + (f" ({focal:.0f} mm)" if focal and str(int(focal)) not in tel else "")
if camera:
    optics += f"; {T['cam']}: {camera}" + (f", {T['sensor']} {rng(min(temps), max(temps), ' °C')}" if temps else "")
lines.append(optics)
if n_lights:
    groups = sorted({e: exps.count(e) for e in exps if e}.items())
    if len(groups) > 1:
        t = f"- {T['lights']}: {n_lights} (" + " + ".join(f"{c} × {e:g} {T['s']}" for e, c in groups) + ")"
    else:
        t = f"- {T['lights']}: {n_lights}" + (f" × {exp:g} {T['s']}" if exp else "")
    if skip:
        t += ", " + T["left_out"].format(n=len(skip)) + (T["no_dark_exp"] if len(groups) > 1 else "")
    if gain is not None:
        t += f", {T['gain']} {gain:g}" + (f" ({T['offset']} {offset:g})" if offset is not None else "")
    if n_stacked:
        t += "; " + T["stacked"].format(n=n_stacked) + (f" = {n_stacked * exp / 60:.0f} {T['minute']}" if exp else "")
    lines.append(t)
cal = []
if env("DARKS"):
    nd_ = len(raw_files(env("DARKS")))
    cal.append(T["darks"].format(n=nd_, w=pl(nd_, T["w_dark"])))
elif env("MASTER_DARK"):
    cal.append(T["mdark"])
elif env("LIGHT_BIAS"):
    cal.append(T["nodark"].format(hot=T["hot"] if env("HOT_DARK") else ""))
if env("FLATS"):
    nf_, ndf_ = len(raw_files(env("FLATS"))), len(raw_files(env("DARKFLATS")))
    cal.append(T["flats"].format(n=nf_, w=pl(nf_, T["w_flat"]))
               + (", " + T["dflats"].format(n=ndf_, w=pl(ndf_, T["w_dflat"])) if env("DARKFLATS") else ""))
elif env("MASTER_FLAT"):
    cal.append(T["mflat"])
else:
    cal.append(T["noflat"])
if n_lights:
    lines.append(f"- {T['cal']}: " + ", ".join(cal))
lines.append("")

lines.append(f"**{T['quality']}**")
if scale:
    lines.append("- " + T["scale"].format(s=num(scale, 2), w=num(field[0]), h=num(field[1])))
if len(frame_fwhm) and scale:
    lines.append("- " + T["fwhm1"].format(p=num(np.median(frame_fwhm)), a=num(np.median(frame_fwhm) * scale),
                                          lo=num(frame_fwhm.min() * scale), hi=num(frame_fwhm.max() * scale)))
if len(stars) and scale:
    f = np.sqrt(stars[:, 0] * stars[:, 1])
    rnd = np.minimum(stars[:, 0], stars[:, 1]) / np.maximum(stars[:, 0], stars[:, 1])
    lines.append("- " + T["fwhm2"].format(n=len(stars), w=pl(len(stars), T["w_star"]), p=num(np.median(f)), a=num(np.median(f) * scale),
                                          r=num(np.median(rnd), 2)))
lines.append("")

lines.append(f"**{T['proc']}**" + (T["tools"].format(v=siril_ver.group(1)) if siril_ver else ""))
if n_lights:
    lines.append("- " + T["stack"] + (T["keep"].format(k=env("KEEP_BEST")) if env("KEEP_BEST", "100%") != "100%" else ""))
bg = env("BG", "poly")
bgt = {"poly": T["bg_poly"], "radial": T["bg_radial"], "none": T["bg_none"]}.get(bg, bg)
if env("BG_RESIDUAL", "0") != "0":
    bgt += T["bg_res"].format(p=env("BG_RESIDUAL"))
if bg != "none":
    bgt += T["bg_mask"].format(x=T["bg_mask_x"] if env("GALAXIES_EXTRA") else "")
lines.append(f"- {T['bg']}: {bgt}")
solver = "astrometry.net" if ("solved offline" in read(os.path.join(siril_dir, "3_solve.log"))
                              or env("SOLVE") == "offline") else "Siril"
lines.append(f"- {T['solve']}: {solver}")
if env("DECONV", "0") == "1":
    lines.append("- " + T["deconv"].format(n=env("DECONV_ITERS", "20"), w=pl(env("DECONV_ITERS", "20"), T["w_iter"])))
if n_gaia:
    lines.append("- " + T["colour"].format(n=n_gaia, w=pl(n_gaia, T["w_star_by"])))
lines.append("- " + T["dn"].format(p=f"{float(env('DENOISE_MOD', '0.7')) * 100:.0f}"))
so = env("STRETCH_OPTS") or ""
stars_txt = "" if "--highlight=hue" in so else T["stars_white"] + ("" if "--clipped=off" in so else T["clipped"])
lines.append("- " + T["stretch"].format(x=stars_txt))
lines.append("")
out_imgs = [f"{name}_{tag}_field.png", f"{name}_{tag}_crop_2x.png"]
if env("NODENOISE", "1") == "1":
    out_imgs += [f"{name}_{tag}_nodn_field.png ({T['nodn']})"]
lines.append(f"{T['images']}: " + ", ".join(out_imgs))

out = os.path.join(work, "final", f"{name}_{tag}_report.md")
os.makedirs(os.path.dirname(out), exist_ok=True)
with open(out, "w") as f:
    f.write("\n".join(lines) + "\n")
print(f"saved {out}")

# ---- the same as a web page, with the images (they sit next to it in final/)
import html as H


def inline(t):
    return re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", H.escape(t))


body, in_list = [], False
for ln in lines:
    if ln.startswith("- "):
        if not in_list:
            body.append("<ul>")
            in_list = True
        body.append(f"<li>{inline(ln[2:])}</li>")
        continue
    if in_list:
        body.append("</ul>")
        in_list = False
    if ln.startswith("# ") or ln.startswith(f"{T['images']}:") or not ln.strip():
        continue                                                 # title and image list are shown separately
    m_ = re.match(r"\*\*(.+?)\*\*(.*)", ln)
    body.append(f"<h2>{H.escape(m_.group(1))}<span>{H.escape(m_.group(2))}</span></h2>" if m_ else f"<p class='about'>{inline(ln)}</p>")
if in_list:
    body.append("</ul>")
about_html = [b_ for b_ in body if b_.startswith("<p class='about'>")]
info_html = [b_ for b_ in body if not b_.startswith("<p class='about'>")]
# sections side by side on wide screens: split at each <h2>
secs, cur = [], []
for b_ in info_html:
    if b_.startswith("<h2>") and cur:
        secs.append(cur)
        cur = []
    cur.append(b_)
if cur:
    secs.append(cur)
f_dn, f_nodn, f_crop = f"{name}_{tag}_field.png", f"{name}_{tag}_nodn_field.png", f"{name}_{tag}_crop_2x.png"
have_nodn = os.path.exists(os.path.join(work, "final", f_nodn))
toggle = ("""<div class="toggle" role="group" aria-label="%s">
  <button aria-pressed="true" data-src="%s">%s</button>
  <button aria-pressed="false" data-src="%s">%s</button></div>""" % (T["h_dn_group"], f_dn, T["h_dn"], f_nodn, T["h_nodn"])) if have_nodn else ""
page = f"""<!doctype html>
<html lang="{LANG}"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{H.escape(title)}</title>
<style>
:root {{ --bg:#f6f7f9; --panel:#ffffff; --text:#1b1f27; --muted:#5d6675; --line:#dde1e7; --accent:#2c5fd6;
         --chip:#e9ecf1; color-scheme:light; }}
:root[data-theme="dark"] {{ --bg:#0b0d12; --panel:#141821; --text:#e6e8ee; --muted:#9aa3b2; --line:#252b37;
         --accent:#8fb4ff; --chip:#1d2230; color-scheme:dark; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; background:var(--bg); color:var(--text); font:16px/1.55 system-ui,-apple-system,"Segoe UI",sans-serif; }}
main {{ max-width:1200px; margin:0 auto; padding:24px 16px 48px; }}
h1 {{ font-size:clamp(1.5rem,3.5vw,2.2rem); margin:0 0 4px; font-weight:650; letter-spacing:-0.01em; }}
.date {{ color:var(--muted); margin:0 0 14px; }}
.about {{ max-width:70ch; margin:0 0 20px; }}
figure {{ margin:0 0 12px; }}
figure img {{ width:100%; height:auto; display:block; border-radius:6px; background:#000; }}
figcaption {{ color:var(--muted); font-size:.9rem; margin-top:6px; }}
.bar {{ display:flex; justify-content:space-between; align-items:center; gap:12px; flex-wrap:wrap; margin:0 0 10px; }}
.toggle {{ display:inline-flex; flex:none; border:1px solid var(--line); border-radius:999px; overflow:hidden; }}
.toggle button {{ background:transparent; color:var(--muted); border:0; padding:6px 14px; font:inherit; cursor:pointer; }}
.toggle button[aria-pressed="true"] {{ background:var(--chip); color:var(--text); }}
.top {{ display:flex; justify-content:space-between; align-items:flex-start; gap:12px; }}
.toggle button:focus-visible {{ outline:2px solid var(--accent); outline-offset:-2px; }}
.grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(300px,1fr)); gap:16px; margin-top:28px; }}
section {{ background:var(--panel); border:1px solid var(--line); border-radius:8px; padding:14px 18px; }}
h2 {{ font-size:1rem; margin:0 0 6px; }} h2 span {{ font-weight:400; color:var(--muted); font-size:.85rem; }}
ul {{ margin:0; padding-left:1.1em; }} li {{ margin:3px 0; }}
a {{ color:var(--accent); }}
</style></head>
<body><main>
<div class="top"><h1>{H.escape(title)}</h1>
<div class="toggle theme" role="group" aria-label="{T['h_theme']}">
  <button aria-pressed="true" data-theme="light">{T['h_light']}</button><button aria-pressed="false" data-theme="dark">{T['h_dark']}</button></div></div>
<p class="date">{H.escape(date_txt)}</p>
{"".join(about_html)}
<div class="bar"><span class="date" style="margin:0">{T['h_orient']}</span>{toggle}</div>
<figure><a id="fieldlink" href="{f_dn}"><img id="field" src="{f_dn}" alt="{H.escape(title)}, {T['h_field']}"></a></figure>
<figure><a href="{f_crop}"><img src="{f_crop}" alt="{H.escape(title)}, {T['h_crop_alt']}"></a>
<figcaption>{T['h_crop']}</figcaption></figure>
<div class="grid">{"".join("<section>" + "".join(sc) + "</section>" for sc in secs)}</div>
</main>
<script>
function pick(group, b) {{
  group.querySelectorAll("button").forEach(function (o) {{ o.setAttribute("aria-pressed", o === b); }});
}}
document.querySelectorAll(".toggle:not(.theme) button").forEach(function (b) {{
  b.addEventListener("click", function () {{
    pick(b.parentNode, b);
    document.getElementById("field").src = b.dataset.src;
    document.getElementById("fieldlink").href = b.dataset.src;
  }});
}});
function setTheme(t) {{
  document.documentElement.setAttribute("data-theme", t);
  pick(document.querySelector(".theme"), document.querySelector('.theme [data-theme="' + t + '"]'));
  try {{ localStorage.setItem("report-theme", t); }} catch (e) {{}}
}}
document.querySelectorAll(".theme button").forEach(function (b) {{
  b.addEventListener("click", function () {{ setTheme(b.dataset.theme); }});
}});
try {{ var t0 = localStorage.getItem("report-theme"); if (t0 === "dark") setTheme("dark"); }} catch (e) {{}}
</script>
</body></html>
"""
out_html = out[:-3] + ".html"
with open(out_html, "w") as f:
    f.write(page)
print(f"saved {out_html}")

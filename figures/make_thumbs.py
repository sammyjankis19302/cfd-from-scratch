"""Build the module thumbnails used by the table on the front page.

    python3 figures/make_thumbs.py

One representative figure per module, letterboxed onto a canvas of identical
aspect so the six cells of that table line up instead of stepping up and down
with whatever each source figure happens to measure. Chosen for legibility at
about 280 px wide, which is the column width GitHub gives a three-column table:
a figure whose point is made by thin lines or small text does not survive that,
so those modules are represented by one that does.
"""
import os
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, "thumbs")
os.makedirs(OUT, exist_ok=True)

W, H, PAD = 560, 374, 10          # 3:2, the shape of most of the sources

# (module, source, crop). The crop is a fraction of the source, and is there
# for the two figures that are unreadable whole at this size: a two-panel
# wavenumber plot and a seven-panel grid of cavities. Taking part of each keeps
# the point legible -- for 06, four Reynolds numbers in which the red
# counter-rotating region visibly takes over.
SOURCES = [
    ("01", "01-linear-advection/figures/report/points-per-wavelength.png", None),
    ("02", "02-spectral-analysis/figures/report/numerical-wavenumbers.png", (0, 0, 0.52, 1)),
    ("03", "03-advection-2d/figures/report/advection-2d.png", None),
    ("04", "04-burgers/figures/report/burgers-inviscid.png", None),
    ("05", "05-poisson/figures/report/laplace-and-poisson.png", None),
    ("06", "06-lid-driven-cavity/figures/report/streamlines-all-re.png", (0, 0, 0.50, 1)),
]

for name, rel, crop in SOURCES:
    src = Image.open(os.path.join(ROOT, rel)).convert("RGB")
    if crop:
        x0, y0, x1, y1 = crop
        src = src.crop((int(x0 * src.width), int(y0 * src.height),
                        int(x1 * src.width), int(y1 * src.height)))
    box = (W - 2 * PAD, H - 2 * PAD)
    s = min(box[0] / src.width, box[1] / src.height)
    im = src.resize((max(int(src.width * s), 1), max(int(src.height * s), 1)),
                    Image.LANCZOS)
    canvas = Image.new("RGB", (W, H), (255, 255, 255))
    canvas.paste(im, ((W - im.width) // 2, (H - im.height) // 2))
    path = os.path.join(OUT, name + ".png")
    canvas.save(path, optimize=True)
    print(f"  thumbs/{name}.png  {os.path.getsize(path)/1024:5.0f} KB  <- {rel}")

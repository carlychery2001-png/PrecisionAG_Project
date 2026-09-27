"""Shared setup for the Etapa 1 maps (EARTH University lot): data, interpolated surface, imagery and map layout.

Imported by publication_map.py (contour and slope maps) and precision_maps.py (wetness, zones, sampling).
Language is taken from the command line of the calling script: python <script>.py [en|es]
"""
import csv
import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import BoundaryNorm, ListedColormap, LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon, Rectangle
from matplotlib.path import Path
from matplotlib.ticker import FuncFormatter, MultipleLocator
from scipy.interpolate import RBFInterpolator

import imagery

# Printed captions contain non-ASCII characters; keep the console from failing on Windows code pages
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

POINTS_CSV = os.path.join("data", "levantamiento_rtk_2026-09-26.csv")   # RTK points: survey grid Pt*, soil samples M*
AREA_CSV = os.path.join("data", "vertices_lote_2026-09-26.csv")         # the four boundary vertices of the lot
BASE_CSV = os.path.join("data", "estacion_base_2026-09-26.csv")         # RTK base station position
OUT_DIR = "outputs"                                                     # maps, tables and captions
LANG = sys.argv[1] if len(sys.argv) > 1 else "en"
DEC = "," if LANG == "es" else "."   # decimal comma for Spanish
MINOR, MAJOR = 0.25, 0.50            # contour interval / index contour (m); 0.25 m meets NMAS for the ~8.0 cm LOOCV RMSE
SLOPE_BREAKS = [0, 0.5, 1, 1.5, 2, 3, 5]   # slope classes (%); 0.5/1/2/5 follow FAO (2006)
GRID = 0.25                          # interpolation grid (m)
ARROW_SPACING = 9                    # downslope arrow spacing (m)
FILL_ALPHA = 1.0                     # data layers are opaque inside the field so map colours match their keys;
ELEV_ALPHA = 1.0                     # the imagery shows around the field
GEOID_N = 11.42                      # EGM2008 geoid height at the field (GeographicLib GeoidEval, 10.2085 N 83.5986 W)
IMG_RES = 0.15                       # resampled imagery pixel size (m)
MM = 1 / 25.4                        # mm -> inch
WIDTH_MM = 180                       # double-column journal width

TXT = {
    "en": dict(east="Easting (m)", north="Northing (m)", elev="Approx. elevation (m a.s.l.)", slope="Slope (%)",
               boundary="Field boundary", index="Index contour", inter="Intermediate contour",
               contour="Contours", survey="RTK survey point", sample="Soil sample",
               flow="Downslope direction", img="Imagery",
               twi="Topographic wetness index (TWI)", drier="drier", wetter="wetter",
               depression="Depression (fill depth ≥ 1 cm)", zone="Zone", zones="Management zones",
               existing="Existing soil sample", proposed="Proposed soil sample"),
    "es": dict(east="Este (m)", north="Norte (m)", elev="Elevación aprox. (m s. n. m.)", slope="Pendiente (%)",
               boundary="Límite del lote", index="Curva índice", inter="Curva intermedia",
               contour="Curvas de nivel", survey="Punto de levantamiento RTK", sample="Muestra de suelo",
               flow="Dirección de escurrimiento", img="Imagen",
               twi="Índice topográfico de humedad (TWI)", drier="más seco", wetter="más húmedo",
               depression="Depresión (profundidad ≥ 1 cm)", zone="Zona", zones="Zonas de manejo",
               existing="Muestra de suelo existente", proposed="Muestra de suelo propuesta"),
}[LANG]

def _register_inter():
    """Use Inter (the report typeface, from the MiKTeX font tree) for all map text; fall back to Arial."""
    import glob
    import os
    import subprocess
    from matplotlib import font_manager
    try:
        path = subprocess.run(["kpsewhich", "Inter-Regular.otf"], capture_output=True, text=True, timeout=60).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        path = ""
    if not path:
        return "Arial"
    for f in glob.glob(os.path.join(os.path.dirname(path), "Inter-*.otf")):
        font_manager.fontManager.addfont(f)
    return "Inter"


mpl.rcParams.update({
    "font.family": _register_inter(),
    "font.size": 8,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
    "xtick.direction": "out", "ytick.direction": "out",
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
})


def num(v, d=1):
    """Format a number with the language's decimal separator."""
    return f"{v:.{d}f}".replace(".", DEC)


def level_label(v):
    """Contour label: one decimal on half-metre levels, two on quarter levels."""
    return num(v, 1) if abs(v * 2 - round(v * 2)) < 1e-6 else num(v, 2)


def short(v):
    """0.5 -> '0.5', 2.0 -> '2' (class limits)."""
    return num(v).replace(DEC + "0", "")


def thousands(v):
    return f"{v:,.0f}".replace(",", " ")


def read(path):
    """Points with approximate orthometric height: the controller's 'Elevacion' is a WGS 84 ellipsoidal
    height (no geoid model applied), so the EGM2008 geoid height is subtracted."""
    with open(path, encoding="utf-8", newline="") as f:
        return [(r["PT Nom"], float(r["Este"]), float(r["Norte"]), float(r["Elevacion"]) - GEOID_N)
                for r in csv.DictReader(f)]


pts = read(POINTS_CSV)
area = read(AREA_CSV)
bxy = np.array([p[1:3] for p in area])
boundary = Path(np.vstack([bxy, bxy[:1]]), closed=True)
area_m2 = 0.5 * abs(np.dot(bxy[:, 0], np.roll(bxy[:, 1], -1)) - np.dot(bxy[:, 1], np.roll(bxy[:, 0], -1)))
perim = np.sum(np.hypot(*np.diff(np.vstack([bxy, bxy[:1]]), axis=0).T))
samples = sorted([p for p in pts if p[0].startswith("M")], key=lambda p: int(p[0][1:]))
grid_pts = np.array([p[1:3] for p in pts if p[0].startswith("Pt")])


def edge_dist(px, py):
    """Distance (m) from a point to the nearest boundary edge."""
    a, b = bxy, np.roll(bxy, -1, axis=0)
    d = b - a
    t = np.clip(((px - a[:, 0]) * d[:, 0] + (py - a[:, 1]) * d[:, 1]) / (d ** 2).sum(1), 0, 1)
    return np.min(np.hypot(px - a[:, 0] - t * d[:, 0], py - a[:, 1] - t * d[:, 1]))


xyz = np.array([p[1:] for p in pts])
origin = xyz[:, :2].mean(axis=0)
xy, z = xyz[:, :2] - origin, xyz[:, 2]


# ------------------------------------------------ interpolation + LOOCV
def loocv(smoothing):
    err = np.empty(len(z))
    for i in range(len(z)):
        m = np.arange(len(z)) != i
        f = RBFInterpolator(xy[m], z[m], kernel="thin_plate_spline", smoothing=smoothing)
        err[i] = f(xy[i:i + 1])[0] - z[i]
    return err


candidates = [0.0, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0]
cv = {s: loocv(s) for s in candidates}
best = min(cv, key=lambda s: np.sqrt(np.mean(cv[s] ** 2)))
err = cv[best]
rmse, mae, bias = np.sqrt(np.mean(err ** 2)), np.mean(np.abs(err)), np.mean(err)
rbf = RBFInterpolator(xy, z, kernel="thin_plate_spline", smoothing=best)

pad = 4
gx = np.arange(bxy[:, 0].min() - pad, bxy[:, 0].max() + pad, GRID)
gy = np.arange(bxy[:, 1].min() - pad, bxy[:, 1].max() + pad, GRID)
GX, GY = np.meshgrid(gx, gy)
Zfull = rbf(np.c_[GX.ravel(), GY.ravel()] - origin).reshape(GX.shape)

# Keep a skirt outside the field so fills and contours reach the clipped edge cleanly
flat = np.c_[GX.ravel(), GY.ravel()]
inside = boundary.contains_points(flat).reshape(GX.shape)
skirt = (boundary.contains_points(flat, radius=6) | boundary.contains_points(flat, radius=-6)).reshape(GX.shape)
Z = np.ma.masked_where(~skirt, Zfull)
Zin = np.ma.masked_where(~inside, Zfull)

# Slope (%) and downslope direction from the unmasked surface (no edge artefacts)
dzdy, dzdx = np.gradient(Zfull, GRID)
S = 100 * np.hypot(dzdx, dzdy)
Sskirt = np.ma.masked_where(~skirt, S)
Sin = np.ma.masked_where(~inside, S)

A = np.c_[np.ones(len(xy)), xy]
_, px, py = np.linalg.lstsq(A, z, rcond=None)[0]
plane_slope = 100 * np.hypot(px, py)
plane_az = np.degrees(np.arctan2(-px, -py)) % 360

zlo = np.floor(Zin.min() / MINOR) * MINOR
zhi = np.ceil(Zin.max() / MINOR) * MINOR
levels = np.round(np.arange(zlo, zhi + MINOR / 2, MINOR), 2)
is_major = np.isclose(np.mod(levels + 1e-9, MAJOR), 0, atol=1e-6)

halo = [pe.withStroke(linewidth=1.8, foreground="white")]
dark_halo = [pe.withStroke(linewidth=2.8, foreground="black")]

# Basemap imagery on exactly the map extent
rgb, img_zoom = imagery.mosaic((gx.min(), gx.max(), gy.min(), gy.max()), res=IMG_RES)
img_extent = (gx.min(), gx.min() + rgb.shape[1] * IMG_RES, gy.max() - rgb.shape[0] * IMG_RES, gy.max())
meta = imagery.metadata(*origin)
img_date = meta["DATE (YYYYMMDD)"]
MONTHS = {"en": "Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec",
          "es": "ene feb mar abr may jun jul ago sep oct nov dic"}[LANG].split()
img_date_txt = f"{int(img_date[6:])} {MONTHS[int(img_date[4:6]) - 1]} {img_date[:4]}"
img_credit = f"{TXT['img']}: Esri, {meta['SOURCE']} ({meta['SOURCE_INFO']}, {img_date_txt}), Earthstar Geographics"


# ------------------------------------------------------------- shared layout
def new_map():
    fig = plt.figure(figsize=(WIDTH_MM * MM, 137 * MM))
    ax = fig.add_axes([0.085, 0.105, 0.80, 0.88])
    cax = fig.add_axes([0.905, 0.25, 0.018, 0.55])
    ax.imshow(rgb, extent=img_extent, origin="upper", interpolation="bilinear", zorder=0)
    clip = Polygon(bxy, closed=True, transform=ax.transData, fc="none", ec="none")
    ax.add_patch(clip)
    return fig, ax, cax, clip


def label_contours(ax, cs, lv, fs, wt):
    texts = ax.clabel(cs, levels=lv, fmt=level_label, fontsize=fs, inline=True, inline_spacing=2)
    for t in texts:
        # Drop labels that would be cut by the field edge
        if not (boundary.contains_point(t.get_position()) and edge_dist(*t.get_position()) > 3.5):
            t.set_visible(False)
            continue
        t.set_fontweight(wt)
        t.set_path_effects(halo)


def area_label():
    """Boundary legend entry with the lot area (grid area, CRTM05)."""
    return f"{TXT['boundary']} ({thousands(area_m2)} m²; {num(area_m2 / 1e4, 2)} ha)"


def decorate(ax, extra_handles, show_survey=True, sample_label=None, show_samples=True):
    """Boundary, points, samples, axes, north arrow, scale bar and legend."""
    ax.add_patch(Polygon(bxy, closed=True, fill=False, ec="white", lw=1.3, joinstyle="miter", zorder=6,
                         path_effects=dark_halo))
    if show_survey:
        ax.scatter(grid_pts[:, 0], grid_pts[:, 1], s=3, c="#111111", ec="white", lw=0.25, zorder=7)
    if show_samples:
        ax.scatter([p[1] for p in samples], [p[2] for p in samples], s=34, marker="^",
                   fc="white", ec="black", lw=0.8, zorder=8)
        for name, x, y, _ in samples:
            ax.annotate(name, (x, y), xytext=(4, 3), textcoords="offset points", fontsize=7,
                        fontweight="bold", path_effects=halo, zorder=9)

    ax.set_xlim(gx.min(), gx.max())
    ax.set_ylim(gy.min(), gy.max())
    ax.set_aspect("equal")
    ax.xaxis.set_major_locator(MultipleLocator(20))
    ax.yaxis.set_major_locator(MultipleLocator(20))
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(MultipleLocator(10))
    fmt = FuncFormatter(lambda v, _: thousands(v))
    ax.xaxis.set_major_formatter(fmt)
    ax.yaxis.set_major_formatter(fmt)
    ax.tick_params(which="both", top=True, right=True, labelsize=7)
    plt.setp(ax.get_yticklabels(), rotation=90, va="center")
    ax.set_xlabel(TXT["east"])
    ax.set_ylabel(TXT["north"])
    # CRS and imagery credit below the frame, right-aligned with it
    ax.figure.text(ax.get_position().x1, 0.008, f"CRTM05 (EPSG:5367) · {img_credit}",
                   ha="right", va="bottom", fontsize=5.5, color="#333")

    nx, ny = 0.955, 0.86
    ax.annotate("", xy=(nx, ny + 0.09), xytext=(nx, ny), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>,head_width=0.35,head_length=0.7", lw=1.1, color="black",
                                path_effects=[pe.withStroke(linewidth=2.6, foreground="white")]), zorder=10)
    ax.text(nx, ny + 0.10, "N", transform=ax.transAxes, ha="center", va="bottom", fontsize=10, fontweight="bold",
            path_effects=halo, zorder=10)

    sx, sy, h = gx.min() + 6, gy.min() + 7, 1.3
    ax.add_patch(Rectangle((sx - 2.5, sy - 2), 49, h + 7, fc="white", ec="none", alpha=0.8, zorder=8))
    for x0, x1, fc in ((0, 10, "black"), (10, 20, "white"), (20, 40, "black")):
        ax.add_patch(Rectangle((sx + x0, sy), x1 - x0, h, fc=fc, ec="black", lw=0.5, zorder=9))
    for v in (0, 10, 20, 40):
        ax.text(sx + v, sy + h + 1.0, f"{v}", ha="center", va="bottom", fontsize=6.5, zorder=9)
    ax.text(sx + 42, sy + h / 2, "m", ha="left", va="center", fontsize=6.5, zorder=9)

    handles = [Line2D([], [], color="white", lw=1.3, path_effects=dark_halo, label=area_label())] + extra_handles
    if show_survey:
        handles.append(Line2D([], [], ls="", marker="o", ms=2, mfc="#111111", mec="white", mew=0.25,
                              label=f"{TXT['survey']} (n = {len(grid_pts)})"))
    if show_samples:
        handles.append(Line2D([], [], ls="", marker="^", ms=5.5, mfc="white", mec="black", mew=0.8,
                              label=f"{sample_label or TXT['sample']} (n = {len(samples)})"))
    leg = ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.005, 0.085), fontsize=6.5,
                    frameon=True, fancybox=False, edgecolor="#555", framealpha=0.95, borderpad=0.6,
                    handlelength=2.2)
    leg.get_frame().set_linewidth(0.5)


def style_colorbar(cb, label):
    cb.ax.tick_params(labelsize=7, width=0.5, length=2.5)
    cb.outline.set_linewidth(0.5)
    cb.set_label(label, fontsize=8)


def out(name):
    """Path of a generated file in outputs/."""
    os.makedirs(OUT_DIR, exist_ok=True)
    return os.path.join(OUT_DIR, name)


def save(fig, name):
    for ext, kw in (("pdf", {"metadata": {"CreationDate": None}}), ("svg", {"metadata": {"Date": None}}), ("png", {"dpi": 600}),
                    ("tif", {"dpi": 600, "pil_kwargs": {"compression": "tiff_lzw"}})):
        try:
            fig.savefig(out(f"{name}.{ext}"), **kw)
        except OSError:
            print(f"WARNING: could not write {name}.{ext} (is it open in another program?)")
    plt.close(fig)

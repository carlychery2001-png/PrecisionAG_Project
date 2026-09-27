"""Generative artwork for the Cotalia report, drawn from the lot's own elevation surface.

  Informe_Etapa1/arte/portada_curvas.pdf  - cover: 2 cm contour lines of the lot (rotated so the lot lies
                                            horizontally), white on transparent, sampling points in copper
  Informe_Etapa1/arte/sello_curvas.pdf    - small square crop used on the back cover

Usage: python report_art.py
"""
import os

from matplotlib.transforms import Affine2D

from mapcore import *  # noqa: F401,F403  (surface Zfull on GX/GY, boundary bxy)

ART = os.path.join("Informe_Etapa1", "arte")
os.makedirs(ART, exist_ok=True)
STEP, INDEX = 0.02, 0.10          # contour spacing and emphasised lines (m)
COPPER = "#D07A45"

# Rotate the lot so its long side (V2 -> V1) runs left to right
d = bxy[0] - bxy[1]
theta = -np.degrees(np.arctan2(d[1], d[0]))
cx, cy = bxy.mean(axis=0)
rot = Affine2D().rotate_deg_around(cx, cy, theta)
bx_r = rot.transform(bxy)

Zc = np.ma.masked_where(~inside, Zfull)
lv = np.arange(np.floor(Zc.min() / STEP) * STEP, Zc.max() + STEP, STEP)
is_idx = np.isclose(np.mod(np.round(lv / STEP), round(INDEX / STEP)), 0)

plan = list(csv.DictReader(open(out("Sampling_Plan.csv"), encoding="utf-8-sig")))
sxy = rot.transform(np.array([[float(r["Easting_CRTM05"]), float(r["Northing_CRTM05"])] for r in plan]))


def draw(fig, ax, lw_scale=1.0, dots=True):
    tr = rot + ax.transData
    ax.contour(GX, GY, Zc, levels=lv[~is_idx], colors="white", linewidths=0.28 * lw_scale, alpha=0.42, transform=tr)
    ax.contour(GX, GY, Zc, levels=lv[is_idx], colors="white", linewidths=0.55 * lw_scale, alpha=0.9, transform=tr)
    ax.add_patch(Polygon(bx_r, closed=True, fill=False, ec="white", lw=0.9 * lw_scale, alpha=0.95))
    if dots:
        ax.scatter(sxy[:, 0], sxy[:, 1], s=9 * lw_scale, c=COPPER, lw=0, zorder=5)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.patch.set_alpha(0)


# Cover: wide panel
pad = 4.0
w_m, h_m = np.ptp(bx_r[:, 0]) + 2 * pad, np.ptp(bx_r[:, 1]) + 2 * pad
W = 190
fig = plt.figure(figsize=(W * MM, W * h_m / w_m * MM))
ax = fig.add_axes([0, 0, 1, 1])
draw(fig, ax)
ax.set_xlim(bx_r[:, 0].min() - pad, bx_r[:, 0].max() + pad)
ax.set_ylim(bx_r[:, 1].min() - pad, bx_r[:, 1].max() + pad)
fig.savefig(os.path.join(ART, "portada_curvas.pdf"), transparent=True, metadata={"CreationDate": None})
plt.close(fig)

# Back-cover seal: circular medallion of contour lines around the lowest (wettest) part of the lot
from matplotlib.patches import Circle

fig = plt.figure(figsize=(60 * MM, 60 * MM))
ax = fig.add_axes([0, 0, 1, 1])
i = np.unravel_index(np.argmin(Zc), Zc.shape)
lx, ly = rot.transform([[GX[i], GY[i]]])[0]
lx, ly = lx + 14, ly - 12                         # centre the medallion inside the lot
R = 20.0
ring = Circle((lx, ly), R, transform=ax.transData, fc="none", ec="none")
ax.add_patch(ring)
tr = rot + ax.transData
for lev, lw, al in ((lv[~is_idx], 0.26, 0.45), (lv[is_idx], 0.5, 0.95)):
    cs = ax.contour(GX, GY, Zc, levels=lev, colors="white", linewidths=lw, alpha=al, transform=tr)
    cs.set_clip_path(ring)
ax.add_patch(Circle((lx, ly), R, fc="none", ec=COPPER, lw=0.8))
ax.set_xlim(lx - R - 1, lx + R + 1)
ax.set_ylim(ly - R - 1, ly + R + 1)
ax.set_aspect("equal")
ax.axis("off")
fig.patch.set_alpha(0)
fig.savefig(os.path.join(ART, "sello_curvas.pdf"), transparent=True, metadata={"CreationDate": None})
plt.close(fig)
print("art written:", sorted(os.listdir(ART)), f"rotation {theta:.1f} deg, {len(lv)} contour levels")

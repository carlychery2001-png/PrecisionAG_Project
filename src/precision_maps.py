"""Precision-agriculture maps derived from the RTK elevation surface:

  1. Wetness map        - topographic wetness index (TWI) + closed depressions (ponding)
  2. Management zones   - fuzzy c-means on elevation, slope and TWI; number of zones by FPI/NCE
  3. Sampling map       - soil samples taken (M1-M8) by zone; optional extra points (stratified, maximin) if PROPOSE_NEW
  S1. Zone validity     - FPI and NCE vs number of zones (supplementary figure)

Outputs (English): Figure_Wetness_Map, Figure_Management_Zones, Figure_Sampling_Plan, Figure_S1_Zone_Validity
                   (.pdf/.svg/.png/.tif) + Sampling_Plan.csv, Management_Zones_Statistics.csv,
                   Figure_Precision_Maps_caption.txt
Outputs (Spanish): Figura_Humedad, Figura_Zonas_de_Manejo, Figura_Plan_de_Muestreo, Figura_S1_Validez_Zonas
                   + Plan_de_Muestreo.csv, Zonas_de_Manejo_Estadisticas.csv, Figura_Mapas_Precision_leyenda.txt

Usage  : python precision_maps.py [en|es]
"""
import heapq

from matplotlib.patches import Patch
from scipy import ndimage

from mapcore import *  # noqa: F401,F403  (data, surface, imagery, layout helpers)

OUT = {
    "en": dict(wet="Figure_Wetness_Map", zones="Figure_Management_Zones", plan="Figure_Sampling_Plan",
               valid="Figure_S1_Zone_Validity", plan_csv="Sampling_Plan.csv",
               stats_csv="Management_Zones_Statistics.csv", cap="Figure_Precision_Maps_caption.txt"),
    "es": dict(wet="Figura_Humedad", zones="Figura_Zonas_de_Manejo", plan="Figura_Plan_de_Muestreo",
               valid="Figura_S1_Validez_Zonas", plan_csv="Plan_de_Muestreo.csv",
               stats_csv="Zonas_de_Manejo_Estadisticas.csv", cap="Figura_Mapas_Precision_leyenda.txt"),
}[LANG]

HG = 1.0                  # hydrology / zoning grid (m)
MIN_DEPTH = 0.01          # depression threshold: fill depth (m)
MIN_DEP_AREA = 5.0        # smallest depression reported (m²)
MFD_P = 1.1               # Freeman (1991) flow-partition exponent
MIN_TANB = 0.001          # slope floor for TWI in flat cells
FEATURE_SIGMA = 2.0       # smoothing of slope and TWI before zoning (m)
FCM_M = 1.3               # fuzziness exponent (MZA default; Fridgen et al., 2004)
C_RANGE = range(2, 7)     # numbers of zones evaluated
MAJORITY_SIZE = 5         # majority filter window (cells)
MIN_PATCH = 100.0         # smallest zone patch kept (m²), about 1% of the field
PROPOSE_NEW = False       # the soil campaign sampled M1-M8 only (26 Sep 2026); True designs extra points
N_TOTAL = 20              # target number of soil samples if extra points are proposed
MIN_PER_ZONE = 4          # minimum samples per zone
CORE_ZONE, CORE_FIELD = 4.0, 5.0   # proposed samples keep this far from zone / field edges (m)
MIN_MEMBERSHIP = 0.6      # proposed samples only where the zone membership is at least this
TWI_ALPHA, ZONE_ALPHA = 1.0, 1.0   # opaque inside the field so map colours match their keys
SAMPLE_TINT = 0.55        # sampling map: zone colours mixed with white (opaque) as a light background
HYD_BUFFER = 1.5          # hydrology domain extends this far beyond the boundary, where survey points exist (m)
# Okabe-Ito colour-blind-safe palettes, ordered from lowest (wettest) to highest zone; blue vs orange
# keeps neighbouring zones distinct over the green imagery
ZONE_PALETTES = {
    2: ["#0072B2", "#E69F00"],
    3: ["#0072B2", "#56B4E9", "#E69F00"],
    4: ["#0072B2", "#56B4E9", "#E69F00", "#D55E00"],
    5: ["#0072B2", "#56B4E9", "#F0E442", "#E69F00", "#D55E00"],
    6: ["#0072B2", "#56B4E9", "#009E73", "#F0E442", "#E69F00", "#D55E00"],
}
LABEL_CLEARANCE = 8.0     # zone labels keep this far from sample points (m)
DISPLAY_UPSAMPLE = 4      # zone edges drawn on a 4x finer grid for smooth outlines
PROPOSED_COLOR = "#FFFFFF"

pc = " %" if LANG == "es" else "%"
NB = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
NB_DIST = [HG * (2 ** 0.5 if di and dj else 1.0) for di, dj in NB]

# ------------------------------------------------------------- grid
hx = np.arange(np.floor(bxy[:, 0].min()) - 3, np.ceil(bxy[:, 0].max()) + 3, HG) + HG / 2
hy = np.arange(np.floor(bxy[:, 1].min()) - 3, np.ceil(bxy[:, 1].max()) + 3, HG) + HG / 2
HX, HY = np.meshgrid(hx, hy)                      # rows run south -> north
H = rbf(np.c_[HX.ravel(), HY.ravel()] - origin).reshape(HX.shape)
fin = boundary.contains_points(np.c_[HX.ravel(), HY.ravel()]).reshape(HX.shape)
h_extent = (hx[0] - HG / 2, hx[-1] + HG / 2, hy[0] - HG / 2, hy[-1] + HG / 2)
cell_area = HG * HG


def poly_dist(x, y):
    """Exact distance (m) from points to the field boundary polygon (vectorized)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    best = np.full(x.shape, np.inf)
    for a, b in zip(bxy, np.roll(bxy, -1, axis=0)):
        d = b - a
        t = np.clip(((x - a[0]) * d[0] + (y - a[1]) * d[1]) / (d @ d), 0, 1)
        best = np.minimum(best, np.hypot(x - a[0] - t * d[0], y - a[1] - t * d[1]))
    return best


# Hydrology runs on the field plus a narrow buffer (survey points reach up to ~1.5 m outside the polygon),
# so the field edge is not an artificial outlet; results are reported for the field only.
hmask = fin | (poly_dist(HX, HY) <= HYD_BUFFER)


def extend(a, mask):
    """Fill cells outside `mask` with the nearest inside value (display only, clipped later)."""
    idx = ndimage.distance_transform_edt(~mask, return_distances=False, return_indices=True)
    return a[tuple(idx)]


def cell_of(x, y):
    return int(np.clip(round((y - hy[0]) / HG), 0, len(hy) - 1)), int(np.clip(round((x - hx[0]) / HG), 0, len(hx) - 1))


# ------------------------------------------------------------- 1. hydrology
def priority_flood(dem, mask, eps=1e-5):
    """Depression filling with an epsilon gradient so flats drain (Barnes et al., 2014)."""
    filled = dem.copy()
    closed = ~mask
    ny, nx = dem.shape
    edge = mask & ~ndimage.binary_erosion(mask, structure=np.ones((3, 3)), border_value=0)
    heap = []
    for i, j in zip(*np.nonzero(edge)):
        heap.append((filled[i, j], i, j))
        closed[i, j] = True
    heapq.heapify(heap)
    while heap:
        e, i, j = heapq.heappop(heap)
        for di, dj in NB:
            a, b = i + di, j + dj
            if 0 <= a < ny and 0 <= b < nx and not closed[a, b]:
                closed[a, b] = True
                if filled[a, b] <= e:
                    filled[a, b] = e + eps
                heapq.heappush(heap, (filled[a, b], a, b))
    return filled


def mfd_accumulation(filled, mask):
    """Multiple-flow-direction contributing area, m² (Freeman, 1991)."""
    ny, nx = filled.shape
    acc = np.where(mask, cell_area, 0.0)
    cells = np.argwhere(mask)
    for i, j in cells[np.argsort(-filled[mask], kind="stable")]:
        zc = filled[i, j]
        w, tgt = [], []
        for (di, dj), d in zip(NB, NB_DIST):
            a, b = i + di, j + dj
            if 0 <= a < ny and 0 <= b < nx and mask[a, b] and filled[a, b] < zc:
                w.append(((zc - filled[a, b]) / d) ** MFD_P)
                tgt.append((a, b))
        if w:
            tot = sum(w)
            for wk, (a, b) in zip(w, tgt):
                acc[a, b] += acc[i, j] * wk / tot
    return acc


filled = priority_flood(H, hmask)
depth = np.where(fin, filled - H, 0.0)
acc = mfd_accumulation(filled, hmask)
gy_, gx_ = np.gradient(H, HG)
tanb = np.hypot(gx_, gy_)
with np.errstate(divide="ignore"):   # acc is 0 outside the hydrology domain; those cells are discarded
    twi = np.where(fin, np.log((acc / HG) / np.maximum(tanb, MIN_TANB)), np.nan)

dep_lab, n_dep = ndimage.label((depth >= MIN_DEPTH) & fin, structure=np.ones((3, 3)))
depressions = []
for k in range(1, n_dep + 1):
    m = dep_lab == k
    a = m.sum() * cell_area
    if a < MIN_DEP_AREA:
        continue
    i, j = np.unravel_index(np.argmax(np.where(m, depth, -1)), depth.shape)
    depressions.append(dict(area=a, max_depth=depth[m].max(), volume=depth[m].sum() * cell_area,
                            x=hx[j], y=hy[i], mask=m))
depressions.sort(key=lambda d: -d["volume"])


# ------------------------------------------------------------- 2. management zones
def smooth_masked(a, mask, sigma):
    num_ = ndimage.gaussian_filter(np.where(mask, a, 0.0), sigma)
    den = ndimage.gaussian_filter(mask.astype(float), sigma)
    return np.where(mask, num_ / np.maximum(den, 1e-9), np.nan)


slope_pct = 100 * tanb
f_elev = H
f_slope = smooth_masked(slope_pct, fin, FEATURE_SIGMA / HG)
f_twi = smooth_masked(np.nan_to_num(twi), fin, FEATURE_SIGMA / HG)
X = np.c_[f_elev[fin], f_slope[fin], f_twi[fin]]
X = (X - X.mean(0)) / X.std(0)


def fcm(X, c, n_init=10, tol=1e-6, max_iter=500):
    """Fuzzy c-means (Bezdek, 1981), Euclidean distance, best of n_init seeded starts."""
    rng = np.random.default_rng(42)
    best = None
    for _ in range(n_init):
        U = rng.random((len(X), c))
        U /= U.sum(1, keepdims=True)
        for _ in range(max_iter):
            Um = U ** FCM_M
            V = (Um.T @ X) / Um.sum(0)[:, None]
            D = np.fmax(np.linalg.norm(X[:, None, :] - V[None], axis=2), 1e-12)
            inv = D ** (-2 / (FCM_M - 1))
            U_new = inv / inv.sum(1, keepdims=True)
            done = np.abs(U_new - U).max() < tol
            U = U_new
            if done:
                break
        J = ((U ** FCM_M) * D ** 2).sum()
        if best is None or J < best[0]:
            best = (J, U, V)
    return best[1], best[2]


def validity(U):
    """Fuzziness performance index and normalized classification entropy (Odeh et al., 1992)."""
    n, c = U.shape
    F = (U ** 2).sum() / n
    Hn = -(U * np.log(np.fmax(U, 1e-300))).sum() / n
    return 1 - (c * F - 1) / (c - 1), Hn / (1 - c / n)


runs = {c: fcm(X, c) for c in C_RANGE}
fpi = {c: validity(runs[c][0])[0] for c in C_RANGE}
nce = {c: validity(runs[c][0])[1] for c in C_RANGE}
norm_ = lambda d: {k: (v - min(d.values())) / (max(d.values()) - min(d.values()) + 1e-12) for k, v in d.items()}
fn, nn = norm_(fpi), norm_(nce)
n_zones = min(C_RANGE, key=lambda c: fn[c] + nn[c])     # both indices as low as possible
U = runs[n_zones][0]

lab = np.zeros(H.shape, int)
lab[fin] = U.argmax(1) + 1
memb = np.zeros(H.shape + (n_zones,))
memb[fin] = U


def majority(lab, mask, size):
    counts = np.stack([ndimage.uniform_filter(((lab == k) & mask).astype(float), size)
                       for k in range(1, n_zones + 1)])
    out = counts.argmax(0) + 1
    return np.where(mask, out, 0)


def remove_small_patches(lab, mask, min_cells):
    changed = True
    while changed:
        changed = False
        for k in range(1, n_zones + 1):
            comp, n = ndimage.label(lab == k)
            for c in range(1, n + 1):
                m = comp == c
                if m.sum() >= min_cells:
                    continue
                ring = ndimage.binary_dilation(m, structure=np.ones((3, 3))) & ~m & mask
                if ring.any():
                    lab[m] = np.bincount(lab[ring]).argmax()
                    changed = True
    return lab


for _ in range(2):
    lab = majority(lab, fin, MAJORITY_SIZE)
lab = remove_small_patches(lab, fin, int(MIN_PATCH / cell_area))

# Renumber zones from lowest to highest mean elevation (Zone 1 = lowest, wettest)
order = sorted(range(1, n_zones + 1), key=lambda k: H[lab == k].mean() if (lab == k).any() else np.inf)
present = [k for k in order if (lab == k).any()]
remap = {old: new for new, old in enumerate(present, 1)}
lab = np.vectorize(lambda v: remap.get(v, 0))(lab)
memb = memb[..., [k - 1 for k in present]]
n_zones = len(present)
ZONE_COLORS = ZONE_PALETTES[n_zones]

zone_stats = []
for k in range(1, n_zones + 1):
    m = lab == k
    zone_stats.append(dict(zone=k, area=m.sum() * cell_area, pct=100 * m.sum() / fin.sum(),
                           elev=(H[m].mean(), H[m].std()), slope=(slope_pct[m].mean(), slope_pct[m].std()),
                           twi=(twi[m].mean(), twi[m].std()), memb=memb[m, k - 1].mean()))

# ------------------------------------------------------------- 3. sampling plan
def smooth_zone_grid():
    """Zone labels on a finer grid from smoothed, bilinearly upsampled zone indicators (display only)."""
    lab_ext = extend(lab, fin)
    ind = [ndimage.zoom(ndimage.gaussian_filter((lab_ext == k).astype(float), 1.0), DISPLAY_UPSAMPLE, order=1)
           for k in range(1, n_zones + 1)]
    fine = np.argmax(np.stack(ind), axis=0) + 1
    fx = np.linspace(hx[0], hx[-1], fine.shape[1])
    fy = np.linspace(hy[0], hy[-1], fine.shape[0])
    FXg, FYg = np.meshgrid(fx, fy)
    return FXg, FYg, fine


FX, FY, LAB_FINE = smooth_zone_grid()
fdx, fdy = FX[0, 1] - FX[0, 0], FY[1, 0] - FY[0, 0]
# Distance (m) from each fine cell to the drawn zone boundary (conservative by half a fine cell)
ZONE_EDGE = {k: ndimage.distance_transform_edt(LAB_FINE == k, sampling=(fdy, fdx)) - max(fdx, fdy) / 2
             for k in range(1, n_zones + 1)}


def fine_of(x, y):
    i = np.clip(np.round((np.asarray(y) - FY[0, 0]) / fdy).astype(int), 0, FY.shape[0] - 1)
    j = np.clip(np.round((np.asarray(x) - FX[0, 0]) / fdx).astype(int), 0, FX.shape[1] - 1)
    return i, j


def zone_at(x, y):
    return int(LAB_FINE[fine_of(x, y)])


def membership_at(x, y, k):
    return float(memb[cell_of(x, y)][k - 1])


def edge_to_zone(x, y, k):
    return float(ZONE_EDGE[k][fine_of(x, y)])


existing = [(n, x, y) for n, x, y, _ in samples]
existing_zone = {n: zone_at(x, y) for n, x, y in existing}
existing_memb = {n: membership_at(x, y, existing_zone[n]) for n, x, y in existing}
existing_edge = {n: edge_to_zone(x, y, existing_zone[n]) for n, x, y in existing}
# An existing sample counts toward its zone's quota only if the terrain there is typical of the zone
counted = {n: existing_memb[n] >= MIN_MEMBERSHIP for n, _, _ in existing}

field_d = poly_dist(HX, HY)
new_points = []
placed = [(x, y) for _, x, y in existing]
targets = {}
for s in zone_stats:
    k = s["zone"]
    have = sum(1 for n, z in existing_zone.items() if z == k and (counted[n] or not PROPOSE_NEW))
    targets[k] = (have, max(MIN_PER_ZONE, round(N_TOTAL * s["pct"] / 100)) if PROPOSE_NEW else have)
fi, fj = fine_of(HX, HY)
for k in (sorted(targets, key=lambda k: -(targets[k][1] - targets[k][0])) if PROPOSE_NEW else []):
    have, want = targets[k]
    cand = ((lab == k) & (LAB_FINE[fi, fj] == k) & (ZONE_EDGE[k][fi, fj] >= CORE_ZONE)
            & (field_d >= CORE_FIELD) & fin & (memb[..., k - 1] >= MIN_MEMBERSHIP))
    if cand.sum() < max(want - have, 1):
        raise RuntimeError(f"zone {k}: not enough candidate cells that meet the sampling constraints")
    cy, cx = np.nonzero(cand)
    cxy = np.c_[hx[cx], hy[cy]]
    for _ in range(max(0, want - have)):
        d = np.min(np.hypot(cxy[:, None, 0] - np.array(placed)[None, :, 0],
                            cxy[:, None, 1] - np.array(placed)[None, :, 1]), axis=1)
        pick = int(np.argmax(d))
        placed.append(tuple(cxy[pick]))
        new_points.append((k, *cxy[pick]))

# Number the new samples as a walking route: nearest-neighbour path from the south corner (next to the base)
route, rest, cur = [], list(new_points), bxy[0]
while rest:
    nxt = min(rest, key=lambda p: np.hypot(p[1] - cur[0], p[2] - cur[1]))
    route.append(nxt)
    rest.remove(nxt)
    cur = np.array(nxt[1:])
first = max(int(n[1:]) for n, _, _ in existing) + 1
proposed = [(f"M{first + i}", k, x, y) for i, (k, x, y) in enumerate(route)]
proposed_check = [(n, float(poly_dist(x, y)), edge_to_zone(x, y, k), membership_at(x, y, k)) for n, k, x, y in proposed]
assert all(fd >= CORE_FIELD and ze >= CORE_ZONE and mb >= MIN_MEMBERSHIP for _, fd, ze, mb in proposed_check)


# ------------------------------------------------------------- figures
def zone_label_points():
    """Most interior cell of each zone's largest patch, kept clear of the sample labels."""
    near_sample = np.zeros(lab.shape, bool)
    for _, x, y, _ in samples:
        near_sample |= np.hypot(HX - x, HY - y) < LABEL_CLEARANCE
    out = {}
    for k in range(1, n_zones + 1):
        comp, n = ndimage.label(lab == k)
        big = np.argmax(np.bincount(comp.ravel())[1:]) + 1
        d = ndimage.distance_transform_edt(comp == big) * ~near_sample
        i, j = np.unravel_index(np.argmax(d), d.shape)
        out[k] = (hx[j], hy[i])
    return out


def tint(color, t):
    rgb = np.array(mpl.colors.to_rgb(color))
    return tuple(rgb + (1 - rgb) * t)


def draw_zones(ax, clip, colors):
    zlev = np.arange(0.5, n_zones + 1)
    zcmap = ListedColormap(colors)
    zf = ax.contourf(FX, FY, LAB_FINE, levels=zlev, cmap=zcmap, norm=BoundaryNorm(zlev, zcmap.N))
    zl = ax.contour(FX, FY, LAB_FINE, levels=zlev[1:-1], colors="white", linewidths=1.0)
    for cs in (zf, zl):
        cs.set_clip_path(clip)


def zone_name(k):
    return f"{TXT['zone']} {k}"


es = LANG == "es"
contour_levels = levels  # contours from mapcore (MINOR interval)

# 1. Wetness map --------------------------------------------------------------
fig, ax, cax, clip = new_map()
twi_ext = extend(np.nan_to_num(twi), fin)
vmin, vmax = np.nanpercentile(twi[fin], [2, 98])
im = ax.imshow(twi_ext, extent=h_extent, origin="lower", cmap="YlGnBu", vmin=vmin, vmax=vmax,
               interpolation="bilinear", alpha=TWI_ALPHA, zorder=1)
im.set_clip_path(clip)
cl = ax.contour(GX, GY, Z, levels=contour_levels, colors="#3a3a3a", linewidths=0.35, alpha=0.7)
cl.set_clip_path(clip)
dl = ax.contour(HX, HY, extend(depth, fin) * fin, levels=[MIN_DEPTH], colors="#D55E00", linewidths=1.1,
                linestyles="--")
dl.set_clip_path(clip)
for d in depressions:
    ax.annotate((f"Depresión potencial\nmáx. ~{num(d['max_depth'] * 100, 0)} cm · ~{num(d['volume'], 0)} m³" if es else
                 f"Potential depression\nmax. ~{d['max_depth'] * 100:.0f} cm · ~{d['volume']:.0f} m³"),
                (d["x"], d["y"]), xytext=(14, -18), textcoords="offset points",
                fontsize=6.5, color="#7a2a00", fontweight="bold", path_effects=halo, zorder=9,
                arrowprops=dict(arrowstyle="-", lw=0.6, color="#7a2a00"))
decorate(ax, [
    Line2D([], [], color="#3a3a3a", lw=0.5, label=f"{TXT['contour']} ({MINOR * 100:.0f} cm)"),
    Line2D([], [], color="#D55E00", lw=1.1, ls="--",
           label="Depresión potencial (≥ 1 cm)" if es else "Potential depression (≥ 1 cm)"),
])
cb = fig.colorbar(im, cax=cax, extend="both")
cb.ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: num(v, 0)))
style_colorbar(cb, TXT["twi"])
cb.ax.text(0.5, -0.06, TXT["drier"], transform=cb.ax.transAxes, ha="center", va="top", fontsize=6.5, style="italic")
cb.ax.text(0.5, 1.06, TXT["wetter"], transform=cb.ax.transAxes, ha="center", va="bottom", fontsize=6.5, style="italic")
save(fig, OUT["wet"])

# 2. Management zones ---------------------------------------------------------
fig, ax, cax, clip = new_map()
cax.remove()
draw_zones(ax, clip, ZONE_COLORS[:n_zones])
for k, (x, y) in zone_label_points().items():
    ax.text(x, y, zone_name(k), ha="center", va="center", fontsize=8, fontweight="bold",
            path_effects=halo, zorder=9)
zone_handles = [Patch(fc=ZONE_COLORS[s["zone"] - 1], ec="white",
                      label=f"{zone_name(s['zone'])}: {thousands(s['area'])} m² ({num(s['pct'], 0)}{pc})")
                for s in zone_stats]
decorate(ax, zone_handles, show_survey=False)
save(fig, OUT["zones"])

# 3. Sampling plan -----------------------------------------------------------
fig, ax, cax, clip = new_map()
cax.remove()
light = [tint(c, SAMPLE_TINT) for c in ZONE_COLORS[:n_zones]]
draw_zones(ax, clip, light)
ax.scatter([p[2] for p in proposed], [p[3] for p in proposed], s=40, marker="o", fc=PROPOSED_COLOR,
           ec="black", lw=0.8, zorder=8)
for name, _, x, y in proposed:
    ax.annotate(name, (x, y), xytext=(4, 3), textcoords="offset points", fontsize=6.5,
                fontweight="bold", path_effects=halo, zorder=9)
if proposed:
    plan_handles = [Patch(fc=light[k - 1], ec="white",
                          label=(f"{zone_name(k)}: {targets[k][0]} existentes + "
                                 f"{len([p for p in proposed if p[1] == k])} propuestas" if es else
                                 f"{zone_name(k)}: {targets[k][0]} existing + "
                                 f"{len([p for p in proposed if p[1] == k])} proposed"))
                    for k in range(1, n_zones + 1)]
    plan_handles.append(Line2D([], [], ls="", marker="o", ms=6, mfc=PROPOSED_COLOR, mec="black", mew=0.8,
                               label=f"{TXT['proposed']} (n = {len(proposed)})"))
    taken_label = TXT["existing"]
else:
    plan_handles = [Patch(fc=light[k - 1], ec="white",
                          label=(f"{zone_name(k)}: {targets[k][0]} muestras" if es else
                                 f"{zone_name(k)}: {targets[k][0]} samples"))
                    for k in range(1, n_zones + 1)]
    taken_label = "Muestra de suelo tomada el 26/09/2026" if es else "Soil sample taken on 26 Sep 2026"
decorate(ax, plan_handles, show_survey=False, sample_label=taken_label)
save(fig, OUT["plan"])

# S1. Zone validity -----------------------------------------------------------
fig = plt.figure(figsize=(88 * MM, 62 * MM))
axv = fig.add_axes([0.16, 0.18, 0.80, 0.76])
cs_ = list(C_RANGE)
axv.plot(cs_, [fpi[c] for c in cs_], "o-", color="#0072B2", lw=1, ms=4, label="FPI")
axv.plot(cs_, [nce[c] for c in cs_], "s--", color="#D55E00", lw=1, ms=4, label="NCE")
axv.axvline(n_zones, color="#555", lw=0.6, ls=":")
axv.set_xticks(cs_)
axv.set_xlabel("Número de zonas" if es else "Number of zones")
axv.set_ylabel("Valor del índice" if es else "Index value")
axv.yaxis.set_major_locator(MultipleLocator(0.05))
axv.yaxis.set_major_formatter(FuncFormatter(lambda v, _: num(v, 2)))
axv.tick_params(labelsize=7)
axv.legend(fontsize=7, frameon=False)
for sp in ("top", "right"):
    axv.spines[sp].set_visible(False)
save(fig, OUT["valid"])

# ------------------------------------------------------------- tables
measured = {n: e for n, _, _, e in samples}   # measured heights (approx. a.s.l.)


def fmt_n(v, d):
    return num(v, d) if es else f"{v:.{d}f}"


with open(out(OUT["plan_csv"]), "w", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f, delimiter=";" if es else ",", lineterminator="\n")
    w.writerow(["Muestra", "Estado", "Zona", "Membresia_zona", "Cuenta_para_zona", "Este_CRTM05", "Norte_CRTM05",
                "Latitud_WGS84", "Longitud_WGS84", "Elevacion_aprox_msnm", "Fuente_elevacion"] if es else
               ["Sample", "Status", "Zone", "Zone_membership", "Counts_for_zone", "Easting_CRTM05", "Northing_CRTM05",
                "Latitude_WGS84", "Longitude_WGS84", "Elevation_approx_masl", "Elevation_source"])
    yes, no = ("sí", "no") if es else ("yes", "no")
    taken = ("tomada 26/09/2026" if es else "taken 2026-09-26") if not PROPOSE_NEW else ("existente" if es else "existing")
    rows = [(n, taken, existing_zone[n], x, y, measured[n],
             "medida RTK" if es else "RTK measured", yes if counted[n] else no) for n, x, y in existing]
    rows += [(n, "propuesta" if es else "proposed", k, x, y, float(rbf(np.array([[x, y]]) - origin)[0]),
              "superficie interpolada" if es else "interpolated surface", yes) for n, k, x, y in proposed]
    for n, status, k, x, y, elev, source, cnt in rows:
        lat, lon = imagery.crtm05_inverse(x, y)
        w.writerow([n, status, k, fmt_n(membership_at(x, y, k), 2), cnt, fmt_n(x, 3), fmt_n(y, 3),
                    fmt_n(float(lat), 7), fmt_n(float(lon), 7), fmt_n(elev, 2), source])

with open(out(OUT["stats_csv"]), "w", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f, delimiter=";" if es else ",", lineterminator="\n")
    w.writerow(["Zona", "Area_m2", "Area_pct", "Elevacion_media_msnm", "Elevacion_DE_m", "Pendiente_media_pct",
                "Pendiente_DE_pct", "TWI_medio", "TWI_DE", "Membresia_media", "Muestras_existentes_contadas",
                "Muestras_propuestas"] if es else
               ["Zone", "Area_m2", "Area_pct", "Elevation_mean_masl", "Elevation_SD_m", "Slope_mean_pct",
                "Slope_SD_pct", "TWI_mean", "TWI_SD", "Membership_mean", "Existing_samples_counted", "Proposed_samples"])
    for s in zone_stats:
        k = s["zone"]
        vals = [s["area"], s["pct"], *s["elev"], *s["slope"], *s["twi"], s["memb"]]
        digs = [0, 1, 2, 2, 2, 2, 2, 2, 2]
        w.writerow([k] + [fmt_n(v, d) for v, d in zip(vals, digs)] +
                   [targets[k][0], len([p for p in proposed if p[1] == k])])

# ------------------------------------------------------------- captions
z_lines = "\n".join(
    f"  {zone_name(s['zone'])}: {thousands(s['area'])} m² ({num(s['pct'], 0)}{pc}); "
    + (f"elevación {num(s['elev'][0], 2)} ± {num(s['elev'][1], 2)} m s. n. m. aprox.; pendiente {num(s['slope'][0])} ± {num(s['slope'][1])}{pc}; "
       f"TWI {num(s['twi'][0])} ± {num(s['twi'][1])}; membresía media {num(s['memb'], 2)}"
       if es else
       f"elevation {s['elev'][0]:.2f} ± {s['elev'][1]:.2f} m a.s.l. (approx.); slope {s['slope'][0]:.1f} ± {s['slope'][1]:.1f}%; "
       f"TWI {s['twi'][0]:.1f} ± {s['twi'][1]:.1f}; mean membership {s['memb']:.2f}")
    for s in zone_stats)
ex_by_zone = ", ".join(f"{zone_name(k)}: {', '.join(n for n, z in existing_zone.items() if z == k) or ('ninguna' if es else 'none')}"
                       for k in range(1, n_zones + 1))
atypical = [n for n in counted if not counted[n]]
border = [n for n in existing_edge if existing_edge[n] < 2.0 and counted[n]]
atyp_txt = ", ".join(f"{n} ({num(existing_memb[n], 2)})" for n in atypical) or ("ninguna" if es else "none")
border_txt = ", ".join(border) or ("ninguna" if es else "none")
dep_txt = "; ".join(
    (f"~{num(d['area'], 0)} m², profundidad máx. ~{num(d['max_depth'] * 100, 0)} cm, volumen ~{num(d['volume'], 0)} m³"
     if es else f"~{d['area']:.0f} m², max. depth ~{d['max_depth'] * 100:.0f} cm, volume ~{d['volume']:.0f} m³")
    for d in depressions) or ("ninguna" if es else "none")
fpi_txt = ", ".join(f"c = {c}: FPI {num(fpi[c], 3)}, NCE {num(nce[c], 3)}" for c in C_RANGE)
edge_share = 100 * np.mean(poly_dist(HX[fin], HY[fin]) <= 3.0)
credit = (f"Imagen de fondo: Esri World Imagery ({meta['SOURCE']} {meta['SOURCE_INFO']}, {img_date_txt}); fuente: Esri, {meta['SOURCE']}, "
          f"Earthstar Geographics y la comunidad de usuarios SIG." if es else
          f"Background: Esri World Imagery ({meta['SOURCE']} {meta['SOURCE_INFO']}, {img_date_txt}); source: Esri, {meta['SOURCE']}, "
          f"Earthstar Geographics, and the GIS User Community.")
prop_fd = min((fd for _, fd, _, _ in proposed_check), default=float("nan"))
prop_ze = min((ze for _, _, ze, _ in proposed_check), default=float("nan"))
min_memb = min(existing_memb.values())
if not proposed:
    z3_es = (f"Figura Z3. Mapa de muestreo de suelos: {len(existing)} muestras tomadas el 26 de septiembre de 2026 en M1–M{len(existing)} "
             f"({ex_by_zone}), con coordenadas RTK centimétricas. Todas tienen una membresía de al menos {num(min_memb, 2)} en su zona; "
             f"muestras a menos de 2 m del límite entre zonas: {border_txt}. Las muestras están en análisis de laboratorio. "
             f"Coordenadas en Plan_de_Muestreo.csv. {credit}")
    z3_en = (f"Figure Z3. Soil sampling map: {len(existing)} samples taken on 26 September 2026 at M1–M{len(existing)} "
             f"({ex_by_zone}), with centimetre RTK coordinates. All have a membership of at least {min_memb:.2f} in their zone; "
             f"samples less than 2 m from a zone border: {border_txt}. The samples are being analysed in the laboratory. "
             f"Coordinates in Sampling_Plan.csv. {credit}")
    note_es = ("Con 8 muestras no es posible ajustar un semivariograma confiable; los mapas de las propiedades del suelo se "
               "apoyarán en las zonas y en las covariables del terreno (elevación y TWI). Para volver a los puntos con RTK, la "
               "base debe instalarse sobre la misma marca con las mismas coordenadas del levantamiento; una nueva posición "
               "autónoma desplazaría todos los puntos varios metros.")
    note_en = ("With 8 samples a reliable semivariogram cannot be fitted; soil property maps will rely on the zones and on the "
               "terrain covariates (elevation and TWI). To return to the points with RTK, set the base up on the same mark with "
               "the same coordinates used in the survey; a new autonomous base position would shift every point by metres.")
else:
    z3_es = (f"Figura Z3. Plan de muestreo por zona: muestras existentes ({ex_by_zone}) y {len(proposed)} muestras propuestas "
             f"({proposed[0][0]}–{proposed[-1][0]}), a ≥ {num(CORE_FIELD, 0)} m del límite (mínimo {num(prop_fd)} m) y ≥ "
             f"{num(CORE_ZONE, 0)} m del borde de zona (mínimo {num(prop_ze)} m). Coordenadas en Plan_de_Muestreo.csv. {credit}")
    z3_en = (f"Figure Z3. Sampling plan by zone: existing samples ({ex_by_zone}) and {len(proposed)} proposed samples "
             f"({proposed[0][0]}–{proposed[-1][0]}), ≥ {CORE_FIELD:g} m from the boundary (minimum {prop_fd:.1f} m) and ≥ "
             f"{CORE_ZONE:g} m from the zone border (minimum {prop_ze:.1f} m). Coordinates in Sampling_Plan.csv. {credit}")
    note_es = "Las muestras propuestas maximizan la distancia mínima a todas las muestras (diseño maximin)."
    note_en = "Proposed samples maximize the minimum distance to all samples (maximin design)."

if es:
    caption = f"""Figura Z1. Índice topográfico de humedad (TWI) del lote. Valores altos (azul) indican donde el agua superficial \
tiende a concentrarse; valores bajos (amarillo), zonas que drenan. Escala de color: percentiles 2–98 del lote. La línea discontinua \
delimita una depresión cerrada potencial con profundidad de llenado ≥ {num(MIN_DEPTH * 100, 0)} cm ({dep_txt}); su profundidad es \
comparable a la incertidumbre de la superficie (RMSE {num(rmse * 100)} cm), por lo que es indicativa y conviene verificarla en campo \
tras una lluvia. El TWI dentro de ~3 m de los bordes pendiente arriba ({num(edge_share, 0)} % del lote está a menos de 3 m del borde) \
está sesgado a valores bajos, porque no se incluye el área de aporte exterior al lote. Curvas de nivel cada {MINOR * 100:.0f} cm. {credit}

Figura Z2. Zonas de manejo provisionales, delimitadas solo con variables del terreno (elevación, pendiente y TWI) mediante \
c-medias difuso; deben validarse con los análisis de suelo. Se seleccionaron {n_zones} zonas (mínimos de FPI y NCE; Figura S1). \
La Zona 1 es, en promedio, la más baja y topográficamente más húmeda.
{z_lines}
(Pendiente y TWI de las estadísticas: valores sin suavizar en malla de {num(HG, 0)} m.) {credit}

{z3_es}

Figura S1. Índice de desempeño de la difusividad (FPI) y entropía de clasificación normalizada (NCE) según el número de zonas \
({fpi_txt}). La línea punteada marca el número seleccionado; ambos índices son mínimos con {n_zones} zonas (la NCE tiende a \
favorecer pocas zonas).

Nota metodológica. La superficie de elevación (spline de placa delgada; ver Figura X) se evaluó en una malla de {num(HG, 0)} m. \
El análisis hidrológico se realizó sobre el lote más una franja de {num(HYD_BUFFER)} m alrededor (donde aún hay puntos medidos), \
de modo que el borde del lote no actúe como salida artificial; fuera de esa franja el agua sale libremente, por lo que el tamaño \
de la depresión es un límite inferior. Las depresiones se rellenaron con Priority-Flood con gradiente épsilon (Barnes et al., 2014); \
la profundidad de llenado indica posibles zonas de encharcamiento (almacenamiento en depresiones). El área de contribución se \
calculó con dirección de flujo múltiple (Freeman, 1991; p = {num(MFD_P)}) sobre la superficie rellenada. TWI = ln(a / tan β) \
(Beven y Kirkby, 1979), donde a es el área de contribución específica (área / ancho de celda) y β la pendiente local de la superficie \
sin rellenar (tan β ≥ {num(MIN_TANB, 3)}). Para la zonificación, la pendiente y el TWI se suavizaron con un filtro gaussiano \
(σ = {num(FEATURE_SIGMA, 0)} m); las tres variables se estandarizaron y se agruparon con c-medias difuso (Bezdek, 1981; distancia \
euclidiana; m = {num(FCM_M)}), evaluando de {min(C_RANGE)} a {max(C_RANGE)} zonas con FPI y NCE (Odeh et al., 1992; Fridgen et al., 2004). \
El mapa de zonas se depuró aplicando dos veces un filtro de mayoría de {MAJORITY_SIZE} × {MAJORITY_SIZE} m y fusionando parches \
menores a {MIN_PATCH:.0f} m². {note_es}

Referencias
Barnes, R., Lehman, C., Mulla, D., 2014. Priority-flood: An optimal depression-filling and watershed-labeling algorithm for digital elevation models. Computers & Geosciences 62, 117–127.
Beven, K.J., Kirkby, M.J., 1979. A physically based, variable contributing area model of basin hydrology. Hydrological Sciences Bulletin 24, 43–69.
Bezdek, J.C., 1981. Pattern Recognition with Fuzzy Objective Function Algorithms. Plenum Press, Nueva York.
Freeman, T.G., 1991. Calculating catchment area with divergent flow based on a regular grid. Computers & Geosciences 17, 413–422.
Fridgen, J.J., Kitchen, N.R., Sudduth, K.A., Drummond, S.T., Wiebold, W.J., Fraisse, C.W., 2004. Management Zone Analyst (MZA): Software for subfield management zone delineation. Agronomy Journal 96, 100–108.
Odeh, I.O.A., McBratney, A.B., Chittleborough, D.J., 1992. Soil pattern recognition with fuzzy-c-means: Application to classification and soil-landform interrelationships. Soil Science Society of America Journal 56, 505–516.
"""
else:
    caption = f"""Figure Z1. Topographic wetness index (TWI) of the field. High values (blue) mark where surface water tends to \
concentrate; low values (yellow) mark areas that drain. Colour stretch: 2nd–98th percentile of the field. The dashed line outlines \
a potential closed depression with fill depth ≥ {MIN_DEPTH * 100:.0f} cm ({dep_txt}); its depth is comparable to the surface \
uncertainty (RMSE {rmse * 100:.1f} cm), so it is indicative and should be checked in the field after rain. TWI within ~3 m of the \
upslope edges ({edge_share:.0f}% of the field lies within 3 m of the edge) is biased low because contributing area outside the field \
is not included. {MINOR * 100:.0f} cm contours are overlaid. {credit}

Figure Z2. Provisional management zones delineated from terrain variables only (elevation, slope and TWI) by fuzzy c-means; \
they must be validated with the soil analyses. {n_zones} zones were selected (minima of FPI and NCE; Figure S1). \
Zone 1 is, on average, the lower and topographically wetter zone.
{z_lines}
(Slope and TWI in the statistics are unsmoothed values on the {HG:g} m grid.) {credit}

{z3_en}

Figure S1. Fuzziness performance index (FPI) and normalized classification entropy (NCE) versus number of zones ({fpi_txt}). \
The dotted line marks the selected number; both indices are lowest at {n_zones} zones (NCE tends to favour few zones).

Methods note. The elevation surface (thin-plate spline; see Figure X) was evaluated on a {HG:g} m grid. Hydrology was computed on the \
field plus a {HYD_BUFFER:g} m strip around it (where measured points still exist), so the field edge is not an artificial outlet; \
beyond that strip water leaves freely, so the depression size is a lower bound. Depressions were filled with the Priority-Flood \
algorithm with an epsilon gradient (Barnes et al., 2014); the fill depth indicates potential ponding (depression storage). \
Contributing area was computed with multiple flow direction routing (Freeman, 1991; p = {MFD_P}) on the filled surface. \
TWI = ln(a / tan β) (Beven and Kirkby, 1979), where a is the specific contributing area (area / cell width) and β the local slope of \
the unfilled surface (tan β ≥ {MIN_TANB}). For zoning, slope and TWI were smoothed with a Gaussian filter (σ = {FEATURE_SIGMA:g} m); \
the three variables were standardized and clustered with fuzzy c-means (Bezdek, 1981; Euclidean distance; m = {FCM_M}), evaluating \
{min(C_RANGE)} to {max(C_RANGE)} zones with FPI and NCE (Odeh et al., 1992; Fridgen et al., 2004). The zone map was cleaned by \
applying a {MAJORITY_SIZE} × {MAJORITY_SIZE} m majority filter twice and merging patches smaller than {MIN_PATCH:.0f} m². {note_en}

References
Barnes, R., Lehman, C., Mulla, D., 2014. Priority-flood: An optimal depression-filling and watershed-labeling algorithm for digital elevation models. Computers & Geosciences 62, 117–127.
Beven, K.J., Kirkby, M.J., 1979. A physically based, variable contributing area model of basin hydrology. Hydrological Sciences Bulletin 24, 43–69.
Bezdek, J.C., 1981. Pattern Recognition with Fuzzy Objective Function Algorithms. Plenum Press, New York.
Freeman, T.G., 1991. Calculating catchment area with divergent flow based on a regular grid. Computers & Geosciences 17, 413–422.
Fridgen, J.J., Kitchen, N.R., Sudduth, K.A., Drummond, S.T., Wiebold, W.J., Fraisse, C.W., 2004. Management Zone Analyst (MZA): Software for subfield management zone delineation. Agronomy Journal 96, 100–108.
Odeh, I.O.A., McBratney, A.B., Chittleborough, D.J., 1992. Soil pattern recognition with fuzzy-c-means: Application to classification and soil-landform interrelationships. Soil Science Society of America Journal 56, 505–516.
"""

with open(out(OUT["cap"]), "w", encoding="utf-8") as f:
    f.write(caption)

print("zones:", n_zones, {c: (round(fpi[c], 3), round(nce[c], 3)) for c in C_RANGE})
print("depressions:", [(round(d["area"]), round(d["max_depth"] * 100, 1), round(d["volume"], 2)) for d in depressions])
for s in zone_stats:
    print("zone", s["zone"], round(s["area"]), round(s["pct"], 1), [round(v, 2) for v in s["elev"]],
          [round(v, 2) for v in s["slope"]], [round(v, 2) for v in s["twi"]], round(s["memb"], 2))
print("existing:", {n: (existing_zone[n], round(existing_memb[n], 2), round(existing_edge[n], 1), counted[n]) for n in existing_zone})
print("proposed check (name, field dist, zone-edge dist, membership):",
      [(n, round(fd, 2), round(ze, 2), round(mb, 2)) for n, fd, ze, mb in proposed_check])

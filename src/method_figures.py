"""Methodology figures and numbers for the Etapa 1 report (Spanish):

  - Figura_Levantamiento     : how the field was surveyed (RTK base, rover path in time order, vertices, samples)
  - Figura_Validacion_MDE    : (a) leave-one-out RMSE versus the spline smoothing parameter, (b) LOOCV residuals
  - Informe_Etapa1/datos/metodo.tex          : survey statistics, LOOCV and zoning numbers as LaTeX macros
  - Informe_Etapa1/datos/interpolacion.tex   : rows of the interpolator comparison table (same LOOCV points)
  - Informe_Etapa1/datos/lambda.tex          : rows of the smoothing-parameter table
  - Informe_Etapa1/datos/validez.tex         : rows of the FPI / NCE table (from the precision_maps.py caption)

Usage  : python method_figures.py es
"""
import os
import re
from datetime import datetime, timedelta

from scipy.interpolate import CloughTocher2DInterpolator, LinearNDInterpolator
from scipy.spatial import cKDTree

from mapcore import *  # noqa: F401,F403  (data, surface, LOOCV, imagery, layout helpers)

REPORT = "Informe_Etapa1"
DATA = os.path.join(REPORT, "datos")
INK, COPPER, SLATE, RULE = "#0C1A17", "#C4622D", "#717B77", "#D8DBD6"
IDW_P, IDW_K = 2, 12


def dms_to_deg(s):
    """'10d12m28.5793sN' -> decimal degrees."""
    d, rest = s.split("d")
    m, rest = rest.split("m")
    val = float(d) + float(m) / 60 + float(rest[:-2]) / 3600
    return -val if rest[-1] in "SWO" else val


def read_rows(path):
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def tstamp(s):
    return datetime.strptime(s, "%Y-%m-%d %H:%M:%S.%f")


rows_pts, rows_area = read_rows(POINTS_CSV), read_rows(AREA_CSV)
reg = rows_area + rows_pts                                  # every RTK point taken (143)
fnum = lambda c: np.array([float(r[c]) for r in reg])

# ------------------------------------------------------------- 1. survey statistics
times = sorted(tstamp(r["Tiempo local"]) for r in reg)
occ = np.array([(tstamp(r["Hora de finalizacion"]) - tstamp(r["Hora de inicio"])).total_seconds() for r in reg])
epochs = np.array([int(r["Cuenta de epocas"]) for r in reg])
tilt = np.array([dms_to_deg(r["Angulo inclinado"] + "N") for r in reg if r["Angulo inclinado"]])
base_row = read_rows(BASE_CSV)[0]                          # RTK base station position
base_lat, base_lon = dms_to_deg(base_row["Latitud"]), dms_to_deg(base_row["Longitud"])
base_e, base_n = imagery.crtm05_forward(base_lat, base_lon)
base_H = float(base_row["Altitud"]) - GEOID_N
nn_grid = cKDTree(grid_pts).query(grid_pts, k=2)[0][:, 1]
t_grid = sorted(tstamp(r["Tiempo local"]) for r in rows_pts if r["PT Nom"].startswith("Pt"))
t_samp = sorted(tstamp(r["Tiempo local"]) for r in rows_pts if r["PT Nom"].startswith("M"))
t_vert = sorted(tstamp(r["Tiempo local"]) for r in rows_area)
dur_min = (times[-1].replace(second=0, microsecond=0) - times[0].replace(second=0, microsecond=0)).total_seconds() / 60
lat_c, _ = imagery.crtm05_inverse(*origin)
tile_m = 156543.03392 * np.cos(np.radians(float(lat_c))) / 2 ** img_zoom

# ------------------------------------------------------------- 2. interpolation: LOOCV of alternatives
n = len(z)


def loo(fun):
    return np.array([fun(np.arange(n) != i, i) for i in range(n)])


def idw(m, i):
    d, ii = cKDTree(xy[m]).query(xy[i], k=IDW_K)
    w = 1 / d ** IDW_P
    return (w @ z[m][ii]) / w.sum() - z[i]


e_tin = loo(lambda m, i: LinearNDInterpolator(xy[m], z[m])(xy[i:i + 1])[0] - z[i])
e_ct = loo(lambda m, i: CloughTocher2DInterpolator(xy[m], z[m])(xy[i:i + 1])[0] - z[i])
e_idw = loo(idw)
hull = ~np.isnan(e_tin)                     # points inside the convex hull of the others: all methods defined


def stats(e, m=None):
    e = e[m] if m is not None else e
    return dict(rmse=100 * np.sqrt(np.mean(e ** 2)), mae=100 * np.mean(np.abs(e)), bias=100 * np.mean(e),
                p90=100 * np.percentile(np.abs(e), 90, method="inverted_cdf"), nmas=100 * np.mean(np.abs(e) <= MINOR / 2))


methods = [("Spline de placa delgada, $\\lambda = %s$ (seleccionado)" % num(best, 0), cv[best], True),
           ("Spline de placa delgada, $\\lambda = 0$ (interpolación exacta)", cv[0.0], True),
           ("Red irregular de triángulos (TIN lineal)", e_tin, False),
           ("Clough-Tocher (cúbica por triángulos)", e_ct, False),
           (f"Distancia inversa ponderada (IDW, $p = {IDW_P}$, {IDW_K} vecinos)", e_idw, True)]

# ------------------------------------------------------------- 3. zoning validity (from the precision_maps.py caption)
cap = open(out("Figura_Mapas_Precision_leyenda.txt"), encoding="utf-8").read()
valid = [(int(c), f, h) for c, f, h in re.findall(r"c = (\d+): FPI (\d+,\d+), NCE (\d+,\d+)", cap)]

# ------------------------------------------------------------- 4. LaTeX data
sign = lambda v, d=2: ("+" if v >= 0 else "−") + num(abs(v), d)
m = dict(
    HoraInicio=t_vert[0].strftime("%H:%M"), HoraFin=times[-1].strftime("%H:%M"),
    HoraVertices=f"{t_vert[0]:%H:%M}–{t_vert[-1]:%H:%M}", HoraGrilla=f"{t_grid[0]:%H:%M}–{t_grid[-1]:%H:%M}",
    HoraMuestras=f"{t_samp[0]:%H:%M}–{t_samp[-1]:%H:%M}",
    DuracionLev=f"{int(dur_min // 60)}~h~{int(round(dur_min % 60)):02d}~min",
    Epocas=str(int(np.median(epochs))), NEpocasCinco=str(int((epochs == np.median(epochs)).sum())),
    Ocupacion=num(np.median(occ), 0),
    SatUsoMin=str(int(fnum("Calculo de satelites").min())), SatUsoMax=str(int(fnum("Calculo de satelites").max())),
    SatRastMin=str(int(fnum("Seguimiento de satelites").min())), SatRastMax=str(int(fnum("Seguimiento de satelites").max())),
    PdopMin=num(fnum("PDOP").min(), 2), HdopMax=num(fnum("HDOP").max(), 2), VdopMax=num(fnum("VDOP").max(), 2),
    HrmsMed=num(100 * np.median(fnum("HRMS")), 1), VrmsMed=num(100 * np.median(fnum("VRMS")), 1),
    BaseMin=num(fnum("Distancia a Ref").min(), 1), BaseMax=num(fnum("Distancia a Ref").max(), 0),
    BaseMed=num(np.median(fnum("Distancia a Ref")), 0),
    InclMed=num(np.median(tilt), 1), InclMax=num(tilt.max(), 1),
    Baston=num(float(reg[0]["Medicion de altura"]), 2), AlturaAnt=num(float(reg[0]["ANT"]), 3),
    EspMediana=num(np.median(nn_grid), 1), EspMin=num(nn_grid.min(), 1), EspMax=num(nn_grid.max(), 1),
    Densidad=num(len(grid_pts) / (area_m2 / 1e4), 0), DensidadTotal=num(len(pts) / (area_m2 / 1e4), 0),
    BaseE=num(float(base_e), 1), BaseN=num(float(base_n), 1),
    ZMedMin=num(z.min(), 2), ZMedMax=num(z.max(), 2),
    Suavizado=num(best, 0), MAE=num(mae * 100, 1), Sesgo=sign(bias * 100, 1),
    PNoventa=num(np.percentile(np.abs(err), 90, method="inverted_cdf") * 100, 1), NMASpct=num(100 * np.mean(np.abs(err) <= MINOR / 2), 0),
    NSSDA=num(1.96 * rmse * 100, 1), RMSEdos=num(rmse * 100, 2), RMSEajuste=num(100 * np.sqrt(np.mean((rbf(xy) - z) ** 2)), 1),
    NHull=str(int(hull.sum())), NFuera=str(int((~hull).sum())),
    MallaCeldas=thousands(inside.sum()), MallaFilas=str(GX.shape[0]), MallaCols=str(GX.shape[1]),
    PlanoBetaUno=sign(px * 100, 2), PlanoBetaDos=sign(py * 100, 2),
    ZoomImagen=str(img_zoom), PixelImagen=num(tile_m, 2), FechaImagen=f"{int(img_date[6:])} de {'enero febrero marzo abril mayo junio julio agosto septiembre octubre noviembre diciembre'.split()[int(img_date[4:6]) - 1]} de {img_date[:4]}",
    RMSELambdaCero=num(stats(cv[0.0])["rmse"], 1), RMSELambdaMil=num(stats(cv[1000.0])["rmse"], 1),
    RMSETPSHull=num(stats(cv[best], hull)["rmse"], 1), RMSETINHull=num(stats(e_tin, hull)["rmse"], 1),
    RMSEIDWHull=num(stats(e_idw, hull)["rmse"], 1), RMSECTHull=num(stats(e_ct, hull)["rmse"], 1),
)
os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "metodo.tex"), "w", encoding="utf-8") as f:
    f.write("% Generado por method_figures.py a partir de los datos de campo; no editar a mano.\n")
    for k, v in m.items():
        f.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")

def write_rows(name, rows):
    """Table body with hairlines between rows; the closing rule sits inside the file, because \\input cannot be
    followed by \\noalign in an alignment."""
    with open(os.path.join(DATA, name), "w", encoding="utf-8") as f:
        f.write(" \\\\ \\filete\n".join(rows) + " \\\\\n\\bottomrule\n")


rows = []
for name, e, all_ok in methods:
    h, a = stats(e, hull), (stats(e, ~np.isnan(e)) if all_ok else None)
    star = "\\semi " if name.endswith("(seleccionado)") else ""
    rows.append(f"{star}{name} & {num(h['rmse'], 1)} & {num(h['mae'], 1)} & {sign(h['bias'], 1)} & "
                f"{num(h['nmas'], 0)}~\\% & {num(a['rmse'], 1) if a else '--'}")
write_rows("interpolacion.tex", rows)
rows = []
for s in candidates:
    st = stats(cv[s])
    lab = f"\\semi {num(s, 0)}" if s == best else num(s, 0)
    rows.append(f"{lab} & {num(st['rmse'], 2)} & {num(st['mae'], 2)} & {num(st['p90'], 1)} & {num(st['nmas'], 1)}~\\%")
write_rows("lambda.tex", rows)
SEMI = "\\semi "
write_rows("validez.tex", [f"{SEMI if c == 2 else ''}{c} & {fpi_} & {nce_}" for c, fpi_, nce_ in valid])

# ------------------------------------------------------------- 5. survey design map
fig, ax, cax, clip = new_map()
order = sorted(reg, key=lambda r: tstamp(r["Tiempo local"]))
E = np.array([float(r["Este"]) for r in order])
N_ = np.array([float(r["Norte"]) for r in order])
tmin = np.array([(tstamp(r["Tiempo local"]) - times[0]).total_seconds() / 60 for r in order])
is_pt = np.array([r in rows_pts and r["PT Nom"].startswith("Pt") for r in order])
is_v = np.array([r in rows_area for r in order])
for rr in (50, 100, 150):
    ax.add_patch(mpl.patches.Circle((base_e, base_n), rr, fill=False, ec="white", lw=0.7, ls=(0, (4, 3)), zorder=3,
                                    alpha=0.9))
    ang = np.radians(121)                     # along the lot's long axis, toward V3, so every label is on the map
    ax.text(base_e + rr * np.cos(ang), base_n + rr * np.sin(ang), f"{rr} m", fontsize=6.3, color="black",
            ha="center", va="center", zorder=8, clip_on=True,
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
gpath = is_pt
ax.plot(E[gpath], N_[gpath], color="white", lw=0.9, zorder=4, alpha=0.95, solid_joinstyle="round")
ax.plot(E[gpath], N_[gpath], color=INK, lw=0.45, zorder=4, alpha=0.9)
cmap_t = LinearSegmentedColormap.from_list("t", ["#FFE3CC", COPPER, "#5A2410"])
sc = ax.scatter(E[gpath], N_[gpath], c=tmin[gpath], cmap=cmap_t, s=16, ec="black", lw=0.4, zorder=6,
                vmin=0, vmax=tmin.max())
ax.scatter(E[is_v], N_[is_v], marker="s", s=34, fc="white", ec="black", lw=0.9, zorder=7)
for i, (x, y) in enumerate(bxy):
    ax.annotate(f"V{i + 1}", (x, y), xytext={0: (6, -9), 1: (-14, 3), 2: (-14, 3), 3: (5, 4)}[i],
                textcoords="offset points", fontsize=7,
                fontweight="bold", path_effects=halo, zorder=9)
ax.scatter([base_e], [base_n], marker="*", s=150, fc=COPPER, ec="black", lw=0.8, zorder=9)
ax.annotate("Base RTK", (base_e, base_n), xytext=(-50, 1), textcoords="offset points", fontsize=7.2,
            fontweight="bold", path_effects=halo, zorder=9)
cb = fig.colorbar(sc, cax=cax)
tick_min = np.arange(0, tmin.max() + 1, 30)
cb.set_ticks(tick_min)
cb.set_ticklabels([(times[0] + timedelta(minutes=float(t_))).strftime("%H:%M") for t_ in tick_min])
style_colorbar(cb, "Hora de medición (26/09/2026)")
decorate(ax, [
    Line2D([], [], color=INK, lw=0.8, label="Recorrido del móvil (orden de medición)"),
    Line2D([], [], ls="", marker="o", ms=4.5, mfc=COPPER, mec="black", mew=0.4,
           label=f"Punto de levantamiento (n = {int(is_pt.sum())})"),
    Line2D([], [], ls="", marker="s", ms=5, mfc="white", mec="black", mew=0.9, label="Vértice del lote (n = 4)"),
    Line2D([], [], ls="", marker="*", ms=9, mfc=COPPER, mec="black", mew=0.8, label="Estación base RTK"),
    Line2D([], [], color="white", lw=0.8, ls=(0, (4, 3)), path_effects=dark_halo,
           label="Distancia a la base: 50, 100 y 150 m"),
], show_survey=False, sample_label="Muestra de suelo")
leg = ax.get_legend()
leg.set_loc("upper right")
leg.set_bbox_to_anchor((0.93, 0.995), transform=ax.transAxes)
save(fig, "Figura_Levantamiento")

# ------------------------------------------------------------- 6. validation figure
plt.rcParams.update({"axes.edgecolor": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK})
fig = plt.figure(figsize=(WIDTH_MM * MM, 64 * MM))
axA = fig.add_axes([0.075, 0.2, 0.39, 0.72])
axB = fig.add_axes([0.575, 0.2, 0.40, 0.72])
xs = np.arange(len(candidates))
r_ = [100 * np.sqrt(np.mean(cv[s] ** 2)) for s in candidates]
a_ = [100 * np.mean(np.abs(cv[s])) for s in candidates]
axA.plot(xs, r_, color=INK, lw=1.2, marker="o", ms=4, mfc="white", mec=INK, mew=1.0, label="RMSE", zorder=3)
axA.plot(xs, a_, color=SLATE, lw=0.9, ls=(0, (3, 2)), marker="o", ms=3, mfc=SLATE, mec=SLATE, label="MAE", zorder=3)
kb = candidates.index(best)
axA.scatter([kb], [r_[kb]], s=70, fc=COPPER, ec=INK, lw=0.8, zorder=5)
axA.annotate(f"λ* = {num(best, 0)}\nRMSE = {num(r_[kb], 2)} cm", (kb, r_[kb]), xytext=(-8, -30),
             textcoords="offset points", fontsize=7, color=INK, ha="center",
             arrowprops=dict(arrowstyle="-", color=SLATE, lw=0.6))
axA.set_xticks(xs)
axA.set_xticklabels([num(s, 0) for s in candidates])
axA.set_xlabel("Parámetro de suavizado λ")
axA.set_ylabel("Error de validación cruzada (cm)")
axA.yaxis.set_major_formatter(FuncFormatter(lambda v, _: num(v, 1)))
axA.legend(frameon=False, fontsize=7, loc="upper center", ncol=2)
axA.set_ylim(5.5, 9.3)
axA.set_title("(a) Selección de λ por validación cruzada", fontsize=8, loc="left", color=INK)

e_cm = 100 * err
bins = np.arange(-36, 37, 3)
axB.axvspan(-100 * MINOR / 2, 100 * MINOR / 2, color=COPPER, alpha=0.10, lw=0, zorder=0)
for s_ in (-1, 1):
    axB.axvline(s_ * 100 * MINOR / 2, color=COPPER, lw=0.8, ls=(0, (3, 2)), zorder=1)
axB.hist(e_cm, bins=bins, color=INK, alpha=0.85, ec="white", lw=0.5, zorder=2)
gx_ = np.linspace(-36, 36, 300)
axB.plot(gx_, n * 3 * np.exp(-0.5 * ((gx_ - e_cm.mean()) / e_cm.std()) ** 2) / (e_cm.std() * np.sqrt(2 * np.pi)),
         color=COPPER, lw=1.1, zorder=3)
axB.text(0.015, 0.97, f"n = {n}\nRMSE {num(rmse * 100, 1)} cm\nMAE {num(mae * 100, 1)} cm\n"
         f"Sesgo {sign(bias * 100, 1)} cm\n{num(100 * np.mean(np.abs(err) <= MINOR / 2), 0)} % ≤ 12,5 cm",
         transform=axB.transAxes, fontsize=6.6, va="top", color=INK, linespacing=1.35)
axB.set_xlim(-42, 42)
axB.yaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
axB.text(100 * MINOR / 2 + 0.8, axB.get_ylim()[1] * 0.93, "½ equidistancia\n(NMAS)", fontsize=6.3, color=COPPER,
         va="top")
axB.set_xlabel("Residuo de validación cruzada: estimado − medido (cm)")
axB.set_ylabel("Número de puntos")
axB.set_title("(b) Distribución de los residuos, λ = " + num(best, 0), fontsize=8, loc="left", color=INK)
for a in (axA, axB):
    a.spines[["top", "right"]].set_visible(False)
    a.tick_params(labelsize=7)
    a.grid(axis="y", color=RULE, lw=0.5, zorder=0)
    a.set_axisbelow(True)
save(fig, "Figura_Validacion_MDE")

print("method outputs written:", sorted(os.listdir(DATA)))
print({k: m[k] for k in ("DuracionLev", "Epocas", "NEpocasCinco", "BaseMed", "EspMediana", "Densidad", "ZoomImagen",
                         "PixelImagen", "RMSETPSHull", "RMSETINHull", "RMSEIDWHull", "RMSECTHull", "NHull")})

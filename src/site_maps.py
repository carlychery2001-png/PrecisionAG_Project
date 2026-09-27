"""Site-description maps for the Etapa 1 report:

  - Location map: Costa Rica (panel a) and the Guácimo / EARTH surroundings (panel b), in CRTM05
  - Lot map: field polygon with vertices V1-V4, side lengths and azimuths, centroid and area

Also writes LaTeX snippets (Spanish run only) with the vertex/side tables and key site numbers, so the
report quotes values computed from the data.

Outputs (English): Figure_Location_Map.*, Figure_Lot_Polygon.*
Outputs (Spanish): Figura_Ubicacion.*, Figura_Lote_Poligono.*  + Informe_Etapa1/datos/*.tex

Usage  : python site_maps.py [en|es]
"""
import os
import re

from mapcore import *  # noqa: F401,F403  (data, surface, imagery, layout helpers)

OUT = {"en": {"loc": "Figure_Location_Map", "lot": "Figure_Lot_Polygon"},
       "es": {"loc": "Figura_Ubicacion", "lot": "Figura_Lote_Poligono"}}[LANG]
REPORT_DATA = os.path.join("Informe_Etapa1", "datos")
es = LANG == "es"

COUNTRY_EXTENT = (270_000, 690_000, 880_000, 1_250_000)   # CRTM05 (m), Costa Rica
COUNTRY_ZOOM, COUNTRY_RES = 7, 600.0
AREA_HALF = (1_200.0, 1_050.0)                            # half width / height of panel b (m)
AREA_ZOOM, AREA_RES = 15, 4.0
LOT_COLOR = "#FF3B30"

L = {
    "en": dict(east_km="Easting (km)", north_km="Northing (km)", site="Study site", lot="Lot",
               vertex="Boundary vertex", centroid="Centroid", side="Side length · grid azimuth",
               panel_a="(a) Costa Rica", panel_b="(b) Guácimo and EARTH University surroundings",
               ref="Esri World Boundaries and Places"),
    "es": dict(east_km="Este (km)", north_km="Norte (km)", site="Área de estudio", lot="Lote",
               vertex="Vértice del lote", centroid="Centroide", side="Longitud del lado · azimut de cuadrícula",
               panel_a="(a) Costa Rica", panel_b="(b) Alrededores de Guácimo y la Universidad EARTH",
               ref="Esri World Boundaries and Places"),
}[LANG]

# ------------------------------------------------------------- polygon geometry
V = bxy                                        # V1..V4 in survey order
cen = np.array([np.mean(V[:, 0]), np.mean(V[:, 1])])
# Area-weighted polygon centroid
xs_, ys_ = V[:, 0], V[:, 1]
cross = xs_ * np.roll(ys_, -1) - np.roll(xs_, -1) * ys_
A_signed = cross.sum() / 2
cen = np.array([((xs_ + np.roll(xs_, -1)) * cross).sum(), ((ys_ + np.roll(ys_, -1)) * cross).sum()]) / (6 * A_signed)
sides = []
for i in range(len(V)):
    a, b = V[i], V[(i + 1) % len(V)]
    d = b - a
    sides.append(dict(frm=i + 1, to=(i + 1) % len(V) + 1, length=float(np.hypot(*d)),
                      az=float(np.degrees(np.arctan2(d[0], d[1])) % 360), mid=(a + b) / 2, vec=d))


def bearing(az):
    """Quadrant bearing, e.g. N 56° 12' O (Spanish) / N 56° 12' W (English)."""
    west = "O" if es else "W"
    if az <= 90:
        ns, ew, ang = "N", "E", az
    elif az <= 180:
        ns, ew, ang = "S", "E", 180 - az
    elif az <= 270:
        ns, ew, ang = "S", west, az - 180
    else:
        ns, ew, ang = "N", west, 360 - az
    dd = int(ang)
    mm = int(round((ang - dd) * 60))
    if mm == 60:
        dd, mm = dd + 1, 0
    return f"{ns} {dd}° {mm:02d}' {ew}"


def dec(s):
    """Parse a Spanish-CSV number ('46,7')."""
    return float(s.replace(",", "."))


def coord(v, d=1):
    """Coordinate with thin-space thousands and the language's decimal separator, e.g. 543 974,9."""
    return f"{v:,.{d}f}".replace(",", " ").replace(".", DEC)


def dms(v, pos, neg, sd=2):
    h = pos if v >= 0 else neg
    v = round(abs(v) * 3600, sd)                 # whole value in rounded seconds, so 59.999" carries over
    d = int(v // 3600)
    m = int((v - d * 3600) // 60)
    s = v - d * 3600 - m * 60
    return f"{d}° {m:02d}' {num(s, sd)}\" {h}"


# Interior angles (degrees)
angles = []
for i in range(len(V)):
    p_, c_, n_ = V[i - 1], V[i], V[(i + 1) % len(V)]
    u1, u2 = p_ - c_, n_ - c_
    angles.append(float(np.degrees(np.arccos(np.dot(u1, u2) / np.linalg.norm(u1) / np.linalg.norm(u2)))))

# ------------------------------------------------------------- 1. location map
def geo_panel(ax, extent, zoom, res, unit_km, ticks, bar, bar_unit):
    rgb, _ = imagery.mosaic(extent, res=res, zoom=zoom)
    ref, _ = imagery.mosaic(extent, res=res, zoom=zoom, layer="reference")
    ext = (extent[0], extent[0] + rgb.shape[1] * res, extent[3] - rgb.shape[0] * res, extent[3])
    ax.imshow(rgb, extent=ext, origin="upper", interpolation="bilinear", zorder=0)
    ax.imshow(ref, extent=ext, origin="upper", interpolation="bilinear", zorder=1)
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_aspect("equal")
    k = 1000.0 if unit_km else 1.0
    fmt = FuncFormatter(lambda v, _: thousands(v / k))
    ax.xaxis.set_major_locator(MultipleLocator(ticks))
    ax.yaxis.set_major_locator(MultipleLocator(ticks))
    ax.xaxis.set_major_formatter(fmt)
    ax.yaxis.set_major_formatter(fmt)
    ax.tick_params(which="both", top=True, right=True, labelsize=6.5)
    plt.setp(ax.get_yticklabels(), rotation=90, va="center")
    ax.set_xlabel(L["east_km"] if unit_km else TXT["east"], fontsize=7.5)
    ax.set_ylabel(L["north_km"] if unit_km else TXT["north"], fontsize=7.5)
    # north arrow
    nx, ny = 0.92, 0.80
    ax.annotate("", xy=(nx, ny + 0.11), xytext=(nx, ny), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>,head_width=0.35,head_length=0.7", lw=1.1, color="black",
                                path_effects=[pe.withStroke(linewidth=2.6, foreground="white")]), zorder=10)
    ax.text(nx, ny + 0.12, "N", transform=ax.transAxes, ha="center", va="bottom", fontsize=9,
            fontweight="bold", path_effects=halo, zorder=10)
    # scale bar (lower left)
    w, hgt = extent[1] - extent[0], extent[3] - extent[2]
    sx, sy, h = extent[0] + 0.05 * w, extent[2] + 0.06 * hgt, 0.012 * hgt
    ax.add_patch(Rectangle((sx - 0.02 * w, sy - 0.025 * hgt), bar[-1] + 0.12 * w, h + 0.085 * hgt,
                           fc="white", ec="none", alpha=0.85, zorder=8))
    for i in range(len(bar) - 1):
        ax.add_patch(Rectangle((sx + bar[i], sy), bar[i + 1] - bar[i], h, fc="black" if i % 2 == 0 else "white",
                               ec="black", lw=0.5, zorder=9))
    for v in bar:
        ax.text(sx + v, sy + h + 0.012 * hgt, f"{v / (1000 if bar_unit == 'km' else 1):g}", ha="center",
                va="bottom", fontsize=6, zorder=9)
    ax.text(sx + bar[-1] + 0.02 * w, sy + h / 2, bar_unit, ha="left", va="center", fontsize=6, zorder=9)


fig = plt.figure(figsize=(WIDTH_MM * MM, 88 * MM))
axa = fig.add_axes([0.07, 0.13, 0.42, 0.80])
axb = fig.add_axes([0.575, 0.13, 0.42, 0.80])
geo_panel(axa, COUNTRY_EXTENT, COUNTRY_ZOOM, COUNTRY_RES, True, 100_000, [0, 50_000, 100_000], "km")
axa.scatter([cen[0]], [cen[1]], s=70, marker="*", fc=LOT_COLOR, ec="white", lw=0.8, zorder=11)
axa.annotate(L["site"], (cen[0], cen[1]), xytext=(7, 7), textcoords="offset points", fontsize=7,
             fontweight="bold", color="white", path_effects=[pe.withStroke(linewidth=2, foreground="black")], zorder=11)
axa.set_title(L["panel_a"], fontsize=8, loc="left", fontweight="bold")

b_ext = (cen[0] - AREA_HALF[0], cen[0] + AREA_HALF[0], cen[1] - AREA_HALF[1], cen[1] + AREA_HALF[1])
geo_panel(axb, b_ext, AREA_ZOOM, AREA_RES, False, 500, [0, 250, 500], "m")
axb.add_patch(Polygon(V, closed=True, fill=False, ec=LOT_COLOR, lw=1.4, zorder=11,
                      path_effects=[pe.withStroke(linewidth=2.6, foreground="white")]))
axb.annotate(f"{L['lot']} ({num(area_m2 / 1e4, 2)} ha)", (V[:, 0].max(), V[:, 1].max()), xytext=(6, 4),
             textcoords="offset points", fontsize=7, fontweight="bold", color="white",
             path_effects=[pe.withStroke(linewidth=2, foreground="black")], zorder=11)
axb.set_title(L["panel_b"], fontsize=8, loc="left", fontweight="bold")
fmt_m = FuncFormatter(lambda v, _: thousands(v))
axb.xaxis.set_major_formatter(fmt_m)
axb.yaxis.set_major_formatter(fmt_m)
# Frame of panel b drawn on panel a would be too small to see; the star marks it.
fig.text(0.995, 0.01, f"CRTM05 (EPSG:5367) · {TXT['img']}: Esri World Imagery; {L['ref']} · Esri, {meta['SOURCE']}, "
         f"Earthstar Geographics", ha="right", va="bottom", fontsize=5.5, color="#333")
save(fig, OUT["loc"])

# ------------------------------------------------------------- 2. lot polygon map
fig, ax, cax, clip = new_map()
cax.remove()
ax.add_patch(Polygon(V, closed=True, fc=LOT_COLOR, alpha=0.12, ec="none", zorder=2))
for i, (x, y) in enumerate(V):
    ax.scatter([x], [y], s=55, marker="o", fc="white", ec="black", lw=1.0, zorder=10)
    side_ = 1 if x > cen[0] else -1          # label beside the vertex (stays inside the map frame)
    if x - gx.min() < 10:
        side_ = 1
    elif gx.max() - x < 10:
        side_ = -1
    ax.annotate(f"V{i + 1}", (x, y), xytext=(side_ * 9, 0), textcoords="offset points",
                ha="left" if side_ > 0 else "right", va="center", fontsize=8.5, fontweight="bold", path_effects=halo, zorder=11)
for s in sides:
    ang = np.degrees(np.arctan2(s["vec"][1], s["vec"][0]))
    if ang > 90:
        ang -= 180
    if ang < -90:
        ang += 180
    normal = np.array([s["vec"][1], -s["vec"][0]]) / s["length"]
    if np.dot(s["mid"] - cen, normal) < 0:
        normal = -normal
    pos = s["mid"] - normal * 5.0                              # just inside the lot
    ax.text(pos[0], pos[1], f"{num(s['length'], 2)} m · Az {num(s['az'], 1)}°", rotation=ang,
            rotation_mode="anchor", ha="center", va="center", fontsize=6.8, fontweight="bold",
            path_effects=halo, zorder=11)
ax.scatter([cen[0]], [cen[1]], s=60, marker="+", c="black", lw=1.4, zorder=11)
ax.annotate(f"{L['centroid']}\nE {coord(cen[0])}\nN {coord(cen[1])}", (cen[0], cen[1]), xytext=(6, 6),
            textcoords="offset points", fontsize=6.5, fontweight="bold", path_effects=halo, zorder=11)
decorate(ax, [
    Line2D([], [], ls="", marker="o", ms=5.5, mfc="white", mec="black", mew=1.0, label=L["vertex"]),
    Line2D([], [], ls="", marker="+", ms=7, color="black", mew=1.4, label=L["centroid"]),
], show_survey=False, show_samples=False)
save(fig, OUT["lot"])

# ------------------------------------------------------------- 3. LaTeX data for the report (Spanish run)
if es:
    os.makedirs(REPORT_DATA, exist_ok=True)
    ell = {r["PT Nom"]: float(r["Elevacion"]) for r in csv.DictReader(open(AREA_CSV, encoding="utf-8"))}
    names = [p[0] for p in area]
    lat, lon = imagery.crtm05_inverse(V[:, 0], V[:, 1])
    with open(os.path.join(REPORT_DATA, "vertices.tex"), "w", encoding="utf-8") as f:
        for i, nm in enumerate(names):
            f.write(f"V{i + 1} & {num(V[i, 0], 3)} & {num(V[i, 1], 3)} & {dms(float(lat[i]), 'N', 'S')} & "
                    f"{dms(float(lon[i]), 'E', 'O')} & {num(ell[nm] - GEOID_N, 2)} & {num(angles[i], 1)}\\textdegree \\\\\n")
        f.write("\\bottomrule\n")   # rules must be inside the file: \input cannot be followed by \noalign
    with open(os.path.join(REPORT_DATA, "lados.tex"), "w", encoding="utf-8") as f:
        for s in sides:
            f.write(f"V{s['frm']}--V{s['to']} & {num(s['length'], 2)} & {num(s['az'], 2)}\\textdegree & "
                    f"{bearing(s['az'])} \\\\\n")
        f.write(f"\\midrule\n\\textbf{{Perímetro}} & \\textbf{{{num(perim, 2)}}} & & \\\\\n\\bottomrule\n")
    clat, clon = imagery.crtm05_inverse(cen[0], cen[1])

    # Register of every RTK point taken (annex B): boundary vertices, soil samples, survey points
    def rtk_rows(path):
        with open(path, encoding="utf-8", newline="") as f:
            return {r["PT Nom"]: r for r in csv.DictReader(f)}
    raw_pts, raw_area = rtk_rows(POINTS_CSV), rtk_rows(AREA_CSV)
    groups = [("Vértices del lote", [(f"V{i + 1}", raw_area[nm]) for i, nm in enumerate(names)]),
              ("Muestras de suelo", [(p[0], raw_pts[p[0]]) for p in samples]),
              ("Puntos de levantamiento", sorted(((k, r) for k, r in raw_pts.items() if k.startswith("Pt")),
                                                 key=lambda kr: int(kr[0][2:])))]
    reg = [r for _, g in groups for _, r in g]
    assert len(reg) == len(raw_pts) + len(raw_area) and all(r["Estado de solucion"] == "Fijo" for r in reg)
    with open(os.path.join(REPORT_DATA, "puntos.tex"), "w", encoding="utf-8") as f:
        f.write("% Generado por site_maps.py; no editar a mano. Se lee con \\CatchFileDef antes del longtable.\n")
        for gi, (title, g) in enumerate(groups):
            f.write(f"\\addlinespace[{7 if gi else 4}pt]\n")
            f.write(f"\\multicolumn{{8}}{{@{{}}l}}{{\\eyebrow{{{title} · {len(g)}}}}} \\\\*[3pt]\n")
            for j, (nm, r) in enumerate(g):
                e, n = float(r["Este"]), float(r["Norte"])
                la, lo = imagery.crtm05_inverse(e, n)
                # hairline every 5 points: by point number for Pt (after Pt5, Pt10, ...), by row for V and M
                block_end = (int(nm[2:]) % 5 == 0) if nm.startswith("Pt") else ((j + 1) % 5 == 0)
                last = j + 1 == len(g)
                f.write(f"{nm} & {num(e, 3)} & {num(n, 3)} & {dms(float(la), 'N', 'S', 3)} & "
                        f"{dms(float(lo), 'E', 'O', 3)} & {num(float(r['Elevacion']) - GEOID_N, 2)} & "
                        f"{num(float(r['HRMS']) * 100, 1)} & {num(float(r['VRMS']) * 100, 1)} "
                        + ("\\\\\n\\filete\n" if block_end and not last else "\\\\\n" if last else "\\\\*\n"))

    zst = list(csv.DictReader(open(out("Zonas_de_Manejo_Estadisticas.csv"), encoding="utf-8-sig"), delimiter=";"))
    plan = list(csv.DictReader(open(out("Plan_de_Muestreo.csv"), encoding="utf-8-sig"), delimiter=";"))
    cap = open(out("Figura_Mapas_Precision_leyenda.txt"), encoding="utf-8").read()
    dep = re.search(r"~(\d+) m², profundidad máx\. ~(\d+) cm, volumen ~(\d+) m³", cap)
    zin, sin = Zin.compressed(), Sin.compressed()
    macros = {
        "AreaLote": thousands(area_m2), "AreaHa": num(area_m2 / 1e4, 2), "Perimetro": num(perim, 1),
        "CentroE": num(cen[0], 1), "CentroN": num(cen[1], 1),
        "CentroLat": dms(float(clat), "N", "S"), "CentroLon": dms(float(clon), "E", "O"),
        "CentroLatDec": f"{float(clat):.5f}", "CentroLonDec": f"{float(clon):.5f}",
        "CentroLatDecCom": num(float(clat), 5).replace("-", "−"), "CentroLonDecCom": num(float(clon), 5).replace("-", "−"),
        "ElevMin": num(zin.min(), 2), "ElevMax": num(zin.max(), 2), "Relieve": num(zin.max() - zin.min(), 2),
        "GradPlano": num(plane_slope, 1), "AzPlano": f"{plane_az:.0f}",
        "PendMedia": num(sin.mean(), 1), "PendMax": num(sin.max(), 1), "PendPnoventaycinco": num(np.percentile(sin, 95), 1),
        "NPuntos": str(len(pts)), "NElev": str(len(grid_pts)), "NMuestras": str(len(samples)),
        "RMSE": num(rmse * 100, 1), "Geoide": num(GEOID_N, 2),
        "ZonaUnoArea": thousands(float(zst[0]["Area_m2"])), "ZonaUnoPct": num(dec(zst[0]["Area_pct"]), 0),
        "ZonaDosArea": thousands(float(zst[1]["Area_m2"])), "ZonaDosPct": num(dec(zst[1]["Area_pct"]), 0),
        "ZonaUnoPend": num(dec(zst[0]["Pendiente_media_pct"]), 1), "ZonaDosPend": num(dec(zst[1]["Pendiente_media_pct"]), 1),
        "ZonaUnoElev": zst[0]["Elevacion_media_msnm"], "ZonaDosElev": zst[1]["Elevacion_media_msnm"],
        "NPropuestas": str(sum(1 for r in plan if r["Estado"] == "propuesta")),
        "NTotalMuestras": str(len(plan)),
        "DepArea": dep.group(1) if dep else "--", "DepProf": dep.group(2) if dep else "--",
        "DepVol": dep.group(3) if dep else "--",
        "NPuntosTotal": str(len(reg)), "NVertices": str(len(groups[0][1])),
        "HrmsMin": num(min(float(r["HRMS"]) for r in reg) * 100, 1),
        "HrmsMax": num(max(float(r["HRMS"]) for r in reg) * 100, 1),
        "VrmsMin": num(min(float(r["VRMS"]) for r in reg) * 100, 1),
        "VrmsMax": num(max(float(r["VRMS"]) for r in reg) * 100, 1),
        "PdopMax": num(max(float(r["PDOP"]) for r in reg), 2),
    }
    with open(os.path.join(REPORT_DATA, "datos_sitio.tex"), "w", encoding="utf-8") as f:
        f.write("% Generado por site_maps.py a partir de los datos de campo; no editar a mano.\n")
        for k, v in macros.items():
            f.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
    print("LaTeX data written:", sorted(os.listdir(REPORT_DATA)))

print("centroid", cen, "sides", [(s["frm"], round(s["length"], 2), round(s["az"], 2)) for s in sides],
      "angles", [round(a, 1) for a in angles], "angle sum", round(sum(angles), 2))

"""Publication-quality contour and slope maps for the Etapa 1 lot (EARTH University).

Inputs : data/levantamiento_rtk_2026-09-26.csv (RTK points), data/vertices_lote_2026-09-26.csv (boundary)
Basemap: Esri World Imagery (see imagery.py), reprojected to CRTM05.
Outputs (English): Figure_Contour_Map.*, Figure_Slope_Map.*  (+ Figure_Contour_Map_caption.txt)
Outputs (Spanish): Figura_Curvas_de_Nivel.*, Figura_Pendientes.*  (+ Figura_Curvas_de_Nivel_leyenda.txt)
         each as .pdf, .svg, 600 dpi .png and .tif

Usage  : python publication_map.py [en|es]
"""
from mapcore import *  # noqa: F401,F403  (data, surface, imagery, layout helpers)

OUT = {
    "en": {"contour": "Figure_Contour_Map", "slope": "Figure_Slope_Map", "cap": "caption"},
    "es": {"contour": "Figura_Curvas_de_Nivel", "slope": "Figura_Pendientes", "cap": "leyenda"},
}[LANG]
SENS_SMOOTHING = (30.0, 300.0)   # LOOCV-equivalent smoothing values used to report slope sensitivity

# ------------------------------------------------------------- 1. contour map
fig, ax, cax, clip = new_map()
# Perceptually uniform, colour-blind-safe ramp (trimmed plasma, light = high); no greens, so it
# stands out against the vegetation in the imagery
cmap = LinearSegmentedColormap.from_list("elev", plt.cm.plasma(np.linspace(0.12, 0.97, 256)))
cf = ax.contourf(GX, GY, Z, levels=levels, cmap=cmap, norm=BoundaryNorm(levels, cmap.N), alpha=ELEV_ALPHA)
c_min = ax.contour(GX, GY, Z, levels=levels[~is_major], colors="#1a1a1a", linewidths=0.35, alpha=0.8)
c_maj = ax.contour(GX, GY, Z, levels=levels[is_major], colors="#1a1a1a", linewidths=0.9)
for cs in (cf, c_min, c_maj):
    cs.set_clip_path(clip)
label_contours(ax, c_maj, levels[is_major], 7, "bold")
label_contours(ax, c_min, levels[~is_major], 5.5, "normal")
decorate(ax, [
    Line2D([], [], color="#1a1a1a", lw=0.9, label=f"{TXT['index']} ({num(MAJOR)} m)"),
    Line2D([], [], color="#1a1a1a", lw=0.35, label=f"{TXT['inter']} ({MINOR * 100:.0f} cm)"),
])
cb = fig.colorbar(cf, cax=cax)
cb.set_ticks(levels)
cb.ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: level_label(v)))
style_colorbar(cb, TXT["elev"])
save(fig, OUT["contour"])

# ------------------------------------------------------------- 2. slope map
fig, ax, cax, clip = new_map()
# Colour-blind-safe sequential classes (ColorBrewer YlOrBr)
scmap = ListedColormap(plt.cm.YlOrBr(np.linspace(0.08, 0.85, len(SLOPE_BREAKS) - 1)))
snorm = BoundaryNorm(SLOPE_BREAKS, scmap.N)
sf = ax.contourf(GX, GY, np.ma.clip(Sskirt, 0, SLOPE_BREAKS[-1] - 1e-6), levels=SLOPE_BREAKS,
                 cmap=scmap, norm=snorm, alpha=FILL_ALPHA)
sc = ax.contour(GX, GY, Z, levels=levels[~is_major], colors="#3a3a3a", linewidths=0.35, alpha=0.7)
sc_maj = ax.contour(GX, GY, Z, levels=levels[is_major], colors="#3a3a3a", linewidths=0.8)
for cs in (sf, sc, sc_maj):
    cs.set_clip_path(clip)
label_contours(ax, sc_maj, levels[is_major], 6.5, "bold")

# Downslope arrows on a regular grid inside the field
ag = np.arange(-200, 200, ARROW_SPACING)
AX_, AY_ = np.meshgrid(ag + origin[0], ag + origin[1])
apts = np.c_[AX_.ravel(), AY_.ravel()]
apts = apts[[boundary.contains_point(p) and edge_dist(*p) > 3 for p in apts]]
ix = np.clip(np.round((apts[:, 0] - gx[0]) / GRID).astype(int), 0, len(gx) - 1)
iy = np.clip(np.round((apts[:, 1] - gy[0]) / GRID).astype(int), 0, len(gy) - 1)
u, v = -dzdx[iy, ix], -dzdy[iy, ix]
mag = np.hypot(u, v)
ax.quiver(apts[:, 0], apts[:, 1], u / mag, v / mag, color="#1f3b73", pivot="mid",
          scale=1 / 5.5, scale_units="xy", angles="xy", width=0.0022, headwidth=4, headlength=4.5,
          headaxislength=4, zorder=5)

decorate(ax, [
    Line2D([], [], color="#3a3a3a", lw=0.8, label=f"{TXT['contour']} ({MINOR * 100:.0f} cm; {num(MAJOR)} m)"),
    Line2D([], [], color="#1f3b73", lw=0.9, marker=">", ms=3.5, markevery=[1], label=TXT["flow"]),
])
cb = fig.colorbar(sf, cax=cax, ticks=SLOPE_BREAKS, spacing="uniform")
cb.ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: short(v)))
style_colorbar(cb, TXT["slope"])
save(fig, OUT["slope"])

# ------------------------------------------------------------- statistics for the captions
zin, sin = Zin.compressed(), Sin.compressed()
zmeas = z                                      # measured heights (approx. a.s.l.)
class_pct = [100 * np.mean((sin >= lo) & (sin < hi)) for lo, hi in zip(SLOPE_BREAKS[:-1], SLOPE_BREAKS[1:])]
class_pct[-1] += 100 * np.mean(sin >= SLOPE_BREAKS[-1])
pc = " %" if LANG == "es" else "%"   # Spanish spaces the percent sign
classes = "; ".join(f"{short(lo)}–{short(hi)}{pc}: {num(p, 0)}{pc}"
                    for lo, hi, p in zip(SLOPE_BREAKS[:-1], SLOPE_BREAKS[1:], class_pct))

abs_err = np.abs(err)
p90 = np.percentile(abs_err, 90)
nmas_ok = p90 <= MINOR / 2
nssda95 = 1.96 * rmse

# Slope sensitivity to the smoothing parameter (the LOOCV curve is flat around the optimum)
sens = {}
flat_in = np.c_[GX[inside], GY[inside]] - origin
for s_ in SENS_SMOOTHING:
    f_ = RBFInterpolator(xy, z, kernel="thin_plate_spline", smoothing=s_)
    Zs = np.full(GX.shape, np.nan)
    Zs[inside] = f_(flat_in)
    gy_s, gx_s = np.gradient(Zs, GRID)
    Ss = (100 * np.hypot(gx_s, gy_s))[inside]
    Ss = Ss[np.isfinite(Ss)]
    sens[s_] = (100 * np.sqrt(np.mean(loocv(s_) ** 2)), Ss.mean(), 100 * np.mean(Ss >= 3))
share3 = 100 * np.mean(sin >= 3)

meta_res = float(meta["RESOLUTION (M)"])
disp_res = 156543.03392 * np.cos(np.radians(10.2085)) / 2 ** img_zoom

if LANG == "es":
    caption = f"""Figura X. Mapa topográfico de curvas de nivel del lote de estudio (área en cuadrícula CRTM05 {thousands(area_m2)} m², \
{num(area_m2 / 1e4, 2)} ha; perímetro {perim:.0f} m). Equidistancia de {MINOR * 100:.0f} cm; curvas índice cada {num(MAJOR)} m. \
La superficie suavizada varía entre {num(zin.min(), 2)} y {num(zin.max(), 2)} m s. n. m. aprox. (puntos medidos: {num(zmeas.min(), 2)}–{num(zmeas.max(), 2)} m). \
El gradiente general del terreno (plano ajustado) es {num(plane_slope)} % hacia {plane_az:.0f}° (norte de cuadrícula). \
Los triángulos indican los {len(samples)} puntos de muestreo de suelo (M1–M{len(samples)}); los puntos negros, los {len(grid_pts)} puntos \
de elevación RTK-GNSS. Coordenadas en CRTM05 (EPSG:5367).

Figura Y. Mapa de pendientes locales, calculadas a partir de la superficie suavizada de spline de placa delgada (representan el \
terreno a la escala del espaciamiento del levantamiento, ~8 m). Las flechas indican la dirección de máxima pendiente (escurrimiento \
superficial); se superponen las curvas de nivel cada {MINOR * 100:.0f} cm. Pendiente local media {num(sin.mean())} % \
(mediana {num(np.median(sin))} %; percentil 95: {num(np.percentile(sin, 95))} %; máxima {num(sin.max())} %). \
Porcentaje del área por clase: {classes}. Los límites de 0,5; 1; 2 y 5 % corresponden a las clases de pendiente de FAO (2006, cuadro 7); \
las clases 1–1,5 %, 1,5–2 % y 2–3 % subdividen las clases de FAO. Con un suavizado de {SENS_SMOOTHING[0]:g} o {SENS_SMOOTHING[1]:g} \
(error de validación cruzada prácticamente igual: {num(sens[SENS_SMOOTHING[0]][0])} y {num(sens[SENS_SMOOTHING[1]][0])} cm) la pendiente media \
sería {num(sens[SENS_SMOOTHING[0]][1])} y {num(sens[SENS_SMOOTHING[1]][1])} %, y el área con pendiente ≥ 3 % sería {num(sens[SENS_SMOOTHING[0]][2], 0)} \
y {num(sens[SENS_SMOOTHING[1]][2], 0)} % (valor publicado: {num(share3, 0)} %).

En ambas figuras, la imagen de fondo es Esri World Imagery ({meta['SOURCE']} {meta['SOURCE_INFO']}, adquirida el {img_date_txt}; \
resolución nativa {num(meta_res, 2)} m, mostrada con teselas de ~{num(disp_res, 1)} m; exactitud horizontal declarada {meta['ACCURACY (M)']} m), \
reproyectada de Web Mercator a CRTM05 por remuestreo bilineal a {num(IMG_RES, 2)} m. Dentro del lote las capas de datos son opacas \
para que los colores coincidan con la leyenda. Fuente: Esri, {meta['SOURCE']}, Earthstar Geographics y la comunidad de usuarios SIG.

Nota metodológica. Las elevaciones se levantaron el 26 de septiembre de 2026 con receptores RTK-GNSS SingularXYZ X1 (base y \
móvil, enlace por radio UHF) en modo de solución fija (RMS horizontal ≈ 1,1–1,5 cm; RMS vertical ≈ 1,9–2,8 cm; PDOP ≤ 1,12), \
con bastón de 1,71 m y compensación de inclinación; el límite del lote se levantó por separado en sus cuatro vértices. La \
superficie se construyó con n = {len(pts)} puntos ({len(grid_pts)} puntos de levantamiento y {len(samples)} muestras). Se ajustó un \
spline de placa delgada (scipy RBFInterpolator, parámetro de suavizado = {best:g}, seleccionado por validación cruzada dejando uno \
fuera) y se evaluó en una malla de {num(GRID, 2)} m. La validación cruzada dio RMSE = {num(rmse * 100)} cm, MAE = {num(mae * 100)} cm \
y sesgo = {'+' if bias >= 0 else '−'}{num(abs(bias) * 100)} cm (NSSDA, 95 %: ±{num(nssda95 * 100, 0)} cm). El \
{num(100 * np.mean(abs_err <= MINOR / 2), 0)} % de los errores de validación (percentil 90: {num(p90 * 100)} cm) está dentro de \
la mitad de la equidistancia, {'lo que cumple' if nmas_ok else 'lo que NO cumple'} el criterio NMAS para curvas de \
{MINOR * 100:.0f} cm. El mayor residuo de validación es de {num(abs_err.max() * 100, 0)} cm. La pendiente es la magnitud del \
gradiente de la superficie suavizada (diferencias finitas centradas). Alturas: el colector entrega alturas elipsoidales WGS 84 \
(no se aplicó modelo geoidal); se convirtieron a elevación aproximada sobre el nivel del mar restando la ondulación geoidal \
EGM2008 N = {num(GEOID_N, 2)} m. Como la base se posicionó de forma autónoma, la elevación absoluta y las coordenadas absolutas \
tienen una incertidumbre de algunos metros; las diferencias de altura y las distancias dentro del lote no se ven afectadas.

Referencias
FAO, 2006. Guidelines for Soil Description, 4.ª ed. Organización de las Naciones Unidas para la Alimentación y la Agricultura, Roma.
Pavlis, N.K., Holmes, S.A., Kenyon, S.C., Factor, J.K., 2012. The development and evaluation of the Earth Gravitational Model 2008 (EGM2008). Journal of Geophysical Research: Solid Earth 117, B04406.
"""
else:
    caption = f"""Figure X. Topographic contour map of the study field (CRTM05 grid area {area_m2:,.0f} m², {area_m2 / 1e4:.2f} ha; \
perimeter {perim:.0f} m). Contour interval {MINOR * 100:.0f} cm; index contours every {MAJOR:.1f} m. The smoothed surface ranges \
from {zin.min():.2f} to {zin.max():.2f} m a.s.l. (approx.; measured points {zmeas.min():.2f}–{zmeas.max():.2f} m). The overall \
(plane-fit) gradient is {plane_slope:.1f}% toward {plane_az:.0f}° (grid north). Triangles mark the {len(samples)} soil-sampling \
locations (M1–M{len(samples)}); dots mark the {len(grid_pts)} RTK-GNSS elevation points. Coordinates in CRTM05 (EPSG:5367).

Figure Y. Map of local slope, computed from the smoothed thin-plate-spline surface (slopes represent terrain at about the ~8 m \
survey spacing). Arrows show the direction of steepest descent (surface runoff direction); {MINOR * 100:.0f} cm contours are \
overlaid. Mean local slope {sin.mean():.1f}% (median {np.median(sin):.1f}%; 95th percentile {np.percentile(sin, 95):.1f}%; \
maximum {sin.max():.1f}%). Share of field area by class: {classes}. Class limits at 0.5, 1, 2 and 5% follow FAO (2006, Table 7); \
the 1–1.5, 1.5–2 and 2–3% classes subdivide the FAO classes. With smoothing {SENS_SMOOTHING[0]:g} or {SENS_SMOOTHING[1]:g} \
(nearly identical cross-validation error: {sens[SENS_SMOOTHING[0]][0]:.1f} and {sens[SENS_SMOOTHING[1]][0]:.1f} cm) the mean slope \
would be {sens[SENS_SMOOTHING[0]][1]:.1f} and {sens[SENS_SMOOTHING[1]][1]:.1f}%, and the area at ≥ 3% slope \
{sens[SENS_SMOOTHING[0]][2]:.0f} and {sens[SENS_SMOOTHING[1]][2]:.0f}% (published value: {share3:.0f}%).

In both figures the background is Esri World Imagery ({meta['SOURCE']} {meta['SOURCE_INFO']}, acquired {img_date_txt}; native \
resolution {meta_res:.2f} m, displayed from ~{disp_res:.1f} m tiles; stated horizontal accuracy {meta['ACCURACY (M)']} m), \
reprojected from Web Mercator to CRTM05 by bilinear resampling to {IMG_RES:.2f} m. Inside the field the data layers are opaque so \
that map colours match the legend. Source: Esri, {meta['SOURCE']}, Earthstar Geographics, and the GIS User Community.

Methods note. Elevations were collected on 26 September 2026 with SingularXYZ X1 RTK-GNSS receivers (base and rover, UHF \
radio link) in fixed-solution mode (horizontal RMS ≈ 1.1–1.5 cm, vertical RMS ≈ 1.9–2.8 cm; PDOP ≤ 1.12) using a 1.71 m pole \
with tilt compensation; the field boundary was surveyed separately at its four vertices. The surface was built from \
n = {len(pts)} points ({len(grid_pts)} survey points and {len(samples)} soil samples). A thin-plate spline (scipy \
RBFInterpolator, smoothing parameter = {best:g}, selected by leave-one-out cross-validation) was fitted and evaluated on a \
{GRID:g} m grid. Cross-validation gave RMSE = {rmse * 100:.1f} cm, MAE = {mae * 100:.1f} cm and bias = {bias * 100:+.1f} cm \
(NSSDA 95%: ±{nssda95 * 100:.0f} cm). {100 * np.mean(abs_err <= MINOR / 2):.0f}% of the cross-validation errors (90th \
percentile {p90 * 100:.1f} cm) lie within half the contour interval, which {'meets' if nmas_ok else 'does NOT meet'} the NMAS \
criterion for {MINOR * 100:.0f} cm contours. The largest cross-validation residual is {abs_err.max() * 100:.0f} cm. Slope is \
the gradient magnitude of the smoothed surface (central finite differences). Heights: the controller reports WGS 84 \
ellipsoidal heights (no geoid model was applied); they were converted to approximate heights above sea level by subtracting \
the EGM2008 geoid height N = {GEOID_N:.2f} m. Because the base was positioned autonomously, absolute elevations and absolute \
coordinates are uncertain by a few metres; height differences and distances within the field are not affected.

References
FAO, 2006. Guidelines for Soil Description, 4th ed. Food and Agriculture Organization of the United Nations, Rome.
Pavlis, N.K., Holmes, S.A., Kenyon, S.C., Factor, J.K., 2012. The development and evaluation of the Earth Gravitational Model 2008 (EGM2008). Journal of Geophysical Research: Solid Earth 117, B04406.
"""

try:
    with open(out(f"{OUT['contour']}_{OUT['cap']}.txt"), "w", encoding="utf-8") as f:
        f.write(caption)
except OSError:
    print("WARNING: could not write caption file (is it open in another program?)")

print({s: round(float(np.sqrt(np.mean(e ** 2))) * 100, 2) for s, e in cv.items()}, "best", best)
print(caption)

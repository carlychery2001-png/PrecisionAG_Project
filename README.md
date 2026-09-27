# Agricultura de precisión en la Finca Académica de Cultivos, Universidad EARTH

Datos de campo, código y resultados del proyecto de la **Etapa 1** del curso IGA-407 Agricultura de Precisión
(Universidad EARTH, Costa Rica), preparado por el equipo **Cotalia**: Tyrone Leslie, Jessie Anastasia Mvula,
Carly Chery y Felix Wanjiku.

El 26 de septiembre de 2026 levantamos con RTK-GNSS un lote de 0,98 ha en Las Mercedes de Guácimo, Limón, y tomamos
ocho muestras de suelo. Con esos datos construimos el modelo de elevación del lote, las curvas de nivel, las
pendientes, la humedad topográfica, dos zonas de manejo y el mapa de muestreo. Todo lo que aparece en el informe se
puede regenerar desde los datos con un solo comando.

**Informe final:** [`Informe_Etapa1/Cotalia_Etapa1.pdf`](Informe_Etapa1/Cotalia_Etapa1.pdf)

*English summary: RTK-GNSS survey (143 points) and soil sampling of a 0.98 ha field at EARTH University, Costa Rica,
with the Python code that builds the elevation model (thin-plate spline, cross-validated), contour, slope, wetness,
management-zone and sampling maps, and the report. Run `python run_all.py` to regenerate everything.*

## Contenido

```
data/                         datos de campo
  levantamiento_rtk_2026-09-26.csv   131 puntos de levantamiento y 8 muestras de suelo (M1–M8)
  vertices_lote_2026-09-26.csv       4 vértices del límite del lote
  estacion_base_2026-09-26.csv       posición de la estación base RTK
src/                          código
  mapcore.py          lectura de datos, modelo de elevación y validación cruzada, diseño común de los mapas
  imagery.py          proyección CRTM05 e imagen satelital Esri reproyectada
  publication_map.py  curvas de nivel y pendientes
  precision_maps.py   humedad topográfica, zonas de manejo y muestreo de suelos
  site_maps.py        ubicación, polígono del lote y datos del informe
  method_figures.py   diseño del levantamiento, validación del modelo y tablas de la metodología
  report_art.py       ilustraciones de portada y contraportada
outputs/                      mapas (PDF y PNG, en español e inglés), tablas CSV y leyendas
Informe_Etapa1/               informe (fuente, datos generados, logos, fotos y PDF)
run_all.py                    regenera todo en orden
```

## Datos

Coordenadas en **CRTM05** (EPSG:5367), la proyección oficial de Costa Rica. Receptores SingularXYZ X1 (base y móvil,
enlace por radio UHF), solo soluciones fijas. La precisión centimétrica es relativa, entre los puntos del lote: la base
se posicionó de forma autónoma, así que las coordenadas y elevaciones absolutas tienen una incertidumbre de algunos
metros. Para volver a los puntos con RTK hay que instalar la base en la misma marca y con las mismas coordenadas.

Los dos archivos de puntos tienen las mismas columnas. Las principales son:

| Columna | Significado |
|---|---|
| `PT Nom` | Identificador del punto. En el levantamiento: `Pt…` (131 puntos) y `M1`…`M8` (muestras de suelo). En los vértices: `Pt1`…`Pt4`, que el código y el informe llaman V1–V4, en ese orden |
| `Este`, `Norte` | Coordenadas CRTM05 (m) |
| `Elevacion` | Altura elipsoidal WGS 84 (m); el código resta la ondulación del geoide EGM2008 (11,42 m) para obtener m s. n. m. |
| `Latitud`, `Longitud` | Coordenadas geográficas WGS 84 (grados, minutos y segundos) |
| `Tiempo local`, `Hora de inicio`, `Hora de finalizacion` | Hora de medición |
| `Estado de solucion` | Solución RTK (todas `Fijo`) |
| `HRMS`, `VRMS` | Precisión horizontal y vertical del punto según el receptor (m) |
| `PDOP`, `HDOP`, `VDOP` | Geometría satelital |
| `Calculo de satelites`, `Seguimiento de satelites` | Satélites usados y rastreados |
| `Cuenta de epocas` | Épocas promediadas por punto |
| `Medicion de altura`, `ANT` | Altura del bastón y del centro de fase de la antena (m) |
| `Angulo inclinado`, `Pitch`, `Roll`, `Yaw` | Inclinación del bastón (compensación de inclinación) |
| `Distancia a Ref` | Distancia a la estación base (m) |

El archivo de la estación base tiene su latitud, longitud y altura elipsoidal (m).

Las muestras de suelo (M1–M8) son compuestas: tres submuestras de 0 a 30 cm por punto, tomadas con barreno de tubo.
Los resultados de laboratorio se agregarán cuando estén disponibles.

## Cómo regenerar todo

1. Python 3.11 o superior (probado con 3.13):
   ```
   pip install -r requirements.txt
   ```
2. Para el informe en PDF: una distribución de LaTeX con XeLaTeX y `latexmk` (MiKTeX o TeX Live). Antes de la
   primera corrida instale las fuentes Inter, IBM Plex Mono y Fira Math, porque los mapas también usan Inter:
   `tlmgr install inter plex firamath` en TeX Live, o los paquetes `inter`, `plex` y `firamath` desde la consola de
   MiKTeX. Los mapas buscan Inter en la distribución de LaTeX; si no la encuentran, usan Arial.
3. Desde la carpeta del repositorio:
   ```
   python run_all.py
   ```
   La primera vez descarga las teselas de la imagen satelital del lote (se guardan en `imagery_cache/`). El proceso
   tarda unos minutos, sobre todo por la validación cruzada del modelo de elevación. Con `--no-report` se generan los
   mapas y tablas sin compilar el informe.

Cada script también puede correrse por separado, siempre desde la carpeta del repositorio. Sin argumento el idioma es
inglés: `python src/publication_map.py es` genera los mapas en español y `python src/publication_map.py en` en inglés.
`site_maps.py` y `method_figures.py` deben correrse con `es` para escribir los datos del informe, y `report_art.py` no
lleva argumento. El orden correcto está en `run_all.py`.

La imagen satelital se descarga de la versión vigente de Esri World Imagery (el informe usó la captura del 21 de
noviembre de 2025). Si Esri actualiza la imagen del lugar, una nueva corrida puede mostrar otro fondo y otra fecha.

## Resultados principales

- Lote de 9 845 m² (0,98 ha) y 426 m de perímetro; 143 puntos RTK-GNSS con precisión relativa centimétrica.
- Modelo de elevación con spline de placa delgada; suavizado elegido por validación cruzada (error 8,0 cm), que
  permite curvas cada 25 cm según el estándar NMAS.
- Relieve de 1,46 m y gradiente general de 1,3 % hacia el norte; una depresión donde puede acumularse agua.
- Dos zonas de manejo provisionales (c-medias difuso sobre elevación, pendiente y humedad topográfica) y ocho muestras
  de suelo, cuatro por zona.

## Créditos

- Imagen satelital: Esri World Imagery (fuente: Esri, Vantor, Earthstar Geographics y la comunidad de usuarios SIG).
  Las teselas no se incluyen en el repositorio; el código las descarga.
- Geoide: EGM2008 (Pavlis et al., 2012), valor obtenido con GeoidEval de GeographicLib.
- Los logotipos de `Informe_Etapa1/logo/stack/` pertenecen a sus dueños (Python Software Foundation, NumPy, SciPy,
  Matplotlib, OpenStreetMap Foundation, Esri y SingularXYZ) y solo identifican las herramientas usadas. El logotipo
  de OpenStreetMap se distribuye bajo CC BY-SA.
- Las referencias completas están en el informe.

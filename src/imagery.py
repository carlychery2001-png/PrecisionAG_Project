"""Esri World Imagery basemap (and the transparent Esri World Boundaries and Places reference layer),
reprojected to CRTM05 (EPSG:5367).

Tiles are cached in imagery_cache/ (reference tiles in imagery_cache/reference/) so they are downloaded only once.
Attribution required on any map: "Esri, Vantor, Earthstar Geographics, and the GIS User Community".
"""
import io
import json
import math
import os
import urllib.request

import numpy as np
from PIL import Image

TILE_URL = "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}"
REF_URL = ("https://server.arcgisonline.com/ArcGIS/rest/services/Reference/World_Boundaries_and_Places/"
           "MapServer/tile/{z}/{y}/{x}")
META_URL = ("https://services.arcgisonline.com/arcgis/rest/services/World_Imagery/MapServer/identify"
            "?geometry={lon},{lat}&geometryType=esriGeometryPoint&sr=4326&layers=all&tolerance=1"
            "&mapExtent={lon0},{lat0},{lon1},{lat1}&imageDisplay=400,400,96&returnGeometry=false&f=json")
CACHE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "imagery_cache")
HEADERS = {"User-Agent": "PrecisionAG-research-map/1.0"}

# ---------------------------------------------------------------- CRTM05 (WGS 84 ellipsoid, TM lon0 = -84, k0 = 0.9999)
A_ = 6378137.0
F_ = 1 / 298.257223563      # WGS 84 ellipsoid (EPSG:5367 CR05 datum)
E2 = F_ * (2 - F_)
EP2 = E2 / (1 - E2)
K0, LON0, FE, FN = 0.9999, math.radians(-84.0), 500000.0, 0.0


def _M(phi):
    return A_ * ((1 - E2 / 4 - 3 * E2 ** 2 / 64 - 5 * E2 ** 3 / 256) * phi
                 - (3 * E2 / 8 + 3 * E2 ** 2 / 32 + 45 * E2 ** 3 / 1024) * np.sin(2 * phi)
                 + (15 * E2 ** 2 / 256 + 45 * E2 ** 3 / 1024) * np.sin(4 * phi)
                 - (35 * E2 ** 3 / 3072) * np.sin(6 * phi))


def crtm05_forward(lat_deg, lon_deg):
    lat, lon = np.radians(lat_deg), np.radians(lon_deg)
    N = A_ / np.sqrt(1 - E2 * np.sin(lat) ** 2)
    T, C, A = np.tan(lat) ** 2, EP2 * np.cos(lat) ** 2, np.cos(lat) * (lon - LON0)
    x = FE + K0 * N * (A + (1 - T + C) * A ** 3 / 6 + (5 - 18 * T + T * T + 72 * C - 58 * EP2) * A ** 5 / 120)
    y = FN + K0 * (_M(lat) + N * np.tan(lat) * (A * A / 2 + (5 - T + 9 * C + 4 * C * C) * A ** 4 / 24
                                                  + (61 - 58 * T + T * T + 600 * C - 330 * EP2) * A ** 6 / 720))
    return x, y


def crtm05_inverse(x, y):
    """CRTM05 easting/northing -> latitude/longitude (degrees). Snyder (1987) eqs. 8-12 to 8-25."""
    e1 = (1 - math.sqrt(1 - E2)) / (1 + math.sqrt(1 - E2))
    mu = ((np.asarray(y) - FN) / K0) / (A_ * (1 - E2 / 4 - 3 * E2 ** 2 / 64 - 5 * E2 ** 3 / 256))
    phi1 = (mu + (3 * e1 / 2 - 27 * e1 ** 3 / 32) * np.sin(2 * mu) + (21 * e1 ** 2 / 16 - 55 * e1 ** 4 / 32) * np.sin(4 * mu)
            + (151 * e1 ** 3 / 96) * np.sin(6 * mu) + (1097 * e1 ** 4 / 512) * np.sin(8 * mu))
    C1, T1 = EP2 * np.cos(phi1) ** 2, np.tan(phi1) ** 2
    N1 = A_ / np.sqrt(1 - E2 * np.sin(phi1) ** 2)
    R1 = A_ * (1 - E2) / (1 - E2 * np.sin(phi1) ** 2) ** 1.5
    D = (np.asarray(x) - FE) / (N1 * K0)
    lat = phi1 - (N1 * np.tan(phi1) / R1) * (D * D / 2 - (5 + 3 * T1 + 10 * C1 - 4 * C1 * C1 - 9 * EP2) * D ** 4 / 24
                                           + (61 + 90 * T1 + 298 * C1 + 45 * T1 * T1 - 252 * EP2 - 3 * C1 * C1) * D ** 6 / 720)
    lon = LON0 + (D - (1 + 2 * T1 + C1) * D ** 3 / 6
                  + (5 - 2 * C1 + 28 * T1 - 3 * C1 * C1 + 8 * EP2 + 24 * T1 * T1) * D ** 5 / 120) / np.cos(phi1)
    return np.degrees(lat), np.degrees(lon)


# ---------------------------------------------------------------- Web Mercator tiles
def _tile_xy(lat, lon, z):
    """Fractional global pixel coordinates (256-px tiles) at zoom z."""
    n = 256 * 2 ** z
    lat_r = np.radians(lat)
    return (lon + 180) / 360 * n, (1 - np.log(np.tan(lat_r) + 1 / np.cos(lat_r)) / math.pi) / 2 * n


def _get(url):
    req = urllib.request.Request(url, headers=HEADERS)
    with urllib.request.urlopen(req, timeout=30) as r:
        return r.read()


def _tile(z, x, y, layer="imagery"):
    if layer == "imagery":
        path, url, mode = os.path.join(CACHE, str(z), str(x), f"{y}.jpg"), TILE_URL, "RGB"
    else:
        path, url, mode = os.path.join(CACHE, "reference", str(z), str(x), f"{y}.png"), REF_URL, "RGBA"
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        data = _get(url.format(z=z, x=x, y=y))           # download first: a failed request leaves no broken file
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path + ".part", "wb") as f:
            f.write(data)
        os.replace(path + ".part", path)
    return Image.open(path).convert(mode)


def _is_placeholder(img):
    """Esri serves a flat grey 'Map data not yet available' tile beyond its native zoom."""
    a = np.asarray(img, dtype=float)
    grey = np.abs(a[..., 0] - a[..., 1]).mean() < 2 and np.abs(a[..., 1] - a[..., 2]).mean() < 2
    return grey and a.std() < 12


def mosaic(extent, res=0.15, zoom=19, layer="imagery"):
    """Return (rgb/rgba array, zoom used) covering extent=(xmin, xmax, ymin, ymax) in CRTM05 at `res` m/pixel.
    For imagery, falls back to coarser zooms where Esri serves 'no data' placeholder tiles."""
    xmin, xmax, ymin, ymax = extent
    lat, lon = crtm05_inverse(np.array([xmin, xmax, xmin, xmax, (xmin + xmax) / 2, (xmin + xmax) / 2]),
                              np.array([ymin, ymin, ymax, ymax, ymin, ymax]))
    for z in range(zoom, max(zoom - 4, 0), -1):
        px, py = _tile_xy(lat, lon, z)
        tx0, tx1 = int(px.min() // 256), int(px.max() // 256)
        ty0, ty1 = int(py.min() // 256), int(py.max() // 256)
        tiles = {(tx, ty): _tile(z, tx, ty, layer) for tx in range(tx0, tx1 + 1) for ty in range(ty0, ty1 + 1)}
        if layer != "imagery" or not any(_is_placeholder(t) for t in tiles.values()):
            break
    big = Image.new("RGB" if layer == "imagery" else "RGBA", ((tx1 - tx0 + 1) * 256, (ty1 - ty0 + 1) * 256))
    for (tx, ty), t in tiles.items():
        big.paste(t, ((tx - tx0) * 256, (ty - ty0) * 256))
    src = np.asarray(big, dtype=float)

    # Inverse-map every output pixel (CRTM05) into the Mercator mosaic, bilinear sampling
    xs = np.arange(xmin + res / 2, xmax, res)
    ys = np.arange(ymax - res / 2, ymin, -res)
    XX, YY = np.meshgrid(xs, ys)
    la, lo = crtm05_inverse(XX, YY)
    sx, sy = _tile_xy(la, lo, z)
    sx, sy = sx - tx0 * 256 - 0.5, sy - ty0 * 256 - 0.5
    x0, y0 = np.floor(sx).astype(int), np.floor(sy).astype(int)
    fx, fy = (sx - x0)[..., None], (sy - y0)[..., None]
    x0 = np.clip(x0, 0, src.shape[1] - 2)
    y0 = np.clip(y0, 0, src.shape[0] - 2)
    out = (src[y0, x0] * (1 - fx) * (1 - fy) + src[y0, x0 + 1] * fx * (1 - fy)
           + src[y0 + 1, x0] * (1 - fx) * fy + src[y0 + 1, x0 + 1] * fx * fy)
    return np.clip(out, 0, 255).astype(np.uint8), z


def metadata(x, y):
    """Capture date, resolution, accuracy and source of the imagery at CRTM05 point (x, y); cached."""
    path = os.path.join(CACHE, "metadata.json")
    if os.path.exists(path):
        with open(path, encoding="utf-8") as f:
            cached = json.load(f)
        if cached:
            return cached
    lat, lon = (float(v) for v in crtm05_inverse(x, y))
    d = json.loads(_get(META_URL.format(lat=lat, lon=lon, lat0=lat - 1e-3, lon0=lon - 1e-3,
                                        lat1=lat + 1e-3, lon1=lon + 1e-3)))
    # The finest-resolution layer that has a capture date
    hits = [r["attributes"] for r in d.get("results", []) if r["attributes"].get("DATE (YYYYMMDD)")]
    found = min(hits, key=lambda a: float(a["RESOLUTION (M)"])) if hits else None
    if found is None:
        raise RuntimeError("Esri imagery metadata unavailable at the lot; try again later")
    os.makedirs(CACHE, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(found, f, indent=1)
    return found

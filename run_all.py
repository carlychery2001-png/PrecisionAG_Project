"""Regenerate every map, table and figure of the Etapa 1 report from the field data, then build the report PDF.

Usage (from the repository folder):
    python run_all.py              # maps (Spanish and English), tables, report data and the report PDF
    python run_all.py --no-report  # everything except the PDF

The first run downloads the Esri World Imagery tiles for the lot into imagery_cache/ (internet needed once).
"""
import os
import shutil
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.abspath(__file__))

# Order matters: site_maps.py and method_figures.py read the tables and captions written by precision_maps.py,
# and report_art.py reads outputs/Sampling_Plan.csv.
STEPS = [
    ("publication_map.py", "es"), ("publication_map.py", "en"),   # contour and slope maps
    ("precision_maps.py", "es"), ("precision_maps.py", "en"),     # wetness, management zones, sampling
    ("site_maps.py", "en"), ("site_maps.py", "es"),               # location and lot maps; es also writes report data
    ("method_figures.py", "es"),                                  # survey design, DEM validation, method tables
    ("report_art.py", None),                                      # cover and back-cover artwork
]


def run(cmd, cwd):
    print("\n>>", " ".join(cmd), flush=True)
    t0 = time.time()
    subprocess.run(cmd, cwd=cwd, check=True, env={**os.environ, "PYTHONIOENCODING": "utf-8"})
    print(f"   done in {time.time() - t0:.0f} s", flush=True)


def main():
    for script, lang in STEPS:
        run([sys.executable, os.path.join("src", script)] + ([lang] if lang else []), ROOT)
    if "--no-report" in sys.argv:
        return
    report = os.path.join(ROOT, "Informe_Etapa1")
    if not shutil.which("latexmk"):
        print("\nlatexmk not found: install MiKTeX or TeX Live, then run in Informe_Etapa1/:"
              "\n    latexmk -xelatex Cotalia_Etapa1.tex")
        return
    run(["latexmk", "-xelatex", "-interaction=nonstopmode", "-halt-on-error", "Cotalia_Etapa1.tex"], report)
    run(["latexmk", "-c", "Cotalia_Etapa1.tex"], report)
    print("\nReport:", os.path.join(report, "Cotalia_Etapa1.pdf"))


if __name__ == "__main__":
    main()

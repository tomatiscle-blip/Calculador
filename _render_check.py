"""Render de la vista de planta para inspección visual."""
import os
import sys
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
RAIZ = Path(__file__).resolve().parent
sys.path.insert(0, str(RAIZ))

from PySide6.QtWidgets import QApplication
from PySide6.QtGui import QImage, QColor

app = QApplication.instance() or QApplication([])

from app.visualizacion import VistaPlantaPorticos
from tests.test_vista_planta import EJES, NIVELES, ESTRUCTURA, LOSAS

vista = VistaPlantaPorticos()
vista.resize(1100, 700)
vista.set_datos(EJES, NIVELES, ESTRUCTURA)
vista.set_losas(LOSAS)
vista.set_portico_activo("Pórtico 1")
vista.set_nivel_activo(3.0)
vista.update()

img = QImage(1100, 700, QImage.Format.Format_ARGB32)
img.fill(QColor("#ffffff"))
vista.render(img)
img.save(str(RAIZ / "_render_planta.png"))
print("render OK")

# Verificación analítica: cada categoría debe tener píxeles en el render.
from collections import Counter

conteo = Counter()
for y in range(0, 700, 2):
    for x in range(0, 1100, 2):
        conteo[img.pixelColor(x, y).name()] += 1

categorias = {
    "eje": vista.COLOR_EJE,
    "portico_activo": vista.COLOR_PORTICO_ACTIVO,
    "portico_inactivo": vista.COLOR_PORTICO,
    "viga": vista.COLOR_VIGA,
    "muro": vista.COLOR_MURO,
    "columna": vista.COLOR_COLUMNA,
    "base": vista.COLOR_BASE,
    "apoyo": vista.COLOR_APOYO,
    "losa": "#059669",
}
print("--- cobertura por categoría ---")
for nombre, hexc in categorias.items():
    c = QColor(hexc).name()
    total = sum(v for k, v in conteo.items() if k.startswith(c[:7]))
    print(f"{nombre:18s} {hexc}  pixeles~{total}")

print("--- colores mas frecuentes ---")
for k, v in conteo.most_common(12):
    print(f"  {k}: {v}")

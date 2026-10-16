import os
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent
if str(RAIZ) not in sys.path:
    sys.path.insert(0, str(RAIZ))

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
from PySide6.QtWidgets import QApplication  # noqa: E402
from app.ejes import EditorMuroPorNiveles  # noqa: E402

app = QApplication(sys.argv)
ejes = {"x": [{"id": "x:a", "nombre": "A", "coord": 0.0}, {"id": "x:b", "nombre": "B", "coord": 5.0}], "y": []}
niveles = [{"nombre": "PB", "cota_m": 0.0}, {"nombre": "PA", "cota_m": 3.0}]
editor = EditorMuroPorNiveles(ejes, niveles)
editor.show()
sys.exit(app.exec())

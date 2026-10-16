import os
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent
if str(RAIZ) not in sys.path:
    sys.path.insert(0, str(RAIZ))

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication  # noqa: E402

app = QApplication.instance() or QApplication([])
from app.ejes import EditorPorticoPorNiveles, EditorMuroPorNiveles  # noqa: E402

ejes = {
    "x": [
        {"id": "x:1", "nombre": "X1", "coordenada_m": 0.0},
        {"id": "x:2", "nombre": "X2", "coordenada_m": 5.0},
    ],
    "y": [
        {"id": "y:1", "nombre": "Y1", "coordenada_m": 0.0},
        {"id": "y:5", "nombre": "Y5", "coordenada_m": 8.0},
    ],
}
niveles = [{"nombre": "PB", "cota_m": 0.0}, {"nombre": "PA", "cota_m": 3.0}]
muro = {
    "tipo": "muro",
    "referencia_planta": {"direccion": "x", "eje_id": "y:5"},
    "posicion_planta_m": 8.0,
    "x_inicio": 0.0,
    "x_fin": 5.0,
}

# Editor de PORTICO con un muro presente: debe dibujarlo como contexto.
ed = EditorPorticoPorNiveles(ejes, niveles, {"Muro 1": muro})
print("portico: muros en visor =", list(ed.vista.muros.keys()))
assert "Muro 1" in ed.vista.muros, "el muro debe verse en el visor de portico"
assert "Muro 1" in ed.estructura_existente, "el muro no debe borrarse"
ed.close()

# Editor de MURO con un portico presente: debe mostrarlo como contexto.
portico = {
    "tipo": "portico",
    "referencia_planta": {"direccion": "x", "eje_id": "y:5", "origen_longitudinal_m": 0.0},
    "posicion_planta_m": 8.0,
    "columnas": {"C0-1": {"x": 0.0}, "C0-2": {"x": 5.0}},
    "vigas": {"V0-1": {"tramos": [{"x_inicio": 0.0, "x_fin": 5.0}]}},
}
edm = EditorMuroPorNiveles(ejes, niveles, {"Portico 5": portico})
print("muro: porticos en visor =", list(edm.vista.porticos.keys()))
assert "Portico 5" in edm.vista.porticos, "el portico debe verse en el visor de muro"
edm.direccion.setCurrentIndex(edm.direccion.findData("y"))
edm.eje_transversal.setCurrentIndex(edm.eje_transversal.findData("y:5"))
cp = edm._portico_coplanar()
print("muro coplanar detectado =", cp[0] if cp else None)
assert cp and cp[0] == "Portico 5"
assert not edm.aviso_coplanar.isHidden(), "el aviso mixto debe estar visible"
assert "Portico 5" in edm.estructura, "el portico NO debe borrarse"
edm.close()
print("HUMO OK")

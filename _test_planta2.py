import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from calc import rutas
import importlib
import app.ejes as ejesModule
from PySide6.QtWidgets import QApplication
rutas.inicializar_obras()
rutas.asegurar_directorios()
app = QApplication.instance() or QApplication(sys.argv)
from app.visualizacion import VistaPlantaPorticos
v = VistaPlantaPorticos()
# load data like app
try:
    ejes = ejesModule.cargar_ejes()
    niveles = ejesModule.cargar_niveles()
    import json
    est = json.loads(rutas.cargar_estructura().__repr__()) if False else None
except Exception as e:
    print("err load", e)
# use rutas helpers
try:
    from calc import rutas as R
    est = R.cargar_estructura()
    print("estructura keys", list(est.keys()) if isinstance(est, dict) else type(est))
except Exception as e:
    print("err est", e)
try:
    v.set_datos(ejes, niveles, est)
    print("set_datos ok")
except Exception as e:
    print("set_datos err", repr(e))
v.resize(1100,750)
v.show()
print("shown", flush=True)
app.exec()
print("done", flush=True)

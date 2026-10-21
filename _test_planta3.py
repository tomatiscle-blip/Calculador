import faulthandler, sys
faulthandler.enable()
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from calc import rutas
import importlib
import app.ejes as E
from PySide6.QtWidgets import QApplication
rutas.inicializar_obras()
rutas.asegurar_directorios()
app = QApplication.instance() or QApplication(sys.argv)
from app.visualizacion import VistaPlantaPorticos
v = VistaPlantaPorticos()
ejes=E.cargar_ejes()
niveles=E.cargar_niveles()
est=rutas.cargar_estructura()
v.set_datos(ejes,niveles,est)
v.resize(1100,750)
v.show()
print("shown", flush=True)
app.exec()
print("done", flush=True)

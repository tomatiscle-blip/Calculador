import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from calc import rutas, pipeline
import importlib
from PySide6.QtWidgets import QApplication
rutas.inicializar_obras()
importlib.reload(pipeline)
rutas.asegurar_directorios()
app = QApplication.instance() or QApplication(sys.argv)
from app.principal import VentanaPrincipal
print("before VP", flush=True)
v = VentanaPrincipal()
print("VP created", flush=True)
v.resize(1000,700)
v.move(50,50)
v.show()
print("shown", flush=True)
app.exec()
print("exec done", flush=True)

#!/usrvenv .venv\Scripts\python.exe
import sys
sys.path.insert(0, '.')
from calc import rutas
rutas.inicializar_obras()
rutas.asegurar_directorios()
from PySide6.QtWidgets import QApplication
app = QApplication([])
from app.visualizacion import VistaPlantaPorticos
v = VistaPlantaPorticos()
v.show()
print('WINDOW_SHOWN', flush=True)
app.exec()
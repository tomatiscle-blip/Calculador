from PySide6.QtWidgets import QApplication
import sys
app = QApplication(sys.argv)
from app.visualizacion import VistaPlantaPorticos
w = VistaPlantaPorticos()
w.show()
print("shown")
app.exec()
print("done")

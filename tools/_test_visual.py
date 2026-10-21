import os
import sys
sys.path.insert(0, '.')
os.environ['QT_QPA_PLATFORM'] = 'offscreen'

from PySide6.QtWidgets import QApplication
app = QApplication.instance() or QApplication([])

from app.visualizacion import VistaPortico, VistaPlantaPorticos, VistaDiagramaMomentos

v = VistaPortico()
print("VistaPortico OK")

vp = VistaPlantaPorticos()
print("VistaPlantaPorticos OK")

vdm = VistaDiagramaMomentos()
print("VistaDiagramaMomentos OK")

print("All classes instantiated successfully")
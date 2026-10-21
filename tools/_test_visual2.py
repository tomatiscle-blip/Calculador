import os
import sys
sys.path.insert(0, '.')
os.environ['QT_QPA_PLATFORM'] = 'offscreen'

from PySide6.QtWidgets import QApplication
app = QApplication.instance() or QApplication([])

from app.visualizacion import VistaPortico

v = VistaPortico()
print("VistaPortico instantiated OK")

# Test actualizar with a non-existent portico
try:
    v.actualizar("test_portico")
    print("actualizar('test_portico') OK")
    print(f"tipo={v.tipo}, nombre={v.nombre!r}, columnas={len(v.columnas)}, vigas={len(v.vigas)}")
except Exception as e:
    print(f"actualizar failed: {type(e).__name__}: {e}")
    import traceback
    traceback.print_exc()

# Test with a muro
try:
    v.actualizar("muro_test")
    print("actualizar('muro_test') OK")
except Exception as e:
    print(f"actualizar('muro_test') failed: {type(e).__name__}: {e}")
    import traceback
    traceback.print_exc()

print("All tests done")
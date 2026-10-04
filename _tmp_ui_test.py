"""
Prueba temporal de la interfaz: levanta la ventana SIN mostrarla, refresca y
muestra en texto lo que quedó en la pestaña Inicio. No abre ninguna ventana.

    py _tmp_ui_test.py
"""

import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parent
sys.path.insert(0, str(RAIZ))

from PySide6.QtWidgets import QApplication, QTabWidget  # noqa: E402
from app.principal import VentanaPrincipal  # noqa: E402

app = QApplication.instance() or QApplication([])
v = VentanaPrincipal()
v.refrescar()

tabs = v.findChild(QTabWidget)
print("UI_OK")
print("tabs:", [tabs.tabText(i) for i in range(tabs.count())])
print("inicio_etiqueta:", v.etiqueta_inicio.text().replace("\n", " | "))
print("cargas:", v.lbl_cargas.text().replace("\n", " | "))
print("motor:", v.lbl_motor.text().replace("\n", " | "))
print("filas_envolvente:", v.tabla_inicio.rowCount())
print("encabezados:", [v.tabla_inicio.horizontalHeaderItem(j).text() for j in range(v.tabla_inicio.columnCount())])
for i in range(v.tabla_inicio.rowCount()):
    print("   fila", i, [v.tabla_inicio.item(i, j).text() for j in range(v.tabla_inicio.columnCount())])
v.close()

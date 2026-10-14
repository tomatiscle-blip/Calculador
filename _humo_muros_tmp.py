import os
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
from PySide6.QtWidgets import QApplication, QTabWidget
app = QApplication.instance() or QApplication([])

from app.principal import VentanaPrincipal
from app.muros_proyecto import EditorCapasMuro, PaginaMuros

v = VentanaPrincipal()
v.refrescar()

# 1) La pestaña Muros existe y está registrada.
tabs = v.centralWidget().findChild(QTabWidget)
titulos = [tabs.tabText(i) for i in range(tabs.count())]
print("pestanas:", titulos)
assert "6 · Muros" in titulos, "falta la pestaña 6 · Muros"
assert v.pagina_muros is not None
# La página de muros refresca sin romperse aunque no haya muros.
v.pagina_muros.refrescar()
print("pagina muros OK, filas:", v.pagina_muros.tabla.rowCount())

# 2) Diálogo de composición de capas funciona y calcula peso propio.
capas = [
    {"clave": "cemento_arena", "grupo": "Morteros", "espesor_m": 0.02},
    {"clave": "ladrillo_hueco_portante", "grupo": "Mamposteria", "espesor_m": 0.18},
    {"clave": "cal_arena", "grupo": "Morteros", "espesor_m": 0.02},
]
d = EditorCapasMuro("Muro 1", capas, 3.0, 7.9)
n_capas, altura, longitud, resultado = d.datos_guardar()
print("capas:", n_capas, "altura:", altura, "long:", longitud)
print("espesor_total:", round(resultado.get("espesor_total_m", 0), 3))
print("peso_m2:", round(resultado.get("peso_propio_kN_m2", 0), 3))
print("lineal:", round(resultado.get("peso_lineal_kN_m", 0), 3))
assert abs(resultado["espesor_total_m"] - 0.22) < 1e-6
assert abs(resultado["peso_propio_kN_m2"] - 2.92) < 1e-6
d.close()
print("HUMO MUROS OK")

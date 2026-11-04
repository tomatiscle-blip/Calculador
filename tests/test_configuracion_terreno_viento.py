import json
import os
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import (  # noqa: E402
    QApplication,
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QLineEdit,
    QPlainTextEdit,
    QTableWidget,
    QTableWidgetItem,
)

from app.cargas_proyecto import PaginaCargas  # noqa: E402
from app.principal import VentanaPrincipal  # noqa: E402
from calc import rutas  # noqa: E402

APP = QApplication.instance() or QApplication([])


class ConfiguracionVientoTests(unittest.TestCase):
    def test_cambiar_direccion_no_cambia_seleccion_de_porticos(self):
        direccion = QComboBox()
        direccion.addItem("X", "x")
        tabla = QTableWidget(2, 4)
        checks = []
        for fila, nombre in enumerate(("Pórtico X", "Pórtico Y")):
            tabla.setItem(fila, 0, QTableWidgetItem(nombre))
            check = QCheckBox()
            check.setChecked(True)
            tabla.setCellWidget(fila, 2, check)
            checks.append(check)
        pagina = SimpleNamespace(direccion_viento=direccion, tabla_viento_porticos=tabla)
        estructura = {
            "Pórtico X": {"referencia_planta": {"direccion": "x"}},
            "Pórtico Y": {"referencia_planta": {"direccion": "y"}},
        }
        with patch("app.cargas_proyecto.rutas.cargar_estructura", return_value=estructura):
            PaginaCargas._actualizar_porticos_compatibles_viento(pagina)

        self.assertTrue(checks[0].isEnabled())
        self.assertFalse(checks[1].isEnabled())
        self.assertTrue(checks[0].isChecked())
        self.assertTrue(checks[1].isChecked())

    def test_cambiar_direccion_no_guarda_sin_accion_explicita(self):
        llamadas = []
        pagina = SimpleNamespace(
            _actualizar_porticos_compatibles_viento=lambda: llamadas.append("compatibles"),
            _marcar_viento_pendiente=lambda: llamadas.append("pendiente"),
            _actualizar_resumen_viento=lambda: llamadas.append("resumen"),
            al_guardar=lambda: llamadas.append("guardar"),
        )

        PaginaCargas._al_cambiar_direccion_viento(pagina, 0)

        self.assertEqual(llamadas, ["compatibles", "pendiente", "resumen"])


class SelectorTerrenoTests(unittest.TestCase):
    def test_lista_obras_anteriores_y_carga_datos_editables(self):
        with tempfile.TemporaryDirectory() as temporal:
            raiz = Path(temporal)
            obra_actual = raiz / "108"
            obra_anterior = raiz / "107"
            (obra_actual / "datos").mkdir(parents=True)
            datos_anteriores = {
                "profundidad_fundacion_m": 1.3,
                "q_adm_kPa": 133.37,
                "fuente_geotecnica": "Informe previo",
                "nivel_freatico_m": 0.0,
                "modulo_balasto_kN_m3": 15985.0,
                "perfil_suelo": "Capa de referencia",
                "recomendaciones_fundacion": [{
                    "tipo_fundacion": "aislada",
                    "profundidad_m": 1.3,
                    "q_adm_kPa": 133.37,
                    "referencia": "Tabla 2",
                }],
            }
            archivo = obra_anterior / "datos" / "terreno.json"
            archivo.parent.mkdir(parents=True)
            archivo.write_text(json.dumps(datos_anteriores), encoding="utf-8")

            selector = QComboBox()
            with (
                patch.object(rutas, "CARPETA_OBRA", obra_actual),
                patch.object(
                    rutas, "listar_obras", return_value=[obra_actual, obra_anterior]
                ),
            ):
                VentanaPrincipal._actualizar_selector_terreno(
                    SimpleNamespace(selector_terreno_anterior=selector)
                )

            self.assertEqual(selector.count(), 2)
            self.assertEqual(selector.itemData(1), str(archivo))

            pagina = SimpleNamespace(
                selector_terreno_anterior=selector,
                spin_profundidad_fundacion=QDoubleSpinBox(),
                txt_q_adm=QLineEdit(),
                txt_fuente_terreno=QLineEdit(),
                txt_napa=QLineEdit(),
                txt_balasto=QLineEdit(),
                txt_perfil_suelo=QPlainTextEdit(),
                recomendaciones_terreno=[],
                tabla_recomendaciones_terreno=QTableWidget(),
            )
            pagina.spin_profundidad_fundacion.setRange(0.1, 20.0)
            pagina._mostrar_recomendaciones_terreno = lambda recomendaciones: (
                setattr(pagina, "recomendaciones_terreno", recomendaciones)
            )

            VentanaPrincipal._aplicar_terreno_anterior(pagina, 1)

            self.assertAlmostEqual(pagina.spin_profundidad_fundacion.value(), 1.3)
            self.assertEqual(pagina.txt_q_adm.text(), "133.37")
            self.assertEqual(pagina.txt_perfil_suelo.toPlainText(), "Capa de referencia")
            self.assertFalse(pagina.txt_q_adm.isReadOnly())
            self.assertFalse(pagina.txt_perfil_suelo.isReadOnly())
            self.assertEqual(len(pagina.recomendaciones_terreno), 1)


if __name__ == "__main__":
    unittest.main()

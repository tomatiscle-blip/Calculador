import os
import unittest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication

from app.ejes import EditorMuroPorNiveles


def _ejes():
    return {
        "x": [
            {"id": "x:1", "nombre": "X1", "coordenada_m": 0.0},
            {"id": "x:2", "nombre": "X2", "coordenada_m": 5.0},
            {"id": "x:3", "nombre": "X3", "coordenada_m": 10.0},
        ],
        "y": [
            {"id": "y:1", "nombre": "Y1", "coordenada_m": 0.0},
            {"id": "y:2", "nombre": "Y2", "coordenada_m": 6.0},
        ],
    }


def _niveles():
    return [
        {"nombre": "PB", "cota_m": 0.0},
        {"nombre": "PA", "cota_m": 3.0},
    ]


def _portico_coplanar():
    """Pórtico que corre en X (direccion='x') ubicado en el eje Y1 (posicion 0)."""
    return {
        "Pórtico 1": {
            "tipo": "portico",
            "referencia_planta": {
                "direccion": "x",
                "eje_id": "y:1",
                "origen_longitudinal_m": 0.0,
            },
            "posicion_planta_m": 0.0,
            "columnas": {
                "C0-1": {"x": 0.0},
                "C0-2": {"x": 5.0},
                "C0-3": {"x": 10.0},
            },
            "vigas": {
                "V0-1": {
                    "tramos": [
                        {"x_inicio": 0.0, "x_fin": 5.0},
                        {"x_inicio": 5.0, "x_fin": 10.0},
                    ],
                },
            },
        },
    }


class EditorMuroMixtoTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def _elegir_direccion(self, editor, familia):
        editor.direccion.setCurrentIndex(editor.direccion.findData(familia))

    def test_detecta_portico_coplanar_sin_borrarlo(self):
        estructura = _portico_coplanar()
        editor = EditorMuroPorNiveles(_ejes(), _niveles(), estructura)
        # «Paralelo a Y»: el muro se aloja en un eje Y (Y1, pos 0) y corre en X,
        # resultando coplanar con el Pórtico 1 (direccion='x', posicion=0).
        self._elegir_direccion(editor, "y")
        editor.eje_transversal.setCurrentIndex(
            editor.eje_transversal.findData("y:1")
        )

        coplanar = editor._portico_coplanar()
        self.assertIsNotNone(coplanar)
        nombre, datos = coplanar
        self.assertEqual(nombre, "Pórtico 1")
        # El pórtico NO se modifica ni se elimina del diccionario de estructura.
        self.assertIn("Pórtico 1", editor.estructura)
        self.assertEqual(editor.estructura["Pórtico 1"], datos)
        # El aviso se muestra para el usuario (isHidden refleja el setVisible,
        # ya que el diálogo no se exhibe con exec() durante el test).
        self.assertFalse(editor.aviso_coplanar.isHidden())
        self.assertFalse(editor.boton_usar_portico.isHidden())
        editor.close()

    def test_usar_longitud_del_portico_ajusta_extremos(self):
        estructura = _portico_coplanar()
        editor = EditorMuroPorNiveles(_ejes(), _niveles(), estructura)
        self._elegir_direccion(editor, "y")
        editor.eje_transversal.setCurrentIndex(
            editor.eje_transversal.findData("y:1")
        )

        editor._usar_longitud_portico()
        datos = editor.datos_muro(mostrar_errores=False)

        # El muro pasa a cubrir exactamente el rango longitudinal del pórtico.
        self.assertAlmostEqual(datos["x_inicio"], 0.0)
        self.assertAlmostEqual(datos["x_fin"], 10.0)
        self.assertAlmostEqual(datos["longitud_m"], 10.0)
        # El pórtico sigue intacto.
        self.assertIn("Pórtico 1", editor.estructura)
        editor.close()

    def test_sin_portico_coplanar_oculta_aviso(self):
        editor = EditorMuroPorNiveles(_ejes(), _niveles(), {})
        self._elegir_direccion(editor, "y")

        self.assertIsNone(editor._portico_coplanar())
        self.assertTrue(editor.aviso_coplanar.isHidden())
        self.assertTrue(editor.boton_usar_portico.isHidden())
        editor.close()


if __name__ == "__main__":
    unittest.main()

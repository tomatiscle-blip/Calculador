import os
import unittest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication

from app.ejes import EditorPorticoPorNiveles


class EditorPorticoPorNivelesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def test_vista_al_desmarcar_primera_columna_en_nivel_superior(self):
        ejes = {
            "x": [
                {"id": "x:a", "nombre": "A", "coordenada_m": 0.0},
                {"id": "x:b", "nombre": "B", "coordenada_m": 5.0},
                {"id": "x:c", "nombre": "C", "coordenada_m": 10.0},
            ],
            "y": [{"id": "y:1", "nombre": "1", "coordenada_m": 0.0}],
        }
        niveles = [
            {"nombre": "PB", "cota_m": 0.0},
            {"nombre": "PA", "cota_m": 3.0},
            {"nombre": "Azotea", "cota_m": 6.0},
        ]
        editor = EditorPorticoPorNiveles(ejes, niveles)
        editor.pestanas.setCurrentIndex(1)

        check = editor.tablas[1].item(0, 0)
        check.setCheckState(Qt.CheckState.Unchecked)

        columnas_en_planta = editor.vista.porticos["Pórtico nuevo"]["columnas"]
        self.assertEqual(columnas_en_planta, [5.0, 10.0])
        geometria = editor.geometria()
        self.assertEqual(
            geometria["referencia_planta"]["eje_longitudinal_id"], "x:a"
        )
        self.assertEqual(
            [
                columna["x_local_m"]
                for columna in geometria["tramos_nivel"][1]["ejes_columna"]
            ],
            [5.0, 10.0],
        )
        editor.close()

    def test_limita_geometria_al_nivel_final_elegido(self):
        ejes = {
            "x": [
                {"id": "x:a", "nombre": "A", "coordenada_m": 0.0},
                {"id": "x:b", "nombre": "B", "coordenada_m": 5.0},
            ],
            "y": [{"id": "y:1", "nombre": "1", "coordenada_m": 0.0}],
        }
        niveles = [
            {"nombre": "PB", "cota_m": 0.0},
            {"nombre": "PA", "cota_m": 3.0},
            {"nombre": "Azotea", "cota_m": 6.0},
        ]
        editor = EditorPorticoPorNiveles(ejes, niveles)
        editor.nivel_final.setCurrentIndex(0)

        geometria = editor.geometria()

        self.assertEqual([nivel["nombre"] for nivel in geometria["niveles"]], ["PB", "PA"])
        self.assertEqual(len(geometria["tramos_nivel"]), 1)
        self.assertTrue(editor.pestanas.isTabVisible(0))
        self.assertFalse(editor.pestanas.isTabVisible(1))
        editor.close()


if __name__ == "__main__":
    unittest.main()

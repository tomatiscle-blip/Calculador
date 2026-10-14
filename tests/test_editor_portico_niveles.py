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

    @staticmethod
    def _elegir_direccion(editor, familia):
        editor.direccion.setCurrentIndex(editor.direccion.findData(familia))

    def test_vista_al_desmarcar_primera_columna_en_nivel_superior(self):
        ejes = {
            "x": [
                {"id": "x:a", "nombre": "A", "coordenada_m": 0.0},
                {"id": "x:b", "nombre": "B", "coordenada_m": 5.0},
                {"id": "x:c", "nombre": "C", "coordenada_m": 10.0},
            ],
            "y": [
                {"id": "y:1", "nombre": "1", "coordenada_m": 0.0},
                {"id": "y:2", "nombre": "2", "coordenada_m": 6.0},
            ],
        }
        niveles = [
            {"nombre": "PB", "cota_m": 0.0},
            {"nombre": "PA", "cota_m": 3.0},
            {"nombre": "Azotea", "cota_m": 6.0},
        ]
        editor = EditorPorticoPorNiveles(ejes, niveles)
        # «Paralelo a Y» ubica el plano sobre un eje Y y corre las columnas en X.
        self._elegir_direccion(editor, "y")
        editor.pestanas.setCurrentIndex(1)

        check = editor.tablas[1].item(0, 0)
        check.setCheckState(Qt.CheckState.Unchecked)

        columnas_en_planta = editor.vista.porticos["Pórtico nuevo"]["columnas"]
        self.assertEqual(columnas_en_planta, [5.0, 10.0])
        geometria = editor.geometria()
        self.assertEqual(
            geometria["referencia_planta"]["eje_longitudinal_id"], "x:a"
        )
        # La referencia guardada usa la familia longitudinal (X), coherente con muros.
        self.assertEqual(geometria["referencia_planta"]["direccion"], "x")
        self.assertEqual(
            [
                columna["x_local_m"]
                for columna in geometria["tramos_nivel"][1]["ejes_columna"]
            ],
            [5.0, 10.0],
        )
        editor.close()

    def test_direccion_x_ubica_plano_en_x_y_columnas_en_y(self):
        ejes = {
            "x": [
                {"id": "x:a", "nombre": "A", "coordenada_m": 0.0},
                {"id": "x:b", "nombre": "B", "coordenada_m": 5.0},
            ],
            "y": [
                {"id": "y:1", "nombre": "1", "coordenada_m": 0.0},
                {"id": "y:2", "nombre": "2", "coordenada_m": 4.0},
            ],
        }
        niveles = [
            {"nombre": "PB", "cota_m": 0.0},
            {"nombre": "PA", "cota_m": 3.0},
        ]
        editor = EditorPorticoPorNiveles(ejes, niveles)
        # «Paralelo a X»: el plano corre sobre un eje X y las columnas sobre ejes Y.
        self._elegir_direccion(editor, "x")
        geometria = editor.geometria()
        self.assertEqual(geometria["referencia_planta"]["direccion"], "y")
        self.assertEqual(
            [col["x_local_m"] for col in geometria["tramos_nivel"][0]["ejes_columna"]],
            [0.0, 4.0],
        )
        editor.close()

    def test_limita_geometria_al_nivel_final_elegido(self):
        ejes = {
            "x": [
                {"id": "x:a", "nombre": "A", "coordenada_m": 0.0},
                {"id": "x:b", "nombre": "B", "coordenada_m": 5.0},
            ],
            "y": [
                {"id": "y:1", "nombre": "1", "coordenada_m": 0.0},
                {"id": "y:2", "nombre": "2", "coordenada_m": 4.0},
            ],
        }
        niveles = [
            {"nombre": "PB", "cota_m": 0.0},
            {"nombre": "PA", "cota_m": 3.0},
            {"nombre": "Azotea", "cota_m": 6.0},
        ]
        editor = EditorPorticoPorNiveles(ejes, niveles)
        self._elegir_direccion(editor, "y")
        editor.nivel_final.setCurrentIndex(0)

        geometria = editor.geometria()

        self.assertEqual([nivel["nombre"] for nivel in geometria["niveles"]], ["PB", "PA"])
        self.assertEqual(len(geometria["tramos_nivel"]), 1)
        self.assertTrue(editor.pestanas.isTabVisible(0))
        self.assertFalse(editor.pestanas.isTabVisible(1))
        editor.close()

    # ------------------------------------------------------------------
    # Ocupación por RANGO longitudinal (aprovechar lugar libre) y reglas
    # verticales muro/pórtico.
    # ------------------------------------------------------------------
    def _ejes_x4(self):
        # Ejes X longitudinales (donde nacen columnas cuando el plano es Y).
        return {
            "x": [
                {"id": "x:1", "nombre": "X1", "coordenada_m": 0.0},
                {"id": "x:2", "nombre": "X2", "coordenada_m": 4.0},
                {"id": "x:3", "nombre": "X3", "coordenada_m": 8.0},
                {"id": "x:4", "nombre": "X4", "coordenada_m": 12.0},
            ],
            "y": [
                {"id": "y:1", "nombre": "Y1", "coordenada_m": 0.0},
                {"id": "y:2", "nombre": "Y2", "coordenada_m": 6.0},
            ],
        }

    @staticmethod
    def _portico_plano_y(nombre_rango, eje_fin_id, origen, hasta):
        return {
            nombre_rango: {
                "tipo": "portico",
                "referencia_planta": {"direccion": "x", "eje_id": "y:1",
                                       "origen_longitudinal_m": origen},
                "posicion_planta_m": 0.0,
                "columnas": {"C0-1": {"x": 0.0}, "C0-2": {"x": hasta}},
                "vigas": {"V0-1": {"tramos": [{"x_inicio": 0.0, "x_fin": hasta}]}},
                "niveles_geometria": [
                    {"nombre": "PB", "cota_m": 0.0, "cota_origen_m": 0.0},
                    {"nombre": "PA", "cota_m": 3.0, "cota_origen_m": 3.0},
                ],
            }
        }

    def test_permite_portico_en_lugar_libre_del_mismo_plano(self):
        # Hay un pórtico en X1 ocupando longitudinalmente 0–4 m (eje X1→X2).
        # Un nuevo pórtico en el mismo plano Y1 pero desde X2 en adelante
        # (4–12 m) NO se solapa: debe permitirse.
        ejes = self._ejes_x4()
        niveles = [
            {"nombre": "PB", "cota_m": 0.0},
            {"nombre": "PA", "cota_m": 3.0},
        ]
        existente = self._portico_plano_y("Pórtico 1", "x:2", 0.0, 4.0)
        editor = EditorPorticoPorNiveles(ejes, niveles, existente)
        self._elegir_direccion(editor, "y")  # plano en Y1, columnas en X
        editor.eje_transversal.setCurrentIndex(editor.eje_transversal.findData("y:1"))
        # Marcar columnas solo en X2, X3, X4 (4–12 m) en TODOS los tramos.
        for tabla in editor.tablas:
            for fila, eje_id in enumerate(["x:1", "x:2", "x:3", "x:4"]):
                estado = Qt.CheckState.Unchecked if eje_id == "x:1" else Qt.CheckState.Checked
                tabla.item(fila, 0).setCheckState(estado)
        # No debe lanzar: hay lugar libre (4–12) que no solapa 0–4.
        geometria = editor.geometria()
        self.assertEqual(geometria["referencia_planta"]["direccion"], "x")
        editor.close()

    def test_rechaza_portico_que_solapa_longitudinal_y_nivel(self):
        # Pórtico existente ocupa 0–4 m en PB–PA. Intentar otro 0–8 m mismo plano
        # y mismos niveles debe fallar por solape.
        ejes = self._ejes_x4()
        niveles = [
            {"nombre": "PB", "cota_m": 0.0},
            {"nombre": "PA", "cota_m": 3.0},
        ]
        existente = self._portico_plano_y("Pórtico 1", "x:2", 0.0, 4.0)
        editor = EditorPorticoPorNiveles(ejes, niveles, existente)
        self._elegir_direccion(editor, "y")
        editor.eje_transversal.setCurrentIndex(editor.eje_transversal.findData("y:1"))
        # Dejar columnas en X1..X3 (0–8 m): solapa con 0–4.
        with self.assertRaises(ValueError):
            editor.geometria()
        editor.close()

    def test_rechaza_portico_sobre_muro_de_planta_baja(self):
        # Un muro arranca en la base (PB) en el plano Y1. No se puede apilar un
        # pórtico encima aunque no solape en niveles.
        ejes = self._ejes_x4()
        niveles = [
            {"nombre": "PB", "cota_m": 0.0},
            {"nombre": "PA", "cota_m": 3.0},
        ]
        muro = {
            "Muro 1": {
                "tipo": "muro",
                "referencia_planta": {"direccion": "x", "eje_id": "y:1"},
                "posicion_planta_m": 0.0,
                "x_inicio": 0.0,
                "x_fin": 4.0,
                "nivel_inicio": 0.0,
                "nivel_fin": 3.0,
            }
        }
        editor = EditorPorticoPorNiveles(ejes, niveles, muro)
        self._elegir_direccion(editor, "y")
        editor.eje_transversal.setCurrentIndex(editor.eje_transversal.findData("y:1"))
        # El nuevo pórtico arranca en PA (nivel_final=PA deja solo PB→PA); para
        # que no solape niveles con el muro (0–3) arrancaría en 3–? Usamos los
        # dos niveles disponibles (0–3): solapa niveles y además hay muro en base.
        with self.assertRaises(ValueError):
            editor.geometria()
        editor.close()


if __name__ == "__main__":
    unittest.main()

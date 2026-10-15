import os
import unittest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication

from calc import muros
from app.muros_proyecto import EditorCapasMuro


class MuroCompuestoTests(unittest.TestCase):
    def _capas_ejemplo(self):
        return [
            {"clave": "cemento_arena", "grupo": "Morteros", "espesor_m": 0.02},
            {"clave": "ladrillo_hueco_portante", "grupo": "Mamposteria", "espesor_m": 0.18},
            {"clave": "cal_arena", "grupo": "Morteros", "espesor_m": 0.02},
        ]

    def test_peso_propio_suma_gamma_por_espesor(self):
        # γ mortero cemento=21, hueco portante=12, mortero cal=17.
        # peso/m2 = 21*0.02 + 12*0.18 + 17*0.02 = 0.42 + 2.16 + 0.34 = 2.92
        resultado = muros.peso_propio_capas(self._capas_ejemplo())
        self.assertAlmostEqual(resultado["espesor_total_m"], 0.22, places=6)
        self.assertAlmostEqual(resultado["peso_propio_kN_m2"], 2.92, places=6)
        self.assertEqual(len(resultado["capas"]), 3)

    def test_calcular_muro_compuesto_geometria(self):
        resultado = muros.calcular_muro_compuesto(
            "Muro 1", self._capas_ejemplo(), altura_m=3.0, longitud_m=7.9,
        )
        self.assertAlmostEqual(resultado["espesor_total_m"], 0.22, places=6)
        self.assertAlmostEqual(resultado["peso_propio_kN_m2"], 2.92, places=6)
        # lineal = peso/m2 * altura = 2.92 * 3 = 8.76 kN/m
        self.assertAlmostEqual(resultado["peso_lineal_kN_m"], 8.76, places=6)
        # total = lineal * longitud = 8.76 * 7.9 = 69.204 kN
        self.assertAlmostEqual(resultado["peso_total_kN"], 69.204, places=3)

    def test_capas_sin_espesor_se_ignoran(self):
        capas = self._capas_ejemplo() + [
            {"clave": "ladrillo_comun", "grupo": "Mamposteria", "espesor_m": 0.0}
        ]
        resultado = muros.peso_propio_capas(capas)
        self.assertEqual(len(resultado["capas"]), 3)
        self.assertAlmostEqual(resultado["espesor_total_m"], 0.22, places=6)

    def test_clave_inexistente_falla_claro(self):
        with self.assertRaisesRegex(ValueError, "No existe la clave"):
            muros.peso_propio_capas(
                [{"clave": "no_existe", "grupo": "Mamposteria", "espesor_m": 0.1}]
            )

    def test_altura_invalida_falla(self):
        with self.assertRaisesRegex(ValueError, "altura"):
            muros.calcular_muro_compuesto("M", self._capas_ejemplo(), 0.0, 5.0)

    def test_pendientes_diseno_y_cimiento_marcadas_como_none(self):
        resultado = muros.calcular_muro_compuesto(
            "Muro 1", self._capas_ejemplo(), altura_m=3.0, longitud_m=5.0,
        )
        # La resistencia térmica AHORA se calcula (todas las capas tienen λ).
        self.assertIsNotNone(resultado["resistencia_termica"])
        # El diseño estructural y el cimiento siguen pendientes.
        self.assertIsNone(resultado["diseno_estructural"])
        self.assertIsNone(resultado["cimiento"])
        self.assertIn("mampostería", resultado["pendiente"])

    def test_resistencia_termica_se_calcula(self):
        # γ: mortero cemento=21, hueco portante=12, mortero cal=17.
        # λ: cemento_arena=1.4, ladrillo_hueco_portante=0.6, cal_arena=0.87.
        # R capas = 0.02/1.4 + 0.18/0.6 + 0.02/0.87
        #          = 0.0142857 + 0.3 + 0.0229885 = 0.3372742 m²·K/W
        resultado = muros.resistencia_termica_capas(self._capas_ejemplo())
        R_capas = 0.02 / 1.4 + 0.18 / 0.6 + 0.02 / 0.87
        self.assertAlmostEqual(resultado["R_capas_m2K_W"], R_capas, places=6)
        # R_total = Rsi(0.13) + R_capas + Rse(0.04)
        self.assertAlmostEqual(resultado["R_total_m2K_W"], 0.13 + R_capas + 0.04, places=6)
        # U = 1/R_total
        self.assertAlmostEqual(
            resultado["U_W_m2K"], 1.0 / (0.13 + R_capas + 0.04), places=6,
        )
        self.assertEqual(len(resultado["capas"]), 3)

    def test_resistencia_termica_en_muro_compuesto(self):
        resultado = muros.calcular_muro_compuesto(
            "Muro 1", self._capas_ejemplo(), altura_m=3.0, longitud_m=5.0,
        )
        termica = resultado["resistencia_termica"]
        self.assertIsNotNone(termica)
        self.assertGreater(termica["R_total_m2K_W"], 0.0)
        self.assertGreater(termica["U_W_m2K"], 0.0)

    def test_camara_de_aire_usa_tabla_no_lambda(self):
        # Tabla ISO 6946 (flujo horizontal): 5 mm → R=0.11 ; 50 mm → R=0.18.
        capas = [{"clave": "camara_aire", "grupo": "Aislantes", "espesor_m": 0.005}]
        resultado = muros.resistencia_termica_capas(capas)
        self.assertAlmostEqual(resultado["R_capas_m2K_W"], 0.11, places=6)
        # Una cámara más ancha (5 cm) ya satura en 0.18 m²·K/W.
        capas = [{"clave": "camara_aire", "grupo": "Aislantes", "espesor_m": 0.05}]
        resultado = muros.resistencia_termica_capas(capas)
        self.assertAlmostEqual(resultado["R_capas_m2K_W"], 0.18, places=6)
        # La cámara no aporta peso propio.
        peso = muros.peso_propio_capas(capas)
        self.assertAlmostEqual(peso["peso_propio_kN_m2"], 0.0, places=6)

    def test_sin_lambda_deja_termica_en_none_pero_calcula_peso(self):
        # Hormigón 'sin_armar' tiene γ (peso) pero no λ cargada: la térmica no
        # se puede calcular, pero el peso propio sí debe seguir funcionando.
        capas = [{"clave": "sin_armar", "grupo": "Hormigon", "espesor_m": 0.1}]
        resultado = muros.calcular_muro_compuesto("M", capas, 3.0, 5.0)
        self.assertIsNone(resultado["resistencia_termica"])
        self.assertGreater(resultado["peso_propio_kN_m2"], 0.0)
        self.assertIn("Resistencia térmica incompleta", resultado["pendiente"])


class EditorCapasMuroTests(unittest.TestCase):
    """Pruebas del diálogo de composición por capas (necesita Qt)."""

    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def _capas(self):
        return [
            {"clave": "cemento_arena", "grupo": "Morteros", "espesor_m": 0.02},
            {"clave": "ladrillo_hueco_portante", "grupo": "Mamposteria", "espesor_m": 0.18},
            {"clave": "cal_arena", "grupo": "Morteros", "espesor_m": 0.02},
        ]

    def test_cargar_capas_selecciona_el_material_correcto(self):
        # Regresión: findData con tupla no matcheaba y todos los combos
        # quedaban en el primer material (peso propio mal calculado).
        d = EditorCapasMuro("Muro 1", self._capas(), altura_m=3.0, longitud_m=7.9)
        claves = [d.tabla.cellWidget(f, 0).currentData()[1]
                  for f in range(d.tabla.rowCount())]
        self.assertEqual(claves, ["cemento_arena", "ladrillo_hueco_portante", "cal_arena"])
        d.close()

    def test_peso_propio_coincide_con_calculo(self):
        # Con los materiales bien seleccionados, el peso debe ser 2,92 kN/m².
        d = EditorCapasMuro("Muro 1", self._capas(), altura_m=3.0, longitud_m=7.9)
        _capas_out, _alt, _lon, res = d.datos_guardar()
        self.assertAlmostEqual(res["espesor_total_m"], 0.22, places=6)
        self.assertAlmostEqual(res["peso_propio_kN_m2"], 2.92, places=6)
        d.close()

    def test_espesor_en_cm_se_guarda_en_metros(self):
        # La UI carga/muestra cm; el cálculo y el guardado usan metros.
        d = EditorCapasMuro("Muro 1", self._capas(), altura_m=3.0, longitud_m=7.9)
        spin = d.tabla.cellWidget(1, 1)  # 0,18 m -> 18 cm en pantalla
        self.assertAlmostEqual(spin.value(), 18.0, places=3)
        capas_out, _alt, _lon, _res = d.datos_guardar()
        self.assertAlmostEqual(capas_out[1]["espesor_m"], 0.18, places=6)
        d.close()

    def test_altura_de_niveles_es_solo_lectura(self):
        # Cuando la altura sale de los niveles, el campo no se edita a mano.
        d = EditorCapasMuro("Muro 1", self._capas(), altura_m=3.0, longitud_m=7.9,
                            altura_de_niveles=True)
        self.assertTrue(d.altura.isReadOnly())
        self.assertAlmostEqual(d.altura.value(), 3.0, places=3)
        d.close()

    def test_recalcular_no_recursa(self):
        # Regresión del crash: _recalcular escribía setItem -> itemChanged ->
        # _recalcular en bucle hasta desbordar la pila y cerrar la app.
        d = EditorCapasMuro("Muro 1", [], altura_m=3.0, longitud_m=7.9)
        d._agregar_capa()
        d._agregar_capa()
        _capas_out, _alt, _lon, res = d.datos_guardar()
        self.assertIn("peso_propio_kN_m2", res)
        d.close()


if __name__ == "__main__":
    unittest.main()

import os
import unittest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from calc import muros


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

    def test_pendientes_marcadas_como_none(self):
        resultado = muros.calcular_muro_compuesto(
            "Muro 1", self._capas_ejemplo(), altura_m=3.0, longitud_m=5.0,
        )
        self.assertIsNone(resultado["resistencia_termica"])
        self.assertIsNone(resultado["diseno_estructural"])
        self.assertIsNone(resultado["cimiento"])
        self.assertIn("conductividades", resultado["pendiente"])


if __name__ == "__main__":
    unittest.main()

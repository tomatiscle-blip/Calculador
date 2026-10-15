import unittest

from calc import mamposteria as mamp


def _muro(e_hoja=0.20):
    return [
        {"clave": "cemento_arena", "grupo": "Morteros", "espesor_m": 0.02},
        {"clave": "ladrillo_hueco_portante", "grupo": "Mamposteria", "espesor_m": e_hoja},
        {"clave": "cal_arena", "grupo": "Morteros", "espesor_m": 0.02},
    ]


class MamposteriaTests(unittest.TestCase):
    def test_fm_ladrillo_hueco_portante_mortero_normal(self):
        # Tabla 6.3: bloque hueco portante cerámico, mortero N = 4.55 MPa.
        self.assertAlmostEqual(mamp.fm_de("ladrillo_hueco_portante", "mortero_normal_n"), 4.55)

    def test_fm_material_no_mapeado_falla_claro(self):
        with self.assertRaisesRegex(ValueError, "no está mapeado"):
            mamp.fm_de("porcelanato")

    def test_esbeltez_placa_0_20m_3m_alto(self):
        # r = t/sqrt(12); h/r = 3000/(200/3.4646) = 51.96
        hr = mamp.esbeltez(3000.0, 200.0)
        self.assertAlmostEqual(hr, 51.9615, places=3)
        self.assertFalse(hr > 99)

    def test_Pn_compresion_valores_esperados(self):
        # Franja 1 m: An = 200mm * 1000mm = 200000 mm2.
        # fm=4.55, As=0 -> Pn = 0.80*[0.80*4.55*200000]*(1-(52/140)^2)
        res = mamp.resistencia_axial_Pn(4.55, 200000.0, 420.0, 0.0, 3000.0, 200.0)
        esperado = 0.80 * (0.80 * 4.55 * 200000.0) * (1 - (res["esbeltez_hr"] / 140.0) ** 2)
        self.assertAlmostEqual(res["Pn_N"], esperado, places=3)
        self.assertEqual(res["rama"], "h/r <= 99")

    def test_phi_simple_vs_armada(self):
        self.assertAlmostEqual(mamp.phi(armada=False), 0.60)
        self.assertAlmostEqual(mamp.phi(armada=True), 0.90)

    def test_verificar_muro_completo(self):
        res = mamp.verificar_muro(_muro(), altura_m=3.0, armada=False)
        self.assertTrue(res["ok"])
        self.assertEqual(res["hoja_portante"]["clave"], "ladrillo_hueco_portante")
        self.assertAlmostEqual(res["hoja_portante"]["fm_mpa"], 4.55)
        self.assertFalse(res["es_esbelta"])
        self.assertGreater(res["compresion"]["phi_Pn_kN"], 0.0)
        self.assertGreater(res["corte"]["phi_Vn_kN"], 0.0)
        self.assertGreaterEqual(len(res["recomendaciones"]), 3)

    def test_verificar_muro_pasa_con_carga_moderada(self):
        res = mamp.verificar_muro(_muro(), altura_m=3.0, carga_axial_kN_m=50.0)
        self.assertTrue(res["compresion"]["pasa"])

    def test_verificar_muro_falla_sin_mamposteria(self):
        # Solo revoque/aislante: no hay hoja portante.
        capas = [{"clave": "cemento_arena", "grupo": "Morteros", "espesor_m": 0.02}]
        res = mamp.verificar_muro(capas, altura_m=3.0)
        self.assertFalse(res["ok"])
        self.assertIn("Mampostería", res["pendiente"])


if __name__ == "__main__":
    unittest.main()

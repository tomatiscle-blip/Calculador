import unittest

from calc.portico import _huella_entradas, _peso_propio_viga_kN_m


class PesoPropioVigaTests(unittest.TestCase):
    def test_predimensiona_peralte_como_dimensionador_y_calcula_carga(self):
        peso, peralte, predimensionado = _peso_propio_viga_kN_m({
            "b_cm": 20,
            "tramos": [{"longitud_m": 7.0}],
        })

        self.assertEqual(peralte, 56)
        self.assertTrue(predimensionado)
        self.assertAlmostEqual(peso, 2.8)

    def test_prefiere_peralte_explicito(self):
        peso, peralte, predimensionado = _peso_propio_viga_kN_m({
            "b_cm": 20,
            "h_cm": 42,
            "tramos": [{"longitud_m": 7.0}],
        })

        self.assertEqual(peralte, 42)
        self.assertFalse(predimensionado)
        self.assertAlmostEqual(peso, 2.1)

    def test_rechaza_seccion_incompleta(self):
        with self.assertRaisesRegex(ValueError, "ancho b_cm"):
            _peso_propio_viga_kN_m({"tramos": [{"longitud_m": 7.0}]})

    def test_cambio_de_seccion_invalida_huella_del_motor(self):
        datos = {"vigas": {"V1-3": {"b_cm": 20, "h_cm": 56}}}
        huella_original = _huella_entradas(datos, {}, {})
        datos["vigas"]["V1-3"]["h_cm"] = 60

        self.assertNotEqual(huella_original, _huella_entradas(datos, {}, {}))


if __name__ == "__main__":
    unittest.main()

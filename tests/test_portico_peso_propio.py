import unittest

from calc.portico import (
    _aviso_viento_no_aplicado,
    _huella_entradas,
    _peso_propio_viga_kN_m,
)


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

    def test_viento_en_otra_direccion_no_es_aviso_del_motor(self):
        aviso = _aviso_viento_no_aplicado({
            "aplicado": False,
            "motivo": "el plano del pórtico no coincide con el eje global del viento",
        })

        self.assertIsNone(aviso)

    def test_portico_desmarcado_no_es_aviso_del_motor(self):
        aviso = _aviso_viento_no_aplicado({
            "aplicado": False,
            "motivo": "pórtico desmarcado",
        })

        self.assertIsNone(aviso)

    def test_viento_con_configuracion_incompleta_sigue_generando_aviso(self):
        aviso = _aviso_viento_no_aplicado({
            "aplicado": False,
            "motivo": "falta ancho tributario",
        })

        self.assertEqual(aviso, "Viento horizontal omitido: falta ancho tributario.")


if __name__ == "__main__":
    unittest.main()

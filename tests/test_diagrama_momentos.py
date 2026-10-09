import unittest

from calc.portico import _diagrama_momento


class Nodo:
    def __init__(self, x, y):
        self.X = x
        self.Y = y


class Miembro:
    def __init__(self, corte_estaciones=None):
        self.i_node = Nodo(0.0, 0.0)
        self.j_node = Nodo(2.0, 0.0)
        self.corte_estaciones = corte_estaciones or [0.0, 1.0, 2.0]

    def moment_array(self, direccion, puntos, combinacion):
        return [0.0, 1.0, 2.0], [0.0, 4.0, 0.0]

    def shear_array(self, direccion, puntos, combinacion):
        return self.corte_estaciones, [-5.0, 2.0, 0.0]


class Modelo:
    members = {"V1": Miembro()}


class DiagramaMomentosTests(unittest.TestCase):
    def test_incluye_momento_corte_y_apoyos(self):
        apoyos = [{"x": 0.0, "y": 0.0, "tipo": "empotramiento"}]
        diagrama = _diagrama_momento(
            Modelo(),
            {
                "columnas": {"C1": "V1"},
                "vigas": {},
                "voladizos": {},
                "apoyos": apoyos,
            },
            "D",
        )

        self.assertEqual(diagrama["max_abs_kNm"], 4.0)
        self.assertEqual(diagrama["max_abs_kN"], 5.0)
        self.assertEqual(diagrama["apoyos"], apoyos)
        self.assertEqual(diagrama["barras"][0]["M_kNm"], [0.0, 4.0, 0.0])
        self.assertEqual(diagrama["barras"][0]["V_kN"], [-5.0, 2.0, 0.0])

    def test_rechaza_series_de_corte_con_otro_numero_de_estaciones(self):
        modelo = Modelo()
        modelo.members = {"V1": Miembro([0.0, 2.0])}

        with self.assertRaisesRegex(ValueError, "diagramas inválidos"):
            _diagrama_momento(
                modelo,
                {
                    "columnas": {"C1": "V1"},
                    "vigas": {},
                    "voladizos": {},
                    "apoyos": [],
                },
                "D",
            )


if __name__ == "__main__":
    unittest.main()

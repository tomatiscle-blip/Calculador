"""
Tests del próximo nombre libre de pórtico (calc/rutas.py).

Fijan la convención del proyecto: los pórticos se llaman "Pórtico N"
(con tilde), como los crea la app y como están en los datos.
"""

import sys
import unittest
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent
if str(RAIZ) not in sys.path:
    sys.path.insert(0, str(RAIZ))

from calc import rutas  # noqa: E402


class ProximoNombrePorticoTests(unittest.TestCase):
    def test_usa_la_convencion_con_tilde(self):
        self.assertTrue(
            rutas.proximo_nombre_portico().startswith("Pórtico "),
            "los pórticos se llaman 'Pórtico N' (con tilde)",
        )

    def test_no_pisa_un_nombre_existente(self):
        self.assertNotIn(
            rutas.proximo_nombre_portico(), rutas.listar_porticos()
        )


if __name__ == "__main__":
    unittest.main()

import unittest
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from calc import rutas

class TestColumnasGlobalesFases(unittest.TestCase):
    def test_fase1_archivos_existen(self):
        SAL = rutas.SAL_COLUMNAS
        archivos = list(SAL.glob("dimensionado_*.json"))
        self.assertGreater(len(archivos), 0, "Falta ejecutar Fase 1")

    def test_fase1_estructura(self):
        SAL = rutas.SAL_COLUMNAS
        archivos = list(SAL.glob("dimensionado_*.json"))
        self.assertGreater(len(archivos), 0)
        d = rutas.leer_json(archivos[0], {})
        self.assertIn("columna_global", d)
        self.assertIn("dimensionados", d)
        self.assertIsInstance(d["dimensionados"], dict)

    def test_fase2_archivos_existen(self):
        SAL = rutas.SAL_COLUMNAS
        archivos = list(SAL.glob("unificado_*.json"))
        self.assertGreater(len(archivos), 0, "Falta ejecutar Fase 2")

    def test_fase4_archivos_existen(self):
        SAL = rutas.SAL_COLUMNAS
        archivos = list(SAL.glob("resultado_final_*.json"))
        self.assertGreater(len(archivos), 0, "Falta ejecutar Fase 4")

    def test_planilla_final_existe(self):
        SAL = rutas.SAL_COLUMNAS
        self.assertTrue((SAL / "planilla_columnas_globales_final.csv").exists())

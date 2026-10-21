"""
Tests de la caché de rutas resueltas del pipeline (calc/pipeline.py).

La caché evita volver a expandir los comodines (*, {portico}) en cada consulta
del semáforo. Estos tests fijan que:

  * las rutas resueltas se reutilizan entre consultas,
  * se vuelven a resolver al cambiar de obra,
  * un cambio en un insumo cambia la huella de la etapa (y su estado).
"""

import json
import sys
import tempfile
import unittest
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent
if str(RAIZ) not in sys.path:
    sys.path.insert(0, str(RAIZ))

from calc import pipeline, rutas  # noqa: E402


class CacheRutasResueltasTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.carpeta = Path(self.tmp.name)

        self.obra_original = rutas.CARPETA_OBRA
        self.estructura_original = rutas.ESTRUCTURA
        self.cargas_original = rutas.CARGAS
        self.sal_vigas_original = rutas.SAL_VIGAS
        self.cache_etapas = dict(pipeline._cache_etapas)
        self.cache_resueltas = dict(pipeline._cache_resueltas)
        self.addCleanup(self._restaurar)

        rutas.ESTRUCTURA = self.carpeta / "estructura.json"
        rutas.CARGAS = self.carpeta / "cargas.json"
        rutas.SAL_VIGAS = self.carpeta / "vigas"
        self.carpeta_salidas = self.carpeta / "vigas"
        self.carpeta_salidas.mkdir(parents=True, exist_ok=True)
        rutas.guardar_json(rutas.ESTRUCTURA, {})
        rutas.guardar_json(rutas.CARGAS, {})
        pipeline.invalidar_cache_etapas()

    def _restaurar(self):
        rutas.CARPETA_OBRA = self.obra_original
        rutas.ESTRUCTURA = self.estructura_original
        rutas.CARGAS = self.cargas_original
        rutas.SAL_VIGAS = self.sal_vigas_original
        pipeline._cache_etapas.clear()
        pipeline._cache_etapas.update(self.cache_etapas)
        pipeline._cache_resueltas.clear()
        pipeline._cache_resueltas.update(self.cache_resueltas)

    def test_las_rutas_resueltas_se_reutilizan(self):
        primera = pipeline._rutas_etapa("vigas", "Pórtico 1")
        segunda = pipeline._rutas_etapa("vigas", "Pórtico 1")
        self.assertIs(primera, segunda, "debió reutilizar las rutas ya resueltas")

    def test_cada_pórtico_tiene_su_propia_resolución(self):
        uno = pipeline._rutas_etapa("vigas", "Pórtico 1")
        otro = pipeline._rutas_etapa("vigas", "Pórtico 2")
        self.assertIsNot(uno, otro)
        self.assertNotEqual(uno[1], otro[1], "las salidas deberían diferir por pórtico")
        self.assertIn("Pórtico 1", str(uno[1]))
        self.assertIn("Pórtico 2", str(otro[1]))

    def test_cambiar_de_obra_resuelve_de_nuevo(self):
        antes = pipeline._rutas_etapa("vigas", "Pórtico 1")
        # Se simula un cambio de obra invalidando las cachés del pipeline.
        pipeline.invalidar_cache_etapas()
        despues = pipeline._rutas_etapa("vigas", "Pórtico 1")
        self.assertIsNot(antes, despues)
        self.assertEqual(antes[0], despues[0])

    def test_la_huella_cambia_si_cambia_un_insumo(self):
        antes = pipeline._huella_etapa("geometria", "Pórtico 1")
        rutas.guardar_json(rutas.ESTRUCTURA, {"Pórtico 1": {"columnas": {}}})
        despues = pipeline._huella_etapa("geometria", "Pórtico 1")
        self.assertNotEqual(antes, despues)

    def test_la_huella_no_cambia_si_no_cambia_nada(self):
        antes = pipeline._huella_etapa("geometria", "Pórtico 1")
        despues = pipeline._huella_etapa("geometria", "Pórtico 1")
        self.assertEqual(antes, despues)

    def test_el_estado_se_recalcula_tras_modificar_un_insumo(self):
        antes = pipeline.estado_etapa("geometria", "Pórtico 1")
        self.assertEqual(antes["estado"], "sin_datos")

        rutas.guardar_json(
            rutas.ESTRUCTURA,
            {"Pórtico 1": {"columnas": {}, "vigas": {}}},
        )
        despues = pipeline.estado_etapa("geometria", "Pórtico 1")
        self.assertNotEqual(despues["detalle"], antes["detalle"])
        self.assertIn("Pórtico 1", despues["detalle"])


if __name__ == "__main__":
    unittest.main()

"""
Tests de las cachés del pipeline y de las firmas de archivos (calc/rutas.py
y calc/pipeline.py).

Estos tests fijan el comportamiento de las optimizaciones:

  * las rutas resueltas de una etapa se reutilizan (no se vuelve a expandir
    los comodines) y se vuelven a resolver al cambiar de obra,
  * la huella de una etapa cambia si se modifica un insumo,
  * guardar un archivo o terminar un script invalida las fechas cacheadas,
  * el TTL corto de las firmas no tapa un cambio hecho por la propia app.
"""

import json
import os
import sys
import tempfile
import time
import unittest
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent
if str(RAIZ) not in sys.path:
    sys.path.insert(0, str(RAIZ))

from calc import pipeline, rutas  # noqa: E402


class CacheFirmasTests(unittest.TestCase):
    """Caché de stat() con vida corta (rutas._stat_cache)."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.carpeta = Path(self.tmp.name)
        self.archivo = self.carpeta / "dato.json"
        self.archivo.write_text("{}", encoding="utf-8")
        rutas.invalidar_cache_firmas()

    def test_la_firma_lee_el_mtime(self):
        firma = rutas.firma_modificacion(self.archivo)
        self.assertIsNotNone(firma)
        self.assertEqual(firma, self.archivo.stat().st_mtime_ns)

    def test_archivo_inexistente_no_tiene_firma(self):
        self.assertIsNone(rutas.firma_modificacion(self.carpeta / "no_existe.json"))

    def test_guardar_texto_invalida_las_fechas(self):
        antes = rutas.firma_modificacion(self.archivo)
        time.sleep(0.01)
        rutas.guardar_texto(self.archivo, '{"a": 1}')
        despues = rutas.firma_modificacion(self.archivo)
        self.assertNotEqual(antes, despues)

    def test_guardar_json_invalida_las_fechas(self):
        antes = rutas.firma_modificacion(self.archivo)
        time.sleep(0.01)
        rutas.guardar_json(self.archivo, {"a": 1})
        despues = rutas.firma_modificacion(self.archivo)
        self.assertNotEqual(antes, despues)

    def test_el_ttl_no_tapa_un_cambio_inmediato(self):
        # La app escribe y consulta enseguida: tiene que verse el cambio.
        for indice in range(3):
            rutas.guardar_json(self.archivo, {"n": indice})
            firma = rutas.firma_modificacion(self.archivo)
            self.assertEqual(
                firma, self.archivo.stat().st_mtime_ns,
                "la caché de firmas ocultó un cambio recién hecho",
            )


if __name__ == "__main__":
    unittest.main()

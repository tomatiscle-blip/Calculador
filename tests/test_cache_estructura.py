"""
Tests de la caché de datos/estructura.json (calc/rutas.py).

La estructura se relee del disco SOLO cuando el archivo cambia (firma
mtime_ns + tamaño). Estos tests fijan ese comportamiento:

  * la caché evita releer el archivo si no cambió,
  * si el archivo cambia (se escribe otra vez), se vuelve a leer,
  * lo que devuelve es una COPIA: modificar lo devuelto no daña la caché,
  * guardar_estructura() refresca la caché con lo recién escrito,
  * cambiar de archivo (obra) invalida la caché.
"""

import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent
if str(RAIZ) not in sys.path:
    sys.path.insert(0, str(RAIZ))

from calc import rutas  # noqa: E402


class CacheEstructuraTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.carpeta = Path(self.tmp.name)
        self.estructura_original = rutas.ESTRUCTURA
        self.cache_original = rutas._cache_estructura
        self.addCleanup(self._restaurar)
        rutas.ESTRUCTURA = self.carpeta / "estructura.json"
        rutas._cache_estructura = None

    def _restaurar(self):
        rutas.ESTRUCTURA = self.estructura_original
        rutas._cache_estructura = self.cache_original

    def _escribir(self, datos):
        (self.carpeta / "estructura.json").write_text(
            json.dumps(datos, ensure_ascii=False), encoding="utf-8"
        )

    def test_sin_archivo_devuelve_vacio(self):
        self.assertEqual(rutas.cargar_estructura(), {})

    def test_no_relee_si_el_archivo_no_cambio(self):
        self._escribir({"Pórtico 1": {"columnas": {}}})
        primero = rutas.cargar_estructura()
        self.assertIn("Pórtico 1", primero)

        lecturas = {"n": 0}
        original = rutas.leer_json

        def _contador(ruta, por_defecto=None):
            if str(ruta) == str(rutas.ESTRUCTURA):
                lecturas["n"] += 1
            return original(ruta, por_defecto)

        rutas.leer_json = _contador
        try:
            segundo = rutas.cargar_estructura()
        finally:
            rutas.leer_json = original

        self.assertEqual(lecturas["n"], 0, "la caché debió evitar releer el archivo")
        self.assertEqual(segundo, primero)

    def test_relee_si_el_archivo_cambio(self):
        self._escribir({"Pórtico 1": {"columnas": {}}})
        self.assertIn("Pórtico 1", rutas.cargar_estructura())

        # Se escribe OTRO contenido (mismo tamaño no importa: cambia el mtime).
        self._escribir({"Pórtico 2": {"columnas": {}, "vigas": {}}})
        datos = rutas.cargar_estructura()
        self.assertIn("Pórtico 2", datos)
        self.assertNotIn("Pórtico 1", datos)

    def test_devuelve_copia_y_no_manda_la_cache(self):
        self._escribir({"Pórtico 1": {"columnas": {"C1": {"x": 0.0}}}})
        copia = rutas.cargar_estructura()
        copia["Pórtico 1"]["columnas"]["C1"]["x"] = 99.0
        copia["Pórtico 1"]["nuevo"] = True

        otra = rutas.cargar_estructura()
        self.assertEqual(otra["Pórtico 1"]["columnas"]["C1"]["x"], 0.0)
        self.assertNotIn("nuevo", otra["Pórtico 1"])

    def test_guardar_estructura_refresca_la_cache(self):
        datos = {"Pórtico 3": {"columnas": {}}}
        rutas.guardar_estructura(datos)
        # Sin tocar el archivo: la caché ya tiene lo recién guardado.
        self.assertEqual(rutas.cargar_estructura(), datos)
        # Y el archivo existe con ese contenido.
        leido = json.loads((self.carpeta / "estructura.json").read_text(encoding="utf-8"))
        self.assertEqual(leido, datos)

    def test_cambiar_de_archivo_invalida_la_cache(self):
        self._escribir({"Pórtico 1": {}})
        self.assertIn("Pórtico 1", rutas.cargar_estructura())

        otro = self.carpeta / "otra_obra.json"
        otro.write_text(json.dumps({"Pórtico 9": {}}), encoding="utf-8")
        rutas.ESTRUCTURA = otro
        datos = rutas.cargar_estructura()
        self.assertIn("Pórtico 9", datos)
        self.assertNotIn("Pórtico 1", datos)


"""
Tests de la caché del semáforo (calc/pipeline.py).

El estado de una etapa se calcula mirando archivos y leyendo JSON. En un
refresco de la app el semáforo se consulta varias veces, así que el estado
se cachea. Estos tests fijan ese comportamiento:

  * dos llamadas seguidas dan el mismo resultado,
  * la segunda llamada NO vuelve a leer los archivos,
  * guardar un JSON invalida la caché (no quedan estados viejos),
  * invalidar_cache_etapas() fuerza el recálculo,
  * lo que se devuelve es una COPIA: modificarla no daña la caché.
"""

import sys
import tempfile
import unittest
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent
if str(RAIZ) not in sys.path:
    sys.path.insert(0, str(RAIZ))

from calc import pipeline, rutas  # noqa: E402


class CacheSemaforoTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.carpeta = Path(self.tmp.name)

        self.estructura_original = rutas.ESTRUCTURA
        self.cargas_original = rutas.CARGAS
        self.cache_etapas_original = dict(pipeline._cache_etapas)
        self.addCleanup(self._restaurar)

        rutas.ESTRUCTURA = self.carpeta / "estructura.json"
        rutas.CARGAS = self.carpeta / "cargas.json"
        rutas.guardar_json(rutas.ESTRUCTURA, {})
        rutas.guardar_json(rutas.CARGAS, {})
        pipeline.invalidar_cache_etapas()

    def _restaurar(self):
        rutas.ESTRUCTURA = self.estructura_original
        rutas.CARGAS = self.cargas_original
        pipeline._cache_etapas.clear()
        pipeline._cache_etapas.update(self.cache_etapas_original)

    def test_dos_llamadas_dan_el_mismo_resultado(self):
        primero = pipeline.estado_etapa("geometria", "Pórtico 1")
        segundo = pipeline.estado_etapa("geometria", "Pórtico 1")
        self.assertEqual(primero, segundo)

    def test_la_segunda_llamada_no_relee_los_archivos(self):
        pipeline.estado_etapa("geometria", "Pórtico 1")

        lecturas = {"n": 0}
        original = rutas.leer_json

        def _contador(ruta, por_defecto=None):
            lecturas["n"] += 1
            return original(ruta, por_defecto)

        rutas.leer_json = _contador
        try:
            pipeline.estado_etapa("geometria", "Pórtico 1")
        finally:
            rutas.leer_json = original

        self.assertEqual(lecturas["n"], 0, "la caché debió evitar releer los archivos")

    def test_guardar_json_invalida_la_cache(self):
        antes = pipeline.estado_etapa("geometria", "Pórtico 1")
        self.assertEqual(antes["estado"], "sin_datos")

        rutas.guardar_json(rutas.ESTRUCTURA, {"Pórtico 1": {"columnas": {}}})
        despues = pipeline.estado_etapa("geometria", "Pórtico 1")

        self.assertNotEqual(
            despues["detalle"], antes["detalle"],
            "el estado debió recalcularse tras guardar",
        )
        self.assertIn("Pórtico 1", despues["detalle"])

    def test_invalidar_cache_etapas_fuerza_el_recalculo(self):
        pipeline.estado_etapa("geometria", "Pórtico 1")
        pipeline.invalidar_cache_etapas()
        self.assertEqual(pipeline._cache_etapas, {})

    def test_devuelve_copia_y_no_manda_la_cache(self):
        primero = pipeline.estado_etapa("geometria", "Pórtico 1")
        primero["estado"] = "inventado"
        primero["detalle"] = "inventado"

        segundo = pipeline.estado_etapa("geometria", "Pórtico 1")
        self.assertNotEqual(segundo["estado"], "inventado")
        self.assertNotEqual(segundo["detalle"], "inventado")

    def test_semaforo_devuelve_lista_coherente(self):
        estados = pipeline.semaforo("Pórtico 1")
        self.assertEqual(len(estados), len(pipeline.ETAPAS))
        for dato in estados:
            self.assertIn("estado", dato)
            self.assertIn("detalle", dato)


if __name__ == "__main__":
    unittest.main()

import os
import unittest
from unittest.mock import patch

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication

from calc import cargas, losas, losas_macizas, muros, portico
from calc.losas_pynite import analizar_paño_macizo
from calc.losas_seccion import espesor_equivalente_hormigon_casetonada
from app.cargas_proyecto import EditorElemento
from app.losas_ejes import EditorLosaDesdeEjes

APP = QApplication.instance() or QApplication([])


BIBLIOTECA = {
    "Forjados": [
        {
            "clave": "Losa_maciza",
            "nombre": "Losa maciza",
            "tipo": "volumetrico",
            "valor": 25.0,
        },
    ],
    "Hormigon": [
        {
            "clave": "armado",
            "nombre": "Hormigón armado",
            "tipo": "volumetrico",
            "valor": 25.0,
        },
    ],
    "Sobrecargas": [
        {"clave": "vivienda", "uso": "Vivienda", "valor_kNm2": 2.0},
    ],
}

SECCION_CASETONADA = {
    "altura_total_m": 0.30,
    "espesor_capa_m": 0.05,
    "ancho_nervio_x_m": 0.10,
    "ancho_nervio_y_m": 0.10,
    "separacion_libre_x_m": 0.60,
    "separacion_libre_y_m": 0.60,
    "kx_mod": 0.50,
    "ky_mod": 0.35,
    "aligeramiento_estructural": True,
}

APOYOS = {
    lado: {"tipo": "portico", "portico": f"Pórtico {lado}", "viga": "V1-1"}
    for lado in ("izq", "der", "inf", "sup")
}


class ModelosDeLosasTests(unittest.TestCase):
    def elemento(self, tipologia):
        elemento = {
            "id": "L1",
            "tipo": "losa",
            "tipologia": tipologia,
            "luz_transversal_m": 4.0,
            "ancho_losa_m": 4.0,
            "ancho_tributario_m": 1.5,
            "sobrecarga": "vivienda",
            "componentes": [],
            "apoya_en": (
                {lado: apoyo for lado, apoyo in APOYOS.items() if lado in ("izq", "der")}
                if tipologia == "maciza"
                else dict(APOYOS)
            ),
        }
        if tipologia == "casetonada":
            elemento["casetonada"] = SECCION_CASETONADA
        else:
            elemento["componentes"] = [{
                "grupo": "Forjados",
                "clave": "Losa_maciza",
                "espesor_m": 0.15,
            }]
        return elemento

    def test_geometria_casetonada_y_limites_cirsoc(self):
        espesor = espesor_equivalente_hormigon_casetonada(SECCION_CASETONADA)
        self.assertAlmostEqual(espesor, 0.1163265306)
        validada = losas_macizas.validar_seccion_casetonada(SECCION_CASETONADA)
        self.assertEqual(validada["referencia_cirsoc"], (
            "Proyecto CIRSOC 201-24, arts. 8.8.1.1–8.8.1.4 y 8.8.2.1.1"
        ))

        geometria_invalida = {**SECCION_CASETONADA, "separacion_libre_x_m": 0.76}
        with self.assertRaisesRegex(ValueError, "750 mm"):
            losas_macizas.validar_seccion_casetonada(geometria_invalida)

    def test_peso_propio_casetonada_se_incluye_en_cargas_globales(self):
        elemento = self.elemento("casetonada")
        elemento["componentes"] = [
            {"grupo": "Forjados", "clave": "Losa_maciza", "espesor_m": 0.15},
            {"grupo": "Pisos", "clave": "porcelanato"},
        ]
        biblioteca = {
            **BIBLIOTECA,
            "Pisos": [{
                "clave": "porcelanato",
                "nombre": "Porcelanato",
                "tipo": "superficial",
                "valor": 0.4,
            }],
        }

        valores = cargas.valores_superficiales(elemento, biblioteca)
        peso_hormigon = 25.0 * espesor_equivalente_hormigon_casetonada(
            SECCION_CASETONADA
        )
        self.assertAlmostEqual(valores["D_kNm2"], peso_hormigon + 0.4)
        items = cargas.items_de_elemento("Losa casetonada", elemento, biblioteca)
        self.assertAlmostEqual(items[0]["valor"], valores["D_kNm2"] * 1.5)

    def test_losas_macizas_y_casetonadas_equilibran_reacciones(self):
        with patch(
            "calc.losas_pynite.configuracion_obra.modulo_elasticidad_hormigon",
            return_value={"Ec_MPa": 25000.0, "articulo": "prueba"},
        ):
            for tipologia in ("maciza", "maciza_cruzada", "casetonada"):
                with self.subTest(tipologia=tipologia):
                    resultado = losas_macizas.calcular_desde_elemento_carga(
                        "Losa test", self.elemento(tipologia), BIBLIOTECA
                    )
                    lados = {reaccion["apoyo"] for reaccion in resultado.reacciones}
                    esperados = (
                        {"izq", "der"}
                        if tipologia == "maciza"
                        else {"izq", "der", "inf", "sup"}
                    )
                    self.assertEqual(lados, esperados)
                    reaccion_total = sum(
                        reaccion["reaccion_total_D_kN"]
                        for reaccion in resultado.reacciones
                    )
                    self.assertAlmostEqual(
                        reaccion_total,
                        resultado.D_kNm2 * resultado.luz_m * resultado.ancho_m,
                        places=6,
                    )

    def test_predimensionado_unidireccional_se_guarda_y_llega_al_calculo(self):
        elemento = {
            "id": "L-auto",
            "tipo": "losa",
            "tipologia": "maciza",
            "luz_transversal_m": 4.0,
            "ancho_losa_m": 4.0,
            "sobrecarga": "vivienda",
            "componentes": [],
            "apoya_en": {
                "izq": APOYOS["izq"],
                "der": APOYOS["der"],
            },
        }
        editor = EditorElemento("Losa auto", elemento, BIBLIOTECA)
        nombre, guardado = editor.resultado()
        componente = next(
            c for c in guardado["componentes"] if c.get("clave") == "Losa_maciza"
        )

        self.assertAlmostEqual(componente["espesor_m"], 0.20)
        self.assertTrue(guardado["espesor_automatico"])
        with patch(
            "calc.losas_pynite.configuracion_obra.modulo_elasticidad_hormigon",
            return_value={"Ec_MPa": 25000.0, "articulo": "prueba"},
        ):
            resultado = losas_macizas.calcular_desde_elemento_carga(
                nombre, guardado, BIBLIOTECA
            )
        self.assertAlmostEqual(resultado.espesor_m, 0.20)

    def test_losa_cruzada_requiere_y_usa_espesor_manual(self):
        elemento = {
            "id": "L-cruzada",
            "tipo": "losa",
            "tipologia": "maciza_cruzada",
            "luz_transversal_m": 4.0,
            "ancho_losa_m": 4.0,
            "sobrecarga": "vivienda",
            "componentes": [],
            "apoya_en": dict(APOYOS),
        }
        editor = EditorElemento("Losa cruzada", elemento, BIBLIOTECA)
        self.assertEqual(editor.espesor_maciza.value(), 0.0)
        editor.espesor_maciza.setValue(0.18)
        nombre, guardado = editor.resultado()
        componente = next(
            c for c in guardado["componentes"] if c.get("clave") == "Losa_maciza"
        )
        self.assertAlmostEqual(componente["espesor_m"], 0.18)
        self.assertFalse(guardado["espesor_automatico"])
        with patch(
            "calc.losas_pynite.configuracion_obra.modulo_elasticidad_hormigon",
            return_value={"Ec_MPa": 25000.0, "articulo": "prueba"},
        ):
            resultado = losas_macizas.calcular_desde_elemento_carga(
                nombre, guardado, BIBLIOTECA
            )
        self.assertAlmostEqual(resultado.espesor_m, 0.18)

    def test_losa_cruzada_de_tres_bordes_equilibra_reacciones(self):
        elemento = self.elemento("maciza_cruzada")
        del elemento["apoya_en"]["sup"]
        with patch(
            "calc.losas_pynite.configuracion_obra.modulo_elasticidad_hormigon",
            return_value={"Ec_MPa": 25000.0, "articulo": "prueba"},
        ):
            resultado = analizar_paño_macizo(
                elemento, {"D_kNm2": 4.0, "L_kNm2": 2.0}, 0.15
            )
        self.assertEqual(
            {reaccion["apoyo"] for reaccion in resultado["reacciones"]},
            {"izq", "der", "inf"},
        )
        self.assertAlmostEqual(
            sum(r["reaccion_total_D_kN"] for r in resultado["reacciones"]),
            4.0 * 4.0 * 4.0,
            places=6,
        )

    def test_placa_cruzada_requiere_al_menos_tres_bordes(self):
        elemento = self.elemento("maciza_cruzada")
        del elemento["apoya_en"]["sup"]
        del elemento["apoya_en"]["inf"]
        with patch(
            "calc.losas_pynite.configuracion_obra.modulo_elasticidad_hormigon",
            return_value={"Ec_MPa": 25000.0, "articulo": "prueba"},
        ):
            with self.assertRaisesRegex(ValueError, "al menos tres bordes"):
                analizar_paño_macizo(elemento, {"D_kNm2": 4.0, "L_kNm2": 2.0}, 0.15)

    def test_maciza_unidireccional_rechaza_apoyos_perpendiculares_adicionales(self):
        elemento = self.elemento("maciza")
        elemento["apoya_en"]["inf"] = APOYOS["inf"]
        with patch(
            "calc.losas_pynite.configuracion_obra.modulo_elasticidad_hormigon",
            return_value={"Ec_MPa": 25000.0, "articulo": "prueba"},
        ):
            with self.assertRaisesRegex(ValueError, "tipología cruzada"):
                analizar_paño_macizo(
                    elemento, {"D_kNm2": 4.0, "L_kNm2": 2.0}, 0.15
                )

    def test_casetonada_unidireccional_y_automatica_por_luz_mayor_a_siete(self):
        elemento = self.elemento("casetonada")
        elemento["apoya_en"] = {
            "izq": APOYOS["izq"],
            "der": APOYOS["der"],
        }
        elemento["panel_ejes"] = {"modelo_placa": "casetonada_unidireccional"}
        with patch(
            "calc.losas_pynite.configuracion_obra.modulo_elasticidad_hormigon",
            return_value={"Ec_MPa": 25000.0, "articulo": "prueba"},
        ):
            resultado = analizar_paño_macizo(
                elemento, {"D_kNm2": 4.0, "L_kNm2": 2.0}, 0.15
            )
        self.assertEqual(
            resultado["modelo"],
            "casetonada equivalente ortótropa unidireccional (dos bordes apoyados)",
        )
        self.assertEqual(
            {reaccion["apoyo"] for reaccion in resultado["reacciones"]},
            {"izq", "der"},
        )
        with patch(
            "calc.losas_pynite.configuracion_obra.modulo_elasticidad_hormigon",
            return_value={"Ec_MPa": 25000.0, "articulo": "prueba"},
        ):
            resultado_losa = losas_macizas.calcular_desde_elemento_carga(
                "Losa casetonada", elemento, BIBLIOTECA
            )
        self.assertIn("LOSA CASETONADA UNIDIRECCIONAL", resultado_losa.memoria)
        self.assertIn("trabajo unidireccional equivalente", resultado_losa.memoria)

    def test_alivianada_no_se_calcula_unidireccional_con_tres_apoyos(self):
        elemento = self.elemento("alivianada")
        with self.assertRaisesRegex(ValueError, "placa cruzada"):
            losas.calcular_desde_elemento_carga("Losa", elemento, BIBLIOTECA)

    def test_tipologia_auto_segun_cantidad_de_apoyos_y_luz(self):
        seleccionar = EditorLosaDesdeEjes._tipologia_segun_apoyos
        self.assertEqual(seleccionar("alivianada", 3, 5.0), "maciza_cruzada")
        self.assertEqual(seleccionar("maciza", 2, 7.0), "maciza")
        self.assertEqual(seleccionar("maciza", 2, 7.01), "casetonada")
        self.assertEqual(seleccionar("maciza", 4, 9.0), "maciza_cruzada")

    def test_viento_solo_considera_porticos_y_estructura_legacy(self):
        estructura = {
            "Pórtico nuevo": {
                "tipo": "portico",
                "referencia_planta": {"direccion": "y"},
            },
            "Pórtico legacy": {
                "referencia_planta": {"direccion": "y"},
            },
            "Muro 1": {
                "tipo": "muro",
                "referencia_planta": {"direccion": "y"},
            },
        }
        self.assertEqual(
            cargas.porticos_afectados_por_viento(estructura, "y"),
            ["Pórtico nuevo", "Pórtico legacy"],
        )

    def test_editor_selecciona_casetonada_y_cruzada_por_apoyos(self):
        ejes = {
            "x": [
                {"id": "x0", "nombre": "X1", "coordenada_m": 0.0},
                {"id": "x1", "nombre": "X2", "coordenada_m": 8.0},
            ],
            "y": [
                {"id": "y0", "nombre": "Y1", "coordenada_m": 0.0},
                {"id": "y1", "nombre": "Y2", "coordenada_m": 4.0},
            ],
        }

        def portico(nombre, direccion, eje_id, eje_longitudinal_id, fin):
            return {
                "tipo": "portico",
                "referencia_planta": {
                    "direccion": direccion,
                    "eje_id": eje_id,
                    "eje_longitudinal_id": eje_longitudinal_id,
                    "origen_longitudinal_m": 0.0,
                },
                "vigas": {
                    "V1": {
                        "cota_m": 3.0,
                        "tramos": [{"x_inicio": 0.0, "x_fin": fin}],
                    }
                },
            }

        estructura = {
            "PY1": portico("PY1", "y", "x0", "y0", 4.0),
            "PY2": portico("PY2", "y", "x1", "y0", 4.0),
            "PX1": portico("PX1", "x", "y0", "x0", 8.0),
            "PX2": portico("PX2", "x", "y1", "x0", 8.0),
        }
        elemento = {
            "tipologia": "maciza",
            "panel_ejes": {
                "direccion_luz": "x",
                "cota_apoyo_m": 3.0,
                "ejes": {
                    "x": {
                        "inicio": {"id": "x0"},
                        "fin": {"id": "x1"},
                    },
                    "y": {
                        "inicio": {"id": "y0"},
                        "fin": {"id": "y1"},
                    },
                },
            },
        }
        editor = EditorLosaDesdeEjes(ejes, [], estructura, elemento=elemento, nombre="Losa")
        self.addCleanup(editor.close)
        editor._actualizar()

        self.assertEqual(editor.tipologia.currentData(), "casetonada")
        _, geometria = editor.resultado()
        self.assertEqual(geometria["tipologia"], "casetonada")
        self.assertEqual(
            geometria["panel_ejes"]["modelo_placa"],
            "casetonada_unidireccional",
        )
        editor_reabierto = EditorLosaDesdeEjes(
            ejes, [], estructura, elemento=geometria, nombre="Losa"
        )
        self.addCleanup(editor_reabierto.close)
        self.assertEqual(editor_reabierto._tipologia_base, "maciza")
        self.assertEqual(editor_reabierto.resultado()[1]["tipologia"], "casetonada")
        editor.selector_vigas["inf"].setCurrentIndex(1)
        self.assertEqual(editor.tipologia.currentData(), "maciza_cruzada")
        _, geometria = editor.resultado()
        self.assertEqual(geometria["tipologia"], "maciza_cruzada")
        self.assertEqual(set(geometria["apoya_en"]), {"izq", "der", "inf"})
        self.assertEqual(
            geometria["panel_ejes"]["tipologia_base_automatica"], "maciza"
        )
        editor_cruzada_reabierta = EditorLosaDesdeEjes(
            ejes, [], estructura, elemento=geometria, nombre="Losa"
        )
        self.addCleanup(editor_cruzada_reabierta.close)
        self.assertEqual(editor_cruzada_reabierta._tipologia_base, "maciza")
        editor_cruzada_reabierta.selector_vigas["inf"].setCurrentIndex(0)
        self.assertEqual(editor_cruzada_reabierta.tipologia.currentData(), "casetonada")

    def test_reacciones_bidireccionales_llegan_a_pórticos_en_ambos_ejes(self):
        elemento = self.elemento("maciza_cruzada")
        elemento["panel_ejes"] = {
            "direccion_luz": "x",
            "coordenada_libre_inicio_m": 0.0,
            "coordenada_libre_fin_m": 4.0,
            "cota_apoyo_m": 3.0,
            "ejes": {
                "x": {
                    "inicio": {"id": "x0", "coordenada_m": 0.0},
                    "fin": {"id": "x1", "coordenada_m": 4.0},
                },
                "y": {
                    "inicio": {"id": "y0", "coordenada_m": 0.0},
                    "fin": {"id": "y1", "coordenada_m": 4.0},
                },
            },
        }
        elemento["apoya_en"] = {
            "izq": {"tipo": "portico", "portico": "PY0", "viga": "V1"},
            "der": {"tipo": "portico", "portico": "PY4", "viga": "V1"},
            "inf": {"tipo": "portico", "portico": "PX0", "viga": "V1"},
            "sup": {"tipo": "portico", "portico": "PX4", "viga": "V1"},
        }
        estructura = {
            "PX0": {
                "tipo": "portico",
                "posicion_planta_m": 0.0,
                "referencia_planta": {
                    "direccion": "x", "eje_longitudinal_id": "x0",
                    "eje_id": "y0", "origen_longitudinal_m": 0.0,
                },
                "vigas": {"V1": {
                    "cota_m": 3.0,
                    "tramos": [{
                        "id": "T1", "longitud_m": 4.0,
                        "x_inicio": 0.0, "x_fin": 4.0,
                    }],
                }},
            },
            "PY0": {
                "tipo": "portico",
                "posicion_planta_m": 0.0,
                "referencia_planta": {
                    "direccion": "y", "eje_longitudinal_id": "y0",
                    "eje_id": "x0", "origen_longitudinal_m": 0.0,
                },
                "vigas": {"V1": {
                    "cota_m": 3.0,
                    "tramos": [{
                        "id": "T1", "longitud_m": 4.0,
                        "x_inicio": 0.0, "x_fin": 4.0,
                    }],
                }},
            },
        }
        ejes = {
            "x": [
                {"id": "x0", "coordenada_m": 0.0},
                {"id": "x1", "coordenada_m": 4.0},
            ],
            "y": [
                {"id": "y0", "coordenada_m": 0.0},
                {"id": "y1", "coordenada_m": 4.0},
            ],
        }
        anchos = {
            "izq": {"ancho_tributario_m": 0.9},
            "der": {"ancho_tributario_m": 1.1},
            "inf": {"ancho_tributario_m": 1.2},
            "sup": {"ancho_tributario_m": 0.8},
        }

        with (
            patch("calc.portico.rutas.leer_json", return_value=ejes),
            patch(
                "calc.losas_macizas.reacciones_y_anchos_por_borde",
                return_value=anchos,
            ),
        ):
            aplicaciones_x, avisos_x = portico._aplicaciones_losas_portico(
                "PX0", {"Losa": elemento}, estructura
            )
            aplicaciones_y, avisos_y = portico._aplicaciones_losas_portico(
                "PY0", {"Losa": elemento}, estructura
            )

        self.assertFalse(avisos_x)
        self.assertFalse(avisos_y)
        self.assertEqual(len(aplicaciones_x), 1)
        self.assertEqual(aplicaciones_x[0]["id"], "losa:L1:inf:T1")
        self.assertEqual(aplicaciones_x[0]["ancho_tributario_m"], 1.2)
        self.assertEqual(len(aplicaciones_y), 1)
        self.assertEqual(aplicaciones_y[0]["id"], "losa:L1:izq:T1")
        self.assertEqual(aplicaciones_y[0]["ancho_tributario_m"], 0.9)

    def test_reacciones_bidireccionales_suman_en_muro(self):
        elemento = self.elemento("maciza_cruzada")
        elemento["panel_ejes"] = {"cota_apoyo_m": 3.0}
        elemento["apoya_en"] = {
            "izq": {"tipo": "muro", "nombre": "Muro 1"},
            "der": {"tipo": "portico", "portico": "P1", "viga": "V1"},
            "inf": {"tipo": "muro", "nombre": "Muro 1"},
            "sup": {"tipo": "portico", "portico": "P2", "viga": "V1"},
        }
        reacciones = {
            lado: {"cargas_kN_m": {"D": valor, "L": valor / 2}}
            for lado, valor in {
                "izq": 1.0, "der": 2.0, "inf": 3.0, "sup": 4.0,
            }.items()
        }
        with (
            patch("calc.muros.materiales.cargar", return_value=BIBLIOTECA),
            patch(
                "calc.losas_macizas.reacciones_y_anchos_por_borde",
                return_value=reacciones,
            ),
        ):
            resultado = muros.cargas_losas_en_muro(
                "Muro 1", {"elementos": {"Losa cruzada": elemento}}
            )

        self.assertAlmostEqual(resultado["D_kN_m"], 4.0)
        self.assertAlmostEqual(resultado["L_kN_m"], 2.0)
        self.assertEqual(len(resultado["losas"]), 1)
        self.assertEqual(resultado["losas"][0]["cota_m"], 3.0)
        self.assertFalse(resultado["avisos"])

    def test_memoria_alivianada_muestra_reacciones_y_combinaciones(self):
        resultado = losas.ResultadoLosa(
            nombre="Losa alivianada",
            luz_libre_m=4.0,
            luz_calculo_m=4.1,
            ancho_losa_m=5.0,
            serie="Serie 4",
            tipo_vigueta="V1",
            combo="EPS",
            bovedilla_altura="12",
            imagen="",
            D1=2.0,
            D2=0.5,
            D=2.5,
            L=2.0,
            Q=5.0,
            combo_critico="1.2D+1.6L",
            uso_sobrecarga="Vivienda",
            momento_req_kgm=100.0,
            momento_busqueda_kgm=111.0,
            momento_disponible_kgm=120.0,
            phi=0.9,
            momento_reducido_kgm=108.0,
            momento_reducido_kNm=1.06,
            area_losa_m2=20.0,
            n_viguetas=10,
            n_bloques=50,
            vol_hormigon_m3=1.0,
            malla={
                "tamanio_pano": "2.4x6",
                "solape_m": 0.2,
                "area_losa_m2": 20.0,
                "area_pano_efectiva_m2": 13.2,
                "cantidad_mallas": 2,
            },
            nervios={
                "n_nervios": 2,
                "barras_por_nervio": 2,
                "diametro_mm": 10,
                "largo_por_nervio_m": 1.8,
                "barras_total": 4,
                "longitud_total_barras_m": 7.2,
            },
            reacciones=[{
                "apoyo": "izq",
                "apoya_en": {"tipo": "muro", "nombre": "Muro 1"},
                "cargas": {"D": 1.25, "L": 1.0, "W": 0.0},
                "combinaciones": {"1.2D+1.6L": 3.1},
            }],
        )

        memoria = losas.generar_memoria(resultado, BIBLIOTECA)
        self.assertIn("D=1.250, L=1.000, W=0.000 kN/m", memoria)
        self.assertIn("1.2D+1.6L: 3.100 kN/m", memoria)
        self.assertIn("muro Muro 1", memoria)


if __name__ == "__main__":
    unittest.main()

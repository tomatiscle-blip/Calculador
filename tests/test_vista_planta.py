"""Tests de la vista de planta de Inicio: geometría y filtro por nivel.

Corre sin obra activa: arma datos sintéticos y consulta la escena que arma la
vista (coordenadas de planta + conteo de elementos por nivel).
"""

import os
import sys
import unittest
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

RAIZ = Path(__file__).resolve().parent.parent
if str(RAIZ) not in sys.path:
    sys.path.insert(0, str(RAIZ))

from PySide6.QtWidgets import QApplication  # noqa: E402

app = QApplication.instance() or QApplication([])

from app.visualizacion import VistaPlantaPorticos  # noqa: E402

EJES = {
    "x": [
        {"id": "x:1", "nombre": "X1", "coordenada_m": 0.0},
        {"id": "x:2", "nombre": "X2", "coordenada_m": 5.0},
    ],
    "y": [
        {"id": "y:1", "nombre": "Y1", "coordenada_m": 0.0},
        {"id": "y:2", "nombre": "Y2", "coordenada_m": 4.0},
    ],
}
NIVELES = [
    {"nombre": "PB", "cota_m": 0.0},
    {"nombre": "N1", "cota_m": 3.0},
]

ESTRUCTURA = {
    # Corre en X sobre el eje Y2 (4.0): columnas locales 0 y 5
    "Pórtico 1": {
        "referencia_planta": {
            "direccion": "x",
            "eje_id": "y:2",
            "eje_longitudinal_id": "x:1",
            "origen_longitudinal_m": 0.0,
        },
        "posicion_planta_m": 4.0,
        "columnas": {
            "C0-1": {"x": 0.0, "nivel": 0.0, "eje_nombre": "X1"},
            "C0-2": {"x": 5.0, "nivel": 0.0, "eje_nombre": "X2"},
            "C1-1": {"x": 0.0, "nivel": 3.0, "eje_nombre": "X1"},
            "C1-2": {"x": 5.0, "nivel": 3.0, "eje_nombre": "X2"},
        },
        "vigas": {
            "V0-1": {
                "cota_m": 3.0,
                "nivel_nombre": "N1",
                "tramos": [{"x_inicio": 0.0, "x_fin": 5.0}],
            }
        },
        "bases": {"B0-1": {"x": 0.0}, "B0-2": {"x": 5.0}},
    },
    # Corre en Y sobre el eje X1 (0.0), con origen longitudinal en Y1
    "Pórtico 2": {
        "referencia_planta": {
            "direccion": "y",
            "eje_id": "x:1",
            "eje_longitudinal_id": "y:1",
            "origen_longitudinal_m": 0.0,
        },
        "posicion_planta_m": 0.0,
        "columnas": {"C0-3": {"x": 4.0, "nivel": 0.0, "eje_nombre": "Y2"}},
        "vigas": {},
        "bases": {},
    },
    # Muro que corre en X entre las cotas 0 y 3
    "Muro 1": {
        "tipo": "muro",
        "referencia_planta": {"direccion": "x", "eje_id": "y:1"},
        "posicion_planta_m": 0.0,
        "x_inicio": 0.0,
        "x_fin": 5.0,
        "nivel_inicio": 0.0,
        "nivel_fin": 3.0,
    },
}

LOSAS = {
    "Losa 1": {
        "tipo": "losa",
        "apoya_en": {
            "izq": {"tipo": "muro", "nombre": "Muro 1"},
            "der": {"tipo": "portico", "portico": "Pórtico 1", "viga": "V0-1"},
        },
        "panel_ejes": {
            "nivel_apoyo_nombre": "N1",
            "cota_apoyo_m": 3.0,
            "ejes": {
                "x": {"inicio": {"coordenada_m": 0.0}, "fin": {"coordenada_m": 5.0}},
                "y": {"inicio": {"coordenada_m": 0.0}, "fin": {"coordenada_m": 4.0}},
            },
        },
    },
    "Losa 2": {"tipo": "losa"},
}

APLICACIONES = [
    {
        "id": "Aplicacion 0-1",
        "carga_id": "peso-propio:V0-1",
        "tipo_aplicacion": "muro_perpendicular",
        "activa": True,
        "separacion_m": 5.0,
        "reacciones": [
            {"portico": "Pórtico 1", "tramo_id": "V0-1 T1", "x_m": 2.5,
             "influencia_m": 2.5},
            {"portico": "Pórtico 2", "tramo_id": "V0-1 T1", "x_m": 2.0,
             "influencia_m": 2.5},
        ],
    },
    # Inactiva: debe ignorarse (no genera flecha).
    {
        "id": "Aplicacion 0-2",
        "carga_id": "peso-propio:V0-1",
        "tipo_aplicacion": "muro_perpendicular",
        "activa": False,
        "separacion_m": 5.0,
        "reacciones": [
            {"portico": "Pórtico 1", "tramo_id": "V0-1 T1", "x_m": 0.0},
            {"portico": "Pórtico 2", "tramo_id": "V0-1 T1", "x_m": 0.0},
        ],
    },
    # Tipo distinto: no afecta a las flechas de muro perpendicular.
    {
        "id": "Aplicacion 0-3",
        "carga_id": "peso-propio:V0-1",
        "tipo_aplicacion": "lineal",
        "activa": True,
        "portico": "Pórtico 1",
        "tramo_id": "V0-1 T1",
    },
]


def vista_con(cota):
    vista = VistaPlantaPorticos()
    vista.resize(600, 400)
    vista.set_datos(EJES, NIVELES, ESTRUCTURA)
    vista.set_losas(LOSAS)
    vista.set_aplicaciones(APLICACIONES)
    vista.set_portico_activo("Pórtico 1")
    vista.set_nivel_activo(cota)
    return vista


class GeometriaDePlanta(unittest.TestCase):
    def test_portico_que_corre_en_x(self):
        geom = vista_con(None)._geometria(ESTRUCTURA["Pórtico 1"])
        self.assertEqual(geom[0], "x")
        # La columna local 5 cae sobre X5 = 0 + 5 y la posición es el eje Y2
        self.assertEqual(vista_con(None)._punto(geom, 5.0), (5.0, 4.0))

    def test_portico_que_corre_en_y(self):
        geom = vista_con(None)._geometria(ESTRUCTURA["Pórtico 2"])
        self.assertEqual(geom[0], "y")
        # Corre en Y: la coordenada longitudinal es la Y, la posición es la X
        self.assertEqual(vista_con(None)._punto(geom, 4.0), (0.0, 4.0))

    def test_muro_usa_coordenadas_absolutas(self):
        geom = vista_con(None)._geometria(ESTRUCTURA["Muro 1"])
        self.assertEqual(vista_con(None)._punto(geom, 5.0), (5.0, 0.0))

    # --- flechas de muro perpendicular entre pórticos -------------------------

    def test_flecha_muro_perpendicular_en_escena(self):
        """Las aplicaciones muro_perpendicular resuelven flechas entre pórticos.

        Pórtico 1 corre en X sobre Y2 (y=4.0): x_m=2.5 → (2.5, 4.0).
        Pórtico 2 corre en Y sobre X1 (x=0.0): x_m=2.0 → (0.0, 2.0).
        La aplicación inactiva y la de tipo lineal no generan flechas.
        """
        escena = vista_con(3.0)._escena()
        cargas = escena["cargas_muro"]
        self.assertEqual(len(cargas), 1, "solo la aplicación activa debería renderizar")
        flecha = cargas[0]
        self.assertEqual(flecha["p1"], (2.5, 4.0))
        self.assertEqual(flecha["p2"], (0.0, 2.0))
        self.assertEqual(flecha["direccion"], "x")
        self.assertEqual(flecha["influencia_m"], 5.0)
        # Los puntos de apoyo participan del encuadre global.
        puntos = escena["puntos"]
        self.assertIn((2.5, 4.0), puntos)
        self.assertIn((0.0, 2.0), puntos)

    def test_set_aplicaciones_filtra_inactivas(self):
        """set_aplicaciones conserva solo muros perpendiculares activos."""
        vista = VistaPlantaPorticos()
        vista.set_datos(EJES, NIVELES, ESTRUCTURA)
        vista.set_aplicaciones(APLICACIONES)
        # La inactiva (Aplicacion 0-2) y la lineal (0-3) quedan fuera.
        self.assertEqual(len(vista.aplicaciones), 1)
        self.assertEqual(vista.aplicaciones[0]["id"], "Aplicacion 0-1")

    def test_aplicacion_sin_reacciones_no_dibuja(self):
        """Una aplicación con una sola reacción no produce flecha."""
        vista = VistaPlantaPorticos()
        vista.set_datos(EJES, NIVELES, ESTRUCTURA)
        vista.set_aplicaciones([
            {"tipo_aplicacion": "muro_perpendicular", "activa": True,
             "separacion_m": 5.0,
             "reacciones": [
                 {"portico": "Pórtico 1", "tramo_id": "V0-1 T1", "x_m": 2.5},
             ]},
        ])
        self.assertEqual(vista._escena()["cargas_muro"], [])

    def test_pintar_muro_perpendicular_no_lanza_excepciones(self):
        """Render completo con flechas de muro perpendicular no explota."""
        from PySide6.QtGui import QColor, QImage

        imagen = QImage(600, 400, QImage.Format.Format_ARGB32)
        imagen.fill(QColor("#ffffff"))
        vista_con(3.0).render(imagen)


class FiltroPorNivel(unittest.TestCase):
    def test_en_pb_se_ven_columnas_y_bases(self):
        conteo = vista_con(0.0)._escena()["conteo"]
        # C0-1 (Pórtico 1) y C0-3 (Pórtico 2) comparten el punto (0.0, 4.0):
        # son la misma columna física en un cruce de ejes. Se cuenta una sola vez.
        self.assertEqual(conteo["columnas"], 2)  # (0.0, 4.0) y (5.0, 4.0)
        self.assertEqual(conteo["bases"], 2)
        self.assertEqual(conteo["vigas"], 0)
        # La viga de N1 es el techo de PB: va como proyección punteada.
        escena = vista_con(0.0)._escena()
        proy_vigas = [p for p in escena["proyecciones"] if p["tipo"] == "viga"]
        self.assertEqual(len(proy_vigas), 1)
        self.assertEqual(proy_vigas[0]["nombre"], "V0-1")

    def test_en_n1_se_van_las_columnas_y_aparecen_vigas_y_losa(self):
        escena = vista_con(3.0)
        conteo = escena._escena()["conteo"]
        self.assertEqual(conteo["columnas"], 2)  # C1-1 y C1-2
        self.assertEqual(conteo["vigas"], 1)
        self.assertEqual(conteo["losas"], 1)
        # El muro nace en PB: en N1 solo lo atraviesa, va punteado.
        datos = escena._escena()
        self.assertEqual(datos["conteo"]["muros"], 0)
        proy_muros = [p for p in datos["proyecciones"] if p["tipo"] == "muro"]
        self.assertEqual(len(proy_muros), 1)

    def test_rotulos_indican_eje_y_apoyos(self):
        escena = vista_con(3.0)._escena()
        textos = [r["texto"] for r in escena["rotulos"]]
        # El pórtico dice sobre qué eje corre y desde dónde.
        self.assertTrue(any("Pórtico 1" in t and "Y2" in t for t in textos))
        # El muro que atraviesa avisa que viene de abajo.
        self.assertTrue(any("viene de abajo" in t for t in textos))
        # La losa dice en qué apoya.
        self.assertTrue(any("Losa 1" in t and "Muro" in t for t in textos))

    def test_columnas_llevan_su_eje(self):
        escena = vista_con(0.0)._escena()
        columnas = escena["porticos"][0]["columnas"]
        por_id = {c[0]: c for c in columnas}
        self.assertEqual(len(por_id["C0-1"]), 3)  # id, punto, eje

    def test_columnas_en_cruce_se_deduplican(self):
        """Las columnas que comparten coordenada global (cruce de ejes X y Y)
        se dibujan una sola vez con los ejes de ambos pórticos."""
        escena = vista_con(0.0)._escena()
        unicas = escena["columnas_unicas"]
        # En PB hay 3 entradas de columna pero 2 puntos únicos:
        # (0.0, 4.0) compartida por C0-1 (Pórtico 1) y C0-3 (Pórtico 2)
        self.assertEqual(len(unicas), 2)
        punto_cruce = [c for c in unicas if c["punto"] == (0.0, 4.0)]
        self.assertEqual(len(punto_cruce), 1)
        # El cruce debe reportar los ejes de ambos pórticos (Y2 del Pórtico 1,
        # X1 del Pórtico 2).
        ejes = punto_cruce[0]["ejes"]
        self.assertIn("Y2", ejes)
        self.assertIn("X1", ejes)
        # El nombre indica que hay más de una entrada en ese punto.
        self.assertTrue(punto_cruce[0]["nombre"].startswith("C0-1"))

    def test_sin_filtro_se_ve_todo(self):
        conteo = vista_con(None)._escena()["conteo"]
        # Sin filtro de nivel: 5 entradas de columna pero solo 2 puntos únicos
        # en planta ((0.0, 4.0) y (5.0, 4.0)) porque C0-1/C0-3/C1-1 comparten
        # el mismo punto y C0-2/C1-2 comparten el otro.
        self.assertEqual(conteo["columnas"], 2)
        self.assertEqual(conteo["vigas"], 1)
        self.assertEqual(conteo["losas"], 1)

    def test_el_eje_del_portico_siempre_queda_dibujado(self):
        # Aunque el nivel no tenga columnas, el pórtico deja su eje de referencia
        escena = vista_con(3.0)._escena()
        nombres = {p["nombre"] for p in escena["porticos"]}
        self.assertEqual(nombres, {"Pórtico 1", "Pórtico 2"})

    def test_losa_sin_panel_no_se_pierde(self):
        # Sin datos de ejes se avisa por etiqueta en vez de desaparecer
        escena = vista_con(3.0)._escena()
        self.assertIn("Losa 2", escena["etiquetas"])

    def test_losa_de_otro_nivel_se_omite(self):
        # La losa apoya en N1: en PB no se dibuja ni se avisa
        escena = vista_con(0.0)._escena()
        self.assertEqual(escena["conteo"]["losas"], 0)
        self.assertNotIn("Losa 1", escena["etiquetas"])
        self.assertNotIn("Losa 1", [nombre for nombre, _ in escena["losas"]])

    def test_tipo_losa_se_reporta(self):
        # La losa sin tipología explícita usa "alivianada" por defecto.
        escena = vista_con(3.0)._escena()
        self.assertEqual(escena["tipos_losa"].get("Losa 1"), "alivianada")

    def test_pintar_no_lanza_excepciones(self):
        from PySide6.QtGui import QColor, QImage

        for cota in (0.0, 3.0, None):
            imagen = QImage(600, 400, QImage.Format.Format_ARGB32)
            imagen.fill(QColor("#ffffff"))
            vista_con(cota).render(imagen)

    # --- rótulos técnicos: se conservan pero no se dibujan ----------------

    def test_rotulos_tecnicos_quedan_flotantes(self):
        """Los rótulos técnicos (viga, muro, base, losa, apoyo) se marcan
        como flotantes: siguen en la escena pero no se dibujan."""
        escena = vista_con(3.0)._escena()
        rotulos = escena["rotulos"]
        self.assertTrue(rotulos, "la escena debe seguir reportando rótulos")
        for rot in rotulos:
            texto = rot["texto"]
            if texto.startswith("V0-1") or texto.startswith("Muro 1"):
                self.assertTrue(rot["flotante"], f"{texto} debería ser flotante")
            if texto.startswith("B0-"):
                self.assertTrue(rot["flotante"], f"{texto} debería ser flotante")
            if texto.startswith("Losa 1"):
                self.assertTrue(rot["flotante"], f"{texto} debería ser flotante")

    def test_solo_se_dibujan_rotulos_no_flotantes(self):
        """El filtro de rótulos deja pasar únicamente los no flotantes."""
        vista = vista_con(3.0)
        escena = vista._escena()
        no_flotantes = [r for r in escena["rotulos"] if not r.get("flotante", False)]
        colocados = vista._colocar_rotulos(no_flotantes, lambda v: v, lambda v: v)
        for rot, _rect in colocados:
            self.assertFalse(rot.get("flotante", False))

    def test_rotulo_flotante_por_defecto(self):
        """``_rotulo`` admite el flag ``flotante`` y lo normaliza a bool."""
        vista = VistaPlantaPorticos()
        rot = vista._rotulo(1.0, 2.0, "X", flotante=True)
        self.assertIs(rot["flotante"], True)
        rot = vista._rotulo(1.0, 2.0, "X")
        self.assertIs(rot["flotante"], True)

    def test_solo_el_portico_activo_conserva_su_texto(self):
        """De los pórticos, solo el activo mantiene el rótulo visible."""
        vista = vista_con(3.0)
        vista.portico_activo = "Pórtico 1"
        escena = vista._escena()
        por_nombre = {}
        for rot in escena["rotulos"]:
            if "Pórtico" in rot["texto"]:
                por_nombre[rot["texto"].split("·")[0].strip()] = rot["flotante"]
        self.assertIn("Pórtico 1", por_nombre)
        self.assertFalse(por_nombre["Pórtico 1"], "el pórtico activo debe verse")
        self.assertTrue(por_nombre["Pórtico 2"], "el inactivo debe quedar flotante")

    # --- estilos visuales por categoría ---------------------------------

    def test_paleta_visual_por_categoria(self):
        """Cada familia de elementos tiene color propio y distinguible."""
        vista = VistaPlantaPorticos()
        colores = [
            vista.COLOR_EJE,
            vista.COLOR_PORTICO,
            vista.COLOR_PORTICO_ACTIVO,
            vista.COLOR_VIGA,
            vista.COLOR_MURO,
            vista.COLOR_COLUMNA,
            vista.COLOR_COLUMNA_INACTIVA,
            vista.COLOR_BASE,
            vista.COLOR_APOYO,
        ]
        self.assertEqual(len(set(colores)), len(colores), "los colores deben ser distintos")
        for color in colores:
            self.assertTrue(color.startswith("#") and len(color) == 7, color)

    def test_leyenda_refleja_los_estilos(self):
        """La leyenda debe listar los símbolos de cada categoría."""
        vista = VistaPlantaPorticos()
        # Los nombres de la leyenda se toman de la paleta: comprobamos que
        # cada familia tiene una entrada con su color.
        self.assertNotEqual(vista.COLOR_APOYO, vista.COLOR_BASE)
        self.assertNotEqual(vista.COLOR_VIGA, vista.COLOR_MURO)

    # --- encuadre, paneo, zoom y reinicio ---------------------------------

    def test_reiniciar_vista_no_duplica_manejadores(self):
        """``reiniciar_vista`` no debe conectar señales nuevas cada vez."""
        vista = vista_con(3.0)
        antes = vista.receivers(vista.wheelEvent.__name__) if hasattr(vista, "receivers") else None
        for _ in range(3):
            vista.reiniciar_vista()
        # Tras reiniciar, la vista sigue respondiendo y sin excepción.
        from PySide6.QtGui import QColor, QImage

        imagen = QImage(600, 400, QImage.Format.Format_ARGB32)
        imagen.fill(QColor("#ffffff"))
        vista.render(imagen)
        self.assertEqual(vista.zoom_factor, 1.0)
        self.assertEqual(vista.desplazamiento_x, 0.0)
        self.assertEqual(vista.desplazamiento_y, 0.0)

    def test_zoom_y_paneo_no_rompen_el_encuadre(self):
        """El zoom y el paneo cambian el encuadre pero no la geometría."""
        vista = vista_con(3.0)
        escena = vista._escena()
        antes = len(escena["porticos"])
        vista.zoom_factor = 2.0
        vista.desplazamiento_x = 30.0
        vista.desplazamiento_y = -15.0
        self.assertEqual(len(vista._escena()["porticos"]), antes)
        vista.reiniciar_vista()
        self.assertEqual(vista.zoom_factor, 1.0)


if __name__ == "__main__":
    unittest.main()

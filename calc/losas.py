"""
calc.losas — Cálculo clásico de losas alivianadas (la ÚNICA cuenta de losa).

Es la versión "de programa" de `L00_Losas_alivianadas.py`: el mismo desarrollo
de siempre (cargas → momento → búsqueda en la tabla de viguetas → cómputo de
materiales), pero como una función `calcular(datos)` que recibe los datos por
parámetro en lugar de pedirlos por teclado. Con eso la puede usar la consola,
la ventana o un test, y se puede calcular UNA losa sola.

De dónde sale cada cosa (un dato, un lugar):
    datos/viguetas.json      la tabla de momentos admisibles por serie y tipo
    datos/materiales.json    Pisos / Contrapiso / Cielorrasos, las sobrecargas y
                             el cómputo por m² (grupo Viguetas_Bovedillas)
    datos/losas.json         los datos de cada losa (materiales, luz, ancho…)

Uso:
    py -m calc.losas --lista          qué losas hay en datos/losas.json
    py -m calc.losas L00              muestra la memoria de la losa L00
    py -m calc.losas L00 --guardar    además la guarda en salidas/losas
    py -m calc.losas                  calcula todas las losas y muestra un resumen
"""

from __future__ import annotations

import csv
import math
import sys
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path

from . import cargas, materiales, rutas

# ---------------------------------------------------------------------------
# Constantes del desarrollo clásico (los mismos valores de siempre)
# ---------------------------------------------------------------------------
PHI = 0.9                        # factor de reducción φ del hormigón
KG_M_POR_KN_M = 101.97           # kN·m → kg·m (así está la tabla de viguetas)
TOLERANCIA_LUZ_M = 0.10          # luz de cálculo = luz libre + 0.10 m
PESO_PROPIO_DEFECTO = 1.881      # kN/m² de peso propio de la losa alivianada
                                 # (REVISAR: debería salir de materiales.json,
                                 #  grupo Forjados/clave Losa_alivianada)
TRAMO_MAX_NERVIO_M = 1.80
PANIO_MALLA_DEFECTO = "2.40x6"
SOLAPE_MALLA_DEFECTO = 0.20
UMBRAL_FRACCION_MALLA = 0.50


# ---------------------------------------------------------------------------
# El resultado de una losa
# ---------------------------------------------------------------------------
@dataclass
class ResultadoLosa:
    """Todo lo que devuelve el cálculo de una losa alivianada."""

    nombre: str
    luz_libre_m: float
    luz_calculo_m: float
    ancho_losa_m: float
    serie: str
    tipo_vigueta: str
    combo: str
    bovedilla_altura: str
    imagen: str

    D1: float
    D2: float
    D: float
    L: float
    Q: float
    combo_critico: str
    uso_sobrecarga: str

    momento_req_kgm: float
    momento_busqueda_kgm: float
    momento_disponible_kgm: float
    phi: float
    momento_reducido_kgm: float
    momento_reducido_kNm: float

    area_losa_m2: float
    n_viguetas: int
    n_bloques: int
    vol_hormigon_m3: float

    seleccion: list = field(default_factory=list)
    malla: dict = field(default_factory=dict)
    nervios: dict = field(default_factory=dict)
    reacciones: list = field(default_factory=list)
    memoria: str = ""
    carga_id: str | None = None


# ---------------------------------------------------------------------------
# Tablas
# ---------------------------------------------------------------------------
def cargar_viguetas() -> dict:
    """Lee datos/viguetas.json (la tabla de momentos admisibles)."""
    datos = rutas.leer_json(rutas.VIGUETAS, por_defecto=None)
    if not datos:
        raise FileNotFoundError(f"No se pudo leer la tabla de viguetas: {rutas.VIGUETAS}")
    return datos


def cargar_losas() -> dict:
    """Las losas definidas en datos/losas.json, por nombre."""
    datos = rutas.leer_json(rutas.LOSAS, por_defecto={}) or {}
    return datos.get("losas", {})


# ---------------------------------------------------------------------------
# Funciones puras (el desarrollo clásico)
# ---------------------------------------------------------------------------
def buscar_serie_mas_cercana(luz_calculo: float, largo_a_serie: dict) -> tuple[str, float]:
    """La serie de vigueta que corresponde a la luz de cálculo (tolera decimales)."""
    claves = [(float(k), k) for k in largo_a_serie.keys()]
    claves.sort(key=lambda x: x[0])
    clave_cercana, clave_original = min(
        claves, key=lambda x: (abs(x[0] - luz_calculo), x[0])
    )
    return largo_a_serie[clave_original], clave_cercana


def seleccionar_vigueta_optima(serie_mom: dict, momento_requerido: float):
    """La vigueta más chica que alcanza el momento requerido (sin pasarse de más)."""
    opciones = []
    for tipo in serie_mom:
        for combo, mom in serie_mom[tipo].items():
            if mom >= momento_requerido:
                opciones.append((mom, tipo, combo))
    if not opciones:
        return None
    return min(opciones, key=lambda x: x[0] - momento_requerido)


def _buscar_material(biblioteca: dict, categoria: str, nombre: str) -> dict | None:
    for mat in biblioteca.get(categoria, []):
        if mat.get("nombre", "").strip() == str(nombre).strip():
            return mat
    return None


def carga_de_material(material: dict) -> float:
    """kN/m² de un material: superficial = su valor; volumétrico = valor × espesor."""
    if material.get("tipo") == "volumetrico":
        espesor = float(material.get("espesor_m", 0.1))
        return float(material["valor"]) * espesor
    return float(material["valor"])


def calcular_D2(seleccion: list, biblioteca: dict) -> float:
    """Suma las cargas accesorias (pisos, contrapiso, cielorraso) en kN/m²."""
    total = 0.0
    for sel in seleccion:
        categoria = sel["categoria"]
        material = _buscar_material(biblioteca, categoria, sel["nombre"])
        if material is None:
            raise KeyError(
                f"El material '{sel['nombre']}' no está en la categoría '{categoria}' "
                f"de {rutas.MATERIALES.name}"
            )
        espesor = sel.get("espesor_m")
        if espesor is not None and material.get("tipo") == "volumetrico":
            total += float(material["valor"]) * float(espesor)
        else:
            total += carga_de_material(material)
    return total


def valor_sobrecarga(seleccion: str, biblioteca: dict) -> tuple[float, str]:
    """Devuelve (valor en kN/m², uso) de la sobrecarga elegida por clave o por uso."""
    for sobrecarga in biblioteca.get("Sobrecargas", []):
        if sobrecarga.get("clave") == seleccion or sobrecarga.get("uso") == seleccion:
            return float(sobrecarga["valor_kNm2"]), sobrecarga.get("uso", str(seleccion))
    raise KeyError(f"No existe la sobrecarga '{seleccion}' en {rutas.MATERIALES.name}")


def computo_malla(
    ancho_losa: float,
    lc: float,
    tamanio_pano: str = PANIO_MALLA_DEFECTO,
    solape_m: float = SOLAPE_MALLA_DEFECTO,
    umbral_fraccion: float = UMBRAL_FRACCION_MALLA,
) -> dict:
    """Cómputo de malla electrosoldada por área (m²)."""
    if tamanio_pano == "2.40x6":
        pano_ancho, pano_largo = 2.40, 6.00
    elif tamanio_pano == "2.40x3":
        pano_ancho, pano_largo = 2.40, 3.00
    else:
        raise ValueError(f"Tamaño de paño no reconocido: {tamanio_pano}")

    area_losa = ancho_losa * lc
    area_pano_efectiva = (pano_ancho - solape_m) * (pano_largo - solape_m)

    n_teorico = area_losa / area_pano_efectiva
    n_enteros = int(n_teorico)
    fraccion = n_teorico - n_enteros

    if n_enteros == 0:
        n_mallas = 1
        nota = None
    elif fraccion <= umbral_fraccion:
        n_mallas = n_enteros
        nota = (
            "Área remanente de malla a completar con barras Ø6 "
            "según criterio habitual de obra."
        )
    else:
        n_mallas = n_enteros + 1
        nota = None

    area_total_mallas = n_mallas * (pano_ancho * pano_largo)

    return {
        "modelo": "Malla electrosoldada Q-131 / R-131 (SIMA)",
        "tamanio_pano": tamanio_pano,
        "solape_m": solape_m,
        "area_losa_m2": round(area_losa, 2),
        "area_pano_efectiva_m2": round(area_pano_efectiva, 2),
        "cantidad_mallas": n_mallas,
        "area_total_mallas_m2": round(area_total_mallas, 2),
        "criterio": "Cómputo por área (m²)",
        "nota": nota,
    }


def computo_nervios_refuerzo(
    luz_calculo: float,
    ancho_losa: float,
    max_tramo_m: float = TRAMO_MAX_NERVIO_M,
    barras_por_nervio: int = 2,
    diametro_mm: int = 8,
) -> dict:
    """Nervios de refuerzo para que no queden tramos mayores al máximo."""
    n_nervios = max(math.ceil(luz_calculo / max_tramo_m) - 1, 0)
    largo_por_nervio_m = ancho_losa
    barras_total = n_nervios * barras_por_nervio
    longitud_total_barras_m = barras_total * largo_por_nervio_m

    return {
        "n_nervios": n_nervios,
        "barras_por_nervio": barras_por_nervio,
        "diametro_mm": diametro_mm,
        "largo_por_nervio_m": largo_por_nervio_m,
        "barras_total": barras_total,
        "longitud_total_barras_m": longitud_total_barras_m,
    }



# ---------------------------------------------------------------------------
# El cálculo de una losa
# ---------------------------------------------------------------------------
def calcular(datos: dict, biblioteca: dict | None = None, tabla_viguetas: dict | None = None) -> ResultadoLosa:
    """
    Calcula una losa alivianada a partir de un diccionario de datos (sin input()).

    `datos` es un renglón de datos/losas.json, con esta forma:
        {
          "nombre": "L00",
          "luz_libre_m": 5.0,
          "ancho_losa_m": 6.0,
          "peso_propio_kNm2": 1.881,          (opcional)
          "materiales": [                      (cargas accesorias D2)
            {"categoria": "Pisos", "nombre": "Ceramica 12mm + pegamento"},
            {"categoria": "Contrapiso", "nombre": "Carpeta nivelacion 2.5cm"},
            {"categoria": "Cielorrasos", "nombre": "Yeso Suspendido"}
          ],
          "sobrecarga": "vivienda",            (clave o uso)
          "malla": {"tamanio_pano": "2.40x6", "solape_m": 0.20}   (opcional)
        }
    """
    bib = biblioteca if biblioteca is not None else materiales.cargar()
    viguetas = tabla_viguetas if tabla_viguetas is not None else cargar_viguetas()

    if "nombre" not in datos:
        raise KeyError("A la losa le falta el nombre en los datos")

    nombre = datos["nombre"]
    luz_libre = float(datos["luz_libre_m"])
    luz_calculo = round(luz_libre + TOLERANCIA_LUZ_M, 1)
    ancho_losa = float(datos["ancho_losa_m"])
    seleccion = datos.get("materiales", [])
    phi = float(datos.get("phi", PHI))
    D1 = float(datos.get("peso_propio_kNm2", PESO_PROPIO_DEFECTO))

    # Cargas permanentes (D) y sobrecarga (L)
    D2 = calcular_D2(seleccion, biblioteca=bib)
    D = D1 + D2
    L, uso_sobrecarga = valor_sobrecarga(datos["sobrecarga"], biblioteca=bib)

    # Combinación crítica de diseño (1.2D+1.6L vs 1.4D), como el cálculo original
    Q1 = 1.2 * D + 1.6 * L
    Q2 = 1.4 * D
    if Q1 >= Q2:
        Q, combo_critico = Q1, "1.2D + 1.6L"
    else:
        Q, combo_critico = Q2, "1.4D"

    # Serie de vigueta según la luz de cálculo
    serie, _luz_usada = buscar_serie_mas_cercana(luz_calculo, viguetas["largo_a_serie"])

    # Momento requerido por vigueta (kN·m → kg·m para coincidir con la tabla)
    momento_req = (Q * luz_calculo ** 2 / 8) * KG_M_POR_KN_M
    momento_busqueda = momento_req / phi

    elegido = seleccionar_vigueta_optima(viguetas[serie], momento_busqueda)
    if elegido is None:
        raise ValueError(
            f"Ninguna vigueta de la serie {serie} cumple con el momento requerido "
            f"({momento_busqueda:.2f} kg·m)"
        )
    mom_disponible, tipo, combo = elegido

    momento_reducido_kgm = mom_disponible * phi
    momento_reducido_kNm = momento_reducido_kgm / KG_M_POR_KN_M

    bovedilla_altura = combo.split("_")[1]
    imagen = f"imagenes/{tipo}_bovedilla_{bovedilla_altura}.png"

    # Cómputo de materiales por m² (grupo Viguetas_Bovedillas de materiales.json)
    tipo_normalizado = tipo.replace("vigueta_", "")
    computo = next(
        (
            item for item in bib.get("Viguetas_Bovedillas", [])
            if item.get("tipo_vigueta") == tipo_normalizado
            and item.get("bovedilla") == int(bovedilla_altura)
        ),
        None,
    )
    if computo is None:
        raise ValueError(
            f"No hay cómputo por m² para '{tipo_normalizado}' bovedilla {bovedilla_altura} "
            f"en {rutas.MATERIALES.name}"
        )

    area_losa = luz_calculo * ancho_losa
    n_viguetas = round(ancho_losa * computo["viguetas_por_m2"])
    n_bloques = round(area_losa * computo["bloques_por_m2"])
    vol_hormigon = area_losa * computo["hormigon_m3_m2"]

    config_malla = datos.get("malla", {})
    malla = computo_malla(
        ancho_losa,
        luz_calculo,
        tamanio_pano=config_malla.get("tamanio_pano", PANIO_MALLA_DEFECTO),
        solape_m=float(config_malla.get("solape_m", SOLAPE_MALLA_DEFECTO)),
    )
    nervios = computo_nervios_refuerzo(luz_calculo, ancho_losa)

    # Reacciones de apoyo (Esquema A: la losa entrega D/L/W, SIN combinar).
    # Losa unidireccional → 2 apoyos; cada uno recibe q · L/2 (kN/m).
    W_sup = float(datos.get("viento_kNm2", 0.0))
    factor = luz_calculo / 2.0
    apoya = datos.get("apoya_en", {})
    reacciones = [
        {
            "apoyo": lado,
            "apoya_en": (apoya.get(lado) if isinstance(apoya, dict) else None),
            "cargas": {"D": D * factor, "L": L * factor, "W": W_sup * factor},
            "combinaciones": cargas.combinaciones_de_componentes(D * factor, L * factor, W_sup * factor),
        }
        for lado in ("izq", "der")
    ]

    resultado = ResultadoLosa(
        nombre=nombre,
        luz_libre_m=luz_libre,
        luz_calculo_m=luz_calculo,
        ancho_losa_m=ancho_losa,
        serie=serie,
        tipo_vigueta=tipo,
        combo=combo,
        bovedilla_altura=bovedilla_altura,
        imagen=imagen,
        D1=D1,
        D2=D2,
        D=D,
        L=L,
        Q=Q,
        combo_critico=combo_critico,
        uso_sobrecarga=uso_sobrecarga,
        momento_req_kgm=momento_req,
        momento_busqueda_kgm=momento_busqueda,
        momento_disponible_kgm=mom_disponible,
        phi=phi,
        momento_reducido_kgm=momento_reducido_kgm,
        momento_reducido_kNm=momento_reducido_kNm,
        area_losa_m2=area_losa,
        n_viguetas=n_viguetas,
        n_bloques=n_bloques,
        vol_hormigon_m3=vol_hormigon,
        seleccion=seleccion,
        malla=malla,
        nervios=nervios,
        reacciones=reacciones,
        carga_id=datos.get("carga_id"),
    )
    resultado.memoria = generar_memoria(resultado, biblioteca=bib)
    return resultado


def calcular_desde_elemento_carga(
    nombre: str,
    elemento: dict,
    biblioteca: dict | None = None,
    tabla_viguetas: dict | None = None,
) -> ResultadoLosa:
    """Dimensiona una losa alivianada usando la composición ya definida en Cargas."""
    if elemento.get("tipo") != "losa" or elemento.get("tipologia") != "alivianada":
        raise ValueError("El cálculo de viguetas está disponible para losas alivianadas.")
    luz = float(elemento.get("luz_transversal_m", 0.0))
    ancho = float(elemento.get("ancho_losa_m", 0.0))
    if luz <= 0 or ancho <= 0:
        raise ValueError("Completá la luz entre apoyos y el ancho del paño en la composición de la losa.")

    bib = biblioteca if biblioteca is not None else materiales.cargar()
    componentes = elemento.get("componentes", [])
    base = [c for c in componentes if c.get("grupo") == "Forjados"]
    if not any(c.get("clave") == "Losa_alivianada" for c in base):
        raise ValueError("Agregá la capa 'Losa alivianada' del grupo Forjados a la composición.")
    accesorios = [c for c in componentes if c.get("grupo") != "Forjados"]
    seleccion = []
    for componente in accesorios:
        material = materiales.buscar(componente["grupo"], componente["clave"], bib)
        sel = {"categoria": componente["grupo"], "nombre": material["nombre"]}
        if "espesor_m" in componente:
            sel["espesor_m"] = componente["espesor_m"]
        seleccion.append(sel)

    datos = {
        "nombre": nombre,
        "carga_id": elemento.get("id"),
        "luz_libre_m": luz,
        "ancho_losa_m": ancho,
        "materiales": seleccion,
        "sobrecarga": elemento.get("sobrecarga", "vivienda"),
        "peso_propio_kNm2": sum(
            cargas.valor_componente_superficial(componente, bib)[1] for componente in base
        ),
    }
    return calcular(datos, biblioteca=bib, tabla_viguetas=tabla_viguetas)



# ---------------------------------------------------------------------------
# La memoria de cálculo (el mismo texto de siempre)
# ---------------------------------------------------------------------------
def generar_memoria(resultado: ResultadoLosa, biblioteca: dict | None = None) -> str:
    """Arma el texto de la memoria, idéntico al que generaba L00."""
    bib = biblioteca if biblioteca is not None else materiales.cargar()
    r = resultado
    m: list[str] = []

    m.append("MEMORIA DE CÁLCULO - LOSA ALIVIANADA\n")
    m.append("=================================\n\n")

    m.append(f"Losa: {r.nombre}\n")
    if r.carga_id:
        m.append(f"Carga del proyecto: {r.carga_id}\n")
    m.append(f"Luz libre: {r.luz_libre_m:.2f} m\n")
    m.append(f"Luz de cálculo (con tolerancia): {r.luz_calculo_m:.2f} m\n\n")

    subtotal_D2 = 0.0
    m.append("Materiales seleccionados:\n")
    for sel in r.seleccion:
        categoria = sel["categoria"]
        nombre = sel["nombre"]
        material = _buscar_material(bib, categoria, nombre)
        if material is None:
            continue
        espesor = sel.get("espesor_m")
        carga = (
            float(material["valor"]) * float(espesor)
            if espesor is not None and material.get("tipo") == "volumetrico"
            else carga_de_material(material)
        )
        subtotal_D2 += carga
        m.append(f" - {categoria}: {nombre} → {carga:.3f} kN/m²\n")
    m.append(f"Subtotal cargas accesorias (D2): {subtotal_D2:.3f} kN/m²\n")
    m.append(f"Total cargas permanentes (D1+D2): {r.D:.3f} kN/m²\n\n")

    m.append("Cargas consideradas:\n")
    m.append(f" - D1 (peso propio): {r.D1:.3f} kN/m²\n")
    m.append(" - D1 se toma del forjado base configurado; no se recalibra automáticamente con la vigueta elegida.\n")
    m.append(f" - D2 (accesoria): {r.D2:.3f} kN/m²\n")
    m.append(f" - D total: {r.D:.3f} kN/m²\n")
    m.append(f" - Sobrecarga seleccionada: {r.uso_sobrecarga} → {r.L:.3f} kN/m²\n")
    m.append(f" - Combinación crítica: {r.combo_critico} → {r.Q:.3f} kN/m²\n\n")

    m.append("Selección de vigueta:\n")
    m.append(f" - Serie de vigueta: {r.serie}\n")
    m.append(f" - Tipo de vigueta seleccionada: {r.tipo_vigueta}\n")
    m.append(f" - Combinación bovedilla/losa: {r.combo}\n")
    m.append(f" - Momento requerido (sin φ): {r.momento_req_kgm:.2f} kg·m\n")
    m.append(f" - Momento requerido para búsqueda (con φ aplicado): {r.momento_busqueda_kgm:.2f} kg·m\n")
    m.append(f" - Momento disponible de la vigueta seleccionada: {r.momento_disponible_kgm:.2f} kg·m\n")
    m.append(f" - Momento reducido aplicado φ={r.phi}: {r.momento_reducido_kgm:.2f} kg·m ({r.momento_reducido_kNm:.2f} kN·m)\n")
    m.append(" - Referencia de capacidad: Tabla 4 – Tensolite\n\n")

    m.append("Cómputo de materiales:\n")
    m.append(f" - Área de losa: {r.area_losa_m2:.2f} m²\n")
    m.append(f" - Viguetas: {r.n_viguetas} unidades de {r.luz_calculo_m:.2f} m\n")
    m.append(f" - Bovedillas: {r.n_bloques} unidades EPS n°{r.bovedilla_altura}\n")
    m.append(f" - Hormigón: {r.vol_hormigon_m3:.3f} m³ (capa de compresión)\n\n")

    malla = r.malla
    m.append("Malla electrosoldada SIMA Q-131 / R-131:\n")
    m.append(f" - Paño comercial: {malla['tamanio_pano']} m (solape {malla['solape_m']:.2f} m)\n")
    m.append(f" - Área losa: {malla['area_losa_m2']:.2f} m²\n")
    m.append(f" - Área efectiva de paño: {malla['area_pano_efectiva_m2']:.2f} m²\n")
    m.append(f" - Cantidad adoptada: {malla['cantidad_mallas']} unidad(es)\n")
    if malla.get("nota"):
        m.append(f" - Nota: {malla['nota']}\n")
    m.append("\n")

    nervios = r.nervios
    m.append("Nervios de refuerzo:\n")
    m.append(f" - Tramo máximo: {TRAMO_MAX_NERVIO_M:.2f} m\n")
    m.append(f" - Nervios requeridos: {nervios['n_nervios']} unidades\n")
    m.append(f" - Especificación por nervio: {nervios['barras_por_nervio']} barras Ø{nervios['diametro_mm']} de {nervios['largo_por_nervio_m']:.2f} m\n")
    m.append(f" - Total barras: {nervios['barras_total']} unidades\n")
    m.append(f" - Longitud total de barras: {nervios['longitud_total_barras_m']:.2f} m\n\n")

    m.append("=================================\n")
    m.append("Fin de la memoria de cálculo\n")
    return "".join(m)



# ---------------------------------------------------------------------------
# Guardar (memoria .txt y fila del cómputo .csv)
# ---------------------------------------------------------------------------
COLUMNAS_COMPUTO = (
    "ID_Losa", "Luz_libre_m", "Luz_calculo_m",
    "Serie", "Tipo_vigueta",
    "N_viguetas", "N_bloques", "Vol_hormigon_m3",
    "Momento_req_kgm", "Momento_reducido_kNm", "Area_losa_m2",
    "Malla_tamanio", "Malla_area_efectiva_m2", "Malla_cantidad",
    "Nervios_n", "Barras_por_nervio", "Diametro_mm",
    "Largo_por_nervio_m", "Barras_total", "Longitud_total_barras_m",
)


def fila_computo(resultado: ResultadoLosa) -> list:
    """La fila del cómputo, en el mismo orden de columnas de siempre."""
    r = resultado
    return [
        r.nombre, r.luz_libre_m, r.luz_calculo_m, r.serie, r.tipo_vigueta,
        r.n_viguetas, r.n_bloques, r.vol_hormigon_m3,
        r.momento_req_kgm, r.momento_reducido_kNm, r.area_losa_m2,
        r.malla["tamanio_pano"], r.malla["area_pano_efectiva_m2"], r.malla["cantidad_mallas"],
        r.nervios["n_nervios"], r.nervios["barras_por_nervio"], r.nervios["diametro_mm"],
        r.nervios["largo_por_nervio_m"], r.nervios["barras_total"], r.nervios["longitud_total_barras_m"],
    ]


def guardar_memoria(resultado: ResultadoLosa, archivo=None) -> Path:
    """Guarda la memoria en salidas/losas (nombre con la fecha, como siempre)."""
    if archivo is None:
        rutas.SAL_LOSAS.mkdir(parents=True, exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        archivo = rutas.SAL_LOSAS / f"memoria_losa_{resultado.nombre}_{timestamp}.txt"
    return rutas.guardar_texto(archivo, resultado.memoria)


def guardar_computo(resultado: ResultadoLosa, archivo=None) -> Path:
    """Agrega una fila al cómputo de losas (crea el archivo con encabezado si no existe)."""
    camino = Path(archivo) if archivo is not None else rutas.COMPUTO_LOSAS
    camino.parent.mkdir(parents=True, exist_ok=True)
    nuevo = not camino.exists()
    with camino.open("a", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        if nuevo:
            writer.writerow(COLUMNAS_COMPUTO)
        writer.writerow(fila_computo(resultado))
    return camino


def ruta_losa(nombre: str) -> Path:
    """La ruta del archivo propio de una losa (un archivo por losa)."""
    return rutas.SAL_LOSAS / f"{rutas.nombre_seguro(nombre)}.json"


def a_dict(resultado: ResultadoLosa) -> dict:
    """El resultado como diccionario para guardar (sin la memoria, que va en .txt)."""
    datos = asdict(resultado)
    datos.pop("memoria", None)
    return datos


def guardar_losa(resultado: ResultadoLosa, archivo=None) -> Path:
    """
    Guarda el resultado de UNA losa en su PROPIO archivo (la fuente de verdad).
    Así ninguna losa pisa a otra: cada una tiene la suya.
    """
    camino = Path(archivo) if archivo is not None else ruta_losa(resultado.nombre)
    return rutas.guardar_json(camino, a_dict(resultado))


def cargar_losa(nombre: str) -> dict:
    """Lee el archivo propio de una losa (diccionario vacío si no está calculada)."""
    return rutas.leer_json(ruta_losa(nombre), por_defecto={}) or {}


def listar_losas_calculadas() -> list[Path]:
    """Los archivos de las losas ya calculadas (salidas/losas/*.json)."""
    return rutas.listar(rutas.SAL_LOSAS, "*.json")


def guardar_resultado(resultado: ResultadoLosa) -> tuple[Path, Path, Path]:
    """Guarda la memoria (.txt), el cómputo (.csv) y el archivo propio (.json) de la losa."""
    return guardar_memoria(resultado), guardar_computo(resultado), guardar_losa(resultado)


# ---------------------------------------------------------------------------
# El conjunto: las reacciones de VARIAS losas (el "análisis de cargas" de losas)
# ---------------------------------------------------------------------------
def _clave_destino(destino: dict) -> str:
    return f"{destino.get('portico', '?')}|{destino.get('viga', '?')}"


def _obra() -> str:
    try:
        return cargas.datos_cargas().get("obra", "")
    except Exception:
        return ""


def _acumular_en_vigas(reacciones: list, por_viga: dict) -> None:
    """Suma las reacciones que declaran destino sobre cada viga (el reparto al pórtico)."""
    for reaccion in reacciones:
        destino = reaccion.get("apoya_en")
        if not destino:
            continue
        clave = _clave_destino(destino)
        acumulado = por_viga.setdefault(
            clave,
            {
                "portico": destino.get("portico"),
                "viga": destino.get("viga"),
                "D": 0.0,
                "L": 0.0,
                "W": 0.0,
            },
        )
        acumulado["D"] += reaccion["cargas"]["D"]
        acumulado["L"] += reaccion["cargas"]["L"]
        acumulado["W"] += reaccion["cargas"]["W"]


def armar_conjunto_desde_archivos() -> dict:
    """
    La VISTA: arma el conjunto de reacciones leyendo los archivos individuales
    (`salidas/losas/*.json`). No es la fuente: se puede rehacer cuando quieras.
    """
    elementos: dict = {}
    por_viga: dict = {}
    for archivo in listar_losas_calculadas():
        datos = rutas.leer_json(archivo, por_defecto={}) or {}
        nombre = datos.get("nombre") or archivo.stem
        reacciones = datos.get("reacciones", [])
        elementos[nombre] = {
            "tipo": "losa_alivianada",
            "direccion": "unidireccional",
            "luz_calculo_m": datos.get("luz_calculo_m"),
            "ancho_m": datos.get("ancho_losa_m"),
            "reacciones": reacciones,
        }
        _acumular_en_vigas(reacciones, por_viga)

    for acumulado in por_viga.values():
        acumulado["combinaciones"] = cargas.combinaciones_de_componentes(
            acumulado["D"], acumulado["L"], acumulado["W"]
        )

    return {
        "obra": _obra(),
        "combinaciones": list(cargas.combinaciones_de_componentes(0.0, 0.0, 0.0)),
        "elementos": elementos,
        "por_viga": por_viga,
    }


def calcular_conjunto(nombres=None, biblioteca=None, tabla_viguetas=None) -> dict:
    """
    Calcula VARIAS losas, guarda el archivo propio de cada una (la fuente) y
    devuelve el conjunto armado desde esos archivos (la vista).
    """
    bib = biblioteca if biblioteca is not None else materiales.cargar()
    viguetas = tabla_viguetas if tabla_viguetas is not None else cargar_viguetas()
    losas = cargar_losas()
    elegidas = list(nombres) if nombres else list(losas)
    for nombre in elegidas:
        datos = {**losas[nombre], "nombre": losas[nombre].get("nombre", nombre)}
        resultado = calcular(datos, biblioteca=bib, tabla_viguetas=viguetas)
        guardar_losa(resultado)
    return armar_conjunto_desde_archivos()


def guardar_reacciones(conjunto: dict, archivo=None) -> Path:
    """Guarda la VISTA del conjunto en salidas/reacciones/_conjunto.json (se regenera)."""
    camino = Path(archivo) if archivo is not None else rutas.SAL_REACCIONES / "_conjunto.json"
    return rutas.guardar_json(camino, conjunto)


def resumen_reacciones(conjunto: dict) -> str:
    lineas = ["REACCIONES DE LOSAS (D / L / W en kN/m)", "=" * 78]
    for nombre, elemento in conjunto["elementos"].items():
        lineas.append(
            f"{nombre}  ({elemento['tipo']}, {elemento['direccion']}, "
            f"L={elemento['luz_calculo_m']:.2f} m)"
        )
        for reaccion in elemento["reacciones"]:
            c = reaccion["cargas"]
            destino = reaccion.get("apoya_en")
            donde = f"  → {destino['portico']} {destino['viga']}" if destino else ""
            lineas.append(
                f"    apoyo {reaccion['apoyo']:<4} D={c['D']:7.2f}  L={c['L']:7.2f}  "
                f"W={c['W']:6.2f}{donde}"
            )
    lineas.append("")
    if conjunto["por_viga"]:
        lineas.append("Reparto por viga (suma de reacciones):")
        for clave, acum in conjunto["por_viga"].items():
            lineas.append(
                f"    {clave:<28} D={acum['D']:7.2f}  L={acum['L']:7.2f}  W={acum['W']:6.2f}"
            )
    else:
        lineas.append("(ninguna losa declara 'apoya_en' todavía: no hay reparto por viga)")
    return "\n".join(lineas)


# ---------------------------------------------------------------------------
# Listado y consola
# ---------------------------------------------------------------------------
def lista_de_losas() -> str:
    losas = cargar_losas()
    lineas = ["LOSAS ALIVIANADAS (datos/losas.json)", "=" * 78]
    if not losas:
        lineas.append("(no hay losas definidas todavía)")
    for nombre, datos in losas.items():
        detalles = []
        if datos.get("luz_libre_m") is not None:
            detalles.append(f"L = {datos['luz_libre_m']} m")
        if datos.get("ancho_losa_m") is not None:
            detalles.append(f"b = {datos['ancho_losa_m']} m")
        if datos.get("sobrecarga"):
            detalles.append(f"uso = {datos['sobrecarga']}")
        lineas.append(f" - {nombre:<22} {' · '.join(detalles)}")
    lineas.append("")
    lineas.append('Para calcular una:   py -m calc.losas "<nombre>"')
    return "\n".join(lineas)


def resumen_texto(resultado: ResultadoLosa) -> str:
    return (
        f"{resultado.nombre}: serie {resultado.serie} · {resultado.tipo_vigueta} "
        f"({resultado.combo}) · L={resultado.luz_calculo_m:.2f} m · "
        f"Q={resultado.Q:.2f} kN/m² ({resultado.combo_critico}) · "
        f"{resultado.n_viguetas} viguetas"
    )


def _main(argv: list[str]) -> int:
    opciones = [a for a in argv if a.startswith("-")]
    nombres = [a for a in argv if not a.startswith("-")]

    if "--lista" in opciones or "-l" in opciones:
        print(lista_de_losas())
        return 0

    if "--conjunto" in opciones:
        conjunto = calcular_conjunto(nombres or None)
        print(resumen_reacciones(conjunto))
        archivo = guardar_reacciones(conjunto)
        print(f"\nReacciones guardadas en: {archivo}")
        return 0

    losas = cargar_losas()
    elegidas = nombres if nombres else list(losas)

    if not elegidas:
        print(lista_de_losas())
        return 0

    faltantes = [n for n in elegidas if n not in losas]
    if faltantes:
        print(f"No existe(n) en datos/losas.json: {', '.join(faltantes)}")
        print("Ver: py -m calc.losas --lista")
        return 1

    for nombre in elegidas:
        datos = {**losas[nombre], "nombre": losas[nombre].get("nombre", nombre)}
        resultado = calcular(datos)
        print(resultado.memoria)
        print("Reacciones de apoyo (D / L / W en kN/m):")
        for reaccion in resultado.reacciones:
            c = reaccion["cargas"]
            print(f"  apoyo {reaccion['apoyo']:<4} D={c['D']:7.2f}  L={c['L']:7.2f}  W={c['W']:6.2f}")
        if "--guardar" in opciones:
            memoria, _computo, json_losa = guardar_resultado(resultado)
            print(f"\nGuardado en: {memoria}")
            print(f"Archivo propio de la losa: {json_losa}")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main(sys.argv[1:]))


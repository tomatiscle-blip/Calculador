"""
calc.cargas — De los elementos que reciben carga a las cargas del proyecto.

Es la ÚNICA cuenta de cargas del programa. Toma los datos de:

    datos/cargas.json      los elementos (losa, muro, techo, encadenado) y el viento
    datos/materiales.json  la biblioteca (pesos específicos, cargas, sobrecargas)

y devuelve las cargas permanentes (D), las sobrecargas (L), el viento (W), las
combinaciones de CIRSOC y el informe de texto que se guarda en
`salidas/analisis_cargas/`.

La misma cuenta sirve para el pórtico, para la losa, para el muro y para el cimiento:
por eso un elemento se puede calcular SOLO, sin armar un pórtico:

    py -m calc.cargas                          todo el conjunto (lo que hacía 00_…)
    py -m calc.cargas "Losa Alivianada L0-1"   ese elemento solo
    py -m calc.cargas --lista                  qué elementos hay y cuáles están activos
"""

from __future__ import annotations

import sys
from datetime import datetime
from math import cos, radians, tan
from pathlib import Path

from . import materiales, rutas

# Orden en que aparecen los elementos en el informe (el mismo de siempre:
# cubiertas, losas/forjados, muros, encadenados y al final el viento general)
ORDEN_TIPOS = ("cubierta", "losa", "muro", "encadenado")


# ---------------------------------------------------------------------------
# Funciones de conversión
# ---------------------------------------------------------------------------
def carga_superficial(gamma: float, espesor: float) -> float:
    """kN/m3 * m = kN/m2"""
    return gamma * espesor


def carga_lineal_muro(gamma: float, espesor: float, altura: float) -> float:
    """kN/m3 * m * m = kN/m"""
    return gamma * espesor * altura


def carga_lineal_desde_superficie(q_superficial: float, ancho_tributario: float) -> float:
    """kN/m2 * m = kN/m"""
    return q_superficial * ancho_tributario


def direccion_portico(datos_portico: dict) -> str:
    """Dirección longitudinal de planta; estructuras anteriores se consideran X."""
    referencia = datos_portico.get("referencia_planta", {}) or {}
    direccion = str(referencia.get("direccion", "x")).strip().lower()
    return direccion if direccion in ("x", "y") else "x"


def ancho_tributario_entre_porticos(estructura: dict, nombre_portico: str) -> float | None:
    """Suma medias separaciones entre pórticos paralelos del mismo plano."""
    portico_objetivo = estructura.get(nombre_portico)
    if not isinstance(portico_objetivo, dict):
        return None
    direccion = direccion_portico(portico_objetivo)
    posiciones = []
    cantidad_porticos = 0
    for nombre, datos in estructura.items():
        if not isinstance(datos, dict) or direccion_portico(datos) != direccion:
            continue
        cantidad_porticos += 1
        try:
            posicion = float(datos["posicion_planta_m"])
        except (KeyError, TypeError, ValueError):
            continue
        posiciones.append((posicion, str(nombre)))
    if (
        not posiciones
        or len(posiciones) != cantidad_porticos
        or len({pos for pos, _ in posiciones}) != len(posiciones)
    ):
        return None
    posiciones.sort()
    indice = next(
        (i for i, (_, nombre) in enumerate(posiciones) if nombre == nombre_portico),
        None,
    )
    if indice is None:
        return None

    posicion = posiciones[indice][0]
    ancho = 0.0
    if indice > 0:
        ancho += (posicion - posiciones[indice - 1][0]) / 2
    if indice + 1 < len(posiciones):
        ancho += (posiciones[indice + 1][0] - posicion) / 2
    return ancho if ancho > 0 else None


def presion_viento(V: float, rho: float = 1.25) -> float:
    """Presión dinámica del viento en kN/m2 (CIRSOC 102)."""
    return 0.5 * rho * (V**2) / 1000


def carga_lineal_viento(V: float, ancho_tributario: float, Cd: float = 1.3, rho: float = 1.25) -> float:
    """Carga lineal de viento (kN/m)."""
    return presion_viento(V, rho) * Cd * ancho_tributario


# ---------------------------------------------------------------------------
# El análisis: una lista de cargas (D / L / W) y su informe
# ---------------------------------------------------------------------------
class AnalisisCargas:
    def __init__(self, nombre: str):
        self.nombre = nombre
        self.id_proyecto: str | None = None
        self.items: list[dict] = []

    def agregar(self, descripcion, tipo, valor, unidad, composicion=None):
        self.items.append({
            "descripcion": descripcion,
            "tipo": tipo,          # "D" permanente, "L" sobrecarga, "W" viento
            "valor": valor,
            "unidad": unidad,
            "composicion": composicion,
        })

    def resumen_texto(self) -> str:
        lineas = []
        D: dict[str, float] = {}
        L: dict[str, float] = {}
        W: dict[str, float] = {}

        lineas.append(f"ANALISIS DE CARGAS: {self.nombre}")
        if self.id_proyecto:
            lineas.append(f"PROYECTO: {self.id_proyecto}")
        lineas.append("-" * 60)

        for i in self.items:
            linea = f"{i['descripcion']:<35} {i['valor']:>7.2f} {i['unidad']} ({i['tipo']})"
            lineas.append(linea)

            if i["composicion"]:
                if isinstance(i["composicion"], list):
                    for linea_comp in i["composicion"]:
                        lineas.append(f"  {linea_comp}")
                else:
                    lineas.append(f"  {i['composicion']}")
            if i["tipo"] == "D":
                D[i["unidad"]] = D.get(i["unidad"], 0) + i["valor"]
            elif i["tipo"] == "L":
                L[i["unidad"]] = L.get(i["unidad"], 0) + i["valor"]
            elif i["tipo"] == "W":
                W[i["unidad"]] = W.get(i["unidad"], 0) + i["valor"]

            lineas.append("")

        lineas.append("-" * 60)

        for u, v in D.items():
            lineas.append(f"D total = {v:.2f} {u}")
        for u, v in L.items():
            lineas.append(f"L total = {v:.2f} {u}")
        for u, v in W.items():
            lineas.append(f"W total = {v:.2f} {u}")

        return "\n".join(lineas)


# ---------------------------------------------------------------------------
# Combinaciones de CIRSOC (una sola cuenta para todo el programa)
# ---------------------------------------------------------------------------
# Factores por combinación: (fD, fL, fW). El motor (calc/portico.py) reusa esta
# MISMA tabla para armar las combinaciones de carga del pórtico: así el análisis
# de cargas y las solicitaciones no pueden dar distinto.
FACTORES_COMBINACIONES: dict[str, tuple[float, float, float]] = {
    "1.4D": (1.4, 0.0, 0.0),
    "1.2D+1.6L": (1.2, 1.6, 0.0),
    "1.2D+0.5L+1.6W": (1.2, 0.5, 1.6),
    "0.9D+1.6W": (0.9, 0.0, 1.6),
}


def combinaciones_de_componentes(d_total: float, l_total: float, w_total: float) -> dict[str, float]:
    """
    Las 4 combinaciones de CIRSOC a partir de D, L y W ya sumados.
    Es la ÚNICA cuenta de combinaciones del programa: la usan el pórtico, la losa
    y las reacciones, así no pueden dar distinto.
    """
    return {
        nombre: fD * d_total + fL * l_total + fW * w_total
        for nombre, (fD, fL, fW) in FACTORES_COMBINACIONES.items()
    }


def combinaciones(items: list[dict]) -> dict[str, float]:
    """
    Combinaciones de carga, en kN/m (o kN/m2, según lo que reciba).
    La usa el pórtico y también la losa: así no pueden dar distinto.
    """
    d_total = sum(i["valor"] for i in items if i["tipo"] == "D")
    l_total = sum(i["valor"] for i in items if i["tipo"] == "L")
    w_total = sum(i["valor"] for i in items if i["tipo"] == "W")

    return combinaciones_de_componentes(d_total, l_total, w_total)


# ---------------------------------------------------------------------------
# Los elementos que reciben carga (datos/cargas.json)
# ---------------------------------------------------------------------------
def datos_cargas() -> dict:
    """Lee datos/cargas.json. Si falta, avisa claro en vez de calcular con nada."""
    datos = rutas.leer_json(rutas.CARGAS, por_defecto=None)
    if not datos:
        raise FileNotFoundError(f"No se pudieron leer los elementos: {rutas.CARGAS}")
    return datos


def elementos(solo_activos: bool = True, datos: dict | None = None) -> dict[str, dict]:
    """
    Los elementos del proyecto en el orden del informe (cubiertas, losas, muros,
    encadenados). Con `solo_activos=False` los devuelve todos, para poder elegir uno.
    """
    fuente = (datos if datos is not None else datos_cargas()).get("elementos", {})
    elegidos: dict[str, dict] = {}
    for tipo in ORDEN_TIPOS:
        for nombre, elemento in fuente.items():
            if elemento.get("tipo") != tipo:
                continue
            if solo_activos and not elemento.get("activo"):
                continue
            elegidos[nombre] = elemento
    return elegidos


def _nombre_componente(componente: dict, material: dict, espesor) -> str:
    plantilla = componente.get("nombre")
    if plantilla:
        if "{" in plantilla:
            metros = float(espesor or 0.0)
            return plantilla.format(espesor=metros, espesor_cm=metros * 100)
        return plantilla
    if espesor is not None and material.get("tipo") == "volumetrico":
        return f"{material['nombre']} {float(espesor):.2f} m"
    return material["nombre"]


def _componente_valor(componente: dict, biblioteca: dict) -> tuple[str, float]:
    """
    Convierte un componente del JSON en (nombre, carga superficial en kN/m2).

    - material 'superficial': el valor es directo.
    - material 'volumetrico' con `espesor_m`: gamma * espesor.
    - `peso_propio_override_kNm2`: ajuste por obra de un forjado superficial.
    - `nombre` puede traer la plantilla `{espesor}` (metros) o `{espesor_cm}`.
    """
    material = materiales.buscar(componente["grupo"], componente["clave"], biblioteca)
    espesor = componente.get("espesor_m")
    valor_override = componente.get("peso_propio_override_kNm2")
    if valor_override is not None:
        valor = float(valor_override)
        if valor <= 0:
            raise ValueError("El peso propio ajustado del forjado debe ser mayor que cero.")
    else:
        valor = float(material["valor"])
    if valor_override is None and espesor is not None and material.get("tipo") == "volumetrico":
        valor = carga_superficial(valor, float(espesor))
    return _nombre_componente(componente, material, espesor), valor


def valor_componente_superficial(
    componente: dict, biblioteca: dict | None = None
) -> tuple[str, float]:
    """Devuelve la descripción y el aporte kN/m² de una capa de losa/cubierta."""
    bib = biblioteca if biblioteca is not None else materiales.cargar()
    return _componente_valor(componente, bib)


def _componentes(elemento: dict, biblioteca: dict) -> list[tuple[str, float]]:
    return [_componente_valor(c, biblioteca) for c in elemento.get("componentes", [])]


def valores_superficiales(elemento: dict, biblioteca: dict | None = None) -> dict:
    """Devuelve D y L de una losa/cubierta por unidad de superficie (kN/m²)."""
    if elemento.get("tipo") not in ("losa", "cubierta"):
        raise ValueError("Los valores superficiales solo aplican a losas y cubiertas.")
    bib = biblioteca if biblioteca is not None else materiales.cargar()
    componentes = _componentes(elemento, bib)
    sobrecargas = materiales.sobrecargas(bib)
    clave_sobrecarga = elemento.get("sobrecarga")
    if clave_sobrecarga not in sobrecargas:
        raise KeyError(f"No existe la sobrecarga '{clave_sobrecarga}' en materiales.json")
    return {
        "D_kNm2": sum(valor for _nombre, valor in componentes),
        "L_kNm2": float(sobrecargas[clave_sobrecarga]),
        "componentes": componentes,
    }


# ---------------------------------------------------------------------------
# Cargas de un elemento, según su tipo
# ---------------------------------------------------------------------------
def items_de_elemento(nombre: str, elemento: dict, biblioteca: dict) -> list[dict]:
    """Devuelve las cargas (D / L / W) de un elemento, según su tipo."""
    tipo = elemento.get("tipo")
    if tipo in ("losa", "cubierta"):
        return _items_horizontal(nombre, elemento, biblioteca)
    if tipo == "muro":
        return _items_muro(nombre, elemento, biblioteca)
    if tipo == "encadenado":
        return _items_encadenado(nombre, elemento, biblioteca)
    raise ValueError(f"El elemento '{nombre}' tiene un tipo que no sé calcular: {tipo!r}")


def _items_horizontal(nombre: str, elemento: dict, biblioteca: dict) -> list[dict]:
    """Losa o cubierta: la carga por m2 pasa a carga por metro con el ancho tributario."""
    valores = valores_superficiales(elemento, biblioteca)
    componentes = valores["componentes"]
    q_d_total = valores["D_kNm2"]
    b = float(elemento["ancho_tributario_m"])

    items = [{
        "descripcion": f"{nombre} – cargas permanentes",
        "tipo": "D",
        "valor": carga_lineal_desde_superficie(q_d_total, b),
        "unidad": "kN/m",
        "composicion": [f"{n:<35} {q:.2f} kN/m2" for n, q in componentes]
                       + [f"q total = {q_d_total:.2f} kN/m2 · b = {b:.2f} m"],
    }]

    q_l = valores["L_kNm2"]
    items.append({
        "descripcion": f"{nombre} – sobrecarga",
        "tipo": "L",
        "valor": carga_lineal_desde_superficie(q_l, b),
        "unidad": "kN/m",
        "composicion": [f"q = {q_l:.2f} kN/m2 · b = {b:.2f} m"],
    })

    if elemento.get("viento_activo"):
        items.append(_item_succion_cubierta(nombre, elemento, biblioteca))
    return items


def _item_succion_cubierta(nombre: str, elemento: dict, biblioteca: dict) -> dict:
    """Succión del viento sobre una cubierta inclinada (CIRSOC 102)."""
    pendiente = float(elemento.get("pendiente_grados", elemento.get("pendiente", 0)) or 0)
    b = float(elemento["ancho_tributario_m"])
    if elemento.get("altura_vertical_m") is None:
        altura_vertical = b * tan(radians(pendiente))
    else:
        altura_vertical = float(elemento["altura_vertical_m"])
    altura_tributaria = altura_vertical * cos(radians(pendiente))

    cd_succion = 0.9
    w_succion = carga_lineal_viento(materiales.viento(biblioteca)["V"], altura_tributaria, cd_succion)
    return {
        "descripcion": f"Viento – succión cubierta {nombre}",
        "tipo": "W",
        "valor": w_succion,
        "unidad": "kN/m",
        "composicion": f"h = {altura_tributaria:.2f} m · Cd_succion = {cd_succion}",
    }


def _items_muro(nombre: str, elemento: dict, biblioteca: dict) -> list[dict]:
    material = materiales.pesos_especificos(biblioteca)["Mamposteria"][elemento["clave"]]
    gamma = material["gamma"]
    espesor = float(elemento["espesor_m"])
    altura = float(elemento["altura_m"])
    return [{
        "descripcion": f"Muro de {material['nombre']} (e={int(espesor*100)} cm · h={altura:.2f} m)",
        "tipo": "D",
        "valor": carga_lineal_muro(gamma, espesor, altura),
        "unidad": "kN/m",
        "composicion": f"γ={gamma} kN/m3 · e={espesor} m · h={altura} m",
    }]


def _items_encadenado(nombre: str, elemento: dict, biblioteca: dict) -> list[dict]:
    hormigon = materiales.pesos_especificos(biblioteca)["Hormigon"]["armado"]
    gamma = hormigon["gamma"]
    base = float(elemento["base_m"])
    altura = float(elemento["altura_m"])
    cantidad = int(elemento.get("cantidad", 1))
    return [{
        "descripcion": f"Encadenado de {hormigon['nombre']} ({int(base*100)}x{int(altura*100)} cm) × {cantidad}",
        "tipo": "D",
        "valor": gamma * base * altura * cantidad,
        "unidad": "kN/m",
        "composicion": f"γ={gamma} kN/m3 · b={base} m · h={altura} m · cantidad={cantidad}",
    }]


def items_viento_general(datos: dict, biblioteca: dict) -> list[dict]:
    """
    El viento que empuja al pórtico (carga horizontal). No es de un elemento:
    es del conjunto, por eso se suma aparte.
    """
    configuracion = datos.get("viento", {})
    if not configuracion.get("activo"):
        return []
    viento = materiales.viento(biblioteca)
    ancho = float(configuracion.get("ancho_tributario_m", 0.0))
    w_viento = carga_lineal_viento(viento["V"], ancho, viento["Cd"], viento["rho"])
    return [{
        "descripcion": f"Viento diseño {viento['ciudad']} (V={viento['V']:.0f} m/s) – general",
        "tipo": "W",
        "valor": w_viento,
        "unidad": "kN/m",
        "composicion": f"b = {ancho:.2f} m · Cd = {viento['Cd']}",
    }]


def _agregar_items(analisis: AnalisisCargas, items: list[dict]) -> None:
    for item in items:
        analisis.agregar(
            item["descripcion"], item["tipo"], item["valor"], item["unidad"], item["composicion"]
        )


# ---------------------------------------------------------------------------
# Las dos formas de pedir las cargas
# ---------------------------------------------------------------------------
def cargas_del_conjunto(biblioteca: dict | None = None) -> AnalisisCargas:
    """Todos los elementos activos + el viento general (lo que hacía 00_Analisis_cargas)."""
    bib = biblioteca if biblioteca is not None else materiales.cargar()
    datos = datos_cargas()
    analisis = AnalisisCargas(datos.get("obra", "Obra"))
    analisis.id_proyecto = datos.get("id_proyecto")
    for nombre, elemento in elementos(solo_activos=True, datos=datos).items():
        identificador = elemento.get("id")
        etiqueta = f"{identificador} · {nombre}" if identificador else nombre
        _agregar_items(analisis, items_de_elemento(etiqueta, elemento, bib))
    _agregar_items(analisis, items_viento_general(datos, bib))
    return analisis


def carga_de_elemento(nombre: str, biblioteca: dict | None = None) -> AnalisisCargas:
    """
    UN elemento solo (una losa, un muro, un techo): la versatilidad del programa.
    No arrastra los demás elementos ni el viento general.
    """
    bib = biblioteca if biblioteca is not None else materiales.cargar()
    disponibles = elementos(solo_activos=False)
    if nombre not in disponibles:
        raise KeyError(f"No existe el elemento '{nombre}'. Ver: py -m calc.cargas --lista")
    datos = datos_cargas()
    elemento = disponibles[nombre]
    analisis = AnalisisCargas(datos.get("obra", nombre))
    analisis.id_proyecto = datos.get("id_proyecto")
    identificador = elemento.get("id")
    etiqueta = f"{identificador} · {nombre}" if identificador else nombre
    _agregar_items(analisis, items_de_elemento(etiqueta, elemento, bib))
    return analisis


# ---------------------------------------------------------------------------
# Informe de texto (el que se guarda en salidas/analisis_cargas)
# ---------------------------------------------------------------------------
def informe_completo(analisis: AnalisisCargas) -> str:
    """Resumen + combinaciones de CIRSOC, con el mismo texto de siempre."""
    texto = analisis.resumen_texto() + "\n\nCOMBINACIONES CIRSOC:\n"
    for nombre, valor in combinaciones(analisis.items).items():
        texto += f"{nombre}: {valor:.2f} kN/m\n"
    return texto


def informe_aplicaciones(datos: dict | None = None, biblioteca: dict | None = None) -> str:
    """Informe de cargas realmente aplicadas, agrupadas por tramo e intervalo."""
    datos = datos if datos is not None else datos_cargas()
    bib = biblioteca if biblioteca is not None else materiales.cargar()
    fuentes = {
        elemento.get("id"): (nombre, elemento)
        for nombre, elemento in datos.get("elementos", {}).items()
        if elemento.get("id")
    }
    estructura = rutas.cargar_estructura()
    por_tramo: dict[tuple[str, str], list[dict]] = {}
    detalle: list[tuple[dict, str, dict, list[dict]]] = []

    for aplicacion in datos.get("aplicaciones", []):
        if not aplicacion.get("activa", True):
            continue
        fuente_info = fuentes.get(aplicacion.get("carga_id"))
        if not fuente_info:
            continue
        nombre, fuente = fuente_info
        if not fuente.get("activo", True):
            continue
        portico = aplicacion.get("portico", "")
        tramo_id = aplicacion.get("tramo_id", "")
        tramo = next((t for v in estructura.get(portico, {}).get("vigas", {}).values()
                      for t in v.get("tramos", []) if t.get("id") == tramo_id), None)
        if not tramo:
            continue
        longitud = float(tramo.get("longitud_m", 0.0))
        x0 = float(aplicacion.get("x_inicio_m", 0.0))
        x1 = float(aplicacion.get("x_fin_m", longitud))
        if longitud <= 0 or x0 < 0 or x1 <= x0 or x1 > longitud + 1e-8:
            continue

        elemento = dict(fuente)
        ancho = None
        if elemento.get("tipo") in ("losa", "cubierta"):
            modo = aplicacion.get("ancho_modo", "manual")
            luz = float(elemento.get("luz_transversal_m", 0.0))
            ancho = (luz / 2 if modo == "media_luz" else luz
                     if modo == "luz_completa" else aplicacion.get("ancho_tributario_m"))
            if ancho is None or float(ancho) <= 0:
                continue
            elemento["ancho_tributario_m"] = float(ancho)

        items = items_de_elemento(f"{fuente.get('id')} · {nombre}", elemento, bib)
        registro = {
            "aplicacion": aplicacion, "nombre": nombre, "fuente": fuente,
            "tramo": tramo, "longitud": longitud, "x0": x0, "x1": x1,
            "ancho": ancho, "items": items,
        }
        detalle.append((aplicacion, nombre, elemento, items))
        por_tramo.setdefault((portico, tramo_id), []).append(registro)

    lineas = [
        f"ANÁLISIS DE CARGAS APLICADAS: {datos.get('obra', 'Obra')}",
        f"PROYECTO: {datos.get('id_proyecto', 'sin ID')}",
        "-" * 72,
        "Las combinaciones se calculan por tramo e intervalo, con D, L y W una sola vez.",
        "Las cargas puntuales existentes de P00 se conservan y el viento horizontal general se informa aparte.",
        "",
        "DETALLE DE APLICACIONES:",
    ]
    if not detalle:
        lineas.append("No hay aplicaciones activas asociadas a tramos.")
    for app, nombre, elemento, items in detalle:
        ancho_txt = ""
        if elemento.get("tipo") in ("losa", "cubierta"):
            ancho_txt = (
                f" · luz completa={float(elemento.get('luz_transversal_m', 0.0)):.2f} m"
                f" · b tributario={float(elemento['ancho_tributario_m']):.2f} m"
            )
        lineas.append(
            f"{app.get('id', 'Aplicación')} · {app.get('carga_id')} {nombre} → "
            f"{app.get('portico')} / {app.get('tramo_id')} · "
            f"x={float(app.get('x_inicio_m', 0)):.2f}–{float(app.get('x_fin_m', 0)):.2f} m{ancho_txt}"
        )
        for item in items:
            lineas.append(f"  {item['tipo']}: {item['valor']:.2f} {item['unidad']} · {item['descripcion']}")
            composicion = item.get("composicion")
            if isinstance(composicion, list):
                lineas.extend(f"    {parte}" for parte in composicion)
            elif composicion:
                lineas.append(f"    {composicion}")

    lineas.extend(["", "DISTRIBUCIÓN Y COMBINACIONES POR TRAMO:"])
    for (portico, tramo_id), registros in sorted(por_tramo.items()):
        longitud = registros[0]["longitud"]
        tramo = registros[0]["tramo"]
        reemplaza = any(r["aplicacion"].get("modo_cargas_previas", "reemplazar") == "reemplazar" for r in registros)
        puntos = {0.0, longitud}
        for r in registros:
            puntos.update((r["x0"], r["x1"]))
        cargas_previas = tramo.get("cargas", {})
        if not reemplaza:
            d_previa = float(cargas_previas.get("D_total", 0.0))
            l_previa = float(cargas_previas.get("L_total", 0.0))
        else:
            d_previa = l_previa = 0.0
        lineas.append(f"\n{portico} · {tramo_id} (L={longitud:.2f} m) · cargas D/L de P00 {'reemplazadas' if reemplaza else 'sumadas'}")
        ordenados = sorted(puntos)
        for a, b in zip(ordenados, ordenados[1:]):
            medio = (a + b) / 2
            totales = {"D": d_previa, "L": l_previa, "W": 0.0}
            for r in registros:
                if r["x0"] - 1e-9 <= medio <= r["x1"] + 1e-9:
                    for item in r["items"]:
                        totales[item["tipo"]] = totales.get(item["tipo"], 0.0) + float(item["valor"])
            combos = combinaciones_de_componentes(totales["D"], totales["L"], totales["W"])
            lineas.append(
                f"  x={a:.2f}–{b:.2f} m: D={totales['D']:.2f}, L={totales['L']:.2f}, "
                f"W={totales['W']:.2f} kN/m"
            )
            lineas.extend(f"    {nombre_combo}: {valor:.2f} kN/m" for nombre_combo, valor in combos.items())

    viento = items_viento_general(datos, bib)
    if viento:
        lineas.extend(["", "VIENTO HORIZONTAL GENERAL (se aplica al pórtico):"])
        for item in viento:
            lineas.append(f"  {item['valor']:.2f} {item['unidad']} · {item['descripcion']}")
            if item.get("composicion"):
                lineas.append(f"    {item['composicion']}")
    return "\n".join(lineas) + "\n"


def guardar_informe(analisis: AnalisisCargas, texto: str | None = None) -> Path:
    """Guarda el informe en salidas/analisis_cargas, con el nombre y la fecha de siempre."""
    rutas.SAL_ANALISIS_CARGAS.mkdir(parents=True, exist_ok=True)
    fecha_hora = datetime.now().strftime("%d-%m-%Y_%H%M")
    archivo = rutas.SAL_ANALISIS_CARGAS / f"{analisis.nombre.replace(' ', '_')}_{fecha_hora}.txt"
    return rutas.guardar_texto(archivo, texto if texto is not None else informe_completo(analisis))


def lista_de_elementos() -> str:
    """Los elementos del proyecto (activos y no), para elegir uno."""
    lineas = ["ELEMENTOS QUE RECIBEN CARGA", "=" * 78]
    for nombre, elemento in elementos(solo_activos=False).items():
        estado = "activo" if elemento.get("activo") else "      "
        detalles = []
        if elemento.get("ancho_tributario_m"):
            detalles.append(f"b = {elemento['ancho_tributario_m']} m")
        if elemento.get("espesor_m"):
            detalles.append(f"e = {elemento['espesor_m']} m")
        if elemento.get("altura_m"):
            detalles.append(f"h = {elemento['altura_m']} m")
        if elemento.get("componentes"):
            detalles.append(f"{len(elemento['componentes'])} componentes")
        lineas.append(
            f"[{estado}] {elemento.get('tipo', '?'):<10} {nombre:<58} {' · '.join(detalles)}"
        )
    lineas.append("")
    lineas.append('Para calcular uno solo:  py -m calc.cargas "<nombre>"')
    return "\n".join(lineas)


def _main(argv: list[str]) -> int:
    opciones = [a for a in argv if a.startswith("-")]
    nombres = [a for a in argv if not a.startswith("-")]

    if "--lista" in opciones or "-l" in opciones:
        print(lista_de_elementos())
        return 0

    analisis = carga_de_elemento(nombres[0]) if nombres else cargas_del_conjunto()
    texto = informe_completo(analisis)
    print(texto)
    if "--guardar" in opciones:
        print(f"\nResumen guardado en: {guardar_informe(analisis, texto)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main(sys.argv[1:]))

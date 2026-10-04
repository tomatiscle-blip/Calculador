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
    - `nombre` puede traer la plantilla `{espesor}` (metros) o `{espesor_cm}`.
    """
    material = materiales.buscar(componente["grupo"], componente["clave"], biblioteca)
    espesor = componente.get("espesor_m")
    valor = float(material["valor"])
    if espesor is not None and material.get("tipo") == "volumetrico":
        valor = carga_superficial(valor, float(espesor))
    return _nombre_componente(componente, material, espesor), valor


def _componentes(elemento: dict, biblioteca: dict) -> list[tuple[str, float]]:
    return [_componente_valor(c, biblioteca) for c in elemento.get("componentes", [])]


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
    componentes = _componentes(elemento, biblioteca)
    q_d_total = sum(q for _, q in componentes)
    b = float(elemento["ancho_tributario_m"])

    items = [{
        "descripcion": f"{nombre} – cargas permanentes",
        "tipo": "D",
        "valor": carga_lineal_desde_superficie(q_d_total, b),
        "unidad": "kN/m",
        "composicion": [f"{n:<35} {q:.2f} kN/m2" for n, q in componentes]
                       + [f"q total = {q_d_total:.2f} kN/m2 · b = {b:.2f} m"],
    }]

    q_l = materiales.sobrecargas(biblioteca)[elemento["sobrecarga"]]
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
    for nombre, elemento in elementos(solo_activos=True, datos=datos).items():
        _agregar_items(analisis, items_de_elemento(nombre, elemento, bib))
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
    analisis = AnalisisCargas(nombre)
    _agregar_items(analisis, items_de_elemento(nombre, disponibles[nombre], bib))
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

"""
calc.materiales — La biblioteca del proyecto (datos/materiales.json).

Es el ÚNICO lugar de donde salen los pesos específicos, las cargas superficiales,
las sobrecargas y el viento. Todos los módulos leen de acá (cargas, losas, vigas…),
así un mismo material no puede valer dos cosas distintas en dos lugares.

En el JSON cada material tiene una "clave" (ej. `armado`, `porcelanato`). Este módulo
arma con esas claves los diccionarios que usa el cálculo, con el nombre que
corresponde a cada uso: `gamma` para los pesos específicos y `q` para las cargas
superficiales.

    py -m calc.materiales        lista los grupos y sus claves
"""

from __future__ import annotations

from . import rutas


def cargar() -> dict:
    """Lee datos/materiales.json. Si no está, avisa claro en vez de calcular con nada."""
    datos = rutas.leer_json(rutas.MATERIALES, por_defecto=None)
    if not datos:
        raise FileNotFoundError(f"No se pudo leer la biblioteca de materiales: {rutas.MATERIALES}")
    return datos


def tabla_por_clave(grupo: str, nombre_valor: str, biblioteca: dict | None = None) -> dict:
    """
    Convierte una lista del JSON en un diccionario por 'clave'.

    Ejemplo: `tabla_por_clave("Hormigon", "gamma")` devuelve
    {"armado": {"nombre": "Hormigón armado", "gamma": 25.0}, ...}

    Los materiales sin clave se saltean: son los que todavía no usa el cálculo.
    """
    datos = biblioteca if biblioteca is not None else cargar()
    tabla: dict[str, dict] = {}
    for elemento in datos.get(grupo, []):
        clave = elemento.get("clave")
        if not clave:
            continue
        tabla[clave] = {"nombre": elemento["nombre"], nombre_valor: float(elemento["valor"])}
    return tabla


def pesos_especificos(biblioteca: dict | None = None) -> dict:
    """Pesos específicos (kN/m3), agrupados y por clave."""
    return {
        "Hormigon": tabla_por_clave("Hormigon", "gamma", biblioteca),
        "Madera": tabla_por_clave("Madera", "gamma", biblioteca),
        "Mamposteria": tabla_por_clave("Mamposteria", "gamma", biblioteca),
        "Morteros": tabla_por_clave("Morteros_Revoques", "gamma", biblioteca),
    }


def sistemas(biblioteca: dict | None = None) -> dict:
    """Cargas superficiales directas (kN/m2), agrupadas y por clave."""
    return {
        "Pisos": tabla_por_clave("Pisos", "q", biblioteca),
        "Cubiertas": tabla_por_clave("Cubiertas", "q", biblioteca),
        "EstructuraCubierta": tabla_por_clave("EstructuraCubierta", "q", biblioteca),
        "Cielorrasos": tabla_por_clave("Cielorrasos", "q", biblioteca),
        "Forjados": tabla_por_clave("Forjados", "q", biblioteca),
    }


def sobrecargas(biblioteca: dict | None = None) -> dict:
    """Sobrecargas de uso (kN/m2) por clave: vivienda, oficina, terraza…"""
    datos = biblioteca if biblioteca is not None else cargar()
    return {
        elemento["clave"]: float(elemento["valor_kNm2"])
        for elemento in datos.get("Sobrecargas", [])
        if elemento.get("clave")
    }


def viento(biblioteca: dict | None = None) -> dict:
    """Datos de viento (CIRSOC 102): ciudad, V (m/s), rho (kg/m3) y Cd."""
    datos = biblioteca if biblioteca is not None else cargar()
    fuente = datos.get("Viento", {})
    return {
        "ciudad": fuente.get("ciudad", ""),
        "V": float(fuente.get("V_ms", 0.0)),
        "rho": float(fuente.get("rho_kg_m3", 1.25)),
        "Cd": float(fuente.get("Cd", 1.3)),
    }


def buscar(grupo: str, clave: str, biblioteca: dict | None = None) -> dict:
    """
    Devuelve el material tal como está escrito en el JSON (nombre, tipo, valor…).
    Sirve para saber si es 'superficial' (kN/m2) o 'volumetrico' (kN/m3).
    """
    datos = biblioteca if biblioteca is not None else cargar()
    for elemento in datos.get(grupo, []):
        if elemento.get("clave") == clave:
            return elemento
    raise KeyError(f"No existe la clave '{clave}' en el grupo '{grupo}' de {rutas.MATERIALES.name}")


def avisos() -> list[str]:
    """Lo que quedó anotado como 'para revisar' dentro de la biblioteca."""
    return list(cargar().get("_revisar", []))


def _main() -> int:
    datos = cargar()
    print(f"BIBLIOTECA DE MATERIALES — {rutas.MATERIALES}")
    print("=" * 78)
    for grupo, elementos in datos.items():
        if not isinstance(elementos, list):
            continue
        con_clave = [e for e in elementos if isinstance(e, dict) and e.get("clave")]
        print(f"{grupo:<22} {len(elementos):>3} materiales, {len(con_clave):>3} con clave")
        for elemento in con_clave:
            nombre = elemento.get("nombre") or elemento.get("uso")
            valor = elemento.get("valor", elemento.get("valor_kNm2"))
            print(f"      {elemento['clave']:<26} {nombre:<48} {valor} {elemento.get('unidad', '')}")
    for aviso in avisos():
        print(f"\n[REVISAR] {aviso}")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())

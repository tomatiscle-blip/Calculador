"""Valores generales de materiales y análisis propios de cada obra."""
from __future__ import annotations

from . import rutas


PARAMETROS_PREDETERMINADOS = {
    "MATERIAL": "Hormigón armado",
    "FC_MPA": 25,
    "FY_MPA": 420,
    "PESO_ESPECIFICO_HORMIGON_KN_M3": 25.0,
    "VERSION_MOTOR": "CIRSOC 201-25",
    "COLUMNA_LADO_X_CM": 25.0,
    "COLUMNA_LADO_Y_CM": 40.0,
    "COLUMNA_DIAMETRO_CM": 30.0,
    "VIGA_ANCHO_PREDIM_CM": 20.0,
}


def cargar() -> dict:
    """Combina los valores de la obra con los valores iniciales del programa."""
    guardados = rutas.leer_json(rutas.PARAMETROS_MOTOR, {}) or {}
    datos = {**PARAMETROS_PREDETERMINADOS, **guardados}
    if str(datos.get("VERSION_MOTOR", "")).strip() == "CIRSOC 201-2024":
        datos["VERSION_MOTOR"] = "CIRSOC 201-25"
    return datos


def fc_predeterminado() -> int:
    try:
        fc = int(cargar().get("FC_MPA", 25))
    except (TypeError, ValueError):
        fc = 25
    return fc if fc in (20, 25, 30) else 25


def peso_especifico_hormigon() -> float:
    try:
        valor = float(cargar().get("PESO_ESPECIFICO_HORMIGON_KN_M3", 25.0))
    except (TypeError, ValueError):
        valor = 25.0
    return valor if valor > 0 else 25.0


def modulo_elasticidad_hormigon_mpa() -> float:
    return float(modulo_elasticidad_hormigon()["Ec_MPa"])


def modulo_elasticidad_hormigon() -> dict:
    """Ec estimado según densidad y f'c, CIRSOC 201-25 art. 19.2.2.1."""
    import math
    fc = float(fc_predeterminado())
    wc = peso_especifico_hormigon()
    densidad_kg_m3 = wc * 1000.0 / 9.80665
    if 1400.0 <= densidad_kg_m3 <= 2500.0:
        ec = 0.043 * densidad_kg_m3**1.5 * math.sqrt(fc)
        formula = "Ec = 0,043·wc^1,5·sqrt(f'c) (MPa), wc en kg/m³"
        articulo = "CIRSOC 201-25, 19.2.2.1(a)"
    else:
        ec = 4700.0 * math.sqrt(fc)
        formula = "Ec = 4700·sqrt(f'c) (MPa), hormigón de peso normal"
        articulo = "CIRSOC 201-25, 19.2.2.1(b)"
    return {
        "Ec_MPa": ec, "formula": formula, "articulo": articulo,
        "fc_MPa": fc, "densidad_kg_m3": densidad_kg_m3,
    }

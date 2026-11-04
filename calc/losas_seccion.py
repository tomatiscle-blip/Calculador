"""Utilidades puras para geometría y peso propio de losas nervuradas."""

from __future__ import annotations

import math


def espesor_equivalente_hormigon_casetonada(seccion: dict) -> float:
    """Calcula el volumen de hormigón por m² como espesor equivalente."""
    if not isinstance(seccion, dict):
        raise ValueError("Definí la geometría de la sección casetonada.")
    claves = (
        "altura_total_m",
        "espesor_capa_m",
        "ancho_nervio_x_m",
        "ancho_nervio_y_m",
        "separacion_libre_x_m",
        "separacion_libre_y_m",
    )
    valores = {}
    for clave in claves:
        try:
            valor = float(seccion.get(clave, 0.0))
        except (TypeError, ValueError):
            raise ValueError(f"El dato de casetonada '{clave}' debe ser numérico.") from None
        if not math.isfinite(valor):
            raise ValueError(f"El dato de casetonada '{clave}' debe ser finito.")
        valores[clave] = valor
    h = valores["altura_total_m"]
    capa = valores["espesor_capa_m"]
    bx = valores["ancho_nervio_x_m"]
    by = valores["ancho_nervio_y_m"]
    sx = valores["separacion_libre_x_m"]
    sy = valores["separacion_libre_y_m"]
    if h <= 0 or capa <= 0 or capa >= h:
        raise ValueError("La casetonada necesita canto total y capa positivos, con capa menor que el canto.")
    if min(bx, by, sx, sy) <= 0:
        raise ValueError("Los anchos y las separaciones libres de los nervios deben ser positivos.")
    paso_x = sx + bx
    paso_y = sy + by
    fraccion_nervios = bx / paso_x + by / paso_y - (bx * by) / (paso_x * paso_y)
    if not math.isfinite(fraccion_nervios) or not 0 < fraccion_nervios < 1:
        raise ValueError("La geometría de los nervios no define una fracción de hormigón válida.")
    return capa + (h - capa) * fraccion_nervios


def componentes_con_hormigon_casetonada(elemento: dict) -> list[dict]:
    """Sustituye capas estructurales ingresadas por el peso geométrico del hormigón."""
    componentes = [
        dict(componente)
        for componente in elemento.get("componentes", [])
        if (componente.get("grupo"), componente.get("clave")) not in {
            ("Forjados", "Losa_casetonada"),
            ("Forjados", "Losa_maciza"),
            ("Hormigon", "armado"),
        }
    ]
    componentes.append({
        "grupo": "Hormigon",
        "clave": "armado",
        "espesor_m": espesor_equivalente_hormigon_casetonada(
            elemento.get("casetonada", {})
        ),
    })
    return componentes

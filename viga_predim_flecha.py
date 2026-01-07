#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
viga_predim_flecha.py
Predimensionado de vigas de hormigón armado por flecha (criterio servicio)
Unidades modernas: m, kN/m, MPa
"""

import math

def predimension_viga(L, q, b=0.25, hormigon="H25", delta_adm_ratio=250):
    """
    L           : Luz de la viga en m
    q           : Carga distribuida lineal en kN/m (sobre ancho b)
    b           : Ancho de la viga en m
    hormigon    : H21 o H25
    delta_adm_ratio : Flecha admisible = L / delta_adm_ratio
    """
    # Módulo de elasticidad en kN/m²
    if hormigon == "H21":
        fc = 21  # MPa
    elif hormigon == "H25":
        fc = 25
    else:
        raise ValueError("Hormigón no reconocido. Use H21 o H25.")

    Ec = 4700 * math.sqrt(fc) * 1e3  # kN/m² (1 MPa = 1e3 kN/m²)

    # Flecha admisible en m
    delta_adm = L / delta_adm_ratio

    # Inercia mínima por flecha
    I_min = (5 * q * L**4) / (384 * Ec * delta_adm)

    # Altura mínima para sección rectangular
    h = (12 * I_min / b) ** (1/3)

    # Redondeo práctico en cm
    h_redondeado = math.ceil(h*100 / 5) * 5 / 100  # múltiplos de 5 cm

    return {
        "L_m": L,
        "q_kN_m": q,
        "b_m": b,
        "h_min_m": round(h, 3),
        "h_redondeado_m": h_redondeado,
        "I_m4": round(I_min, 4),
        "delta_adm_m": round(delta_adm, 3),
        "Ec_kN_m2": round(Ec, 1)
    }

# ==========================
# Bloque principal para prueba
if __name__ == "__main__":
    L = 4.80         # Luz en m
    q = 44.0         # Carga lineal en kN/m
    b = 0.20         # Ancho en m
    hormigon = "H21"
    delta_ratio = 250

    resultado = predimension_viga(L, q, b, hormigon, delta_ratio)

    print("=== Predimensionado de viga por flecha ===")
    for k, v in resultado.items():
        print(f"{k}: {v}")

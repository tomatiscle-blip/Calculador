#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
viga_predim_flecha_ext_print.py
Predimensionado de vigas de hormigón armado por flecha (criterio servicio)
Incluye: cálculo con inercia bruta y efectiva, ratio L/h, flecha corto y largo plazo
Unidades: m, kN/m, MPa
"""

import math

def predimension_viga(L, q, b=0.25, hormigon="H25", delta_adm_ratio=250,
                      alpha_I=0.35, phi=2.0):
    # Diccionario de resistencias
    fc_dict = {"H21": 21, "H25": 25, "H30": 30, "H35": 35}
    fc = fc_dict.get(hormigon, None)
    if fc is None:
        raise ValueError("Hormigón no reconocido. Use H21, H25, H30 o H35.")

    # Módulo de elasticidad en kN/m²
    Ec = 4700 * math.sqrt(fc) * 1e3

    # Flecha admisible en m
    delta_adm = L / delta_adm_ratio

    # Inercia mínima por flecha (bruta)
    I_min = (5 * q * L**4) / (384 * Ec * delta_adm)

    # Altura mínima para sección rectangular
    h = (12 * I_min / b) ** (1/3)

    # Redondeo práctico en cm (múltiplos de 5 cm)
    h_redondeado = math.ceil(h*100 / 5) * 5 / 100

    # Ratio L/h
    L_h_ratio = round(L / h_redondeado, 1)

    # Inercia efectiva (fisuración)
    I_eff = alpha_I * I_min

    # Flecha con I_eff (corto plazo)
    delta_CP = (5 * q * L**4) / (384 * Ec * I_eff)

    # Flecha largo plazo (fluencia)
    delta_LP = delta_CP * (1 + phi)

    return {
        # === Datos de entrada ===
        "L_m": L,
        "q_kN_m": q,
        "b_m": b,
        "h_min_m": round(h, 3),
        "h_redondeado_m": h_redondeado,
        "L/h": L_h_ratio,

        # === Propiedades del material ===
        "Hormigon": hormigon,
        "Ec_kN_m2": round(Ec, 1),

        # === Inercias ===
        "I_bruta_m4": round(I_min, 4),
        "I_bruta_cm4": round(I_min * 1e8, 1),
        "I_efectiva_m4": round(I_eff, 4),
        "I_efectiva_cm4": round(I_eff * 1e8, 1),

        # === Flechas ===
        "delta_adm_m": round(delta_adm, 3),
        "delta_adm_cm": round(delta_adm * 100, 1),
        "delta_CP_m": round(delta_CP, 3),
        "delta_CP_cm": round(delta_CP * 100, 1),
        "delta_LP_m": round(delta_LP, 3),
        "delta_LP_cm": round(delta_LP * 100, 1)
    }

def imprimir_resultado(resultado):
    print("=== Datos de entrada ===")
    print(f"L (m): {resultado['L_m']}")
    print(f"q (kN/m): {resultado['q_kN_m']}")
    print(f"b (m): {resultado['b_m']}")
    print(f"h_min (m): {resultado['h_min_m']}")
    print(f"h_redondeado (m): {resultado['h_redondeado_m']}")
    print(f"L/h: {resultado['L/h']}")

    print("\n=== Propiedades del material ===")
    print(f"Hormigón: {resultado['Hormigon']}")
    print(f"Ec (kN/m²): {resultado['Ec_kN_m2']}")

    print("\n=== Inercias ===")
    print(f"I_bruta: {resultado['I_bruta_m4']} m4   ({resultado['I_bruta_cm4']} cm4)")
    print(f"I_efectiva: {resultado['I_efectiva_m4']} m4   ({resultado['I_efectiva_cm4']} cm4)")

    print("\n=== Flechas ===")
    # δ_adm: límite reglamentario (L/x), no depende de cálculo de rigidez
    print(f"δ_adm (flecha admisible reglamento): {resultado['delta_adm_m']} m   ({resultado['delta_adm_cm']} cm)")

    # δ_CP: flecha calculada con inercia efectiva (fisuración), estado de servicio a corto plazo
    print(f"δ_CP (flecha con I_eff, corto plazo): {resultado['delta_CP_m']} m   ({resultado['delta_CP_cm']} cm)")
    # δ_LP: flecha calculada con I_eff + fluencia, estado de servicio a largo plazo
    print(f"δ_LP (flecha con I_eff + fluencia, largo plazo): {resultado['delta_LP_m']} m   ({resultado['delta_LP_cm']} cm)")

# ==========================
# Bloque principal para prueba
# Criterios de flecha admisible según CIRSOC 201:
# - Vigas y losas en general ............ L/250
# - Vigas con tabiques o acabados frágiles L/360
# - Losas con cielorrasos/terminaciones sensibles L/480
# - Voladizos ........................... L/180
#
# Nota: se evalúa con cargas de servicio (no mayoradas).
#       Considerar fisuración y fluencia para flecha real.
# ==========================

if __name__ == "__main__":
    L = 6.20
    q = 20
    b = 0.20
    hormigon = "H21"
    delta_ratio = 360  # L/360 para vigas con tabiques o acabados frágiles

    resultado = predimension_viga(L, q, b, hormigon, delta_ratio)
    imprimir_resultado(resultado)

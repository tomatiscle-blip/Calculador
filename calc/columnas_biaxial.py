"""Calculo de flexocompresion biaxial para columnas globales (Metodo de Bresler)."""
from __future__ import annotations


def calcular_p0(fck_mpa, fy_mpa, b_cm, h_cm, ast_cm2=0.0):
    ag_m2 = (b_cm / 100.0) * (h_cm / 100.0)
    ast_m2 = ast_cm2 / 10000.0
    pn_n = 0.85 * fck_mpa * 1e6 * (ag_m2 - ast_m2) + fy_mpa * 1e6 * ast_m2
    return pn_n / 1000.0


def bresler_pn(pnx_kn, pny_kn, p0_kn):
    if pnx_kn == 0 or pny_kn == 0 or p0_kn == 0:
        return 0.0
    inv = 1.0 / pnx_kn + 1.0 / pny_kn - 1.0 / p0_kn
    if inv <= 0:
        return float("inf")
    return 1.0 / inv


def verificar_biaxial(pu_kn, pnx_kn, pny_kn, p0_kn):
    pn = bresler_pn(pnx_kn, pny_kn, p0_kn)
    cumple = (pn == float("inf")) or (pu_kn <= pn)
    return {"Pu_kN": pu_kn, "Pnx_kN": pnx_kn, "Pny_kN": pny_kn, "P0_kN": p0_kn,
            "Pn_biaxial_kN": (None if pn == float("inf") else pn),
            "cumple": bool(cumple), "metodo": "Bresler (carga reciproca)"}


def generar_punto_diagrama_interaccion(c, b, h, fck, fy, armadura_capas):
    epsilon_cu = 0.003
    Es = 200000
    beta1 = 0.85 if fck <= 30 else max(0.65, 0.85 - 0.05 * (fck - 30) / 7)
    a = beta1 * c
    Cc = 0.85 * fck * a * b
    F_acero_total = 0
    M_acero_total = 0
    h_baricentro = h / 2
    for As_i, d_i in armadura_capas:
        epsilon_si = epsilon_cu * (c - d_i) / c
        sigma_si = Es * epsilon_si
        if sigma_si > fy:
            sigma_si = fy
        elif sigma_si < -fy:
            sigma_si = -fy
        F_i = As_i * sigma_si
        F_acero_total += F_i
        brazo = h_baricentro - d_i
        M_acero_total += F_i * brazo
    Pn = Cc + F_acero_total
    M_cc = Cc * (h_baricentro - a / 2)
    Mn = M_cc + M_acero_total
    return Pn, Mn


def combinar_momentos_globales(mx, my):
    return {"Mx_kNm": mx, "My_kNm": my,
            "Mu_max_kNm": max(abs(mx), abs(my)),
            "biaxial": abs(mx) > 0 and abs(my) > 0}

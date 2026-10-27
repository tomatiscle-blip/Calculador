"""Fase 4: verificacion biaxial (Bresler) para columnas globales unificadas."""
from __future__ import annotations

import importlib.util
import json
import csv
import csv
import sys
from pathlib import Path
from typing import Any

_RAIZ = Path(__file__).resolve().parent.parent
if str(_RAIZ) not in sys.path:
    sys.path.insert(0, str(_RAIZ))

from calc import rutas, columnas_biaxial




def cargar(seccion_tipo, fc, gamma):
    gamma_str = f"{int(gamma*100):03}"
    carpeta = rutas.DIAGRAMAS_INTERACCION
    for archivo in carpeta.iterdir():
        if archivo.name.startswith(f"diagrama{seccion_tipo}") and f"_fc{fc}_" in archivo.name and f"gamma{gamma_str}" in archivo.name:
            with open(archivo, encoding="utf-8") as f:
                return json.load(f)
    raise FileNotFoundError(f"No se encontro diagrama para fc={fc}, gamma={gamma}")

def interpola_y(points, x_obj):
    for i in range(len(points)-1):
        x1, y1 = points[i]["x"], points[i]["y"]
        x2, y2 = points[i+1]["x"], points[i+1]["y"]
        if (x1 <= x_obj <= x2) or (x2 <= x_obj <= x1):
            t = (x_obj - x1) / (x2 - x1)
            return y1 + t * (y2 - y1)
    return None

def interp_rho(curves, rho_target):
    curves_sorted = sorted(curves, key=lambda c: c["parameterValue"])
    if rho_target < curves_sorted[0]["parameterValue"] or rho_target > curves_sorted[-1]["parameterValue"]:
        raise ValueError(f"rho_target={rho_target} fuera del rango de curvas")
    c_inf = max([c for c in curves_sorted if c["parameterValue"] <= rho_target], key=lambda c: c["parameterValue"])
    c_sup = min([c for c in curves_sorted if c["parameterValue"] >= rho_target], key=lambda c: c["parameterValue"])
    if c_inf["parameterValue"] == c_sup["parameterValue"]:
        return c_inf["points"]
    t_rho = (rho_target - c_inf["parameterValue"]) / (c_sup["parameterValue"] - c_inf["parameterValue"])
    puntos_interp = []
    for p_inf, p_sup in zip(c_inf["points"], c_sup["points"]):
        x = p_inf["x"] + t_rho*(p_sup["x"] - p_inf["x"])
        y = p_inf["y"] + t_rho*(p_sup["y"] - p_inf["y"])
        puntos_interp.append({"x": x, "y": y})
    return puntos_interp


def _gamma_normativo(gamma: float) -> float:
    opciones = [0.70, 0.80, 0.90]
    return min(opciones, key=lambda g: abs(g - gamma))


def _pn_uniaxial(seccion_tipo: str, fck: float, gamma: float, b: float, h: float,
                 ast: float, mu_knm: float) -> float:
    """Capacidad axial nominal Pn para un Mu dado, usando el diagrama de interaccion."""
    if mu_knm == 0:
        # Pure axial not in diagram; approximate with P0
        return columnas_biaxial.calcular_p0(fck, 420.0, b, h, ast)
    rho = ast / (b * h) if b * h > 0 else 0.0
    gamma_n = _gamma_normativo(gamma)
    diag = cargar(seccion_tipo, int(fck), gamma_n)
    curva = interp_rho(diag["curves"], rho)
    x_target = 10.0 * mu_knm * 100.0 / (b * h * h)
    y = interpola_y(curva, x_target)
    if y is None:
        return 0.0
    pn = y * b * h / 10.0
    return float(pn)


def ejecutar() -> list[dict]:
    SAL = rutas.SAL_COLUMNAS
    SAL.mkdir(parents=True, exist_ok=True)

    archivos_unif = sorted(SAL.glob("unificado_*.json"))
    resultados: list[dict] = []

    for archivo in archivos_unif:
        u = rutas.leer_json(archivo, {}) or {}
        gid = u.get("columna_global") or archivo.stem.replace("unificado_", "")

        # tipo_seccion y recubrimiento desde Fase 1
        dim_file = SAL / ("dimensionado_" + gid + ".json")
        dim = rutas.leer_json(dim_file, {}) or {}
        dims = dim.get("dimensionados", {}) or {}
        tipo_seccion = "I"
        recubrimiento = 4.0
        if dims:
            d0 = list(dims.values())[0]
            tipo_seccion = d0.get("tipo_seccion", "I") or "I"
            recubrimiento = float(d0.get("recubrimiento_cm", 4.0) or 4.0)

        b = float(u.get("b_cm", 0) or 0)
        h = float(u.get("h_cm", 0) or 0)
        fck = float(u.get("fck_MPa", 0) or 0)
        ast = float(u.get("Ast_cm2", 0) or 0)
        n_val = float(u.get("N_kN", 0) or 0)
        mx = float(u.get("Mx_kNm", 0) or 0)
        my = float(u.get("My_kNm", 0) or 0)
        biaxial = bool(u.get("biaxial", False))

        gamma = (h - 2 * recubrimiento) / h if h > 0 else 0.8

        if not biaxial:
            resultado = {
                "columna_global": gid,
                "tipo": "uniaxial",
                "biaxial": False,
                "verificacion": "RESUELTA (uniaxial)",
                "N_kN": n_val,
                "Mx_kNm": mx,
                "My_kNm": my,
                "b_cm": b,
                "h_cm": h,
                "fck_MPa": fck,
                "Ast_cm2": ast,
            }
        else:
            pnx = _pn_uniaxial(tipo_seccion, fck, gamma, b, h, ast, mx)
            pny = _pn_uniaxial(tipo_seccion, fck, gamma, b, h, ast, my)
            p0 = columnas_biaxial.calcular_p0(fck, 420.0, b, h, ast)
            verif = columnas_biaxial.verificar_biaxial(n_val, pnx, pny, p0)
            resultado = {
                "columna_global": gid,
                "tipo": "biaxial",
                "biaxial": True,
                "N_kN": n_val,
                "Mx_kNm": mx,
                "My_kNm": my,
                "b_cm": b,
                "h_cm": h,
                "fck_MPa": fck,
                "Ast_cm2": ast,
                "gamma": gamma,
                "Pnx_kN": pnx,
                "Pny_kN": pny,
                "P0_kN": p0,
                "verificacion_biaxial": verif,
            }

        rutas.guardar_json(SAL / ("resultado_final_" + gid + ".json"), resultado)
        resultados.append(resultado)

    return resultados


if __name__ == "__main__":
    res = ejecutar()
    print("Resultados finales generados:", len(res))
    for r in res:
        print(r["columna_global"], r.get("tipo"), r.get("verificacion", r.get("verificacion_biaxial", {}).get("cumple")))

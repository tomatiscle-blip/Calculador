from __future__ import annotations
import json
import sys
from pathlib import Path

_RAIZ = Path(__file__).resolve().parent.parent
if str(_RAIZ) not in sys.path:
    sys.path.insert(0, str(_RAIZ))

from calc import rutas

def gamma_normativo(gamma):
    opciones = [0.70, 0.80, 0.90]
    return min(opciones, key=lambda g: abs(g - gamma))

def cargar(seccion_tipo, fc, gamma):
    gamma_str = f"{int(gamma*100):03}"
    carpeta = rutas.DIAGRAMAS_INTERACCION
    for archivo in carpeta.iterdir():
        if archivo.name.startswith(f"diagrama{seccion_tipo}") and f"_fc{fc}_" in archivo.name and f"gamma{gamma_str}" in archivo.name:
            with open(archivo, encoding="utf-8") as f:
                return json.load(f)
    return None

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
        return None
    c_inf = max([c for c in curves_sorted if c["parameterValue"] <= rho_target], key=lambda c: c["parameterValue"])
    c_sup = min([c for c in curves_sorted if c["parameterValue"] >= rho_target], key=lambda c: c["parameterValue"])
    if c_inf["parameterValue"] == c_sup["parameterValue"]:
        return c_inf["points"]
    t_rho = (rho_target - c_inf["parameterValue"]) / (c_sup["parameterValue"] - c_inf["parameterValue"])
    puntos = []
    for p_inf, p_sup in zip(c_inf["points"], c_sup["points"]):
        x = p_inf["x"] + t_rho*(p_sup["x"] - p_inf["x"])
        y = p_inf["y"] + t_rho*(p_sup["y"] - p_inf["y"])
        puntos.append({"x": x, "y": y})
    return puntos

def pn_uniaxial(seccion_tipo, fck, gamma, b, h, ast, mu_knm):
    if mu_knm == 0 or b*h == 0:
        return 0.0
    rho = ast/(b*h)
    gamma_n = gamma_normativo(gamma)
    diag = cargar(seccion_tipo, int(fck), gamma_n)
    if diag is None:
        return 0.0
    curva = interp_rho(diag["curves"], rho)
    if curva is None:
        return 0.0
    x_target = 10.0*mu_knm*100.0/(b*h*h)
    y = interpola_y(curva, x_target)
    if y is None:
        return 0.0
    return float(y*b*h/10.0)

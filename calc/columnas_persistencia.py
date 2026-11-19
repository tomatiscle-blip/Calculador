"""Persistencia estructurada de dimensionado de columnas por pórtico (Fase 1)."""
from __future__ import annotations
import csv, json, re, sys
from pathlib import Path
from typing import Any

_RAIZ = Path(__file__).resolve().parent.parent
if str(_RAIZ) not in sys.path:
    sys.path.insert(0, str(_RAIZ))

from calc import rutas
from calc import columnas_diagramas as cd
from calc.formato import parse_numero

def _parse_dimensiones(texto):
    m = re.search(r"b\s*=\s*([\d.,]+)\s*,\s*h\s*=\s*([\d.,]+)", str(texto))
    if not m:
        return 0.0, 0.0
    return parse_numero(m.group(1)), parse_numero(m.group(2))

def _parse_nota(texto):
    out = {"gamma": 0.0, "X": 0.0, "Y": 0.0, "rho": 0.0}
    for clave, key in (("γ", "gamma"), ("X", "X"), ("Y", "Y"), ("ρg", "rho")):
        m = re.search(rf"{re.escape(clave)}\s*=\s*([\d.,]+)", str(texto))
        if m:
            out[key] = parse_numero(m.group(1))
    return out

def _leer_csv():
    ruta = rutas.PLANILLA_COLUMNAS
    if not ruta.exists():
        return []
    with ruta.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter=";")
        return [row for row in reader]

def _fila_a_dimensionado(row):
    b, h = _parse_dimensiones(row.get("Dimensiones", ""))
    nota = _parse_nota(row.get("Nota_diagrama", ""))
    n_barras = int(parse_numero(row.get("n_barras"), 0))
    diam_long = parse_numero(row.get("diam_long(cm)"), 0)
    diam_estribo = parse_numero(row.get("diam_estribo(cm)"), 0)
    paso_estribo = parse_numero(row.get("paso_estribo(cm)"), 0)
    rec = parse_numero(row.get("Recubrimiento(cm)"), 4)
    gamma = (h - 2*rec)/h if h > 0 else 0.8
    seccion = row.get("Tipo_seccion", "I") or "I"
    fck = parse_numero(row.get("f'c(MPa)"), 0)
    ast = parse_numero(row.get("Ast(cm²)"), 0)
    mu = parse_numero(row.get("Mu(kNm)"), 0)
    pn = cd.pn_uniaxial(seccion, fck, gamma, b, h, ast, mu)
    return {
        "portico_origen": row.get("Pórtico", ""),
        "columna_local": row.get("Columna", ""),
        "tipo": row.get("Tipo", ""),
        "tipo_seccion": row.get("Tipo_seccion", ""),
        "b_cm": b,
        "h_cm": h,
        "lado_x_cm": h,
        "lado_y_cm": b,
        "altura_libre_m": parse_numero(row.get("Altura_libre(m)"), 0),
        "recubrimiento_cm": parse_numero(row.get("Recubrimiento(cm)"), 0),
        "Pu_kN": parse_numero(row.get("Pu(kN)"), 0),
        "Mu_kNm": parse_numero(row.get("Mu(kNm)"), 0),
        "fck_MPa": parse_numero(row.get("f'c(MPa)"), 0),
        "lambda": parse_numero(row.get("λ"), 0),
        "clasificacion": row.get("Clasificación", ""),
        "Ast_cm2": parse_numero(row.get("Ast(cm²)"), 0),
        "armadura_longitudinal": {
            "n_barras": n_barras,
            "diametro_mm": diam_long,
            "descripcion": f"{n_barras}Ø{int(round(diam_long))}" if n_barras else "",
        },
        "estribos": {
            "diametro_mm": diam_estribo,
            "paso_cm": paso_estribo,
            "cumple": (row.get("cumple_estribo", "") or "").strip(),
        },
        "gamma": nota["gamma"],
        "X": nota["X"],
        "Y": nota["Y"],
        "rho": nota["rho"],
        "nota_diagrama": row.get("Nota_diagrama", ""),
        "Pn_kN": pn,
    }

def generar_dimensionados_json(portico=None):
    estructura = rutas.cargar_estructura()
    catalogo = rutas._catalogo_columnas_globales(estructura) if hasattr(rutas, "_catalogo_columnas_globales") else {}
    if not catalogo:
        return {}
    filas = _leer_csv()
    if portico is not None:
        filas = [r for r in filas if r.get("Pórtico") == portico]
    resultados = {}
    SAL = rutas.SAL_COLUMNAS
    SAL.mkdir(parents=True, exist_ok=True)
    por_gid = {}
    for row in filas:
        p = row.get("Pórtico", "")
        cid = row.get("Columna", "")
        mapeo = rutas._mapear_cids_a_globales(p, estructura, catalogo)
        gid = mapeo.get(cid)
        if gid is None:
            continue
        dim = _fila_a_dimensionado(row)
        datos_portico = rutas.pórtico(estructura, p) or {}
        direccion = str((datos_portico.get("referencia_planta", {}) or {}).get("direccion", "")).lower()
        if direccion not in ("x", "y"):
            # Sin orientación confiable no se debe adivinar el lado global.
            dim["lado_x_cm"] = None
            dim["lado_y_cm"] = None
            dim["orientacion_global_valida"] = False
        else:
            dim["direccion_plano"] = direccion
            dim["lado_x_cm"], dim["lado_y_cm"] = (dim["h_cm"], dim["b_cm"]) if direccion == "x" else (dim["b_cm"], dim["h_cm"])
            dim["orientacion_global_valida"] = True
        por_gid.setdefault(gid, {})[p] = dim
    for gid, dims_por_portico in por_gid.items():
        archivo = SAL / f"dimensionado_{gid}.json"
        previo = rutas.leer_json(archivo, {}) or {}
        dimensionados = previo.get("dimensionados", {}) or {}
        dimensionados.update(dims_por_portico)
        salida = {
            "columna_global": gid,
            "catalogo": catalogo.get(gid, {}),
            "dimensionados": dimensionados,
            "n_porticos": len(dimensionados),
        }
        rutas.guardar_json(archivo, salida)
        resultados[gid] = archivo
    return resultados

if __name__ == "__main__":
    res = generar_dimensionados_json()
    print(f"JSON de dimensionado generados: {len(res)}")
    for gid, ruta in list(res.items())[:3]:
        print(gid, "->", ruta)

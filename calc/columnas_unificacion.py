"""Fase 2: unificacion de seccion/armadura por columna global (sin Bresler)."""
from __future__ import annotations

import csv
import sys
from pathlib import Path
from typing import Any

_RAIZ = Path(__file__).resolve().parent.parent
if str(_RAIZ) not in sys.path:
    sys.path.insert(0, str(_RAIZ))

from calc import rutas

def _cargar_resultados_globales() -> list[dict]:
    archivo = rutas.SAL_COLUMNAS / "resultados_columnas_globales.json"
    return rutas.leer_json(archivo, []) or []

def unificar(criterio: str = "max_ast") -> list[dict]:
    SAL = rutas.SAL_COLUMNAS
    SAL.mkdir(parents=True, exist_ok=True)
    if criterio not in ("max_ast", "max_area", "max_pu"):
        criterio = "max_ast"
    resultados_globales = {r["id"]: r for r in _cargar_resultados_globales()}
    archivos_dim = sorted(SAL.glob("dimensionado_*.json"))
    unificados = []
    for archivo in archivos_dim:
        datos = rutas.leer_json(archivo, {}) or {}
        gid = datos.get("columna_global") or archivo.stem.replace("dimensionado_", "")
        catalogo = datos.get("catalogo", {}) or {}
        dimensionados = datos.get("dimensionados", {}) or {}
        if not dimensionados:
            continue
        mejor_portico = None
        mejor_score = None
        for p, d in dimensionados.items():
            if criterio == "max_area":
                score = (float(d.get("b_cm", 0) or 0) * float(d.get("h_cm", 0) or 0), float(d.get("Ast_cm2", 0) or 0))
            elif criterio == "max_pu":
                score = (float(d.get("Pu_kN", 0) or 0), float(d.get("Ast_cm2", 0) or 0))
            else:
                score = (float(d.get("Ast_cm2", 0) or 0), float(d.get("b_cm", 0) or 0) * float(d.get("h_cm", 0) or 0))
            if mejor_score is None or score > mejor_score:
                mejor_score = score
                mejor_portico = p
        base = dimensionados[mejor_portico]
        # detectar secciones distintas entre porticos
        secciones = {(float(d.get("b_cm",0) or 0), float(d.get("h_cm",0) or 0)) for d in dimensionados.values()}
        secciones_distintas = len(secciones) > 1
        fck = max(float(d.get("fck_MPa", 0) or 0) for d in dimensionados.values())
        ast_max = max(float(d.get("Ast_cm2", 0) or 0) for d in dimensionados.values())
        portico_ast_max = max(dimensionados, key=lambda p: float(dimensionados[p].get("Ast_cm2", 0) or 0))
        armadura_desc = dimensionados[portico_ast_max]["armadura_longitudinal"].get("descripcion", "")
        n_barras = dimensionados[portico_ast_max]["armadura_longitudinal"].get("n_barras", 0)
        diam_long = dimensionados[portico_ast_max]["armadura_longitudinal"].get("diametro_mm", 0)
        diam_estribo = max(float(d.get("estribos", {}).get("diametro_mm", 0) or 0) for d in dimensionados.values())
        pasos = [float(d.get("estribos", {}).get("paso_cm", 0) or 0) for d in dimensionados.values() if float(d.get("estribos", {}).get("paso_cm", 0) or 0) > 0]
        paso_estribo = min(pasos) if pasos else 0.0
        clasificaciones = [str(d.get("clasificacion", "")).upper() for d in dimensionados.values()]
        clasificacion = "ESBELTA" if "ESBELTA" in clasificaciones else (clasificaciones[0] if clasificaciones else "")
        lambda_val = max(float(d.get("lambda", 0) or 0) for d in dimensionados.values())
        rg = resultados_globales.get(gid, {})
        mx = float(rg.get("Mx", 0.0) or 0.0)
        my = float(rg.get("My", 0.0) or 0.0)
        n_val = float(rg.get("N", 0.0) or 0.0)
        porticos_nombres = rg.get("porticos", list(dimensionados.keys()))
        biaxial = (abs(mx) > 0 and abs(my) > 0)
        unificado = {
            "columna_global": gid,
            "eje_nombre": catalogo.get("eje_nombre"),
            "nivel": catalogo.get("nivel"),
            "N_kN": n_val,
            "Mx_kNm": mx,
            "My_kNm": my,
            "Mu_max_kNm": max(abs(mx), abs(my)),
            "biaxial": biaxial,
            "b_cm": float(base.get("b_cm", 0) or 0),
            "h_cm": float(base.get("h_cm", 0) or 0),
            "fck_MPa": fck,
            "Ast_cm2": ast_max,
            "armadura_longitudinal": {"n_barras": n_barras, "diametro_mm": diam_long, "descripcion": armadura_desc},
            "estribos": {"diametro_mm": diam_estribo, "paso_cm": paso_estribo},
            "clasificacion": clasificacion,
            "lambda": lambda_val,
            "porticos": porticos_nombres,
            "n_porticos": len(porticos_nombres),
            "seccion_origen": mejor_portico,
            "armadura_origen": portico_ast_max,
            "secciones_distintas": secciones_distintas,
            "criterio_usado": criterio,
            "criterio_unificacion": "max Ast (desempate mayor area b*h); fck max; estribos max diam/min paso; clasificacion si alguno ESBELTA",
        }
        rutas.guardar_json(SAL / ("unificado_" + gid + ".json"), unificado)
        unificados.append(unificado)


    csv_path = SAL / "planilla_columnas_globales.csv"
    columnas = ["Columna","Eje","Nivel","N(kN)","Mx(kNm)","My(kNm)","Mu_max(kNm)","Biaxial","b(cm)","h(cm)","f'c(MPa)","Ast(cm2)","Armadura","Estribos","Clasificacion","lambda","Porticos","n_porticos","Seccion_origen","Armadura_origen"]
    with csv_path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.writer(f, delimiter=";")
        writer.writerow(columnas)
        for u in unificados:
            writer.writerow([
                u["columna_global"], u["eje_nombre"], u["nivel"],
                round(u["N_kN"],2), round(u["Mx_kNm"],2), round(u["My_kNm"],2), round(u["Mu_max_kNm"],2),
                "SI" if u["biaxial"] else "NO",
                u["b_cm"], u["h_cm"], u["fck_MPa"], u["Ast_cm2"],
                u["armadura_longitudinal"]["descripcion"],
                "O" + str(u["estribos"]["diametro_mm"]) + " c/" + str(u["estribos"]["paso_cm"]),
                u["clasificacion"], u["lambda"], ";".join(u["porticos"]), u["n_porticos"],
                u["seccion_origen"], u["armadura_origen"],
            ])
    return unificados

if __name__ == "__main__":
    res = unificar()
    print("Columnas globales unificadas:", len(res))
    print("Planilla:", rutas.SAL_COLUMNAS / "planilla_columnas_globales.csv")

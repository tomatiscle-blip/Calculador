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
from calc.formato import valor_csv

def _cargar_resultados_globales() -> list[dict]:
    archivo = rutas.SAL_COLUMNAS / "resultados_columnas_globales.json"
    return rutas.leer_json(archivo, []) or []

def unificar(criterio: str = "max_ast") -> list[dict]:
    SAL = rutas.SAL_COLUMNAS
    SAL.mkdir(parents=True, exist_ok=True)
    if criterio not in ("max_ast", "max_area", "max_pu"):
        criterio = "max_ast"
    resultados_globales = {r["id"]: r for r in _cargar_resultados_globales() if r.get("id")}
    estructura = rutas.cargar_estructura()
    catalogo_global = rutas._catalogo_columnas_globales(estructura)
    dimensiones_por_gid = {}
    for archivo in sorted(SAL.glob("dimensionado_*.json")):
        datos = rutas.leer_json(archivo, {}) or {}
        gid = datos.get("columna_global") or archivo.stem.replace("dimensionado_", "")
        if datos.get("dimensionados"):
            dimensiones_por_gid[gid] = datos
    ids = sorted(set(resultados_globales) | set(catalogo_global) | set(dimensiones_por_gid))
    unificados = []
    for gid in ids:
        datos = dimensiones_por_gid.get(gid, {})
        catalogo = datos.get("catalogo", {}) or catalogo_global.get(gid, {})
        dimensionados = datos.get("dimensionados", {}) or {}
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
        base = dimensionados[mejor_portico] if mejor_portico else {}
        # detectar secciones distintas entre porticos
        secciones = {(float(d.get("lado_x_cm") or 0), float(d.get("lado_y_cm") or 0)) for d in dimensionados.values()}
        secciones_distintas = len(secciones) > 1
        fck = max((float(d.get("fck_MPa", 0) or 0) for d in dimensionados.values()), default=0.0)
        ast_max = max((float(d.get("Ast_cm2", 0) or 0) for d in dimensionados.values()), default=0.0)
        portico_ast_max = max(dimensionados, key=lambda p: float(dimensionados[p].get("Ast_cm2", 0) or 0)) if dimensionados else None
        detalle_armadura = dimensionados[portico_ast_max].get("armadura_longitudinal", {}) if portico_ast_max else {}
        armadura_desc = detalle_armadura.get("descripcion", "")
        n_barras = detalle_armadura.get("n_barras", 0)
        diam_long = detalle_armadura.get("diametro_mm", 0)
        diam_estribo = max((float(d.get("estribos", {}).get("diametro_mm", 0) or 0) for d in dimensionados.values()), default=0.0)
        pasos = [float(d.get("estribos", {}).get("paso_cm", 0) or 0) for d in dimensionados.values() if float(d.get("estribos", {}).get("paso_cm", 0) or 0) > 0]
        paso_estribo = min(pasos) if pasos else 0.0
        clasificaciones = [str(d.get("clasificacion", "")).upper() for d in dimensionados.values()]
        clasificacion = "ESBELTA" if "ESBELTA" in clasificaciones else (clasificaciones[0] if clasificaciones else "")
        lambda_val = max((float(d.get("lambda", 0) or 0) for d in dimensionados.values()), default=0.0)
        rg = resultados_globales.get(gid, {})
        mx = float(rg.get("Mx", 0.0) or 0.0)
        my = float(rg.get("My", 0.0) or 0.0)
        n_val = float(rg.get("N", 0.0) or 0.0)
        porticos_nombres = rg.get("porticos", list(dimensionados.keys()))
        biaxial = (abs(mx) > 0 and abs(my) > 0)
        dimensionados_con_seccion = [d for d in dimensionados.values() if d.get("b_cm") and d.get("h_cm")]
        orientacion_valida = all(d.get("orientacion_global_valida") for d in dimensionados_con_seccion)
        lado_x = float(base.get("lado_x_cm") or 0) if base else 0.0
        lado_y = float(base.get("lado_y_cm") or 0) if base else 0.0
        estado_dimensionado = "Dimensionada" if base and lado_x and lado_y else "Pendiente de dimensionado"
        if base and not orientacion_valida:
            estado_dimensionado = "Revisar orientación del pórtico"
        elif base and secciones_distintas:
            estado_dimensionado = "Revisar secciones entre pórticos"
        altura = catalogo.get("altura_m", catalogo.get("altura", ""))
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
            "lado_x_cm": lado_x,
            "lado_y_cm": lado_y,
            "altura_m": altura,
            "fck_MPa": fck,
            "Ast_cm2": ast_max,
            "armadura_longitudinal": {"n_barras": n_barras, "diametro_mm": diam_long, "descripcion": armadura_desc},
            "estribos": {"diametro_mm": diam_estribo, "paso_cm": paso_estribo, "cumple": any((d.get("estribos", {}).get("cumple") or "").lower() in ("si", "sí", "true", "ok") for d in dimensionados.values())},
            "recubrimiento_cm": float(base.get("recubrimiento_cm", 0) or 0) if base else 0.0,
            "tipo_seccion": base.get("tipo_seccion", "") if base else "",
            "orientacion_global_valida": orientacion_valida,
            "estado_dimensionado": estado_dimensionado,
            "clasificacion": clasificacion,
            "lambda": lambda_val,
            "porticos": porticos_nombres,
            "n_porticos": len(porticos_nombres),
            "seccion_origen": mejor_portico or "",
            "armadura_origen": portico_ast_max or "",
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
            writer.writerow([valor_csv(valor) for valor in [
                u["columna_global"], u["eje_nombre"], u["nivel"],
                round(u["N_kN"],2), round(u["Mx_kNm"],2), round(u["My_kNm"],2), round(u["Mu_max_kNm"],2),
                "SI" if u["biaxial"] else "NO",
                u["b_cm"], u["h_cm"], u["fck_MPa"], u["Ast_cm2"],
                u["armadura_longitudinal"]["descripcion"],
                "O" + str(u["estribos"]["diametro_mm"]) + " c/" + str(u["estribos"]["paso_cm"]),
                u["clasificacion"], u["lambda"], ";".join(u["porticos"]), u["n_porticos"],
                u["seccion_origen"], u["armadura_origen"],
            ]])
    return unificados

if __name__ == "__main__":
    res = unificar()
    print("Columnas globales unificadas:", len(res))
    print("Planilla:", rutas.SAL_COLUMNAS / "planilla_columnas_globales.csv")

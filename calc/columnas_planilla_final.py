"""Genera la planilla maestra constructiva de columnas de la obra."""
from __future__ import annotations

import csv
import sys
from pathlib import Path

_RAIZ = Path(__file__).resolve().parent.parent
if str(_RAIZ) not in sys.path:
    sys.path.insert(0, str(_RAIZ))

from calc import rutas
from calc.formato import valor_csv


def generar() -> Path:
    """Una fila por columna física, con sección referida a los ejes globales."""
    sal = rutas.SAL_COLUMNAS
    sal.mkdir(parents=True, exist_ok=True)
    unificados = {
        p.stem.replace("unificado_", ""): rutas.leer_json(p, {}) or {}
        for p in sal.glob("unificado_*.json")
    }
    finales = {
        p.stem.replace("resultado_final_", ""): rutas.leer_json(p, {}) or {}
        for p in sal.glob("resultado_final_*.json")
    }
    ruta_csv = sal / "planilla_columnas_globales_final.csv"
    encabezados = [
        "Columna", "Eje", "X global (m)", "Y global (m)", "Nivel (m)",
        "Altura (m)", "Pórticos", "Lado X (cm)", "Lado Y (cm)",
        "Sección", "Hormigón f'c (MPa)", "Armadura longitudinal",
        "Barras", "Ø longitudinal (mm)", "Ast (cm²)", "Estribos Ø (mm)",
        "Separación estribos (cm)", "Recubrimiento (cm)", "Flexión",
        "Estado dimensionado", "Verificación biaxial",
    ]
    estructura = rutas.cargar_estructura()
    catalogo = rutas._catalogo_columnas_globales(estructura)
    ids = sorted(set(catalogo) | set(unificados) | set(finales))

    with ruta_csv.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.writer(f, delimiter=";")
        writer.writerow(encabezados)
        for gid in ids:
            u = unificados.get(gid, {})
            fin = finales.get(gid, {})
            cat = catalogo.get(gid, {})
            armado = u.get("armadura_longitudinal", {}) or {}
            estribos = u.get("estribos", {}) or {}
            verif = fin.get("verificacion_biaxial", {}) or {}
            biaxial = bool(u.get("biaxial", fin.get("biaxial", False)))
            if fin.get("tipo") == "pendiente" or not fin:
                estado_verif = "PENDIENTE"
            elif biaxial:
                estado_verif = "CUMPLE" if verif.get("cumple") else "NO CUMPLE"
            else:
                estado_verif = "No aplica: uniaxial" if fin else "Pendiente"
            writer.writerow([valor_csv(valor) for valor in [
                gid,
                cat.get("eje_nombre", u.get("eje_nombre", "")),
                cat.get("X", cat.get("x", "")),
                cat.get("Y", cat.get("y", "")),
                cat.get("nivel", u.get("nivel", "")),
                u.get("altura_m", ""),
                ";".join(u.get("porticos", [])),
                u.get("lado_x_cm", ""),
                u.get("lado_y_cm", ""),
                u.get("tipo_seccion", ""),
                u.get("fck_MPa", ""),
                armado.get("descripcion", ""),
                armado.get("n_barras", ""),
                armado.get("diametro_mm", ""),
                u.get("Ast_cm2", ""),
                estribos.get("diametro_mm", ""),
                estribos.get("paso_cm", ""),
                u.get("recubrimiento_cm", ""),
                "Biaxial" if biaxial else "Uniaxial",
                u.get("estado_dimensionado", "Pendiente de dimensionado"),
                estado_verif,
            ]])
    print("Planilla maestra de columnas generada:", ruta_csv)
    return ruta_csv


if __name__ == "__main__":
    generar()

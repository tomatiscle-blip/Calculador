from __future__ import annotations
import csv
import sys
from pathlib import Path

_RAIZ = Path(__file__).resolve().parent.parent
if str(_RAIZ) not in sys.path:
    sys.path.insert(0, str(_RAIZ))

from calc import rutas

def generar():
    SAL = rutas.SAL_COLUMNAS
    SAL.mkdir(parents=True, exist_ok=True)
    unifs = {p.stem.replace("unificado_", ""): rutas.leer_json(p, {}) for p in SAL.glob("unificado_*.json")}
    finals = {p.stem.replace("resultado_final_", ""): rutas.leer_json(p, {}) for p in SAL.glob("resultado_final_*.json")}
    csv_path = SAL / "planilla_columnas_globales_final.csv"
    cols = ["Columna", "Eje", "Nivel", "N(kN)", "Mx(kNm)", "My(kNm)", "Mu_max(kNm)", "Biaxial", "b(cm)", "h(cm)", "f'c(MPa)", "Ast(cm2)", "Armadura", "Clasificacion", "Pnx(kN)", "Pny(kN)", "P0(kN)", "Pn_biaxial(kN)", "Cumple", "Porticos"]
    with csv_path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(cols)
        for gid, u in unifs.items():
            fin = finals.get(gid, {})
            verif = fin.get("verificacion_biaxial", {}) or {}
            biax = bool(u.get("biaxial", False))
            w.writerow([
                gid,
                u.get("eje_nombre", ""),
                u.get("nivel", ""),
                round(u.get("N_kN", 0), 2),
                round(u.get("Mx_kNm", 0), 2),
                round(u.get("My_kNm", 0), 2),
                round(u.get("Mu_max_kNm", 0), 2),
                "SI" if biax else "NO",
                u.get("b_cm", ""),
                u.get("h_cm", ""),
                u.get("fck_MPa", ""),
                u.get("Ast_cm2", ""),
                u.get("armadura_longitudinal", {}).get("descripcion", ""),
                u.get("clasificacion", ""),
                round(fin.get("Pnx_kN", 0), 2) if biax else "",
                round(fin.get("Pny_kN", 0), 2) if biax else "",
                round(fin.get("P0_kN", 0), 2) if biax else "",
                round(verif.get("Pn_biaxial_kN", 0), 2) if biax else "",
                verif.get("cumple", "") if biax else "RESUELTA",
                ";".join(u.get("porticos", [])),
            ])
    print("Planilla final generada:", csv_path)

if __name__ == "__main__":
    generar()

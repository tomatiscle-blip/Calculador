"""
Control temporal: corre anaStruct FRESCO por combinación (con voladizos y
cargas puntuales enteras, igual que P01) y lo compara contra calc.portico.

    py _tmp_validar_motor.py
"""

import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parent
sys.path.insert(0, str(RAIZ))

from calc import portico as motor  # noqa: E402
from calc import rutas  # noqa: E402
from anastruct import SystemElements  # noqa: E402

FACT = motor.combinaciones()


def crear(p, ss):
    eids = {}
    for cid, col in p["columnas"].items():
        x, y0, h = col["x"], col.get("nivel", 0.0), col["altura_m"]
        eids["C-" + cid] = int(ss.add_element(location=[[x, y0], [x, y0 + h]]))
    for base in p.get("bases", {}).values():
        n = ss.find_node_id([base["x"], 0.0])
        if base.get("tipo") == "articulado":
            ss.add_support_hinged(node_id=n)
        else:
            ss.add_support_fixed(node_id=n)
    for viga_id, viga in p["vigas"].items():
        piso = int(viga_id.split("-")[0][1:])
        cols = sorted((c for cid, c in p["columnas"].items() if cid.startswith(f"C{piso}-")),
                      key=lambda c: c["x"])
        if not cols:
            continue
        y = cols[0].get("nivel", 0.0) + cols[0]["altura_m"]
        for i, tramo in enumerate(viga["tramos"]):
            if tramo.get("es_voladizo"):
                if "izq" in tramo["id"]:
                    tip_x, base_x = cols[0]["x"] - tramo["longitud_m"], cols[0]["x"]
                else:
                    tip_x, base_x = cols[-1]["x"] + tramo["longitud_m"], cols[-1]["x"]
                eid = ss.add_element(location=[[tip_x, y], [base_x, y]])
            else:
                if i >= len(cols) - 1:
                    continue
                eid = ss.add_element(location=[[cols[i]["x"], y], [cols[i + 1]["x"], y]])
            eids["V-" + tramo["id"]] = int(eid)
    return eids


def aplicar(p, ss, eids, caso):
    fD, fL, fW, fP = FACT[caso]
    for viga_id, viga in p["vigas"].items():
        piso = int(viga_id.split("-")[0][1:])
        cols = sorted((c for cid, c in p["columnas"].items() if cid.startswith(f"C{piso}-")),
                      key=lambda c: c["x"])
        if not cols:
            continue
        if fW:
            w_total = float(viga["tramos"][0].get("cargas", {}).get("W_total", 0.0))
            por_col = w_total * cols[0]["altura_m"] / len(cols)
            for c in cols:
                n = ss.find_node_id([c["x"], c.get("nivel", 0.0) + c["altura_m"]])
                ss.point_load(node_id=n, Fx=fW * por_col)
        for tramo in viga["tramos"]:
            car = tramo.get("cargas", {})
            q = fD * float(car.get("D_total", 0.0)) + fL * float(car.get("L_total", 0.0))
            eid = eids["V-" + tramo["id"]]
            if abs(q) > 1e-9:
                ss.q_load(element_id=eid, q=-q)
            for cp in tramo.get("cargas_puntuales", []):
                n = ss.find_node_id([cp["coordenadas"]["x"], cp["coordenadas"]["y"]])
                ss.point_load(node_id=n, Fy=-fP * float(cp["valor_kN"]))


est = rutas.cargar_estructura()
p = est["Portico 1"]
res_motor = motor.calcular(p, "Portico 1")

sal = []
for caso in FACT:
    ss = SystemElements(mesh=50)
    eids = crear(p, ss)
    aplicar(p, ss, eids, caso)
    ss.solve()

    ev = ss.get_element_results(element_id=eids["V-V0-1 T1"], verbose=True)
    Mv, Qv = list(ev["M"]), list(ev["Q"])
    ec = ss.get_element_results(element_id=eids["C-C0-a"], verbose=True)
    Mc, Nc = list(ec["M"]), list(ec["N"])
    nb = ss.find_node_id([0.0, 0.0])
    R = ss.reaction_forces[nb]

    m = res_motor["solicitaciones"][caso]
    a = m["vigas"]["V0-1 T1"]
    b = m["columnas"]["C0-a"]
    c = m["bases"]["B0-a"]

    sal.append(f"=== {caso} ===   (anaStruct | Pynite-motor | dif)")
    filas = [
        ("viga M_izq", Mv[0], a["M_izq"]),
        ("viga M_der", Mv[-1], a["M_der"]),
        ("viga M_campo", min(Mv), a["M_campo"]),
        ("viga V_izq", Qv[0], a["V_izq"]),
        ("viga V_der", Qv[-1], a["V_der"]),
        ("col M_inf", Mc[0], b["M_inf"]),
        ("col M_sup", Mc[-1], b["M_sup"]),
        ("col N", sum(Nc) / len(Nc), b["N"]),
        ("base Fx", R.Fx, c["Fx"]),
        ("base Fy", R.Fy, c["Fy"]),
        ("base Mz", R.Tz, c["Mz"]),
    ]
    for etq, ana, mot in filas:
        dif = abs(abs(ana) - abs(mot))
        bandera = "" if dif < 0.05 else "   <-- DIF"
        sal.append(f"   {etq:<14} {ana:>10.2f} | {mot:>10.2f} | {dif:>7.2f}{bandera}")

print("\n".join(sal))

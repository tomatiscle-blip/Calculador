import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parent
sys.path.insert(0, str(RAIZ))
from calc import rutas  # noqa: E402
from anastruct import SystemElements  # noqa: E402

FACT = {
    "servicio": {"D": 1.0, "L": 1.0, "W": 0.0},
    "1.4D": {"D": 1.4, "L": 0.0, "W": 0.0},
    "1.2D+1.6L": {"D": 1.2, "L": 1.6, "W": 0.0},
    "1.2D+0.5L+1.6W": {"D": 1.2, "L": 0.5, "W": 1.6},
    "0.9D+1.6W": {"D": 0.9, "L": 0.0, "W": 1.6},
}


def crear(portico, ss):
    for cid, col in portico["columnas"].items():
        x = col["x"]; y0 = col.get("nivel", 0.0); h = col["altura_m"]
        eid = ss.add_element(location=[[x, y0], [x, y0 + h]])
        col["eid"] = int(eid)
        col["Nodo_superior"] = ss.element_map[eid].node_2.id
    for viga_id, viga in portico["vigas"].items():
        piso = int(viga_id.split("-")[0][1:])
        cols = sorted((c for cid, c in portico["columnas"].items() if cid.startswith(f"C{piso}-")),
                      key=lambda c: c["x"])
        if not cols:
            continue
        y_viga = cols[0].get("nivel", 0.0) + cols[0]["altura_m"]
        # viento
        q_w = float(viga["tramos"][0]["cargas"].get("W_total", 0.0))
        w_col = q_w * cols[0]["altura_m"] / len(cols)
        for c in cols:
            c["carga_viento_kN"] = w_col
        for i, tramo in enumerate(viga["tramos"]):
            if tramo.get("es_voladizo"):
                continue
            if i >= len(cols) - 1:
                continue
            eid = ss.add_element(location=[[cols[i]["x"], y_viga], [cols[i + 1]["x"], y_viga]])
            tramo["subtramos"] = [{"eid": int(eid)}]
    for base in portico["bases"].values():
        nodo = ss.find_node_id([base["x"], 0])
        if nodo:
            ss.add_support_fixed(node_id=nodo)
    return ss


def cargas(portico, ss, caso):
    f = FACT[caso]
    for cid, col in portico["columnas"].items():
        if f["W"] > 0 and col.get("carga_viento_kN"):
            ss.point_load(node_id=col["Nodo_superior"], Fx=f["W"] * col["carga_viento_kN"])
    for viga in portico["vigas"].values():
        for tramo in viga["tramos"]:
            if "cargas" not in tramo:
                continue
            q = f["D"] * float(tramo["cargas"].get("D_total", 0.0)) + f["L"] * float(tramo["cargas"].get("L_total", 0.0))
            if abs(q) < 1e-9:
                continue
            for sub in tramo.get("subtramos", []):
                ss.q_load(element_id=sub["eid"], q=-q)


est = rutas.cargar_estructura()
p = est["Portico 1"]
sal = ["anaStruct FRESCO por combinación — V0-1 T1 (Pynite dio: serv 75.09/55.52/-61.23, 1.2D+1.6L 95.07/71.06/-78.28)"]
for caso in FACT:
    ss = SystemElements(mesh=50)
    crear(p, ss)
    cargas(p, ss, caso)
    ss.solve()
    eid = p["vigas"]["V0-1"]["tramos"][0]["subtramos"][0]["eid"]
    res = ss.get_element_results(element_id=eid, verbose=True)
    M = list(res["M"])
    sal.append(f"  {caso:16} M_izq={M[0]:8.2f} M_der={M[-1]:8.2f} minM={min(M):8.2f}")
print("\n".join(sal))
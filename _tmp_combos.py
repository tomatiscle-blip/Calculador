import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parent
sys.path.insert(0, str(RAIZ))
from calc import rutas  # noqa: E402
from Pynite import FEModel3D  # noqa: E402

E_MOD = 1.0; G_MOD = 1.0; EA = 15000.0; EI = 5000.0
A_SEC = EA / E_MOD; I_SEC = EI / E_MOD; J_SEC = EI / E_MOD

COMBOS = {
    "servicio": {"D": 1.0, "L": 1.0, "W": 0.0},
    "1.4D": {"D": 1.4, "L": 0.0, "W": 0.0},
    "1.2D+1.6L": {"D": 1.2, "L": 1.6, "W": 0.0},
    "1.2D+0.5L+1.6W": {"D": 1.2, "L": 0.5, "W": 1.6},
    "0.9D+1.6W": {"D": 0.9, "L": 0.0, "W": 1.6},
}


def k(x):
    return round(float(x), 4)


def nodo(x, y):
    return f"N{k(x)}_{k(y)}"


def construir(portico):
    m = FEModel3D()
    m.add_material("mat", E_MOD, G_MOD, 0.3, 0.0)
    m.add_section("sec", A_SEC, I_SEC, I_SEC, J_SEC)
    for cid, col in portico["columnas"].items():
        x = k(col["x"]); y0 = k(col.get("nivel", 0.0)); h = k(col["altura_m"])
        nb, nt = nodo(x, y0), nodo(x, y0 + h)
        if nb not in m.nodes:
            m.add_node(nb, x, y0, 0.0)
        if nt not in m.nodes:
            m.add_node(nt, x, y0 + h, 0.0)
        m.def_support(nb, False, False, True, True, True, False)
        m.def_support(nt, False, False, True, True, True, False)
        m.add_member(f"C-{cid}", nb, nt, "mat", "sec")
    for base in portico["bases"].values():
        nn = nodo(base["x"], 0.0)
        if nn in m.nodes:
            m.def_support(nn, True, True, True, True, True, True)
    for viga_id, viga in portico["vigas"].items():
        piso = int(viga_id.split("-")[0][1:])
        cols = sorted((c for cid, c in portico["columnas"].items() if cid.startswith(f"C{piso}-")),
                      key=lambda c: c["x"])
        if not cols:
            continue
        y_viga = k(cols[0].get("nivel", 0.0) + cols[0]["altura_m"])
        # viento: fuerza por columna en la punta
        q_w = float(viga["tramos"][0]["cargas"].get("W_total", 0.0))
        w_col = q_w * cols[0]["altura_m"] / len(cols)
        for c in cols:
            nn = nodo(c["x"], c.get("nivel", 0.0) + c["altura_m"])
            m.add_node_load(nn, "FX", w_col, case="W")
        for i, tramo in enumerate(viga["tramos"]):
            cargas = tramo.get("cargas", {})
            wD = float(cargas.get("D_total", 0.0)); wL = float(cargas.get("L_total", 0.0))
            if tramo.get("es_voladizo"):
                if "izq" in tramo["id"]:
                    tip_x = k(cols[0]["x"] - tramo["longitud_m"]); base_x = k(cols[0]["x"])
                else:
                    tip_x = k(cols[-1]["x"] + tramo["longitud_m"]); base_x = k(cols[-1]["x"])
                mid = f"V-{tramo['id']}"
                ntip, nbase = nodo(tip_x, y_viga), nodo(base_x, y_viga)
                if ntip not in m.nodes:
                    m.add_node(ntip, tip_x, y_viga, 0.0)
                    m.def_support(ntip, False, False, True, True, True, False)
                m.add_member(mid, ntip, nbase, "mat", "sec")
            else:
                if i >= len(cols) - 1:
                    continue
                mid = f"V-{tramo['id']}"
                m.add_member(mid, nodo(cols[i]["x"], y_viga), nodo(cols[i + 1]["x"], y_viga), "mat", "sec")
            if abs(wD) > 1e-9:
                m.add_member_dist_load(mid, "FY", -wD, -wD, case="D")
            if abs(wL) > 1e-9:
                m.add_member_dist_load(mid, "FY", -wL, -wL, case="L")
            for cp in tramo.get("cargas_puntuales", []):
                nn = nodo(cp["coordenadas"]["x"], cp["coordenadas"]["y"])
                if nn in m.nodes:
                    m.add_node_load(nn, "FY", -float(cp["valor_kN"]), case="D")
    for nombre, f in COMBOS.items():
        m.add_load_combo(nombre, {"D": f["D"], "L": f["L"], "W": f["W"]})
    return m


est = rutas.cargar_estructura()
p = est["Portico 1"]
m = construir(p)
m.analyze()

sal = []
viga = m.members["V-V0-1 T1"]
sal.append("VIGA V0-1 T1  (guardado: Mu_kNm_apoyo_izq=-46.52, Mu_kNm_apoyo_der=-92.63, campo=66.36; Vu_izq=-101.82 Vu_der=120.83)")
for nombre in COMBOS:
    sal.append(f"  {nombre:16} M_izq={viga.moment('Mz',0,nombre):8.2f} M_der={viga.moment('Mz',viga.L(),nombre):8.2f} "
               f"minM={viga.min_moment('Mz',nombre):8.2f} V_izq={viga.shear('Fy',0,nombre):8.2f} V_der={viga.shear('Fy',viga.L(),nombre):8.2f}")

col = m.members["C-C0-a"]
sal.append("COL C0-a  (guardado: Mu_kNm_inf=26.46, Mu_kNm_sup=11.59, P_kN=-163.85)")
for nombre in COMBOS:
    sal.append(f"  {nombre:16} M_inf={col.moment('Mz',0,nombre):8.2f} M_sup={col.moment('Mz',col.L(),nombre):8.2f} N={col.axial(0,nombre):8.2f}")

sal.append("BASES (guardado B0-a: Fx_kN=-4.8 Fy_kN=-163.85 Tz_kNm=26.46; Fx_servicio=21.62 Fy_servicio=-167.91 Tz_servicio=-25.91)")
for nombre in COMBOS:
    n = m.nodes[nodo(0.0, 0.0)]
    sal.append(f"  {nombre:16} Fx={float(n.RxnFX[nombre]):8.2f} Fy={float(n.RxnFY[nombre]):8.2f} Tz={float(n.RxnMZ[nombre]):8.2f}")

print("\n".join(sal))
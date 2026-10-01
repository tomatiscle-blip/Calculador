import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parent
sys.path.insert(0, str(RAIZ))
from calc import rutas  # noqa: E402
from Pynite import FEModel3D  # noqa: E402

# Misma rigidez que los valores por defecto de anaStruct en P01
EA, EI, E_MOD = 15000.0, 5000.0, 1.0
G_MOD = 1.0
A_SEC, I_SEC, J_SEC = EA / E_MOD, EI / E_MOD, EI / E_MOD


def k(x):
    return round(float(x), 4)


def construir(portico):
    m = FEModel3D()
    m.add_material("mat", E_MOD, G_MOD, 0.3, 0.0)
    m.add_section("sec", A_SEC, I_SEC, I_SEC, J_SEC)

    col_top = {}
    tops = []
    for cid, col in portico["columnas"].items():
        x = col["x"]; y0 = col.get("nivel", 0.0); h = col["altura_m"]
        nb, nt = f"B@{k(x)}", f"T@{k(x)}"
        m.add_node(nb, x, y0, 0.0)
        m.add_node(nt, x, y0 + h, 0.0)
        m.add_member(f"C-{cid}", nb, nt, "mat", "sec")
        m.def_support(nb, True, True, True, True, True, True)
        m.def_support(nt, False, False, True, True, True, False)
        col_top[k(x)] = nt
        tops.append(round(y0 + h, 4))
    y_beam = max(tops) if tops else 0.0

    for viga in portico["vigas"].values():
        for tramo in viga["tramos"]:
            cargas = tramo.get("cargas", {})
            w = float(cargas.get("D_total", 0.0)) + float(cargas.get("L_total", 0.0))
            xi, xf = k(tramo["x_inicio"]), k(tramo["x_fin"])
            mid = f"M-{tramo['id']}"
            if tramo.get("es_voladizo"):
                if xi in col_top and xf not in col_top:
                    base, tip = xi, xf
                elif xf in col_top and xi not in col_top:
                    base, tip = xf, xi
                else:
                    base, tip = xi, xf
                tn = f"V@{tip}"
                m.add_node(tn, tip, y_beam, 0.0)
                m.def_support(tn, False, False, True, True, True, False)
                m.add_member(mid, col_top[base], tn, "mat", "sec")
            else:
                m.add_member(mid, col_top[xi], col_top[xf], "mat", "sec")
            m.add_member_dist_load(mid, "FY", -w, -w)
            for c in tramo.get("cargas_puntuales", []):
                cx = k(c["coordenadas"]["x"])
                nn = f"V@{cx}"
                if nn in m.nodes:
                    m.add_node_load(nn, "FY", -float(c["valor_kN"]))
    return m


portico = rutas.cargar_estructura()["Portico 1"]
m = construir(portico)
m.analyze()
print("Nodos:", list(m.nodes))
for nb in ["B@0.0", "B@4.85"]:
    n = m.nodes[nb]
    print(f"{nb}: RxnFX={n.RxnFX}  RxnFY={n.RxnFY}  RxnMZ={n.RxnMZ}")
print("anaStruct (guardado): B0-a Fy=-167.91 Mz=-25.91 | B0-b Fy=-100.17 Mz=11.49")

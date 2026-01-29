import json
import os
import numpy as np
import matplotlib.pyplot as plt
# Ruta al archivo
base_dir = os.path.dirname(os.path.abspath(__file__))
archivo = os.path.join(base_dir, "datos", "estructura.json")

with open(archivo, "r") as f:
    estructura = json.load(f)

# -------------------------------
# 1. Análisis estructural de cada pórtico
# -------------------------------
for portico_id, portico in estructura.items():
    print(f"\n=== {portico_id} ===")

    # Columnas
    for col_id, col in portico["columnas"].items():
        print(f"Columna {col_id}: altura = {col['altura_m']} m")

    # Vigas y tramos
    for viga_id, viga in portico["vigas"].items():
        print(f"\nViga {viga_id}:")
        for tramo in viga["tramos"]:
            print(f"  Tramo {tramo['id']} - L = {tramo['longitud_m']} m")
            if "cargas" in tramo:
                print("    Cargas distribuidas:")
                for k,v in tramo["cargas"].items():
                    if k != "combinaciones":
                        print(f"      {k}: {v}")
                print("    Combinaciones:")
                for comb, val in tramo["cargas"]["combinaciones"].items():
                    print(f"      {comb}: {val}")

    # Cargas puntuales
    if "cargas_puntuales" in portico:
        print("\nCargas puntuales:")
        for carga in portico["cargas_puntuales"]:
            print(f"  Piso {carga['piso']} - X={carga['coordenadas']['x']} m, "
                  f"Y={carga['coordenadas']['y']} m, "
                  f"Valor={carga['valor_kN']} kN, Tipo={carga['tipo']}")

# -------------------------------
# 2. Dimensionado de pórtico de hormigón armado
# ------------------------------

from anastruct import SystemElements
import math, json, os


# Funciones auxiliares
def aplicar_cargas(portico, ss, caso="servicio"):
    for viga in portico["vigas"].values():
        for tramo in viga["tramos"]:
            if "cargas" in tramo:
                # obtener valores de carga
                if caso == "servicio":
                    q_d = float(tramo["cargas"].get("D_total", 0.0))
                    q_l = float(tramo["cargas"].get("L_total", 0.0))
                    # NO sumar viento aquí
                    q_total = q_d + q_l
                else:
                    q_total = float(tramo["cargas"]["combinaciones"].get(caso, 0.0))

                # aplicar en cada subtramo
                for sub in tramo.get("subtramos", []):
                    eid = sub.get("eid")
                    if q_total > 0 and eid:
                        ss.q_load(element_id=eid, q=-q_total)

def guardar_resultados_servicio(portico, ss):
    # vigas
    for viga in portico["vigas"].values():
        for tramo in viga["tramos"]:
            resultados_M = []
            resultados_V = []

            # recorrer todos los subtramos si existen
            for sub in tramo.get("subtramos", []):
                eid = sub.get("eid")
                if not eid:
                    continue
                res = ss.get_element_results(element_id=eid, verbose=True)
                M = list(res.get("M")) or []
                V = list(res.get("Q")) or []
                resultados_M.extend(M)
                resultados_V.extend(V)

            # consolidar resultados en el tramo padre
            if resultados_M:
                tramo.update({
                    "Mu_servicio_izq": round(resultados_M[0], 2),
                    "Mu_servicio_der": round(resultados_M[-1], 2),
                    "Mu_servicio_campo": round(min(resultados_M), 2)
                })
            else:
                tramo.update({
                    "Mu_servicio_izq": 0.0,
                    "Mu_servicio_der": 0.0,
                    "Mu_servicio_campo": 0.0
                })

            if resultados_V:
                tramo.update({
                    "Vu_servicio_izq": round(resultados_V[0], 2),
                    "Vu_servicio_der": round(resultados_V[-1], 2)
                })
            else:
                tramo.update({
                    "Vu_servicio_izq": 0.0,
                    "Vu_servicio_der": 0.0
                })

    # columnas
    for col in portico["columnas"].values():
        eid = col.get("eid")
        if not eid:
            continue
        res = ss.get_element_results(element_id=eid, verbose=True)

        M = list(res.get("M")) or []
        N = list(res.get("N")) or []

        col.update({
            "Mu_servicio_inf": round(M[0], 2) if M else 0.0,
            "Mu_servicio_sup": round(M[-1], 2) if M else 0.0,
            "N_servicio": round(np.mean(N), 2) if N else 0.0
        })

    # bases
    for base_id, base in portico["bases"].items():
        node_id = ss.find_node_id([base["x"], 0])
        if node_id in ss.reaction_forces:
            reaction = ss.reaction_forces[node_id]
            Fx = round(reaction.Fx, 2)
            Fy = round(reaction.Fy, 2)
            Tz = round(reaction.Tz, 2)

            base.update({
                "Fx_servicio": Fx,
                "Fy_servicio": Fy,
                "Tz_servicio": Tz
            })

def guardar_resultados(portico, ss, caso):
    # vigas
    for viga in portico["vigas"].values():
        for tramo in viga["tramos"]:
            eid = tramo.get("eid")
            if not eid:
                continue

            res = ss.get_element_results(element_id=eid, verbose=True)

            M = list(res.get("M", [])) if res.get("M") is not None else []
            Mu_izq = round(M[0], 2) if M else 0.0
            Mu_der = round(M[-1], 2) if M else 0.0
            Mu_campo = round(min(M, key=lambda x: abs(x)), 2) if M else 0.0

            V = list(res.get("Q", [])) if res.get("Q") is not None else []
            Vu_izq = round(V[0], 2) if V else 0.0
            Vu_der = round(V[-1], 2) if V else 0.0

            N1 = res.get("axial_1", 0.0)
            N2 = res.get("axial_2", 0.0)

            tramo.update({
                f"Mu_{caso}_izq": Mu_izq,
                f"Mu_{caso}_der": Mu_der,
                f"Mu_{caso}_campo": Mu_campo,
                f"Vu_{caso}_izq": Vu_izq,
                f"Vu_{caso}_der": Vu_der,
                f"N_{caso}_izq": N1,
                f"N_{caso}_der": N2
            })

    # columnas
    for col in portico["columnas"].values():
        eid = col.get("eid")
        if not eid:
            continue
        res = ss.get_element_results(element_id=eid, verbose=True)

        M = list(res.get("M")) or []
        N = list(res.get("N")) or []

        col.update({
            f"Mu_{caso}_inf": round(M[0], 2) if M else 0.0,
            f"Mu_{caso}_sup": round(M[-1], 2) if M else 0.0,
            f"N_{caso}": round(np.mean(N), 2) if N else 0.0
        })

    # bases
    for base_id, base in portico["bases"].items():
        node_id = ss.find_node_id([base["x"], 0])
        if node_id in ss.reaction_forces:
            reaction = ss.reaction_forces[node_id]
            Fx = round(reaction.Fx, 2)
            Fy = round(reaction.Fy, 2)
            Tz = round(reaction.Tz, 2)

            base.update({
                f"Fx_{caso}": Fx,
                f"Fy_{caso}": Fy,
                f"Tz_{caso}": Tz
            })

def seleccionar_maximos(portico):
    # recorrer vigas
    for viga in portico["vigas"].values():
        for tramo in viga["tramos"]:
            mu_izq = [v for k,v in tramo.items() if k.startswith("Mu_") and "_izq" in k and "servicio" not in k and v is not None]
            mu_der = [v for k,v in tramo.items() if k.startswith("Mu_") and "_der" in k and "servicio" not in k and v is not None]
            mu_campo = [v for k,v in tramo.items() if k.startswith("Mu_") and "_campo" in k and "servicio" not in k and v is not None]
            vu_izq = [v for k,v in tramo.items() if k.startswith("Vu_") and "_izq" in k and "servicio" not in k and v is not None]
            vu_der = [v for k,v in tramo.items() if k.startswith("Vu_") and "_der" in k and "servicio" not in k and v is not None]

            tramo["Mu_kNm_apoyo_izq"] = max(mu_izq, key=lambda x: abs(x), default=0)
            tramo["Mu_kNm_apoyo_der"] = max(mu_der, key=lambda x: abs(x), default=0)
            tramo["Mu_kNm_campo"] = max(mu_campo, key=lambda x: abs(x), default=0)
            tramo["Vu_kN_izq"] = max(vu_izq, key=lambda x: abs(x), default=0)
            tramo["Vu_kN_der"] = max(vu_der, key=lambda x: abs(x), default=0)

    # recorrer columnas
    for col in portico["columnas"].values():
        momentos = [v for k,v in col.items() if k.startswith("Mu_") and "servicio" not in k and v is not None]
        axiales = [v for k,v in col.items() if k.startswith("N_") and "servicio" not in k and v is not None]

        col["Mu_kNm_max"] = max(momentos, key=lambda x: abs(x), default=0)
        col["P_kN_max"] = max(axiales, key=lambda x: abs(x), default=0)

def crear_portico(portico, ss):
    """
    Crea el modelo de pórtico en Anastruct:
    - Columnas
    - Vigas entre columnas
    - Voladizos
    - Apoyos
    - Cargas puntuales y viento
    """

    # -------------------------------
    # Columnas
    # -------------------------------
    for col_id, col in portico["columnas"].items():
        x = col["x"]
        y0 = col.get("nivel", 0.0)
        h = col["altura_m"]

        eid = ss.add_element(location=[[x, y0], [x, y0 + h]])
        col["eid"] = int(eid)
        col["nodos"] = [
            ss.element_map[eid].node_1.id,
            ss.element_map[eid].node_2.id
        ]
        col["Nodo_inferior"] = ss.element_map[eid].node_1.id
        col["Nodo_superior"] = ss.element_map[eid].node_2.id

    # -------------------------------
    # Vigas entre columnas
    # -------------------------------
    for viga_id, viga in portico["vigas"].items():
        piso_viga = int(viga_id.split("-")[0][1:])  # ej. V1-1 → piso 1
        columnas_piso = [c for cid, c in portico["columnas"].items() if cid.startswith(f"C{piso_viga}-")]
        columnas_ordenadas = sorted(columnas_piso, key=lambda c: c["x"])

        if not columnas_ordenadas:
            continue

        y_viga = columnas_ordenadas[0]["nivel"] + columnas_ordenadas[0]["altura_m"]

        for i, tramo in enumerate(viga["tramos"]):
            if tramo.get("es_voladizo", False):
                continue

            if i >= len(columnas_ordenadas) - 1:
                print(f"[WARN] No hay suficientes columnas para tramo {tramo['id']}, se ignora")
                continue

            x0 = columnas_ordenadas[i]["x"]
            x1 = columnas_ordenadas[i + 1]["x"]

            if abs(x1 - x0) < 1e-6:
                print(f"[WARN] Tramo con longitud cero entre columnas x0={x0} x1={x1}, tramo {tramo['id']} ignorado")
                continue

            # buscar cargas puntuales en este tramo
            cargas_tramo = tramo.get("cargas_puntuales", [])

            # nodos intermedios = extremos + posiciones de cargas
            nodos_x = [x0] + [c["coordenadas"]["x"] for c in cargas_tramo] + [x1]
            nodos_x = sorted(set(nodos_x))

            # discretizar la viga en subtramos
            tramo["subtramos"] = []
            for j in range(len(nodos_x)-1):
                xi, xj = nodos_x[j], nodos_x[j+1]
                eid = ss.add_element(location=[[xi, y_viga], [xj, y_viga]])
                tramo["subtramos"].append({
                    "eid": int(eid),
                    "nodos": [ss.element_map[eid].node_1.id,
                          ss.element_map[eid].node_2.id]
                })

            # aplicar cargas puntuales en los nodos intermedios
            for carga in cargas_tramo:
                nodo = ss.find_node_id([carga["coordenadas"]["x"], carga["coordenadas"]["y"]])
                if nodo:
                    ss.point_load(node_id=nodo, Fy=-float(carga["valor_kN"]))

    # -------------------------------
    # Voladizos
    # -------------------------------
    for viga_id, viga in portico["vigas"].items():
        piso_viga = int(viga_id.split("-")[0][1:])
        columnas_piso = [c for cid, c in portico["columnas"].items() if cid.startswith(f"C{piso_viga}-")]
        columnas_ordenadas = sorted(columnas_piso, key=lambda c: c["x"])

        if not columnas_ordenadas:
            continue

        y_viga = columnas_ordenadas[0]["nivel"] + columnas_ordenadas[0]["altura_m"]

        for tramo in viga["tramos"]:
            if not tramo.get("es_voladizo", False):
                continue

            if "izq" in tramo["id"]:
                x0 = columnas_ordenadas[0]["x"]
                x1 = x0 - tramo["longitud_m"]
            elif "der" in tramo["id"]:
                x0 = columnas_ordenadas[-1]["x"]
                x1 = x0 + tramo["longitud_m"]
            else:
                print(f"[WARN] Voladizo sin dirección en {tramo['id']}, se ignora")
                continue

            if abs(x1 - x0) < 1e-6:
                print(f"[WARN] Voladizo con longitud cero {tramo['id']}, se ignora")
                continue

            # buscar cargas puntuales en este voladizo
            cargas_tramo = tramo.get("cargas_puntuales", [])

            nodos_x = [x0] + [c["coordenadas"]["x"] for c in cargas_tramo] + [x1]
            nodos_x = sorted(set(nodos_x))

            tramo["subtramos"] = []
            for j in range(len(nodos_x)-1):
                xi, xj = nodos_x[j], nodos_x[j+1]
                eid = ss.add_element(location=[[xi, y_viga], [xj, y_viga]])
                tramo["subtramos"].append({
                    "eid": int(eid),
                    "nodos": [ss.element_map[eid].node_1.id,
                          ss.element_map[eid].node_2.id]
                })

            # aplicar cargas puntuales en los nodos intermedios
            for carga in cargas_tramo:
                nodo = ss.find_node_id([carga["coordenadas"]["x"], carga["coordenadas"]["y"]])
                if nodo:
                    ss.point_load(node_id=nodo, Fy=-float(carga["valor_kN"]))
    # -------------------------------
    # Apoyos
    # -------------------------------
    for base_id, base in portico["bases"].items():
        x = base["x"]
        nodo = ss.find_node_id([x, 0])
        if nodo is None:
            print(f"[WARN] No se encontró nodo para apoyo en x={x}")
            continue

        if base["tipo"] == "empotramiento":
            ss.add_support_fixed(node_id=nodo)
        elif base["tipo"] == "articulado":
            ss.add_support_hinged(node_id=nodo)

    # -------------------------------
    # Viento
    # -------------------------------
    for viga_id, viga in portico["vigas"].items():
        total_w = 0.0
        for tramo in viga.get("tramos", []):
            if "cargas" in tramo and "W_total" in tramo["cargas"]:
                L_tramo = tramo["longitud_m"]
                q_w = float(tramo["cargas"]["W_total"])
                total_w += q_w * L_tramo

        piso_viga = int(viga_id.split("-")[0][1:])
        columnas_piso = [c for cid, c in portico["columnas"].items() if cid.startswith(f"C{piso_viga}-")]

        if columnas_piso and total_w > 0:
            carga_por_columna = total_w / len(columnas_piso)
            for col in columnas_piso:
                nodo_sup = col.get("Nodo_superior")
                if nodo_sup:
                    ss.point_load(node_id=nodo_sup, Fx=carga_por_columna)
                    col["carga_viento_kN"] = carga_por_columna

    return ss

# Funciones auxiliares para retrolimentación de VIGAS.

def obtener_Mx_global(ss, eid_list):
    xg, Mg = [], []
    x_offset = 0.0

    for eid in eid_list:
        res = ss.get_element_results(element_id=eid, verbose=True)
        if not res or res.get("M") is None:
            continue

        M = list(res["M"])
        xloc = res.get("x")

        if xloc is None:
            L = ss.element_map[eid].l
            xloc = np.linspace(0, L, len(M))

        xloc = [xi + x_offset for xi in xloc]

        xg.extend(xloc)
        Mg.extend(M)

        x_offset = xg[-1]

    return xg, Mg

def obtener_Q_global(ss, eid_list, devolver_x=False):
    xg, Qg = [], []
    x_offset = 0.0

    for eid in eid_list:
        res = ss.get_element_results(element_id=eid, verbose=True)
        if not res or res.get("Q") is None:
            continue

        Q = list(res["Q"])
        xloc = res.get("x")

        if xloc is None:
            L = ss.element_map[eid].l
            xloc = np.linspace(0, L, len(Q))

        xloc = [xi + x_offset for xi in xloc]

        xg.extend(xloc)
        Qg.extend(Q)

        x_offset = xg[-1]

    if devolver_x:
        return xg, Qg
    return Qg

def puntos_inflexion(xg, Mg):
    pis = []
    for i in range(len(Mg) - 1):
        if Mg[i] * Mg[i+1] < 0:
            x0, x1 = xg[i], xg[i+1]
            m0, m1 = Mg[i], Mg[i+1]
            xi = x0 - m0 * (x1 - x0) / (m1 - m0)
            pis.append(round(float(xi), 3))
    return pis

def extremos_voladizo(Mg, Qg, tol=1e-6):
    if abs(Mg[0]) > tol:
        return Mg[0], Qg[0]
    return Mg[-1], Qg[-1]

def armar_geometria(portico):
    """
    Devuelve una lista de líneas (x1, y1, x2, y2) de vigas y columnas
    para poder plotear la estructura completa.
    """
    lineas = []

    # Columnas
    for col in portico.get("columnas", {}).values():
        x = col.get("x", 0.0)
        y_inf = 0.0  # base del pórtico
        y_sup = col.get("altura_m", 3.0)
        lineas.append((x, y_inf, x, y_sup))

    # Vigas
    for viga in portico.get("vigas", {}).values():
        for tramo in viga.get("tramos", []):
            x1 = tramo.get("x_inicio", 0.0)
            x2 = tramo.get("x_fin", 0.0)
            # altura de la viga: asumimos top de columnas
            y = max([col.get("altura_m", 3.0) for col in portico.get("columnas", {}).values()] + [3.0])
            lineas.append((x1, y, x2, y))

    return lineas



#-------------------------------
# Menu: Análisis estructural y dimensionado
#-------------------------------
# Mostrar lista de pórticos disponibles
print("Pórticos disponibles en estructura.json:")
for i, pid in enumerate(estructura.keys(), start=1):
    print(f"{i}. {pid}")

# Preguntar cuál procesar
opcion = int(input("¿Qué pórtico querés dimensionar primero? (número): "))
portico_id = list(estructura.keys())[opcion - 1]
portico = estructura[portico_id]


ya_calculado = any("Mu_kNm_inf" in col for col in portico["columnas"].values())
if ya_calculado:
    print(f"\n=== {portico_id} ya está calculado, se omite ===")
else:
    print(f"\n=== Dimensionando {portico_id} ===")
    from anastruct import SystemElements

    # -------------------------------
    # Caso de servicio (D+L + W)
    # -------------------------------
    caso = "servicio"
    ss = SystemElements(mesh=50)
    ss = crear_portico(portico, ss)
    aplicar_cargas(portico, ss, caso)
    ss.solve()
    guardar_resultados_servicio(portico, ss)

    # -------------------------------
    # Casos mayorados
    # -------------------------------
    for comb in list(portico["vigas"].values())[0]["tramos"][0]["cargas"]["combinaciones"].keys():
        caso = comb
        ss = SystemElements(mesh=50)
        ss = crear_portico(portico, ss)
        aplicar_cargas(portico, ss, caso)
        ss.solve()
        guardar_resultados(portico, ss, caso)

    # -------------------------------
    # Selección de máximos
    # -------------------------------
    seleccionar_maximos(portico)

    # -------------------------------
    # Visualización de varias combinaciones en una sola gráfica
    # -------------------------------

    combinaciones = ["1.2D+1.6L", "servicio"] #"servicio", "1.2D+1.6L", "1.4D", "1.2D+1.0L+0.5W", "1.2D+1.0L+1.0W", "0.9D+1.0W"

    for comb in combinaciones:
        ss = SystemElements(mesh=50)
        ss = crear_portico(portico, ss)
        aplicar_cargas(portico, ss, comb)
        ss.solve()
        print(f"=== Diagrama para {comb} ===")
        ss.show_structure()
        ss.show_bending_moment()
        ss.show_shear_force()
        ss.show_axial_force()

# -------------------------------
# Retroalimentar VIGAS (servicio y mayorados)
# -------------------------------
for viga_id, viga in portico["vigas"].items():
    nuevos_tramos = []

    for tramo in viga["tramos"]:

        # =====================================================
        # VOLADIZOS → solo M y Q de empotramiento
        # =====================================================
        if tramo.get("es_voladizo", False):

            eid_list = [
                sub["eid"] for sub in tramo.get("subtramos", [])
                if sub.get("eid") is not None
            ]

            if not eid_list:
                nuevos_tramos.append(tramo)
                continue

            # resultados globales (arrays continuos)
            xg, Mg = obtener_Mx_global(ss, eid_list)
            Vg = obtener_Q_global(ss, eid_list)

            # empotramiento real
            M_emp, V_emp = extremos_voladizo(Mg, Vg)

            if caso == "servicio":
                nuevo_tramo = {
                    **tramo,
                    "Mu_servicio": -abs(float(M_emp)),   # ← FORZADO
                    "Vu_servicio": float(V_emp),
                }
            else:
                nuevo_tramo = {
                    **tramo,
                    "Mu_kNm": -abs(float(M_emp)),         # ← FORZADO
                    "Vu_kN_emp": float(V_emp),
                }

            nuevos_tramos.append(nuevo_tramo)
            continue  # ← CLAVE: no cae al tramo normal

        # =====================================================
        # TRAMOS NORMALES
        # =====================================================
        eid_list = [
            sub["eid"] for sub in tramo.get("subtramos", [])
            if sub.get("eid") is not None
        ]

        if not eid_list:
            nuevos_tramos.append(tramo)
            continue

        # resultados globales
        xg, Mg = obtener_Mx_global(ss, eid_list)
        _, Qg = obtener_Q_global(ss, eid_list, devolver_x=True)

        if caso == "servicio":
            nuevo_tramo = {
                **tramo,
                "Mu_servicio_izq": -abs(float(Mg[0])),
                "Mu_servicio_der":-abs(float(Mg[-1])),
                "Mu_servicio_campo": abs(float(min(Mg))),
                "Vu_servicio_izq": float(Qg[0]),
                "Vu_servicio_der": float(Qg[-1]),
            }
        else:
            nuevo_tramo = {
                **tramo,
                "Mu_kNm_apoyo_izq": -abs(float(Mg[0])),
                "Mu_kNm_apoyo_der": -abs(float(Mg[-1])),
                "Mu_kNm_campo": abs(float(min(Mg))),
                "Vu_kN_izq": float(Qg[0]),
                "Vu_kN_der": float(Qg[-1]),
                "puntos_inflexion_m": puntos_inflexion(xg, Mg),
            }

        nuevos_tramos.append(nuevo_tramo)

    viga["tramos"] = nuevos_tramos

# -------------------------------
# Retroalimentar COLUMNAS (servicio y mayorados)
# -------------------------------
node_results = ss.get_node_displacements()

for col_id, col in portico["columnas"].items():
    eid = col.get("eid")
    elem = ss.get_element_results(element_id=eid, verbose=True) if eid else {}
    L = elem.get("length", col.get("altura_m", 3.0))

    Mu_inf = Mu_sup = 0.0
    distancia_campo_cero = None

    # nodos de la columna (del JSON)
    nodos_col = col.get("nodos", [])
    nodo_inf = col.get("Nodo_inferior")
    nodo_sup = col.get("Nodo_superior")

    # momentos
    if "M" in elem and elem["M"] is not None and len(elem["M"]) > 0:
        M = list(elem["M"])

        y_inf = node_results[nodo_inf]["uy"] * 1000 if nodo_inf < len(node_results) else 0.0
        y_sup = node_results[nodo_sup]["uy"] * 1000 if nodo_sup < len(node_results) else 0.0

        if y_inf <= y_sup:
            Mu_inf, Mu_sup = round(M[0], 2), round(M[-1], 2)
        else:
            Mu_inf, Mu_sup = round(M[-1], 2), round(M[0], 2)

        for i in range(len(M) - 1):
            if M[i] * M[i+1] < 0:
                frac = abs(M[i]) / (abs(M[i]) + abs(M[i+1]))
                distancia_campo_cero = round(L * (i + frac) / (len(M)-1), 3)
                break

    # axiales
    N_val = 0.0
    N1 = N2 = 0.0
    if "N" in elem and elem["N"] is not None and len(elem["N"]) > 0:
        N = list(elem["N"])
        N_val = round(np.mean(N), 2)
        N1, N2 = round(N[0], 2), round(N[-1], 2)


    # deriva horizontal
    altura = col.get("altura_m", L)
    ux_sup = node_results[nodo_sup]["ux"] * 1000 if nodo_sup < len(node_results) else 0.0
    deriva_h = round(ux_sup / (altura * 1000), 6)

    nuevo_paquete = {
        "x": col.get("x", 0.0),
        "altura_m": col.get("altura_m", L),
        "eid": eid,
        "nodos": nodos_col,
        "Nodo_inferior": nodo_inf,
        "Nodo_superior": nodo_sup,
        "Mu_kNm_inf": Mu_inf,
        "Mu_kNm_sup": Mu_sup,
        "P_kN": N_val,
        "Mu_servicio_inf": Mu_inf if caso == "servicio" else col.get("Mu_servicio_inf"),
        "Mu_servicio_sup": Mu_sup if caso == "servicio" else col.get("Mu_servicio_sup"),
        "P_servicio": N_val if caso == "servicio" else col.get("N_servicio"),
        "deriva_h": deriva_h,
        "distancia_campo_cero_m": distancia_campo_cero
    }

    portico["columnas"][col_id] = nuevo_paquete
# -------------------------------
# Retroalimentar BASES
# -------------------------------
for base_id, base in portico["bases"].items():
    node_id = ss.find_node_id([base["x"], 0])
    if node_id in ss.reaction_forces:
        reaction = ss.reaction_forces[node_id]

        Fx = round(reaction.Fx, 2)
        Fy = round(reaction.Fy, 2)
        Tz = round(reaction.Tz, 2)

        nuevo_paquete = {
            "x": base.get("x", 0.0),
            "tipo": base.get("tipo", "empotramiento"),
            # servicio
            "Fx_servicio": Fx if caso == "servicio" else base.get("Fx_servicio", 0.0),
            "Fy_servicio": Fy if caso == "servicio" else base.get("Fy_servicio", 0.0),
            "Tz_servicio": Tz if caso == "servicio" else base.get("Tz_servicio", 0.0),
            # mayorado
            "Fx_kN": Fx if caso != "servicio" else base.get("Fx_kN", 0.0),
            "Fy_kN": Fy if caso != "servicio" else base.get("Fy_kN", 0.0),
            "Tz_kNm": Tz if caso != "servicio" else base.get("Tz_kNm", 0.0)
        }


        portico["bases"][base_id] = nuevo_paquete

# -------------------------------
# Guardar retroalimentado
# -------------------------------
estructura[portico_id] = portico


# Guardar todo el JSON
with open(archivo, "w") as f:
    json.dump(estructura, f, indent=2)
print(f"\nDatos del pórtico '{portico_id}' actualizados en '{archivo}'")
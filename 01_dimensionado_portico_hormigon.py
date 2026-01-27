import json
import os
import numpy as np

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
# -------------------------------

from anastruct import SystemElements
import math, json, os

# Funciones auxiliares
def aplicar_cargas(portico, ss, caso="servicio"):
    for viga in portico["vigas"].values():
        for tramo in viga["tramos"]:
            if "cargas" in tramo:
                eid = tramo.get("eid")
                if eid:
                    if caso == "servicio":
                        D = float(tramo["cargas"].get("D_total", 0.0))
                        L = float(tramo["cargas"].get("L_total", 0.0))
                        ss.q_load(element_id=eid, q=-(D + L))
                    else:
                        q_comb = float(tramo["cargas"]["combinaciones"][caso])
                        ss.q_load(element_id=eid, q=-q_comb)

def guardar_resultados_servicio(portico, ss):
    # vigas
    for viga in portico["vigas"].values():
        for tramo in viga["tramos"]:
            eid = tramo.get("eid")
            if not eid:
                continue
            res = ss.get_element_results(element_id=eid, verbose=True)

            M = list(res.get("M")) or []
            V = list(res.get("Q")) or []

            tramo.update({
                "Mu_servicio_izq": round(M[0],2) if M else 0.0,
                "Mu_servicio_der": round(M[-1],2) if M else 0.0,
                "Mu_servicio_campo": round(min(M),2) if M else 0.0,
                "Vu_servicio_izq": round(V[0],2) if V else 0.0,
                "Vu_servicio_der": round(V[-1],2) if V else 0.0
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
            "Mu_servicio_inf": round(M[0],2) if M else 0.0,
            "Mu_servicio_sup": round(M[-1],2) if M else 0.0,
            "P_servicio": round(np.mean(N),2) if N else 0.0
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
    for viga in portico["vigas"].values():
        for tramo in viga["tramos"]:
            eid = tramo.get("eid")
            if not eid:
                continue

            # pedir resultados detallados (sin combination)
            res = ss.get_element_results(element_id=eid, verbose=True)

            # momentos
            M = list(res.get("M", [])) if res.get("M") is not None else []
            Mu_izq = round(M[0], 2) if M else 0.0
            Mu_der = round(M[-1], 2) if M else 0.0
            Mu_campo = round(min(M, key=lambda x: abs(x)), 2) if M else 0.0

            # cortantes
            V = list(res.get("Q", [])) if res.get("Q") is not None else []
            Vu_izq = round(V[0], 2) if V else 0.0
            Vu_der = round(V[-1], 2) if V else 0.0

            # axiales
            N1 = res.get("axial_1", 0.0)
            N2 = res.get("axial_2", 0.0)

            # guardar con prefijo del caso
            tramo.update({
                f"Mu_{caso}_izq": Mu_izq,
                f"Mu_{caso}_der": Mu_der,
                f"Mu_{caso}_campo": Mu_campo,
                f"Vu_{caso}_izq": Vu_izq,
                f"Vu_{caso}_der": Vu_der,
                f"N_{caso}_izq": N1,
                f"N_{caso}_der": N2
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
    # Crear columnas
    for col_id, col in portico["columnas"].items():
        x = col["x"]
        h = col["altura_m"]
        eid = ss.add_element(location=[[x, 0], [x, h]])
        col["eid"] = int(eid)
        col["nodos"] = [
            ss.element_map[eid].node_1.id,
            ss.element_map[eid].node_2.id
        ]

    # Crear vigas
    for viga_id, viga in portico["vigas"].items():
        y_viga = list(portico["columnas"].values())[0]["altura_m"]
        columnas_ordenadas = sorted(portico["columnas"].values(), key=lambda c: c["x"])
        # tramos entre columnas
        for i in range(len(columnas_ordenadas) - 1):
            x0 = columnas_ordenadas[i]["x"]
            x1 = columnas_ordenadas[i + 1]["x"]
            if i < len(viga["tramos"]):
                tramo = viga["tramos"][i]
                eid = ss.add_element(location=[[x0, y_viga], [x1, y_viga]])
                tramo["eid"] = int(eid)
                tramo["nodos"] = [
                    ss.element_map[eid].node_1.id,
                    ss.element_map[eid].node_2.id
                ]
        # voladizos
        for tramo in viga["tramos"]:
            if tramo.get("es_voladizo", False):
                if "izq" in tramo["id"]:
                    x0 = columnas_ordenadas[0]["x"]
                    x1 = x0 - tramo["longitud_m"]
                elif "der" in tramo["id"]:
                    x0 = columnas_ordenadas[-1]["x"]
                    x1 = x0 + tramo["longitud_m"]
                else:
                    continue
                eid = ss.add_element(location=[[x0, y_viga], [x1, y_viga]])
                tramo["eid"] = int(eid)
                tramo["nodos"] = [
                    ss.element_map[eid].node_1.id,
                    ss.element_map[eid].node_2.id
                ]

    # Crear apoyos
    for base_id, base in portico["bases"].items():
        x = base["x"]
        nodo = ss.find_node_id([x, 0])
        if base["tipo"] == "empotramiento":
            ss.add_support_fixed(node_id=nodo)
        elif base["tipo"] == "articulado":
            ss.add_support_hinged(node_id=nodo)

    # Cargas puntuales
    if "cargas_puntuales" in portico:
        for carga in portico["cargas_puntuales"]:
            x = carga["coordenadas"]["x"]
            y = carga["coordenadas"]["y"]
            nodo = ss.find_node_id([x, y])
            if nodo is not None:
                carga["valor_kN"] = float(carga["valor_kN"])
                ss.point_load(node_id=nodo, Fy=-carga["valor_kN"])

    # Viento
    for viga_id, viga in portico["vigas"].items():
        if "tramos" in viga:
            total_w = 0.0
            for tramo in viga["tramos"]:
                if "cargas" in tramo and "W_total" in tramo["cargas"]:
                    L_tramo = tramo["longitud_m"]
                    q_w = float(tramo["cargas"]["W_total"])
                    total_w += q_w * L_tramo
            columnas_piso = list(portico["columnas"].values())
            if columnas_piso and total_w > 0:
                carga_por_columna = total_w / len(columnas_piso)
                for col in columnas_piso:
                    nodo_sup = col["nodos"][1]
                    ss.point_load(node_id=nodo_sup, Fx=carga_por_columna)
                    col["carga_viento_kN"] = carga_por_columna

    return ss
# Mostrar lista de pórticos disponibles
print("Pórticos disponibles en estructura.json:")
for i, pid in enumerate(estructura.keys(), start=1):
    print(f"{i}. {pid}")

# Preguntar cuál procesar
opcion = int(input("¿Qué pórtico querés dimensionar primero? (número): "))
portico_id = list(estructura.keys())[opcion - 1]
portico = estructura[portico_id]

print(f"\n=== Dimensionando {portico_id} ===")

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
# 3- Resolver el sistema de forma sensata
# -------------------------------

# Elegir la combinacion a mostrar (puede ser "servicio" o cualquier mayorada)
combinacion_ver = "1.2D+1.6L"  # servicio, 1.2D+1.6L <- cambiar según lo que quieras inspeccionar

# Crear el sistema solo una vez
ss = SystemElements(mesh=50)
ss = crear_portico(portico, ss)

# Aplicar la combinacion elegida
aplicar_cargas(portico, ss, combinacion_ver)
ss.solve()

# Mostrar estructura y diagramas para esa combinacion
print(f"\n--- Mostrando combinacion: {combinacion_ver} ---")
ss.show_structure()
ss.show_bending_moment()
#ss.show_shear_force()
#ss.show_displacement()
ss.show_axial_force()

# -------------------------------
# Print general de resultados
# -------------------------------
print("\n--- DEBUG ELEMENT RESULTS ---")
for eid in ss.element_map:
    elem = ss.get_element_results(eid, verbose=True)

    print(f"\nElemento {eid}")
    print("keys:", elem.keys())

    if elem.get("M") is not None:
        # Mostrar todo el array
        print("M =", elem["M"])
        # Mostrar un valor intermedio (ejemplo el segundo punto)
        if len(elem["M"]) > 2:
            print("M[0] =", elem["M"][0])
        # Mostrar último valor
        print("M[-1] =", elem["M"][-1])
    else:
        print("Mmin =", elem.get("Mmin"))
        print("Mmax =", elem.get("Mmax"))

# -------------------------------
# Retroalimentar VIGAS (servicio y mayorados en un solo paso)
# -------------------------------
for viga_id, viga in portico["vigas"].items():
    nuevos_tramos = []
    for tramo in viga["tramos"]:
        eid = tramo.get("eid")

        # si es voladizo, salida simplificada
        if tramo.get("es_voladizo", False):
            elem = ss.get_element_results(element_id=eid, verbose=True)

            M = list(elem.get("M")) if elem.get("M") is not None else []
            V = list(elem.get("Q")) if elem.get("Q") is not None else []

            Mu_raw = round(M[0], 2) if abs(M[0]) > 1e-6 else round(M[-1], 2)
            Vu_raw = round(V[0], 2) if abs(V[0]) > 1e-6 else round(V[-1], 2)
            # aplicar convención de signo
            Mu_emp = -round(Mu_raw, 2)   # 👉 siempre negativo en empotramiento
            Vu_emp = round(Vu_raw, 2)    # 👉 cortante positivo

            nuevo_tramo = {
                "id": tramo["id"],
                "longitud_m": tramo["longitud_m"],
                "es_voladizo": True,
                "cargas": tramo.get("cargas", {}),
                "eid": tramo.get("eid"),
                "nodos": tramo.get("nodos", []),
                # mayorados
                "Mu_kNm": Mu_emp,
                "Vu_kN_emp": Vu_emp,
                # servicio
                "Mu_servicio": -round(tramo.get("Mu_servicio_der", 0), 2),
                "Vu_servicio_emp": round(tramo.get("Vu_servicio_der", 0), 2)
            }

        else:
            # inicializar acumuladores
            max_Mu_izq = max_Mu_der = max_Mu_campo = 0.0
            max_Vu_izq = max_Vu_der = 0.0

            # recorrer TODAS las combinaciones
            for comb in tramo.get("cargas", {}).get("combinaciones", {}).keys():
                elem = ss.get_element_results(element_id=eid, verbose=True)

                M = list(elem.get("M")) if elem.get("M") is not None else []
                V = list(elem.get("Q")) if elem.get("Q") is not None else []

                Mu_izq   = round(M[0], 2) if len(M) > 0 else 0.0
                Mu_der   = round(M[-1], 2) if len(M) > 0 else 0.0
                Mu_campo = round(min(M), 2) if len(M) > 0 else 0.0   # 👉 máximo negativo

                Vu_izq   = round(V[0], 2) if len(V) > 0 else 0.0
                Vu_der   = round(V[-1], 2) if len(V) > 0 else 0.0

                # acumular máximos absolutos
                max_Mu_izq = max(max_Mu_izq, abs(Mu_izq))
                max_Mu_der = max(max_Mu_der, abs(Mu_der))
                max_Mu_campo = max(max_Mu_campo, abs(Mu_campo))
                max_Vu_izq = max(max_Vu_izq, abs(Vu_izq))
                max_Vu_der = max(max_Vu_der, abs(Vu_der))

            # tramo final con mayorados y servicio
            nuevo_tramo = {
                "id": tramo["id"],
                "longitud_m": tramo["longitud_m"],
                "es_voladizo": False,
                "cargas": tramo.get("cargas", {}),
                "eid": tramo.get("eid"),
                "nodos": tramo.get("nodos", []),
                # mayorados
                "Mu_kNm_apoyo_izq": -max_Mu_izq,   # convención: apoyos negativos
                "Mu_kNm_apoyo_der": -max_Mu_der,
                "Mu_kNm_campo": abs(max_Mu_campo), # campo positivo
                "Vu_kN_izq": max_Vu_izq,
                "Vu_kN_der": max_Vu_der,
                # servicio
                "Mu_servicio_izq": -tramo.get("Mu_servicio_izq", 0),
                "Mu_servicio_der": -tramo.get("Mu_servicio_der", 0),
                "Mu_servicio_campo": abs(tramo.get("Mu_servicio_campo", 0)),
                "Vu_servicio_izq": tramo.get("Vu_servicio_izq", 0),
                "Vu_servicio_der": tramo.get("Vu_servicio_der", 0)
            }

        nuevos_tramos.append(nuevo_tramo)

    viga["tramos"] = nuevos_tramos

# -------------------------------
# Retroalimentar COLUMNAS (servicio y mayorados)
# -------------------------------
for col_id, col in portico["columnas"].items():
    eid = col.get("eid")
    elem = ss.get_element_results(element_id=eid, verbose=True) if eid else {}
    L = elem.get("length", col.get("altura_m", 3.0))

    Mu_inf = Mu_sup = 0.0
    distancia_campo_cero = None

    # momentos desde nodos
    if "M" in elem and elem["M"] is not None and len(elem["M"]) > 0:
        M = list(elem["M"])
        y_inf = col["Nodo_inferior"]["uy_mm"]
        y_sup = col["Nodo_superior"]["uy_mm"]

        # corregir orientación según desplazamientos
        if y_inf <= y_sup:
            Mu_inf, Mu_sup = round(M[0], 2), round(M[-1], 2)
        else:
            Mu_inf, Mu_sup = round(M[-1], 2), round(M[0], 2)

        # momento cero en campo (interpolación lineal)
        for i in range(len(M)-1):
            if M[i] * M[i+1] < 0:  # cambio de signo
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
    ux_sup = col["Nodo_superior"].get("ux_mm", 0.0)
    deriva_h = round(ux_sup / (altura * 1000), 6)

    # paquete nuevo para la columna
    nuevo_paquete = {
        "x": col.get("x", 0.0),
        "altura_m": col.get("altura_m", L),
        "eid": eid,
        "nodos": col.get("nodos", []),
        "Nodo_inferior": col.get("Nodo_inferior", {}),
        "Nodo_superior": col.get("Nodo_superior", {}),
        # combinación crítica (mayorado)
        "Mu_kNm_inf": Mu_inf,
        "Mu_kNm_sup": Mu_sup,
        "P_kN": N_val,
        # servicio
        "Mu_servicio_inf": Mu_inf if caso == "servicio" else col.get("Mu_servicio_inf"),
        "Mu_servicio_sup": Mu_sup if caso == "servicio" else col.get("Mu_servicio_sup"),
        "P_servicio": N_val if caso == "servicio" else col.get("N_servicio"),
        # extras
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
        }

        if caso == "servicio":
            nuevo_paquete.update({
                "Fx_servicio": Fx,
                "Fy_servicio": Fy,
                "Tz_servicio": Tz
            })
        else:  # combinación crítica (mayorado)
            nuevo_paquete.update({
                "Fx_kN": Fx,
                "Fy_kN": Fy,
                "Tz_kNm": Tz
            })

        portico["bases"][base_id].update(nuevo_paquete)

# -------------------------------
# Guardar retroalimentado
# -------------------------------
estructura[portico_id] = portico


# Guardar todo el JSON
with open(archivo, "w") as f:
    json.dump(estructura, f, indent=2)
print(f"\nDatos del pórtico '{portico_id}' actualizados en '{archivo}'")
import re

from calc import rutas

def leer_cargas_txt():
    archivos = rutas.listar(rutas.SAL_ANALISIS_CARGAS, "*.txt")
    print("\nArchivos de análisis disponibles:")
    for i, archivo in enumerate(archivos, 1):
        print(f"{i}. {archivo.name}")
    if not archivos:
        raise FileNotFoundError(
            f"No hay análisis de cargas en {rutas.SAL_ANALISIS_CARGAS}"
        )
    opcion = int(input("Seleccione archivo de cargas (número): "))
    if not 1 <= opcion <= len(archivos):
        raise ValueError(f"Opción fuera de rango: {opcion}")
    archivo_elegido = archivos[opcion - 1]

    with open(archivo_elegido, "r", encoding="utf-8") as f:
        contenido = f.read()

    def extraer_valor(patron):
        m = re.search(patron, contenido)
        return float(m.group(1)) if m else None

    return {
        "D_total": extraer_valor(r"D total\s*=\s*([\d.]+)"),
        "L_total": extraer_valor(r"L total\s*=\s*([\d.]+)"),
        "W_total": extraer_valor(r"W total\s*=\s*([\d.]+)"),
        "combinaciones": {
            "1.4D": extraer_valor(r"1\.4D:\s*([\d.]+)"),
            "1.2D+1.6L": extraer_valor(r"1\.2D\+1\.6L:\s*([\d.]+)"),
            "1.2D+0.5L+1.6W": extraer_valor(r"1\.2D\+0\.5L\+1\.6W:\s*([\d.]+)"),
            "0.9D+1.6W": extraer_valor(r"0\.9D\+1\.6W:\s*([\d.]+)")
        }
    }

def ingresar_datos_estructura():
    # 1. Cargar la obra y elegir un nombre libre sin pisar pórticos existentes.
    estructura = rutas.cargar_estructura()
    archivo = rutas.ESTRUCTURA
    portico_id = rutas.proximo_nombre_portico()
    numero_portico = rutas.numero_portico(portico_id, estructura)
    estructura[portico_id] = {"vigas": {}, "columnas": {}, "bases": {}, "cargas_puntuales": []}

    # 2. Preguntar cantidad de pisos
    n_pisos = int(input("Cantidad de pisos (0 = PB, 1 = primer piso, etc.): "))

    longitudes_pb = []   # se guarda solo una vez en PB
    nivel_actual = 0.0   # acumulador de nivel

    for piso in range(n_pisos+1):
        print(f"\nPiso {piso}:")
        if piso == 0:
            n_tramos = int(input("Cantidad de tramos: "))
        else:
            n_tramos = len(longitudes_pb)  # reutiliza los tramos de PB

        altura_col = float(input("Altura columnas [m]: "))

        # Generar columnas automáticas
        x_pos = 0.0
        for j in range(n_tramos+1):
            letra = chr(97+j)
            col_id = f"C{piso}-{letra}"

            estructura[portico_id]["columnas"][col_id] = {
                "x": x_pos,
                "altura_m": altura_col,
                "nivel": nivel_actual
            }

            # bases solo en planta baja
            if piso == 0:
                base_id = f"B{piso}-{letra}"
                estructura[portico_id]["bases"][base_id] = {
                    "x": x_pos,
                    "tipo": "empotramiento"
                }

            if j < n_tramos:
                if piso == 0:
                    L = float(input(f"Longitud tramo {j+1} [m]: "))
                    longitudes_pb.append(L)
                else:
                    L = longitudes_pb[j]  # reutiliza
                x_pos += L

        # Crear viga principal
        viga_id = f"V{piso}-{numero_portico}"
        estructura[portico_id]["vigas"][viga_id] = {"tramos": []}

        # Guardar tramos
        x_pos = 0.0
        for t, L in enumerate(longitudes_pb, 1):
            tramo_id = f"{viga_id} T{t}"
            tramo = {
                "id": tramo_id,
                "longitud_m": L,
                "es_voladizo": False,
                "x_inicio": x_pos,
                "x_fin": x_pos + L
            }
            asignar = input("¿Asignar cargas desde TXT? (s/n): ").lower() == "s"
            if asignar:
                tramo["cargas"] = leer_cargas_txt()
            estructura[portico_id]["vigas"][viga_id]["tramos"].append(tramo)
            x_pos += L

        # Voladizos
        vol_izq = input("¿Voladizo a la izquierda? (s/n): ").lower() == "s"
        if vol_izq:
            Lvi = float(input("Longitud voladizo izq [m]: "))
            tramo_vol = {
                "id": f"{viga_id} tv_izq",
                "longitud_m": Lvi,
                "es_voladizo": True,
                "x_inicio": 0.0 - Lvi,
                "x_fin": 0.0
            }
            asignar = input("¿Asignar cargas al voladizo desde TXT? (s/n): ").lower() == "s"
            if asignar:
                tramo_vol["cargas"] = leer_cargas_txt()
            estructura[portico_id]["vigas"][viga_id]["tramos"].append(tramo_vol)

        vol_der = input("¿Voladizo a la derecha? (s/n): ").lower() == "s"
        if vol_der:
            Lvd = float(input("Longitud voladizo der [m]: "))
            tramo_vol = {
                "id": f"{viga_id} tv_der",
                "longitud_m": Lvd,
                "es_voladizo": True,
                "x_inicio": x_pos,
                "x_fin": x_pos + Lvd
            }
            asignar = input("¿Asignar cargas al voladizo desde TXT? (s/n): ").lower() == "s"
            if asignar:
                tramo_vol["cargas"] = leer_cargas_txt()
            estructura[portico_id]["vigas"][viga_id]["tramos"].append(tramo_vol)

        # Cargas puntuales
        while True:
            agregar = input("¿Ingresar carga puntual en este piso? (s/n): ").lower() == "s"
            if not agregar:
                break

            # calcular rango total de la viga en este piso
            viga_id = f"V{piso}-{numero_portico}"
            x_min = min(tramo["x_inicio"] for tramo in estructura[portico_id]["vigas"][viga_id]["tramos"])
            x_max = max(tramo["x_fin"] for tramo in estructura[portico_id]["vigas"][viga_id]["tramos"])

            print(f"Rango válido de coordenadas X: {x_min} a {x_max} m")
            x = float(input("Coordenada X desde el origen [m]: "))
            if x < x_min or x > x_max:
                print(f"[WARN] La carga puntual en x={x} está fuera del rango ({x_min} a {x_max}).")
                print("Se asignará al tramo más cercano.")

            valor = float(input("Valor de la carga puntual [kN]: "))
            tipo = input("Tipo de carga (ej. maquinaria, tabique, sobrecarga): ")   
            columnas_piso = [c for cid, c in estructura[portico_id]["columnas"].items() if cid.startswith(f"C{piso}-")]
            if columnas_piso:
                y_viga = columnas_piso[0]["nivel"] + columnas_piso[0]["altura_m"]
            else:
                y_viga = nivel_actual

            carga_puntual = {
                "piso": piso,
                "coordenadas": {"x": x, "y": y_viga},
                "valor_kN": valor,
                "tipo": tipo
            }
            estructura[portico_id]["cargas_puntuales"].append(carga_puntual)

        nivel_actual += altura_col


    # 3. Distribuir cargas puntuales en cada tramo (siempre asignar)
    for viga_id, viga in estructura[portico_id]["vigas"].items():
        for tramo in viga["tramos"]:
            tramo.setdefault("cargas_puntuales", [])

        for carga in estructura[portico_id].get("cargas_puntuales", []):
            x_carga = carga["coordenadas"]["x"]
            asignado = False
            for tramo in viga["tramos"]:
                if tramo["x_inicio"] <= x_carga <= tramo["x_fin"]:
                    tramo["cargas_puntuales"].append(carga)
                    asignado = True
                    break
            if not asignado:
                # asignar al tramo más cercano
                if x_carga < viga["tramos"][0]["x_inicio"]:
                    viga["tramos"][0]["cargas_puntuales"].append(carga)
                else:
                    viga["tramos"][-1]["cargas_puntuales"].append(carga)

    estructura[portico_id]["cargas_puntuales"] = []

    # 4. Guardar JSON
    rutas.guardar_estructura(estructura)

    print(f"\nDatos guardados en {archivo}")

if __name__ == "__main__":
    ingresar_datos_estructura()

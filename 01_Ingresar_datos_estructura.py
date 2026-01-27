import json, os, re
from pathlib import Path

def leer_cargas_txt():
    carpeta = Path("salidas") / "analisis_cargas"
    archivos = list(carpeta.glob("*.txt"))
    print("\nArchivos de análisis disponibles:")
    for i, archivo in enumerate(archivos, 1):
        print(f"{i}. {archivo.name}")
    opcion = int(input("Seleccione archivo de cargas (número): "))
    archivo_elegido = archivos[opcion-1]

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
    # 1. Autonumerar Portico
    base_dir = os.path.dirname(os.path.abspath(__file__))
    carpeta = os.path.join(base_dir, "datos")
    os.makedirs(carpeta, exist_ok=True)
    archivo = os.path.join(carpeta, "estructura.json")

    if os.path.exists(archivo):
        with open(archivo, "r") as f:
            estructura = json.load(f)
    else:
        estructura = {}

    nro_portico = len(estructura) + 1
    portico_id = f"Portico {nro_portico}"
    estructura[portico_id] = {"vigas": {}, "columnas": {}, "bases": {}}

    # 2. Preguntar cantidad de pisos
    n_pisos = int(input("Cantidad de pisos (0 = PB, 1 = primer piso, etc.): "))

    for piso in range(n_pisos+1):
        print(f"\nPiso {piso}:")
        n_tramos = int(input("Cantidad de tramos: "))
        altura_col = float(input("Altura columnas [m]: "))

        # Generar columnas automáticas
        x_pos = 0.0
        for j in range(n_tramos+1):
            letra = chr(97+j)
            col_id = f"C{piso}-{letra}"
            estructura[portico_id]["columnas"][col_id] = {
                "x": x_pos,
                "altura_m": altura_col
            }
            base_id = f"B{piso}-{letra}"
            estructura[portico_id]["bases"][base_id] = {
                "x": x_pos,
                "tipo": "empotramiento"
            }
            if j < n_tramos:
                L = float(input(f"Longitud tramo {j+1} [m]: "))
                x_pos += L

        # Crear viga principal
        viga_id = f"V{piso}-1"
        estructura[portico_id]["vigas"][viga_id] = {"tramos": []}

        # Guardar tramos
        x_pos = 0.0
        for t in range(n_tramos):
            L = float(input(f"Longitud tramo {t+1} [m] (confirmar): "))
            tramo_id = f"{viga_id} T{t+1}"
            tramo = {
                "id": tramo_id,
                "longitud_m": L,
                "es_voladizo": False
            }
            asignar = input("¿Asignar cargas desde TXT? (s/n): ").lower() == "s"
            if asignar:
                tramo["cargas"] = leer_cargas_txt()
            estructura[portico_id]["vigas"][viga_id]["tramos"].append(tramo)
            x_pos += L

        # Voladizo izquierdo
        vol_izq = input("¿Voladizo a la izquierda? (s/n): ").lower() == "s"
        if vol_izq:
            Lvi = float(input("Longitud voladizo izq [m]: "))
            tramo_vol = {
                "id": f"{viga_id} tv_izq",
                "longitud_m": Lvi,
                "es_voladizo": True
            }
            asignar = input("¿Asignar cargas al voladizo desde TXT? (s/n): ").lower() == "s"
            if asignar:
                tramo_vol["cargas"] = leer_cargas_txt()
            estructura[portico_id]["vigas"][viga_id]["tramos"].append(tramo_vol)

        # Voladizo derecho
        vol_der = input("¿Voladizo a la derecha? (s/n): ").lower() == "s"
        if vol_der:
            Lvd = float(input("Longitud voladizo der [m]: "))
            tramo_vol = {
                "id": f"{viga_id} tv_der",
                "longitud_m": Lvd,
                "es_voladizo": True
            }
            asignar = input("¿Asignar cargas al voladizo desde TXT? (s/n): ").lower() == "s"
            if asignar:
                tramo_vol["cargas"] = leer_cargas_txt()
            estructura[portico_id]["vigas"][viga_id]["tramos"].append(tramo_vol)

        # Cargas puntuales
        estructura[portico_id].setdefault("cargas_puntuales", [])
        while True:
            agregar = input("¿Ingresar carga puntual en este piso? (s/n): ").lower() == "s"
            if not agregar:
                break

            x = float(input("Coordenada X desde el origen [m]: "))
            valor = float(input("Valor de la carga puntual [kN]: "))
            tipo = input("Tipo de carga (ej. maquinaria, tabique, sobrecarga): ")

            altura_y = sum(
                estructura[portico_id]["columnas"][col]["altura_m"]
                for col in estructura[portico_id]["columnas"]
                if col.startswith("C") and int(col.split("-")[0][1:]) < piso
            )

            carga_puntual = {
                "piso": piso,
                "coordenadas": {"x": x, "y": altura_y},
                "valor_kN": valor,
                "tipo": tipo
            }
            estructura[portico_id]["cargas_puntuales"].append(carga_puntual)

    # 3. Guardar JSON
    with open(archivo, "w") as f:
        json.dump(estructura, f, indent=2)

    print(f"\nDatos guardados en {archivo}")

if __name__ == "__main__":
    ingresar_datos_estructura()
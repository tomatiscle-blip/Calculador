import os
import json
import math
import csv

# ============================================================
# CONFIGURACIÓN GENERAL
# ============================================================

OUTPUT_DIR = "salidas/bases"
CSV_COLUMNAS = "salidas/columnas/planilla_columnas.csv"
JSON_ESTRUCTURA = "datos/estructura.json"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================
# ENTRADAS GENERALES
# ============================================================

fck = float(input("Ingrese fck [MPa] (ej: 20): "))
fy = 420  # acero ADN 420 MPa

prof = float(input("Ingrese profundidad fundación [m] (0.8 / 1.3): "))

if abs(prof - 0.8) < 0.05:
    q_adm_kgcm2 = 0.87
elif abs(prof - 1.3) < 0.05:
    q_adm_kgcm2 = 1.36
else:
    q_adm_kgcm2 = float(input("Ingrese q_adm [kg/cm²]: "))

q_adm_kPa = q_adm_kgcm2 * 98.1
print(f"\nq_adm adoptado = {q_adm_kPa:.2f} kPa")

# ============================================================
# LECTURA COLUMNAS
# ============================================================

def leer_columnas_csv(portico):
    columnas = {}
    with open(CSV_COLUMNAS, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter=';')
        for r in reader:
            if r["Pórtico"] == portico:
                columnas[r["Columna"]] = {
                    "tipo_seccion": r["Tipo_seccion"],
                    "dimensiones": r["Dimensiones"],
                    "recubrimiento_cm": float(r["Recubrimiento(cm)"]),
                    "n_barras": int(r["n_barras"]),
                    "diam_long_cm": float(r["diam_long(cm)"])
                }
    return columnas

# ============================================================
# GEOTECNIA ZAPATA
# ============================================================

def dimensionar_base_geotecnia(Fy, M, q_adm):
    N = abs(Fy)
    area = N / q_adm
    L = math.sqrt(area)
    e = M / N if N > 0 else 0.0

    q_med = N / area
    q_max = q_med * (1 + 6 * e / L)
    q_min = q_med * (1 - 6 * e / L)

    return {
        "lado_m": round(L, 2),
        "area_m2": round(area, 2),
        "excentricidad_m": round(e, 3),
        "q_max_kPa": round(q_max, 2),
        "q_min_kPa": round(q_min, 2),
        "requiere_viga_fundacion": q_min < 0
    }

# ============================================================
# ARMADURA MÍNIMA ZAPATA
# ============================================================

def armadura_minima_zapata(lado_m, h_m):
    """
    Calcula armadura mínima de zapata en ambos sentidos
    y metraje de acero.
    """

    # selección diámetro
    if h_m <= 0.40:
        diam_mm = 8
        area_barra_mm2 = 50.3
    else:
        diam_mm = 10
        area_barra_mm2 = 78.5

    paso_m = 0.20
    recubrimiento_m = 0.07

    # longitud útil
    L_util = lado_m - 2 * recubrimiento_m

    # cantidad de barras por sentido
    n_barras = int(L_util / paso_m) + 1

    # longitud total por sentido
    long_total_sentido = n_barras * L_util

    # total ambos sentidos
    long_total = 2 * long_total_sentido

    peso_kg_m = {
        8: 0.395,
        10: 0.617
    }[diam_mm]

    peso_total_kg = long_total * peso_kg_m

    return {
        "diametro": f"Ø{diam_mm}",
        "paso_cm": int(paso_m * 100),
        "barras_por_sentido": n_barras,
        "longitud_total_m": round(long_total, 2),
        "peso_total_kg": round(peso_total_kg, 1)
    }

# ============================================================
# TRONCO (PEDESTAL)
# ============================================================

def calcular_tronco(col, altura_fundacion, espesor_zapata):
    """
    Tronco siempre CUADRADO
    lado = dimensión máxima de columna + 3 cm por lado
    altura = altura_fundacion - espesor_zapata
    """

    # ----------------------------------
    # Dimensión de la columna
    # ----------------------------------
    if col["tipo_seccion"] == "III":  # circular
        diam_cm = float(col["dimensiones"].replace("Ø", ""))
        dim_col_m = diam_cm / 100
    else:  # rectangular
        partes = col["dimensiones"].replace("b=", "").replace("h=", "").split(",")
        b = float(partes[0]) / 100
        h = float(partes[1]) / 100
        dim_col_m = max(b, h)

    # ----------------------------------
    # Tronco: +3 cm por lado
    # ----------------------------------
    lado = dim_col_m + 0.06  # 3 cm por lado

    # Redondeo constructivo a 5 cm
    lado = math.ceil(lado / 0.05) * 0.05

    # ----------------------------------
    # Altura del tronco
    # ----------------------------------
    altura = altura_fundacion - espesor_zapata
    altura = round(altura, 2)

    volumen = lado * lado * altura

    return {
        "lado_m": round(lado, 2),
        "altura_m": altura,
        "volumen_m3": round(volumen, 3)
    }

# ============================================================
# PELOS DE COLUMNA
# ============================================================

def calcular_pelos(col, altura_tronco):
    """
    Calcula pelos de arranque de columna:
    - anclaje en zapata con gancho
    - subida por tronco
    - empalme con columna
    """

    diam_cm = col["diam_long_cm"] / 10   # mm → cm
    n = col["n_barras"]

    # ----------------------------------
    # Longitudes reglamentarias
    # ----------------------------------
    L_anclaje = 30 * diam_cm / 100     # gancho en zapata
    L_empalme = 50 * diam_cm / 100     # empalme con columna

    # ----------------------------------
    # Longitud total del pelo
    # ----------------------------------
    L_total = L_anclaje + altura_tronco + L_empalme

    return {
        "cantidad": n,
        "diametro_cm": diam_cm,
        "anclaje_m": round(L_anclaje, 2),
        "empalme_m": round(L_empalme, 2),
        "altura_tronco_m": round(altura_tronco, 2),
        "longitud_total_m": round(L_total, 2)
    }


# ============================================================
# VIGA DE FUNDACIÓN
# ============================================================

def calcular_viga_fundacion(base_i, base_j):
    L = abs(base_j["x"] - base_i["x"])

    M = max(abs(base_i["Tz_kNm"]), abs(base_j["Tz_kNm"]))
    T_kN = M / L

    fy = 420
    phi = 0.9

    As_req = (T_kN * 1e3) / (phi * fy)

    # adopción lógica
    if As_req <= 100:
        armadura = "2Ø12"
        As_adop = 226
    else:
        armadura = "2Ø16"
        As_adop = 402

    # estribos
    paso = 0.20
    n_estribos = int(L / paso) + 1
    long_estribo = 2 * (0.20 + 0.40)  # perímetro aproximado
    peso_estribos = n_estribos * long_estribo * 0.395

    return {
        "longitud_m": round(L, 2),
        "seccion_cm": "20x40",
        "armadura_longitudinal": armadura,
        "As_requerida_mm2": round(As_req, 1),
        "estribos": {
            "diam": "Ø8",
            "paso_cm": 20,
            "cantidad": n_estribos,
            "peso_kg": round(peso_estribos, 1)
        }
    }

# ============================================================
# LECTURA ESTRUCTURA
# ============================================================

with open(JSON_ESTRUCTURA, "r", encoding="utf-8") as f:
    estructura = json.load(f)

porticos = list(estructura.keys())
for i, p in enumerate(porticos, 1):
    print(f"{i}. {p}")

portico = porticos[int(input("Seleccione pórtico: ")) - 1]

columnas = leer_columnas_csv(portico)
bases = estructura[portico]["bases"]

# ============================================================
# CÁLCULO GENERAL
# ============================================================

resultados = {}

for nb, b in bases.items():
    geo = dimensionar_base_geotecnia(b["Fy_kN"], b["Tz_kNm"], q_adm_kPa)
    h = max(0.35, geo["lado_m"] / 10)

    col = columnas.get(nb.replace("B", "C"))
    altura_fundacion = prof

    resultados[nb] = {
        "geotecnia": geo,
        "espesor_m": round(h, 2),
        "armadura_zapata": armadura_minima_zapata(
            geo["lado_m"], h
        ),
        "tronco": calcular_tronco(col, altura_fundacion, h) if col else None,
        "pelos_columna": calcular_pelos(col, altura_fundacion - h) if col else None
    }

# ============================================================
# VIGAS DE FUNDACIÓN
# ============================================================

vigas = {}
nombres = list(bases.keys())

for i in range(len(nombres) - 1):
    bi = nombres[i]
    bj = nombres[i + 1]

    # condición geotécnica (resultado)
    necesita_viga = (
        resultados[bi]["geotecnia"]["requiere_viga_fundacion"] or
        resultados[bj]["geotecnia"]["requiere_viga_fundacion"]
    )

    if necesita_viga:
        vigas[f"{bi}-{bj}"] = calcular_viga_fundacion(
            bases[bi],   # tiene x y Tz_kNm
            bases[bj]
        )

# guardar en resultados finales
resultados["vigas_fundacion"] = vigas

# ============================================================
# GUARDAR
# ============================================================

out = os.path.join(OUTPUT_DIR, f"bases_{portico}.json")
with open(out, "w", encoding="utf-8") as f:
    json.dump(resultados, f, indent=4, ensure_ascii=False)

print(f"\n✔ Resultados guardados en {out}")

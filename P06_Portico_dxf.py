import ezdxf
import json
from pathlib import Path
import csv
import re
import math
import datetime
# ======================
# Configuración
# ======================

PORTICO = "Portico 1"

ESTRUCTURA_PATH = Path("datos/estructura.json")
PLANILLAS_VIGAS_DIR = Path("salidas/vigas")

BASES_PATH    = Path(f"salidas/bases/bases_{PORTICO}.json")
COLUMNAS_PATH = Path("salidas/columnas/planilla_columnas.csv")
OUT_PATH      = Path(f"salidas/dxf/{PORTICO}_lateral.dxf")

CAPAS = {
    "hormigon": {"color": 9},
    "armadura": {"color": 1},
    "texto":    {"color": 2},
    "ejes":     {"color": 8, "linetype": "MY_CENTER"},
    "suelo":    {"color": 7, "linetype": "CONTINUOUS"},
    "puntos_inflexion": {"color": 4, "linetype": "MY_DOT"},
    "recubrimiento": {"color": 6, "linetype": "MY_DASHED"},
    "estribos": {"color": 3, "linetype": "CONTINUOUS"},
    "sunchos": {"color": 11, "linetype": "CONTINUOUS"},  # agregado
    "cotas": {"color": 5, "linetype": "CONTINUOUS"} 
}
LINETYPES = {
    "MY_DOT":    [0.05, 0.0, -0.05],                     # punto + espacio
    "MY_CENTER": [0.18, 0.12, -0.02, 0.02, -0.02],       # centro
    "MY_DASHED": [0.07, 0.04, -0.05],                    # trazos
}
TEXTO_ALTURA_TITULO  = 0.07  # para títulos (ej. nombre de viga/tramo)
TEXTO_ALTURA = 0.06
TEXTO_ALTURA_DETALLE = 0.05

CONFIG_BARRAS = {
    # --- Superiores ---
    "sup. voladizo": {
        "base": "sup_vol",      # fila voladizo
        "offset_y": 0.0,
        "capa": "arm_sup"
    },
    "sup. apoyo izq": {
        "base": "sup_ap",       # fila apoyos
        "offset_y": 0.0,
        "capa": "arm_sup"
    },
    "sup. apoyo der": {
        "base": "sup_ap",       # fila apoyos
        "offset_y": 0.0,
        "capa": "arm_sup"
    },
    "auxiliar superior": {
        "base": "sup_aux",      # debajo de apoyos
        "offset_y":0.0,
        "capa": "arm_sup"
    },
    "arm. comprimida": {
        "base": "sup_aux",       # misma fila que auxiliares/compresión
        "offset_y": 0.0,       # un poco más arriba que auxiliar
        "capa": "arm_sup_compresion"
    },

    # --- Inferiores ---
    "inf. tramo": {
        "base": "inf",
        "offset_y": 0.0,
        "capa": "arm_inf"
    },
    "inf. tramo completo": {
        "base": "inf",
        "offset_y": -0.20,
        "capa": "arm_inf"
    },
    "inf. mínima": {
        "base": "inf",
        "offset_y": -0.40,
        "capa": "arm_inf"
    },

    # --- Estribos ---
    "estribo apoyo_izq": {
        "base": "estribo",
        "offset_y": 0.0,
        "capa": "estribos"
    },
    "estribo apoyo_der": {
        "base": "estribo",
        "offset_y": 0.0,
        "capa": "estribos"
    },
    "estribo central": {
        "base": "estribo",
        "offset_y": 0.0,
        "capa": "estribos"
    },
    "estribo empotramiento": {
        "base": "estribo",
        "offset_y": 0.0,
        "capa": "estribos"
    },
    "estribo extremo_libre": {
        "base": "estribo",
        "offset_y": 0.0,
        "capa": "estribos"
    }
}


# ======================
# Utilidades
# ======================

def leer_estructura(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)

def crear_doc():
    doc = ezdxf.new(setup=True)
    msp = doc.modelspace()

    # Crear capas
    for c, p in CAPAS.items():
        if c not in doc.layers:
            doc.layers.add(c, **p)

    # Crear linetypes
    for name, pattern in LINETYPES.items():
        if name not in doc.linetypes:
            doc.linetypes.add(name, pattern)

    doc.set_modelspace_vport(10, center=(0, 0))
    return doc, msp

def leer_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)

def leer_csv(path, portico):
    with open(path, newline='', encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter=";")
        filas = [fila for fila in reader if fila["Pórtico"] == portico]
    return filas

def leer_planilla_tramo(path_txt):
    barras = []
    estribos = []
    pi_izq = None
    pi_der = None
    b = h = d = rec = None
    zonas = {}

    with open(path_txt, "r", encoding="utf-8") as f:
        lines = f.readlines()

    for line in lines:
        # --- Tabla de barras y estribos ---
        if re.match(r"^\s*\d+\s*\|", line):
            parts = [p.strip() for p in line.split("|")]
            if len(parts) >= 6:
                item = {
                    "pos": int(parts[0]),
                    "cant": int(parts[1]),
                    "diam": int(parts[2]),
                    "detalle": parts[3],
                    "longitud_m": float(parts[4]),
                    "peso": float(parts[5])
                }
                if "Estribo" in item["detalle"]:
                    # capturar también spacing si está al final
                    match = re.search(r"c/(\d+)cm", line)
                    if match:
                        item["spacing_m"] = float(match.group(1)) / 100.0
                    estribos.append(item)
                else:
                    barras.append(item)

        # --- Puntos de inflexión (solo interiores) ---
        if "PI izquierdo" in line:
            pi_izq = float(re.findall(r"[\d\.]+", line)[0])
        if "PI derecho" in line:
            pi_der = float(re.findall(r"[\d\.]+", line)[0])

        # --- Caso voladizo: línea compacta con b/h/rec ---
        if "b =" in line and "|" in line and "recubrimiento" in line.lower():
            nums = re.findall(r"[\d\.]+", line)
            if len(nums) >= 3:
                b   = float(nums[0]) / 100.0
                h   = float(nums[1]) / 100.0
                rec = float(nums[2]) / 100.0

        # --- Caso interior: líneas separadas ---
        elif line.strip().startswith("b ="):
            nums = re.findall(r"[\d\.]+", line)
            if nums:
                b = float(nums[0]) / 100.0

        elif line.strip().startswith("h ="):
            nums = re.findall(r"[\d\.]+", line)
            if nums:
                h = float(nums[0]) / 100.0

        elif line.strip().startswith("d ="):
            nums = re.findall(r"[\d\.]+", line)
            if nums:
                d = float(nums[0]) / 100.0

        elif "recubrimiento" in line.lower():
            nums = re.findall(r"[\d\.]+", line)
            if nums:
                rec = float(nums[0]) / 100.0

        # --- Zonificación ---
        if "Zonificación" in line:
            partes = line.split("|")
            for p in partes:
                m = re.search(r"(\w+)=([\d\.]+)-([\d\.]+)", p)
                if m:
                    nombre = m.group(1).lower()
                    ini = float(m.group(2))
                    fin = float(m.group(3))
                    zonas[nombre] = (ini, fin)

    return {
        "barras": barras,
        "estribos": estribos,
        "zonas": zonas,
        "pi_izq": pi_izq,
        "pi_der": pi_der,
        "b": b,
        "h": h,
        "d": d,
        "rec": rec
    }

# ======================
# Auxiliares
# ======================

def adaptar_columna(fila_csv):
    col = {
        "b": 0.20,
        "h": 0.20,
        "forma": "rectangular",
        "recubrimiento_m": 0.0,
        "tipo_armado": "",
        "n_barras": 0,
        "diam_long": 0.0,
        "diam_estribo": 0.0,
        "paso_estribo": 0.0,
    }

    if fila_csv:
        dim = fila_csv.get("Dimensiones", "").strip()
        if "b=" in dim and "h=" in dim:
            partes = [p.strip() for p in dim.split(",")]
            col["b"] = float(partes[0].split("=")[1]) / 100.0
            col["h"] = float(partes[1].split("=")[1]) / 100.0
            col["forma"] = "rectangular"
        elif "Ø" in dim:
            diam = float(dim.replace("Ø", "").strip()) / 100.0
            col["b"] = diam
            col["h"] = diam
            col["forma"] = "circular"

        col["recubrimiento_m"] = float(fila_csv.get("Recubrimiento(cm)", 0) or 0) / 100.0
        col["tipo_armado"] = fila_csv.get("Tipo", "").strip()
        col["n_barras"] = int(fila_csv.get("n_barras", 0) or 0)
        col["diam_long"] = float(fila_csv.get("diam_long(cm)", 0) or 0) / 100.0
        col["diam_estribo"] = float(fila_csv.get("diam_estribo(cm)", 0) or 0) / 100.0
        col["paso_estribo"] = float(fila_csv.get("paso_estribo(cm)", 0) or 0) / 100.0

    return col

def calcular_cotas(portico_data):
    cotas = {}
    y_actual = 0.0
    niveles = sorted(set(n[:2] for n in portico_data["columnas"].keys()))
    for nivel in niveles:
        col = next(c for n,c in portico_data["columnas"].items() if n.startswith(nivel))
        altura = col["altura_m"]
        cotas[nivel] = (y_actual, y_actual + altura)
        y_actual += altura
    return cotas

def nivel_viga_a_columna(nivel_viga: str) -> str:
    # "V0" -> "C0", "V1" -> "C1", etc.
    return "C" + nivel_viga[1:]

def normalizar_id(tramo_id: str) -> str:
    """
    Normaliza los IDs de tramo para que coincidan entre estructura, materiales, flexión, armaduras y estribos.
    Ejemplo: 'V0-1.V0-1 T1' -> 'V0-1 T1'
    """
    if "." in tramo_id:
        return tramo_id.split(".")[-1].strip()
    return tramo_id.strip()

def dibujar_contorno_viga(msp, x1, x2, y_centro, h):
    msp.add_lwpolyline([
        (x1, y_centro - h/2),
        (x2, y_centro - h/2),
        (x2, y_centro + h/2),
        (x1, y_centro + h/2),
        (x1, y_centro - h/2),
    ], dxfattribs={"layer": "hormigon"})

def dibujar_recubrimientos(msp, x1, x2, y_sup, y_inf, rec):
    # línea superior
    msp.add_line(
        (x1, y_sup - rec),
        (x2, y_sup - rec),
        dxfattribs={"layer": "recubrimiento", "linetype": "MY_DASHED"}
    )
    # línea inferior
    msp.add_line(
        (x1, y_inf + rec),
        (x2, y_inf + rec),
        dxfattribs={"layer": "recubrimiento", "linetype": "MY_DASHED"}
    ) 

def posicionar_barra(barra, x_ini, x_fin, y_sup, y_inf, rec, pi_izq=None, pi_der=None):
    detalle = barra["detalle"].lower().strip()
    L = barra["longitud_m"]

    if "sup" in detalle:
        if "izq" in detalle:
            x1, x2 = x_ini, x_ini + L
        elif "der" in detalle:
            x1, x2 = x_fin - L, x_fin
        else:
            x1, x2 = x_ini, x_ini + L

    elif "inf" in detalle:
        if "tramo completo" in detalle:
            x1, x2 = x_ini, x_fin
        elif "tramo" in detalle and pi_izq and pi_der:
            xc = (pi_izq + pi_der) / 2
            x1, x2 = xc - L/2, xc + L/2
        else:
            x1, x2 = x_ini, x_ini + L

    else:
        x1, x2 = x_ini, x_ini + L

    return float(x1), float(x2)

def dibujar_estribo_voladizo(msp, estribo, x1, x2, y_sup, y_inf, rec, lado_voladizo="der", longitud_tramo=None):
    detalle = estribo["detalle"].lower()
    match = re.search(r"c/(\d+)cm", estribo["detalle"])
    spacing = float(match.group(1)) / 100.0 if match else 0.16

    L = longitud_tramo if longitud_tramo else (x2 - x1)

    # definir zona según detalle
    if "empotramiento" in detalle:
        zona_ini, zona_fin = 0.00, 0.40
    elif "extremo_libre" in detalle:
        zona_ini, zona_fin = L-0.40, L
    else:
        zona_ini, zona_fin = 0.00, L

    # dibujar estribos
    x = x1 + zona_ini
    while x <= x1 + zona_fin:
        msp.add_line((x, y_inf + rec), (x, y_sup - rec), dxfattribs={"layer": "estribos"})
        x += spacing

    # texto sintético
    diam = estribo.get("diam", "?")
    etiqueta = f"{estribo['detalle']} | Ø{diam} c/{spacing*100:.0f}cm"

    if "empotramiento" in detalle:
        if lado_voladizo == "izq":
            pos_x = x1
            attach = 1  # Middle Left
        else:
            pos_x = x2-0.15
            attach = 3  # Middle Right
        pos_y = y_sup + 0.25

    elif "extremo_libre" in detalle:
        if lado_voladizo == "izq":
            pos_x = x1
            attach = 1
        else:
            pos_x = x2-0.15
            attach = 3
        pos_y = y_sup + 0.10

    else:
        pos_x = (x1 + x2) / 2
        pos_y = y_sup + 0.25
        attach = 5

    mtext = msp.add_mtext(
        etiqueta,
        dxfattribs={"layer": "texto", "char_height": 0.05, "attachment_point": attach}
    )
    mtext.set_location((pos_x, pos_y))

def dibujar_estribo(msp, estribo, x1, x2, y_sup, y_inf, rec, longitud_tramo=None, zonas=None):
    detalle = estribo["detalle"].lower()
    spacing = estribo.get("spacing_m", 0.16)
    L = longitud_tramo if longitud_tramo else (x2 - x1)

    L = longitud_tramo if longitud_tramo else (x2 - x1)

    # usar zonas reales si están disponibles
    if zonas:
        if "apoyo_izq" in detalle:
            zona_ini, zona_fin = zonas["apoyo_izq"]
        elif "central" in detalle:
            zona_ini, zona_fin = zonas["central"]
        elif "apoyo_der" in detalle:
            zona_ini, zona_fin = zonas["apoyo_der"]
        else:
            zona_ini, zona_fin = 0.00, L
    else:
        # fallback proporcional
        if "apoyo_izq" in detalle:
            zona_ini, zona_fin = 0.00, L*0.25
        elif "central" in detalle:
            zona_ini, zona_fin = L*0.25, L*0.75
        elif "apoyo_der" in detalle:
            zona_ini, zona_fin = L*0.75, L
        else:
            zona_ini, zona_fin = 0.00, L


    # dibujar estribos como líneas verticales
    x = x1 + zona_ini
    while x <= x1 + zona_fin:
        msp.add_line(
            (x, y_inf + rec),
            (x, y_sup - rec),
            dxfattribs={"layer": "estribos"}
        )
        x += spacing

    # --- texto de cota / sonificación ---
    xc = x1 + (zona_ini + zona_fin) / 2
    if "empotramiento" in detalle or "extremo_libre" in detalle:
        # llamar a la auxiliar para voladizo
        dibujar_estribo_voladizo(msp, estribo, x1, x2, y_sup, y_inf, rec, longitud_tramo)
    else:
        # texto arriba para tramos (con Ø incluido)
        yc = y_sup + 0.25
        diam = estribo.get("diam", "?")
        etiqueta = f"{estribo['detalle']} | Ø{diam} c/{spacing*100:.0f}cm"
        mtext = msp.add_mtext(
            etiqueta,
            dxfattribs={"layer": "texto", "char_height": 0.05, "attachment_point": 5}
        )
        mtext.set_location((xc, yc))

        # línea horizontal de cota
        msp.add_line(
            (x1 + zona_ini, y_sup + 0.20),
            (x1 + zona_fin, y_sup + 0.20),
            dxfattribs={"layer": "cotas"}
        )

        # remates verticales en los extremos (simulan flechas)
        msp.add_line(
            (x1 + zona_ini, y_sup + 0.20),
            (x1 + zona_ini, y_sup + 0.15),
            dxfattribs={"layer": "cotas"}
        )
        msp.add_line(
            (x1 + zona_fin, y_sup + 0.20),
            (x1 + zona_fin, y_sup + 0.15),
            dxfattribs={"layer": "cotas"}
        )

        # texto de longitud de zona (ej. "1.10 m") abajo de la cota
        longitud_zona = zona_fin - zona_ini
        etiqueta_cota = f"{longitud_zona:.2f} m"
        yc_cota = y_sup + 0.05
        mtext_cota = msp.add_mtext(
            etiqueta_cota,
            dxfattribs={"layer": "cotas", "char_height": 0.05, "attachment_point": 5}
        )
        mtext_cota.set_location((xc, yc_cota))

def dibujar_puntos_inflexion(msp, x1, y_sup, y_inf, pi_izq=None, pi_der=None):
    EXT_UP = 0.40   # extensión hacia arriba
    EXT_DOWN = 1.60 # extensión hacia abajo

    def dibujar_pi(x_pi):
        # línea de puntos de inflexión
        msp.add_line(
            (x_pi, y_inf - EXT_DOWN),
            (x_pi, y_sup + EXT_UP),
            dxfattribs={"layer": "puntos_inflexion", "linetype": "MY_DOT"}
        )
        # texto "P.I." centrado arriba
        mtext = msp.add_mtext(
            "P.I.",
            dxfattribs={
                "layer": "texto",
                "char_height": TEXTO_ALTURA,
                "attachment_point": 5  # Middle Center
            }
        )
        mtext.set_location((x_pi, y_sup + EXT_UP + 0.10))

    if pi_izq is not None:
        x_pi_izq = x1 + pi_izq
        dibujar_pi(x_pi_izq)

    if pi_der is not None:
        x_pi_der = x1 + pi_der
        dibujar_pi(x_pi_der)

def dibujar_despiece_barra(msp, barra, x1, x2, y_pos):
#- dibujar_despiece_barra(...) → se encarga de una sola barra, aplicando las reglas de ganchos y poniendo la etiqueta.

    diam_mm = barra["diam"]
    gancho = 0.10  # fijo, como lo tenés
    detalle = barra["detalle"].lower().strip()

    # polilínea base de la barra
    puntos = [(x1, y_pos), (x2, y_pos)]

    # reglas de ganchos
    if "sup" in detalle:
        if "voladizo" in detalle:
            # sup voladizo: sin gancho
            pass
        elif "izq" in detalle:
            # sup apoyo izq: gancho en el lado derecho
            puntos.append((x2, y_pos - gancho))
        elif "der" in detalle:
            # sup apoyo der: gancho en el lado izquierdo
            puntos.insert(0, (x1, y_pos - gancho))
        elif "aux" in detalle:
            # sup auxiliar: gancho hacia abajo en ambos extremos
            puntos.insert(0, (x1, y_pos - gancho))
            puntos.append((x2, y_pos - gancho))

    elif "inf" in detalle:
        if "tramo" in detalle:
            # inf tramo: ganchos hacia arriba en diagonal
            puntos.insert(0, (x1 - gancho/1.414, y_pos + gancho/1.414))
            puntos.append((x2 + gancho/1.414, y_pos + gancho/1.414))
        else:
            # inf completas: sin gancho
            pass

    # dibujar la barra con sus ganchos
    msp.add_lwpolyline(
        puntos,
        dxfattribs={"layer": "armadura", "lineweight": diam_mm}
    )

    # etiqueta centrada arriba
    etiqueta = f"{barra['cant']}Ø{diam_mm} | {barra['detalle']} | L={barra['longitud_m']:.2f}m"
    xc = (x1 + x2) / 2
    yc = y_pos + 0.05  # un poco arriba de la barra

    mtext = msp.add_mtext(
        etiqueta,
        dxfattribs={
            "layer": "texto",
            "char_height": 0.05,
            "attachment_point": 5   # Middle Center
        }
    )
    mtext.set_location((xc, yc))

def dibujar_barras_despiece(msp, barras, x1, x2, y_inf, y_sup, rec=None, pi_izq=None, pi_der=None):
    h = y_sup - y_inf

    # bases superiores
    y_base_sup_vol = y_sup - 0.20 - h        # fila voladizo
    y_base_sup_ap  = y_base_sup_vol - 0.20   # fila apoyos
    y_base_sup_aux = y_base_sup_ap - 0.20    # fila auxiliares/comprimida

    # base inferior: debajo de las superiores
    y_base_inf     = y_base_sup_aux - 0.40   # fila inferior

    inf_count = 0

    for barra in barras:
        detalle = barra["detalle"].lower().strip()

        # buscar configuración exacta
        config = None
        for key, val in CONFIG_BARRAS.items():
            if key in detalle:
                config = val
                break

        if config:
            base = config["base"]
            if base == "sup_vol":
                y_pos = y_base_sup_vol + config["offset_y"]
            elif base == "sup_ap":
                y_pos = y_base_sup_ap + config["offset_y"]
            elif base == "sup_aux":
                y_pos = y_base_sup_ap - 0.20 + config["offset_y"]
            elif base == "inf":
                y_pos = y_base_inf + config["offset_y"] - inf_count * 0.20
                inf_count += 1
            elif base == "estribo":
                # acá podés dibujar estribos con otra función
                continue
            else:
                # fallback
                y_pos = y_base_inf - (inf_count + 1) * 0.20
                inf_count += 1

            capa = config["capa"]

        else:
            # sin config → fallback inferior
            y_pos = y_base_inf - (inf_count + 1) * 0.20
            capa = "armadura"
            inf_count += 1

        # ahora posicionar horizontalmente
        x_ini, x_fin = posicionar_barra(barra, x1, x2, y_sup, y_inf, rec, pi_izq, pi_der)
        dibujar_despiece_barra(msp, barra, x_ini, x_fin, y_pos)

def dibujar_barra_longitudinal_L(msp, base, x_pos, ht, esp, lado="izq"):
    rec_zapata = base.get("recubrimiento_m", 0.05)
    pelos = base["pelos_columna"]
    empalme = pelos["empalme_m"]

    # nivel de arranque en la parrilla inferior
    y_arm = -(ht + esp) + rec_zapata + 0.05
    desarrollo = 0.20  # 20 cm de tramo horizontal

    if lado == "izq":
        x_start = x_pos - desarrollo
        x_end   = x_pos
    else:  # lado derecho
        x_start = x_pos + desarrollo
        x_end   = x_pos

    # tramo horizontal
    msp.add_line(
        (x_start, y_arm),
        (x_end, y_arm),
        dxfattribs={"layer": "armadura"}
    )

    # tramo vertical: arranca donde termina el horizontal
    y_top_tronco = 0.0
    y_end = y_top_tronco + empalme
    msp.add_line(
        (x_end, y_arm),
        (x_end, y_end),
        dxfattribs={"layer": "armadura"}
    )

def texto_base(nombre, base):
    arm = base["armadura_zapata"]
    tronco = base["tronco"]
    pelos = base["pelos_columna"]
    geo = base["geotecnia"]

    # convertir diámetro de pelos a mm
    diam_mm = round(pelos["diametro_cm"] * 10)

    detalle = (
        f"Armadura zapata: {arm['diametro']} @{arm['paso_cm']} cm, "
        f"{arm['barras_por_sentido']} barras/sentido\n"
        f"Zapata: lado {geo['lado_m']:.2f} m, espesor {base['espesor_m']:.2f} m\n"
        f"Tronco: {tronco['lado_m']:.2f} x {tronco['altura_m']:.2f} m\n"
        f"Pelos: {pelos['cantidad']} Ø{diam_mm} mm, "
        f"anclaje {pelos['anclaje_m']:.2f} m, empalme {pelos['empalme_m']:.2f} m"
    )
    return detalle

# ======================
# Dibujo lateral
# ======================

def dibujar_base_lateral(msp, base, nombre, x0):
    lado = base["geotecnia"]["lado_m"]
    esp  = base["espesor_m"]

    # tronco (desde terreno hacia abajo)
    t = base["tronco"]
    lt = t["lado_m"]
    ht = t["altura_m"]

    msp.add_lwpolyline([
        (x0 - lt/2, 0),
        (x0 + lt/2, 0),
        (x0 + lt/2, -ht),
        (x0 - lt/2, -ht),
        (x0 - lt/2, 0),
    ], dxfattribs={"layer": "hormigon"})

    # zapata (debajo del tronco)
    msp.add_lwpolyline([
        (x0 - lado/2, -ht),
        (x0 + lado/2, -ht),
        (x0 + lado/2, -(ht+esp)),
        (x0 - lado/2, -(ht+esp)),
        (x0 - lado/2, -ht),
    ], dxfattribs={"layer": "hormigon"})

    # recubrimiento solo abajo
    rec = base.get("recubrimiento_m", 0.05)
    msp.add_line(
        (x0 - lado/2 + rec, -(ht+esp) + rec),
        (x0 + lado/2 - rec, -(ht+esp) + rec),
        dxfattribs={"layer": "recubrimiento", "linetype": "MY_DASHED"}
    )

    # armadura inferior (una sola línea con ganchos)
    y_arm = -(ht+esp) + rec + 0.01  # un poquito arriba del recubrimiento

    # línea horizontal
    msp.add_line(
        (x0 - lado/2 + rec, y_arm),
        (x0 + lado/2 - rec, y_arm),
        dxfattribs={"layer": "armadura"}
    )

    # ganchos levantados de 10 cm
    msp.add_line(
        (x0 - lado/2 + rec, y_arm),
        (x0 - lado/2 + rec, y_arm + 0.10),
        dxfattribs={"layer": "armadura"}
    )
    msp.add_line(
        (x0 + lado/2 - rec, y_arm),
        (x0 + lado/2 - rec, y_arm + 0.10),
        dxfattribs={"layer": "armadura"}
    )
    # barras longitudinales en los dos extremos del tronco
    x_left  = x0 - lt/2 + 0.08
    x_right = x0 + lt/2 - 0.08

    dibujar_barra_longitudinal_L(msp, base, x_left, ht, esp, lado="izq")
    dibujar_barra_longitudinal_L(msp, base, x_right, ht, esp, lado="der")
    # --- Texto grande: nombre de la base ---
    mtext_nombre = msp.add_mtext(
        nombre,
        dxfattribs={"layer": "texto", "char_height": TEXTO_ALTURA_TITULO, "attachment_point": 1}
    )
    mtext_nombre.set_location((x0, -(ht+esp) - 0.20))

    detalle_txt = texto_base(nombre, base)
    mtext_detalle = msp.add_mtext(
        detalle_txt,
        dxfattribs={"layer": "texto", "char_height": TEXTO_ALTURA, "attachment_point": 1}
    )
    mtext_detalle.set_location((x0, -(ht+esp) - 0.40))

    return 0.0

def dibujar_columna_lateral(msp, columna, x0, y_base, y_top):
    b   = columna["b"]   # lado menor (m)
    h   = columna["h"]   # lado mayor (m)
    rec = columna["recubrimiento_m"]

    # --- print de control ---
    #print("Columna:", columna["nombre"],"forma:", columna["forma"],"b=", b, "h=", h,"rec=", rec,"tipo_armado=", columna.get("tipo_armado"),"paso=", columna.get("paso_estribo"),"n_barras=", columna.get("n_barras"),"diam_long=", columna.get("diam_long"))

    if columna["forma"] == "circular":
        diam = b  # aquí b=h=Ø
        # contorno
        msp.add_lwpolyline([
            (x0 - diam/2, y_base),
            (x0 + diam/2, y_base),
            (x0 + diam/2, y_top),
            (x0 - diam/2, y_top),
            (x0 - diam/2, y_base),
        ], dxfattribs={"layer": "hormigon"})
        # recubrimiento
        msp.add_lwpolyline([
            (x0 - diam/2 + rec, y_base + rec),
            (x0 + diam/2 - rec, y_base + rec),
            (x0 + diam/2 - rec, y_top - rec),
            (x0 - diam/2 + rec, y_top - rec),
            (x0 - diam/2 + rec, y_base + rec),
        ], dxfattribs={"layer": "recubrimiento", "linetype": "MY_DASHED"})
        texto = f"{columna['nombre']} (Ø{int(diam*100)})"

    else:  # rectangular
        # contorno
        msp.add_lwpolyline([
            (x0 - h/2, y_base),
            (x0 + h/2, y_base),
            (x0 + h/2, y_top),
            (x0 - h/2, y_top),
            (x0 - h/2, y_base),
        ], dxfattribs={"layer": "hormigon"})
        # recubrimiento
        msp.add_lwpolyline([
            (x0 - h/2 + rec, y_base + rec),
            (x0 + h/2 - rec, y_base + rec),
            (x0 + h/2 - rec, y_top - rec),
            (x0 - h/2 + rec, y_top - rec),
            (x0 - h/2 + rec, y_base + rec),
        ], dxfattribs={"layer": "recubrimiento", "linetype": "MY_DASHED"})
        texto = f"{columna['nombre']} ({int(b*100)}x{int(h*100)})"

    # texto identificador
    mtext = msp.add_mtext(
        texto,
        dxfattribs={"layer": "texto", "char_height": TEXTO_ALTURA_TITULO, "attachment_point": 3}
    )
    desplazamiento_x = -(h/2 + 0.10)
    desplazamiento_y = 1.05
    mtext.set_location((x0 + desplazamiento_x, y_base + desplazamiento_y))

    # estribos / sunchos
    paso = columna.get("paso_estribo", 0)
    if columna.get("tipo_armado") == "estribos" and paso > 0:
        y = y_base + rec
        while y < y_top - rec:
            msp.add_line(
                (x0 - h/2 + rec, y),
                (x0 + h/2 - rec, y),
                dxfattribs={"layer": "estribos"}
            )
            y += paso

    elif columna.get("tipo_armado") == "sunchos" and paso > 0:
        # espiral simplificado
        y = y_base + rec
        toggle = True
        while y < y_top - rec:
            if toggle:
                msp.add_line((x0 - h/2 + rec, y), (x0 + h/2 - rec, y + paso),
                             dxfattribs={"layer": "sunchos"})
            else:
                msp.add_line((x0 + h/2 - rec, y), (x0 - h/2 + rec, y + paso),
                             dxfattribs={"layer": "sunchos"})
            y += paso
            toggle = not toggle

    # texto de armado (longitudinal + estribos)
    diam_long = columna.get("diam_long", 0)  # ya en metros
    n_barras  = columna.get("n_barras", 0)
    armado_txt = ""
    if diam_long > 0 and n_barras > 0:
        armado_txt = f"{n_barras}Ø{int(diam_long*100)}"

    diam_estribo = columna.get("diam_estribo", 0)
    paso_estribo = columna.get("paso_estribo", 0)
    if diam_estribo > 0 and paso_estribo > 0:
        armado_txt += f"  est Ø{int(diam_estribo*100)} c/{int(paso_estribo*100)}"

    if armado_txt:
        mtext2 = msp.add_mtext(
            armado_txt,
            dxfattribs={"layer": "texto", "char_height": TEXTO_ALTURA, "attachment_point": 3}
        )
        mtext2.set_location((x0 + desplazamiento_x, y_base + desplazamiento_y - 0.15))

def dibujar_armaduras_viga(msp, nombre_viga, tramo, cotas, nivel_col, carpeta_txt="salidas/vigas"):
    nombre_archivo = f"planilla_{PORTICO}_{nombre_viga}_{tramo['id']}.txt"
    path_txt = Path(carpeta_txt) / nombre_archivo
    if not path_txt.exists():
        print("No se encontró planilla:", path_txt)
        return

    datos = leer_planilla_tramo(path_txt)
    x1, x2 = tramo["x_inicio"], tramo["x_fin"]
    y_centro = cotas[nivel_col][1]

    if datos["h"] is None or datos["b"] is None:
        print("⚠️ No se pudo leer b/h en:", path_txt)
        return

    y_sup = y_centro + datos["h"]/2
    y_inf = y_centro - datos["h"]/2
    rec   = datos.get("rec", 0.04)

    # ============================
    # Contorno y recubrimientos
    # ============================
    dibujar_contorno_viga(msp, x1, x2, y_centro, datos["h"])
    dibujar_recubrimientos(msp, x1, x2, y_sup, y_inf, rec)

    # ============================
    # Barras en despiece (llamada a la auxiliar)
    # ============================
    dibujar_barras_despiece(
        msp,
        datos["barras"],
        x1,
        x2,
        y_inf,
        y_sup,   # <-- este faltaba
        rec,
        datos.get("pi_izq"),
        datos.get("pi_der")
    )

    # ============================
    # Estribos
    # ============================
    for e in datos["estribos"]:
        dibujar_estribo(msp, e, x1, x2, y_sup, y_inf, rec, tramo["longitud_m"], datos["zonas"])

    # ============================
    # Puntos de inflexión
    # ============================
    dibujar_puntos_inflexion(msp, x1, y_sup, y_inf, datos.get("pi_izq"), datos.get("pi_der"))
    
    # ============================
    # Nombre del tramo con dimensiones
    # ============================
    b_cm = int(datos["b"] * 100)  # pasar a cm
    h_cm = int(datos["h"] * 100)  # pasar a cm
    xc = (x1 + x2) / 2
    texto_viga = f"{tramo['id']} ({b_cm}x{h_cm})"

    mtext_viga = msp.add_mtext(
        texto_viga,
        dxfattribs={
            "layer": "texto",
            "char_height": TEXTO_ALTURA_TITULO,
            "attachment_point": 5  # Middle Center
        }
    )
    mtext_viga.set_location((xc, y_sup + 0.50))  # misma altura que PI

               
# ======================
# MAIN
# ======================

def main():
    # 1. Leer datos base
    bases_json   = leer_json(BASES_PATH)
    filas_csv    = leer_csv(COLUMNAS_PATH, PORTICO)
    columnas_csv = {f["Columna"]: adaptar_columna(f) for f in filas_csv}
    estructura   = leer_estructura(ESTRUCTURA_PATH)
    portico_data = estructura[PORTICO]

    cotas = calcular_cotas(portico_data)
    doc, msp = crear_doc()

    # 2. Dibujar bases
    for nombre, base in bases_json.items():
        if nombre == "vigas_fundacion":
            continue
        x0 = portico_data["bases"][nombre]["x"]
        dibujar_base_lateral(msp, base, nombre, x0)

    # 3. Dibujar columnas
    for nombre, col_data in portico_data["columnas"].items():
        # acá ya tenés la columna adaptada
        col = columnas_csv.get(nombre, {}).copy()

        # completar con lo del JSON
        col["nombre"] = nombre
        col["x"] = col_data["x"]
        col["altura_m"] = col_data["altura_m"]
        #print("Columna:", col)
        y_base, _ = cotas.get(nombre[:2], (0,0))
        y_top = y_base + col["altura_m"]

        dibujar_columna_lateral(msp, col, col["x"], y_base, y_top)


    # 4. Dibujar vigas desde planillas TXT
    for nombre_viga, viga in portico_data["vigas"].items():
        nivel_col = "C" + nombre_viga[1]
        for tramo in viga["tramos"]:
            dibujar_armaduras_viga(msp, nombre_viga, tramo, cotas, nivel_col, carpeta_txt=PLANILLAS_VIGAS_DIR)
    # 5. Extras
    x_min = min(col["x"] for col in portico_data["columnas"].values())
    x_max = max(col["x"] for col in portico_data["columnas"].values())
    msp.add_line((x_min-1, 0), (x_max+1, 0), dxfattribs={"layer": "suelo"})
    
    # 6. Rotulación general
    portico_nombre = f"Pórtico {PORTICO}"
    autor = "Arq. Claudio E. Tomatis"
    fecha = datetime.date.today().strftime("%d/%m/%Y")

    cartucho_txt = (
        f"{portico_nombre}\n"
        f"{autor}\n"
        f"Fecha: {fecha}"
    )

    mtext_cartucho = msp.add_mtext(
        cartucho_txt,
        dxfattribs={"layer": "texto", "char_height": TEXTO_ALTURA_TITULO, "attachment_point": 1}
    )
    # calcular el mínimo Y de los textos de base
    y_min = min(-(base["espesor_m"] + base["tronco"]["altura_m"]) - 0.40
                for base in bases_json.values()
                if isinstance(base, dict) and "tronco" in base)

    # ubicar el cartucho justo debajo del último detalle
    mtext_cartucho.set_location((x_max + 1.0, y_min - 0.5))



    # 7. Guardar DXF
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    doc.saveas(OUT_PATH)
    print(f"DXF lateral generado: {OUT_PATH}")

if __name__ == "__main__":
    main()

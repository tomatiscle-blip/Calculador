import json
from PIL import Image
import xml.etree.ElementTree as ET
import math
import os
from datetime import datetime


# ==========================
# Cargar archivos JSON
# ==========================
with open("datos/materiales.json", "r", encoding="utf-8") as f:
    materiales = json.load(f)

with open("datos/viguetas.json", "r", encoding="utf-8") as f:
    viguetas = json.load(f)

# ==========================
# Ingreso nombre de la losa
# ==========================
losa_nombre = input("Nombre de la losa: ")

# ==========================
# Selección de materiales
# ==========================
seleccion = []
for categoria in ["Pisos", "Contrapiso", "Cielorrasos"]:
    print(f"\n{categoria}:")
    for i, mat in enumerate(materiales[categoria], start=1):
        print(f"{i}. {mat['nombre']}")
    
    elegidos = input(f"Seleccione los números del {categoria} separados por coma (o Enter para ninguno): ")
    if elegidos.strip():
        indices = [int(x.strip())-1 for x in elegidos.split(",")]
        for idx in indices:
            seleccion.append({"categoria": categoria, "nombre": materiales[categoria][idx]['nombre']})

# ==========================
# Función para calcular D2
# ==========================
def calcular_D2(elementos, materiales, categorias_busqueda):
    D2_total = 0
    for elem in elementos:
        encontrado = False
        nombre_elem = elem.get("nombre", "")
        for categoria in categorias_busqueda:
            for mat in materiales[categoria]:
                if mat.get("nombre", "").strip() == nombre_elem.strip():
                    encontrado = True
                    if mat["tipo"] == "volumetrico":
                        espesor = mat.get("espesor_m", 0.1)  # puedes ajustar según necesidad
                        D2_total += mat["valor"] * espesor
                    else:
                        D2_total += mat["valor"]
        if not encontrado:
            print(f"⚠️ Elemento no encontrado en JSON: {nombre_elem}")
    return D2_total
# ==========================

D1 = 1.881  # peso propio estándar kN/m²
categorias_a_buscar = ["Pisos", "Contrapiso", "Cielorrasos"]
D2 = calcular_D2(seleccion, materiales, categorias_a_buscar)
D = D1 + D2

# ==========================
# Selección de sobrecarga
# ==========================
print("\nTipos de sobrecarga disponibles:")
for i, s in enumerate(materiales["Sobrecargas"], start=1):
    print(f"{i}. {s['uso']} → {s['valor_kNm2']} kN/m²")

seleccion_sobrecarga = int(input("Seleccione el número de sobrecarga: ")) - 1
L = materiales["Sobrecargas"][seleccion_sobrecarga]['valor_kNm2']
uso_sobrecarga = materiales["Sobrecargas"][seleccion_sobrecarga]['uso']
print(f"Sobrecarga seleccionada: {uso_sobrecarga}")

# ==========================
# Combinaciones de carga
# ==========================
Q1 = 1.2*D + 1.6*L
Q2 = 1.4*D
if Q1 >= Q2:
    Q = Q1
    combo_critico = "1.2D + 1.6L"
else:
    Q = Q2
    combo_critico = "1.4D"

print(f"Combinación crítica: {combo_critico} → {Q:.3f} kN/m²")

# ==========================
# Luz de la losa/vigueta
# ==========================
luz_libre = float(input("Ingrese luz libre de la losa (m): "))
luz_calculo = round(luz_libre + 0.10, 1)  # metros, con tolerancia

# ==========================
# Selección de serie y momento
# ==========================
def buscar_serie_mas_cercana(luz_calculo, largo_a_serie):
    # Convertir las claves del JSON a float pero guardar también la original
    claves = [(float(k), k) for k in largo_a_serie.keys()]
    # Ordenar por valor numérico
    claves.sort(key=lambda x: x[0])
    # Buscar la más cercana, pero si hay empate, elegir la menor
    clave_cercana_float, clave_original = min(
        claves,
        key=lambda x: (abs(x[0] - luz_calculo), x[0])  # segundo criterio: menor valor
    )
    return largo_a_serie[clave_original], clave_cercana_float

# Uso
serie, luz_usada = buscar_serie_mas_cercana(luz_calculo, viguetas["largo_a_serie"])
print(f"Serie seleccionada: {serie} (usando luz {luz_usada} m)")

# ==========================
# Factor de reducción fijo
# ==========================
phi = 0.9  # factor de reducción del hormigón

# ==========================
# Momento requerido por la vigueta
# ==========================
# Convertimos a kg·m para coincidir con la tabla de viguetas
momento_req = (Q * luz_calculo**2 / 8) * 101.97  # kN·m → kg·m

# Ajuste para la búsqueda: aumentar el momento requerido para tener en cuenta φ
momento_para_busqueda = momento_req / phi

# ==========================
# Función para elegir la vigueta que cumpla sin pasarse
# ==========================
def seleccionar_vigueta_optima(serie_mom, momento_req):
    opciones = []
    for tipo in serie_mom:  # recorre 'vigueta_simple' y 'vigueta_doble'
        for combo, mom in serie_mom[tipo].items():  # recorre combinaciones bovedilla/losa
            if mom >= momento_req:  # solo los que cumplen
                opciones.append((mom, tipo, combo))
    if not opciones:
        return None  # ninguna opción suficiente
    # Elegir la que tenga el menor exceso sobre el momento requerido
    mom_sel, tipo_sel, combo_sel = min(opciones, key=lambda x: x[0] - momento_req)
    return tipo_sel, combo_sel, mom_sel

# ==========================
# Selección de vigueta óptima según momento requerido
# ==========================
resultado = seleccionar_vigueta_optima(viguetas[serie], momento_para_busqueda)
if resultado:
    tipo, combo, mom_disponible = resultado
    # Aplicar φ para reportar los resultados finales
    momento_reducido_kgm = mom_disponible * phi       # en kg·m
    momento_reducido_kNm = momento_reducido_kgm / 101.97  # opcional, en kN·m
else:
    raise ValueError("⚠️ Ninguna vigueta de la serie cumple con el momento requerido")

# Extraer altura de bovedilla del combo
# combo = "bovedilla_12_losa_17"
bovedilla_altura = combo.split("_")[1]

# Generar nombre de imagen
imagen_path = f"imagenes/{tipo}_bovedilla_{bovedilla_altura}.png"

# Abrir imagen
from PIL import Image
try:
    img = Image.open(imagen_path)
    img.show()
    print(f"Imagen seleccionada: {imagen_path}")
except FileNotFoundError:
    print(f"⚠️ No se encontró la imagen: {imagen_path}")

# ==========================
# Cómputo de materiales desde materiales.json
# ==========================

# Buscar datos en JSON
# Extraer altura de bovedilla del combo
bovedilla_altura = combo.split("_")[1]
if "_" in bovedilla_altura:
    bovedilla_altura = bovedilla_altura.split("_")[0]

# Normalizar tipo: "vigueta_doble" → "doble"
tipo_normalizado = tipo.replace("vigueta_", "")

# Buscar datos en JSON
datos = next(
    (item for item in materiales["Viguetas_Bovedillas"]
     if item["tipo_vigueta"] == tipo_normalizado and item["bovedilla"] == int(bovedilla_altura)),
    None
)

if datos:
    ancho_losa = float(input("Ingrese ancho de la losa (m): "))
    area_losa = luz_calculo * ancho_losa

    # Viguetas: por metro de ancho (2.00 simple, 3.17 doble)
    n_viguetas = round(ancho_losa * datos["viguetas_por_m2"])

    # Bovedillas: por m²
    n_bloques = round(area_losa * datos["bloques_por_m2"])

    # Hormigón: por m²
    vol_hormigon = area_losa * datos["hormigon_m3_m2"]

    # Extraer altura de bovedilla desde el combo (ej: "bovedilla_12_losa_17")
    bovedilla_altura = combo.split("_")[1]

    

    print("\n=== CÓMPUTO DE MATERIALES ===")
    print(f"- Área de losa: {area_losa:.2f} m²")
    print(f"- Viguetas: {n_viguetas} unidades de {luz_calculo:.2f} m")
    print(f"- Bovedillas: {n_bloques} unidades EPS n°{bovedilla_altura}")
    print(f"- Hormigón: {vol_hormigon:.3f} m³ (capa de compresión)")
    
else:
    print("⚠️ No se encontraron datos en materiales.json para esa combinación")

# ==========================
# Cómputo de malla electrosoldada (criterio por m²)

def computo_malla(
    ancho_losa,
    lc,
    tamanio_pano="2.40x6",
    solape_m=0.15,
    umbral_fraccion=0.50
):
    # --- Tamaños comerciales ---
    if tamanio_pano == "2.40x6":
        pano_ancho, pano_largo = 2.40, 6.00
    elif tamanio_pano == "2.40x3":
        pano_ancho, pano_largo = 2.40, 3.00
    else:
        raise ValueError(f"Tamaño de paño no reconocido: {tamanio_pano}")

    # --- Áreas ---
    area_losa = ancho_losa * lc

    # Área efectiva del paño (descontando solape)
    area_pano_efectiva = (pano_ancho - solape_m) * (pano_largo - solape_m)

    # --- Cómputo teórico ---
    n_teorico = area_losa / area_pano_efectiva
    n_enteros = int(n_teorico)
    fraccion = n_teorico - n_enteros

    # --- Regla de decisión ---
    if n_enteros == 0:
        n_mallas = 1
        nota = None
    elif fraccion <= umbral_fraccion:
        n_mallas = n_enteros
        nota = (
            "Área remanente de malla a completar con barras Ø6 "
            "según criterio habitual de obra."
        )
    else:
        n_mallas = n_enteros + 1
        nota = None

    area_total_mallas = n_mallas * (pano_ancho * pano_largo)

    return {
        "modelo": "Malla electrosoldada Q-131 / R-131 (SIMA)",
        "tamanio_pano": tamanio_pano,
        "solape_m": solape_m,
        "area_losa_m2": round(area_losa, 2),
        "area_pano_efectiva_m2": round(area_pano_efectiva, 2),
        "cantidad_mallas": n_mallas,
        "area_total_mallas_m2": round(area_total_mallas, 2),
        "criterio": "Cómputo por área (m²)",
        "nota": nota
    }
# ==========================

def computo_nervios_refuerzo(luz_calculo, ancho_losa,
                             max_tramo_m=1.80, barras_por_nervio=2, diametro_mm=8):
    # Nervios para que no haya tramos > 1.80 m
    n_nervios = max(math.ceil(luz_calculo / max_tramo_m) - 1, 0)

    # Cada nervio recorre el ancho de la losa
    largo_por_nervio_m = ancho_losa

    barras_total = n_nervios * barras_por_nervio
    longitud_total_barras_m = barras_total * largo_por_nervio_m

    return {
        "n_nervios": n_nervios,
        "barras_por_nervio": barras_por_nervio,
        "diametro_mm": diametro_mm,
        "largo_por_nervio_m": largo_por_nervio_m,
        "barras_total": barras_total,
        "longitud_total_barras_m": longitud_total_barras_m
    }

def generar_memoria(losa_nombre, luz_libre, luz_calculo,
                    D1, D2, D, L, Q,
                    serie, tipo, combo,
                    momento_req, momento_para_busqueda,
                    mom_disponible, phi,
                    momento_reducido_kgm, momento_reducido_kNm,
                    seleccion, materiales, seleccion_sobrecarga,
                    area_losa, n_viguetas, n_bloques, vol_hormigon,
                    combo_critico, ancho_losa):
    memoria = []
    memoria.append("MEMORIA DE CÁLCULO - LOSA ALIVIANADA\n")
    memoria.append("=================================\n\n")

    # 1) Datos generales
    memoria.append(f"Losa: {losa_nombre}\n")
    memoria.append(f"Luz libre: {luz_libre:.2f} m\n")
    memoria.append(f"Luz de cálculo (con tolerancia): {luz_calculo:.2f} m\n\n")

    # 2) Materiales seleccionados (cargas accesorias)
    subtotal_D2 = 0.0
    memoria.append("Materiales seleccionados:\n")
    for sel in seleccion:
        categoria = sel["categoria"]
        nombre = sel["nombre"]
        for mat in materiales[categoria]:
            if mat["nombre"] == nombre:
                if mat["tipo"] == "volumetrico":
                    espesor = mat.get("espesor_m", 0.1)
                    carga = mat["valor"] * espesor
                else:
                    carga = mat["valor"]
                subtotal_D2 += carga
                memoria.append(f" - {categoria}: {nombre} → {carga:.3f} kN/m²\n")
    memoria.append(f"Subtotal cargas accesorias (D2): {subtotal_D2:.3f} kN/m²\n")
    memoria.append(f"Total cargas permanentes (D1+D2): {D:.3f} kN/m²\n\n")

    # 3) Cargas consideradas
    memoria.append("Cargas consideradas:\n")
    memoria.append(f" - D1 (peso propio): {D1:.3f} kN/m²\n")
    memoria.append(f" - D2 (accesoria): {D2:.3f} kN/m²\n")
    memoria.append(f" - D total: {D:.3f} kN/m²\n")
    uso = materiales["Sobrecargas"][seleccion_sobrecarga]["uso"]
    memoria.append(f" - Sobrecarga seleccionada: {uso} → {L:.3f} kN/m²\n")
    memoria.append(f" - Combinación crítica: {combo_critico} → {Q:.3f} kN/m²\n\n")

    # 4) Selección de vigueta
    memoria.append("Selección de vigueta:\n")
    memoria.append(f" - Serie de vigueta: {serie}\n")
    memoria.append(f" - Tipo de vigueta seleccionada: {tipo}\n")
    memoria.append(f" - Combinación bovedilla/losa: {combo}\n")
    memoria.append(f" - Momento requerido (sin φ): {momento_req:.2f} kg·m\n")
    memoria.append(f" - Momento requerido para búsqueda (con φ aplicado): {momento_para_busqueda:.2f} kg·m\n")
    memoria.append(f" - Momento disponible de la vigueta seleccionada: {mom_disponible:.2f} kg·m\n")
    memoria.append(f" - Momento reducido aplicado φ={phi}: {momento_reducido_kgm:.2f} kg·m ({momento_reducido_kNm:.2f} kN·m)\n")
    memoria.append(" - Referencia de capacidad: Tabla 4 – Tensolite\n\n")

    # 5) Cómputo de materiales
    memoria.append("Cómputo de materiales:\n")
    memoria.append(f" - Área de losa: {area_losa:.2f} m²\n")
    memoria.append(f" - Viguetas: {n_viguetas} unidades de {luz_calculo:.2f} m\n")
    bovedilla_altura = combo.split("_")[1]
    memoria.append(f" - Bovedillas: {n_bloques} unidades EPS n°{bovedilla_altura}\n")
    memoria.append(f" - Hormigón: {vol_hormigon:.3f} m³ (capa de compresión)\n\n")

    # Malla electrosoldada
    malla = computo_malla(ancho_losa, luz_calculo, tamanio_pano="2.40x6", solape_m=0.20)

    memoria.append("Malla electrosoldada SIMA Q-131 / R-131:\n")
    memoria.append(f" - Paño comercial: {malla['tamanio_pano']} m (solape {malla['solape_m']:.2f} m)\n")
    memoria.append(f" - Área losa: {malla['area_losa_m2']:.2f} m²\n")
    memoria.append(f" - Área efectiva de paño: {malla['area_pano_efectiva_m2']:.2f} m²\n")
    memoria.append(f" - Cantidad adoptada: {malla['cantidad_mallas']} unidad(es)\n")

    #  Nota opcional si el área remanente se completa con barras Ø6
    if malla.get("nota"):
        memoria.append(f" - Nota: {malla['nota']}\n")

    memoria.append("\n")
    
    # Nervios de refuerzo
    nervios = computo_nervios_refuerzo(luz_calculo, ancho_losa)
    memoria.append("Nervios de refuerzo:\n")
    memoria.append(f" - Tramo máximo: 1.80 m\n")
    memoria.append(f" - Nervios requeridos: {nervios['n_nervios']} unidades\n")
    memoria.append(f" - Especificación por nervio: {nervios['barras_por_nervio']} barras Ø{nervios['diametro_mm']} "
                   f"de {nervios['largo_por_nervio_m']:.2f} m\n")
    memoria.append(f" - Total barras: {nervios['barras_total']} unidades\n")
    memoria.append(f" - Longitud total de barras: {nervios['longitud_total_barras_m']:.2f} m\n\n")

    memoria.append("=================================\n")
    memoria.append("Fin de la memoria de cálculo\n")
    return memoria
# ==========================
# Generar memoria técnica
# ==========================
memoria = generar_memoria(
    losa_nombre, luz_libre, luz_calculo,
    D1, D2, D, L, Q,
    serie, tipo, combo,
    momento_req, momento_para_busqueda,
    mom_disponible, phi,
    momento_reducido_kgm, momento_reducido_kNm,
    seleccion, materiales, seleccion_sobrecarga,
    area_losa, n_viguetas, n_bloques, vol_hormigon,
    combo_critico, ancho_losa
)

# --------------------------
# Guardado en carpeta salidas (sin sobreescribir)
# --------------------------
SALIDAS_DIR = "salidas"
os.makedirs(SALIDAS_DIR, exist_ok=True)

timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
nombre_archivo = f"memoria_losa_{losa_nombre}_{timestamp}.txt"
ruta_salida = os.path.join(SALIDAS_DIR, nombre_archivo)

with open(ruta_salida, "w", encoding="utf-8") as f:
    f.writelines(memoria)

print(f"📄 Memoria técnica generada en '{ruta_salida}'")
# ==========================

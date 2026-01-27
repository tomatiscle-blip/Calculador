from pathlib import Path
import json
import re
from anastruct import SystemElements

# -------------------------------
# 1. Carpeta de salidas
# -------------------------------
SALIDAS_DIR = Path("salidas") / "analisis_cargas"

# Listar archivos disponibles
archivos = list(SALIDAS_DIR.glob("*.txt"))
print("Archivos disponibles:")
for i, archivo in enumerate(archivos, 1):
    print(f"{i}. {archivo.name}")

# Elegir archivo con input
opcion = int(input("\nIngrese el número del archivo a analizar: "))
archivo_elegido = archivos[opcion - 1]
print("\nUsando:", archivo_elegido.name)

# Leer contenido
with open(archivo_elegido, "r", encoding="utf-8") as f:
    contenido = f.read()


# -------------------------------
# 2. Función para extraer valores
# -------------------------------
def extraer_valor(patron, texto):
    match = re.search(patron, texto)
    return float(match.group(1)) if match else None

D = extraer_valor(r"D total\s*=\s*([\d.]+)", contenido)
L = extraer_valor(r"L total\s*=\s*([\d.]+)", contenido)
W = extraer_valor(r"W total\s*=\s*([\d.]+)", contenido)

comb1 = extraer_valor(r"1.4D:\s*([\d.]+)", contenido)
comb2 = extraer_valor(r"1.2D\+1.6L:\s*([\d.]+)", contenido)
comb3 = extraer_valor(r"1.2D\+0.5L\+1.6W:\s*([\d.]+)", contenido)
comb4 = extraer_valor(r"0.9D\+1.6W:\s*([\d.]+)", contenido)

print("\nValores extraídos:")
print("D =", D, "kN/m")
print("L =", L, "kN/m")
print("W =", W, "kN/m")
print("Combinaciones:", comb1, comb2, comb3, comb4)

# -------------------------------
# 3. Carga de diseño crítica
# -------------------------------
q_design = max(comb1, comb2, comb3, comb4)
# Carga de servicio (rigidez)
q_service = D + L   # kN/m, sin mayorar

print(f"\nCarga crítica (resistencia) = {q_design:.2f} kN/m")
print(f"Carga de servicio (rigidez)   = {q_service:.2f} kN/m")


# -------------------------------
# 4. Modelo estructural
# -------------------------------
# Lista de longitudes de cada vano (pueden ser distintas)
lc_list = [4.05] # 4.05, 3.80, 3.90
n_vanos = len(lc_list)

ss = SystemElements()

# Crear elementos consecutivos con longitudes variables
x = 0.0
for i, L in enumerate(lc_list):
    ss.add_element(location=[[x, 0], [x + L, 0]])
    x += L

# Apoyos
ss.add_support_hinged(node_id=1)          # primer nodo
for i in range(2, n_vanos+1):             # intermedios
    ss.add_support_roll(node_id=i)
ss.add_support_roll(node_id=n_vanos+1)    # último nodo

# Cargas distribuidas en cada vano
for i in range(1, n_vanos+1):
    ss.q_load(element_id=i, q=-q_design)

# Resolver
ss.solve()
ss.show_structure()
ss.show_bending_moment()
ss.show_shear_force()


# Resultados del elemento
res = ss.get_element_results(element_id=1)
Mmax = res["Mmax"]   # kN·m
Mmin = res["Mmin"]   # kN·m
Qmax = res["Qmax"]   # kN
Qmin = res["Qmin"]   # kN

print("\nResultados del modelo:")
print(f"Momento máximo = {Mmax:.2f} kN·m")
print(f"Momento mínimo = {Mmin:.2f} kN·m")
print(f"Cortante máximo = {Qmax:.2f} kN")
print(f"Cortante mínimo = {Qmin:.2f} kN")

# Momento crítico (el mayor en valor absoluto)
Msol = max(abs(Mmax), abs(Mmin))  # kN·m
fy = 235  # MPa, acero típico
# Calcular Wx requerido
Msol_Nmm = Msol * 1e6             # convertir kN·m → N·mm
Wreq_mm3 = Msol_Nmm / fy          # mm³
Wreq_cm3 = Wreq_mm3 / 1000        # cm³

print(f"\nWx requerido = {Wreq_cm3:.2f} cm³")


# -------------------------------
# 5. Leer perfiles metálicos JSON
# -------------------------------
archivo_perfiles = Path("datos") / "perfiles_metalicos.json"
with open(archivo_perfiles, "r", encoding="utf-8") as f:
    data = json.load(f)

fy = 235  # MPa, acero típico

# Conversión de unidades
def convertir_modulo_resistente(valor_cm3):
    return float(valor_cm3) * 1e-6  # cm³ → m³

# -------------------------------
# 6. Iterar solo sobre la primera tabla (PERFILES C)
# -------------------------------
print("\n--- Verificación de perfiles (solo PERFILES C) ---")

# Momento solicitado del modelo (kN·m)
Msol = max(abs(Mmax), abs(Mmin))

# Wx requerido en cm³
Wreq_cm3 = (Msol * 1e6) / fy   # N·mm / (N/mm²) = mm³
Wreq_cm3 = Wreq_cm3 / 1000     # mm³ → cm³


# Buscar la tabla que tenga el título "PERFILES C"
tabla_c = next((t for t in data["tablas"] if "PERFILES C" in t["titulo_tabla"]), None)

if tabla_c:
    print("\nTabla:", tabla_c["titulo_tabla"])
    print("Norma:", tabla_c.get("norma_o_tipo_general", ""))
    print("Notas:", tabla_c.get("notas_adicionales", ""))

    perfil_elegido = None
    ESPESOR_MAX = 2.0  # mm

    for perfil in tabla_c["perfiles"]:
        nombre = perfil["nombre"]
        Wx_cm3 = float(perfil["propiedades"]["modulo_resistente_x"])
        peso = float(perfil["propiedades"]["peso_lineal"])
        espesor_alma = float(perfil["dimensiones"]["e_espesor_alma"])

        # Filtro: solo perfiles con espesor <= 2.0 mm
        if espesor_alma <= ESPESOR_MAX and Wx_cm3 >= Wreq_cm3:
            perfil_elegido = (nombre, Wx_cm3, peso, espesor_alma)
            break  # nos quedamos con el primero que cumple

    if perfil_elegido:
        nombre, Wx_cm3, peso, espesor = perfil_elegido
        
        print(f"\nPerfil elegido: {nombre}")
        print(f"Wx = {Wx_cm3:.2f} cm³, Peso = {peso:.2f} kg/m, Espesor = {espesor:.2f} mm")
    
        #   Buscar el perfil completo en la tabla por nombre
        perfil_data = next(p for p in tabla_c["perfiles"] if p["nombre"] == nombre)
    
        # Ahora tenés acceso a todas las propiedades
        Ix_cm4 = float(perfil_data["propiedades"]["inercia_x"])
        Ix_mm4 = Ix_cm4 * 10000  # cm⁴ → mm⁴
    
        # Parámetros de cálculo
        E = 210000  # MPa = N/mm²
        L_vano = lc_list[0]      # m
        L_vano_mm = L_vano * 1000  # convertir a mm

        # Conversión a N/mm
        q_service_Nmm = q_service * 1000 / 1000 # kN/m → N/mm

        # Flecha máxima con carga de servicio
        delta_mm = (5 * q_service_Nmm * (L_vano_mm**4)) / (384 * E * Ix_mm4)
        # Límite normativo
        delta_adm_mm = L_vano_mm / 300

        # Relación L/δ
        relacion = L_vano_mm / delta_mm
        print("\n--- Verificación de flecha ---")
        print(f"Carga crítica (resistencia) = {q_design:.2f} kN/m")
        print(f"Carga de servicio (rigidez) = {q_service:.2f} kN/m")
        print(f"δ calculada (servicio) = {delta_mm:.2f} mm")
        print(f"δ admisible = {delta_adm_mm:.2f} mm (L/300)")
        print(f"Relación L/δ = L/{relacion:.0f}")
        if delta_mm <= delta_adm_mm:
            print("✅ Cumple rigidez")
        else:
            print("❌ No cumple rigidez")
            Ix_req_mm4 = (5 * q_service_Nmm * (L_vano_mm**4)) / (384 * E * delta_adm_mm)
            Ix_req_cm4 = Ix_req_mm4 / 10000  # mm4 → cm4
            print(f"Inercia mínima requerida = {Ix_req_cm4:.2f} cm⁴")
        # -------------------------------
        # 7. Guardar resultados en el archivo TXT
        # -------------------------------
        # Retroalimentar el archivo TXT con el perfil elegido
        with open(archivo_elegido, "a", encoding="utf-8") as f:
            f.write("\n--- Resultado final ---\n")
            f.write(f"Perfil elegido: {nombre}\n")
            f.write(f"Wx = {Wx_cm3:.2f} cm³, Peso = {peso:.2f} kg/m, Espesor = {espesor:.2f} mm\n")
            f.write(f"δ calculada (servicio) = {delta_mm:.2f} mm\n")
            f.write(f"δ admisible = {delta_adm_mm:.2f} mm (L/300)\n")
            f.write(f"Relación L/δ = L/{relacion:.0f}\n")
            if delta_mm <= delta_adm_mm:
                f.write("✅ Cumple rigidez\n")
            else:
                f.write("❌ No cumple rigidez\n")
                f.write(f"Inercia mínima requerida = {Ix_req_cm4:.2f} cm⁴\n")



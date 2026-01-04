import json
from pathlib import Path
from datetime import datetime

# ==========================
# Cargar catálogo de materiales
# ==========================
def cargar_catalogo(ruta="datos/materiales.json"):
    with open(Path(ruta), "r", encoding="utf-8") as f:
        return json.load(f)

# ==========================
# Utilidades
# ==========================
def buscar_material(catalogo, categoria, nombre):
    for m in catalogo.get(categoria, []):
        if m.get("nombre", "").strip() == nombre.strip():
            return m
    return None

def ensure_dir(path: Path):
    path.mkdir(parents=True, exist_ok=True)

# ==========================
# Funciones de cálculo lineal
# ==========================
def carga_losa(material, ancho):
    return material["valor"] * ancho   # kN/m² × m = kN/m

def carga_muro(material, espesor, altura):
    return material["valor"] * espesor * altura   # kN/m³ × m × m = kN/m

def carga_encadenado(material, ancho, alto, cantidad=1):
    return material["valor"] * (ancho * alto) * cantidad   # kN/m³ × m² × cantidad

def carga_techo(material, ancho):
    return material["valor"] * ancho

def carga_cielorraso(material, ancho):
    return material["valor"] * ancho

# ==========================
# Combinaciones
# ==========================
def combinaciones_cargas(D_kN_m, L_kN_m, Q_kN_m=0.0):
    return {
        "1.4D": 1.4 * D_kN_m,
        "1.2D + 1.6L": 1.2 * D_kN_m + 1.6 * L_kN_m,
        "1.2D + 1.6Q": 1.2 * D_kN_m + 1.6 * Q_kN_m,
        "0.9D + 1.6Q": 0.9 * D_kN_m + 1.6 * Q_kN_m
    }

# ==========================
# Inventario con nombres automáticos
# ==========================
class InventarioCargas:
    def __init__(self):
        self.cargas = {}
        self.contador = {"qlosa":0, "qmuro":0, "qencadenado":0, "qtecho":0, "qcielorraso":0}

    def agregar(self, tipo, valor_kN_m, origen):
        self.contador[tipo] += 1
        nombre = f"{tipo}{self.contador[tipo]}"
        self.cargas[nombre] = {
            "valor_kN_m": round(valor_kN_m, 4),
            "origen": origen
        }
        return nombre

    def total_D(self):
        return round(sum(v["valor_kN_m"] for v in self.cargas.values()), 4)

# ==========================
# Flujo principal (solo cargas y combinaciones)
# ==========================
if __name__ == "__main__":
    catalogo = cargar_catalogo()

    # Nombre del cálculo de cargas (no de la viga)
    nombre_calculo = input("Nombre del cálculo de cargas: ").strip()
    ancho_tributario = float(input("Ancho tributario (m): ").strip())

    inv = InventarioCargas()

    # Ejemplos: reemplazá por tus selecciones reales
    mat_piso = buscar_material(catalogo, "Pisos", "Porcelanato 10mm + pegamento")
    if mat_piso:
        inv.agregar("qlosa", carga_losa(mat_piso, ancho_tributario),
                    {"categoria":"Pisos","nombre":mat_piso["nombre"],"ancho_m":ancho_tributario})

    mat_contra = buscar_material(catalogo, "Contrapiso", "Carpeta nivelacion 3cm")
    if mat_contra:
        inv.agregar("qlosa", carga_losa(mat_contra, ancho_tributario),
                    {"categoria":"Contrapiso","nombre":mat_contra["nombre"],"ancho_m":ancho_tributario})

    mat_cielo = buscar_material(catalogo, "Cielorrasos", "Yeso Suspendido")
    if mat_cielo:
        inv.agregar("qcielorraso", carga_cielorraso(mat_cielo, ancho_tributario),
                    {"categoria":"Cielorrasos","nombre":mat_cielo["nombre"],"ancho_m":ancho_tributario})

    mat_muro = buscar_material(catalogo, "Mamposteria", "Ladrillo hueco (portante)")
    if mat_muro:
        inv.agregar("qmuro", carga_muro(mat_muro, espesor=0.20, altura=2.60),
                    {"categoria":"Mamposteria","nombre":mat_muro["nombre"],"espesor_m":0.20,"altura_m":2.60})

    mat_horm = buscar_material(catalogo, "Hormigon", "Hormigon armado")
    if mat_horm:
        inv.agregar("qencadenado", carga_encadenado(mat_horm, ancho=0.20, alto=0.25, cantidad=1),
                    {"categoria":"Hormigon","nombre":mat_horm["nombre"],"ancho_m":0.20,"alto_m":0.25,"cantidad":1})

    mat_techo = buscar_material(catalogo, "Cubiertas", "Chapa trapezoidal T101")
    if mat_techo:
        inv.agregar("qtecho", carga_techo(mat_techo, ancho_tributario),
                    {"categoria":"Cubiertas","nombre":mat_techo["nombre"],"ancho_m":ancho_tributario})

    # Sobrecarga normativa L
    print("\nSobrecargas disponibles:")
    for i, sob in enumerate(catalogo.get("Sobrecargas", []), start=1):
        print(f"{i}. {sob['uso']} ({sob['valor_kNm2']} kN/m²)")
    idx = int(input("Seleccione sobrecarga: ").strip()) - 1
    sob = catalogo["Sobrecargas"][idx]
    L_lineal = round(sob["valor_kNm2"] * ancho_tributario, 4)

    # Totales y combinaciones
    D_total = inv.total_D()
    combos = combinaciones_cargas(D_total, L_lineal)

    # Persistencia en JSON
    salida_dir = Path("datos/salidas")
    ensure_dir(salida_dir)
    salida_path = salida_dir / f"{nombre_calculo}_cargas.json"

    salida = {
        "metadata": {
            "nombre_calculo": nombre_calculo,
            "timestamp": datetime.now().isoformat(timespec="seconds"),
            "ancho_tributario_m": ancho_tributario
        },
        "cargas_lineales_kN_m": inv.cargas,
        "totales": {
            "D_total_kN_m": D_total,
            "L_lineal_kN_m": L_lineal,
            "D_plus_L_kN_m": round(D_total + L_lineal, 4)
        },
        "combinaciones_kN_m": {k: round(v, 4) for k, v in combos.items()}
    }

    with open(salida_path, "w", encoding="utf-8") as f:
        json.dump(salida, f, ensure_ascii=False, indent=2)

    # Resumen en consola
    print("\nDetalle de cargas lineales (kN/m):")
    for nombre, info in inv.cargas.items():
        print(f"- {nombre}: {info['valor_kN_m']:.4f} kN/m | origen: {info['origen']}")

    print(f"\nD total = {D_total:.4f} kN/m")
    print(f"L lineal ({sob['uso']}) = {L_lineal:.4f} kN/m")

    print("\nCombinaciones:")
    for k, v in combos.items():
        print(f"{k}: {v:.4f} kN/m")

    print(f"\nGuardado en: {salida_path}")
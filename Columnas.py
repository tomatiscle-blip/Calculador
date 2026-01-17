import math

# ==============================
# Parámetros normativos CIRSOC
# ==============================
PHI = {
    "estribos": 0.65,   # columnas con estribos rectos
    "sunchos": 0.70     # columnas zunchadas (espiral)
}
ALPHA = {
    "estribos": 0.80,
    "sunchos": 0.85
}

# ==============================
# 1. Cálculo de sección mínima
# ==============================
def calcular_seccion(Pu_kN, fck, fy=420, tipo="estribos"):
    Pu = Pu_kN * 1000  # pasa a N
    phi = PHI[tipo]
    alpha = ALPHA[tipo]
    rho_g = 0.01

    # Resistencia nominal requerida
    Pn_max = Pu / phi

    # Área bruta mínima teórica
    Ag_calc = Pn_max / (alpha * (0.85*fck*(1-rho_g) + fy*rho_g))  # mm2
    Ast_min_calc = rho_g * Ag_calc

    # Ajuste normativo de armadura mínima
    if tipo == "estribos":
        Ast_min_adopt = max(Ast_min_calc, 4 * (math.pi*(12/10)**2/4)*100)  # 4Ø12
    else:
        Ast_min_adopt = max(Ast_min_calc, 6 * (math.pi*(12/10)**2/4)*100)  # 6Ø12

    return Ag_calc, Ast_min_calc, Ast_min_adopt

# ==============================
# 2. Dimensiones según tipo
# ==============================
def dimensiones(Ag, tipo):
    if tipo == "estribos":
        lado_calc = math.sqrt(Ag) / 10  # cm
        lado_adopt = max(lado_calc, 20) # mínimo 20 cm
        return lado_calc, lado_adopt, f"Sección cuadrada ≈ {lado_calc:.1f} cm (calc), adoptada ≥ {lado_adopt:.1f} cm"
    elif tipo == "sunchos":
        diam_calc = math.sqrt(4*Ag/math.pi) / 10  # cm
        diam_adopt = max(diam_calc, 20) # mínimo 20 cm
        return diam_calc, diam_adopt, f"Sección circular ≈ Ø{diam_calc:.1f} cm (calc), adoptada ≥ Ø{diam_adopt:.1f} cm"
import math

# ==============================
# 3. Condiciones de armado
# ==============================

def armado_estribos(Ast_min, lado_cm, diametros=[12, 16, 20], diam_estribo=6):
    opciones = []
    Ast_min_cm2 = Ast_min/100

    # mínimo normativo: 4Ø12
    area_barra_min = math.pi*(12/10)**2/4
    Ast_min_norma = 4 * area_barra_min
    texto_norma = f"Ast mínimo normativo ≈ {Ast_min_norma:.2f} cm² → 4Ø12"

    # probar distintos diámetros
    for d in diametros:
        area_barra = math.pi*(d/10)**2/4
        n_barras = max(4, math.ceil(Ast_min_cm2 / area_barra))
        Ast_real = n_barras * area_barra
        opciones.append((n_barras, d, Ast_real))

    # elegir la opción más ajustada (Ast_real más cercano al mínimo requerido)
    mejor = min(opciones, key=lambda x: x[2])
    n_barras_calc, diam_long, Ast_real_calc = mejor
    texto_prop = f"Armadura propuesta: {n_barras_calc}Ø{diam_long} → Ast real ≈ {Ast_real_calc:.2f} cm²"

    # Normativa diámetro mínimo de estribo
    if diam_long <= 16:
        diam_min = 6
    elif diam_long <= 25:
        diam_min = 8
    elif diam_long <= 32:
        diam_min = 10
    else:
        diam_min = 12

    cumple_diam = diam_estribo >= diam_min

    # Límites de separación
    limite1 = 12*diam_long        # mm
    limite2 = 48*diam_estribo     # mm
    limite3 = lado_cm*10          # mm
    s_max = min(limite1, limite2, limite3)

    if s_max == limite1:
        gobernante = f"12·db = {limite1/10:.1f} cm"
    elif s_max == limite2:
        gobernante = f"48·Øestribo = {limite2/10:.1f} cm"
    else:
        gobernante = f"lado menor = {limite3/10:.1f} cm"

    return (
        f"{texto_norma}\n"
        f"{texto_prop}\n"
        f"Estribos Ø{diam_estribo} mm (mínimo normativo Ø{diam_min}, cumple={cumple_diam})\n"
        f"Separación máxima ≈ {s_max/10:.1f} cm (gobierna {gobernante})"
    )


def armado_sunchos(Ast_min, diam_cm, fck, fy, diametros=[12, 16, 20], diam_espiral=10, recubrimiento=40):
    opciones = []
    Ast_min_cm2 = Ast_min/100

    # mínimo normativo: 6Ø12
    area_barra_min = math.pi*(12/10)**2/4
    Ast_min_norma = 6 * area_barra_min
    texto_norma = f"Ast mínimo normativo ≈ {Ast_min_norma:.2f} cm² → 6Ø12"

    # probar distintos diámetros
    for d in diametros:
        area_barra = math.pi*(d/10)**2/4
        n_barras = max(6, math.ceil(Ast_min_cm2 / area_barra))
        Ast_real = n_barras * area_barra
        opciones.append((n_barras, d, Ast_real))

    # elegir la opción más ajustada
    mejor = min(opciones, key=lambda x: x[2])
    n_barras_calc, diam_long, Ast_real_calc = mejor
    texto_prop = f"Armadura propuesta: {n_barras_calc}Ø{diam_long} → Ast real ≈ {Ast_real_calc:.2f} cm²"

    # núcleo confinado
    Dc = diam_cm*10 - 2*recubrimiento  # mm
    Ac = math.pi*(Dc/2)**2             # mm2
    Ag = math.pi*(diam_cm*10/2)**2     # mm2

    # cuantía volumétrica mínima
    rho_s = 0.425*((Ag/Ac)-1)*(fck/fy)

    # área de barra del espiral (cm2)
    Asp = math.pi*(diam_espiral/10)**2/4

    # relación Asp/s → cm2/m
    Asp_s = rho_s*(Dc/10/4)*100

    # paso resultante
    s = Asp/Asp_s * 100   # cm
    cumple = 2.5 <= s <= 8.0

    return (
        f"{texto_norma}\n"
        f"{texto_prop}\n"
        f"Espiral Ø{diam_espiral} mm → paso ≈ {s:.1f} cm "
        f"(cumple rango 2.5–8.0 cm: {cumple})\n"
        f"Asp/s requerido ≈ {Asp_s:.2f} cm²/m"
    )

# ==============================
# Ejemplo interactivo
# ==============================
if __name__ == "__main__":
    Pu = float(input("Ingrese Pu (kN): "))
    fck = int(input("Ingrese f'c (20,25,30 MPa): "))
    # mapeo de opciones
    TIPOS = {
        "1": "estribos",
        "2": "sunchos",
        "estribos": "estribos",
        "sunchos": "sunchos"
    }
    # en el main
    tipo = input("Tipo de columna (1=estribos / 2=sunchos): ").lower()
    tipo = TIPOS.get(tipo, "estribos")  # por defecto estribos si no coincide
  # normaliza a minúsculas

    Ag, Ast_min_calc, Ast_min_adopt = calcular_seccion(Pu, fck, tipo=tipo)

    dim_calc, dim_adopt, texto_dim = dimensiones(Ag, tipo)
    print(texto_dim)

    if tipo == "estribos":
        print(armado_estribos(Ast_min_adopt, dim_adopt))
    else:
        print(armado_sunchos(Ast_min_adopt, dim_adopt, fck, fy=420))
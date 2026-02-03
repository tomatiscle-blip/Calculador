import math
import json
import os
import csv, os

# ==============================
# Funciones de predimensión
# ==============================

def definir_fck(nombre_columna, Pu_kN):
    """
    Define el fck de la columna.
    C0 = planta baja
    C1, C2 = pisos superiores
    """

    # planta baja → permitir H25
    if nombre_columna.startswith("C0"):
        if Pu_kN > 600:      # umbral razonable
            return 25
        else:
            return 20

    # pisos superiores → casi siempre H20
    return 20

import math

def predimensionar_seccion(Pu_kN, Mu_kNm, fck_MPa):
    """
    Predimensionado de columna rectangular:
    - b = 20 cm mínimo normativo
    - h = 20 cm base
    - Ajuste por excentricidad si e > 20 cm
    """

    # Convertir unidades
    Pu = abs(Pu_kN) * 1000  # N
    Mu = abs(Mu_kNm) * 1e6  # Nmm

    # Excentricidad en cm
    if Pu > 0:
        e_cm = Mu / Pu / 10.0  # pasar de mm a cm
    else:
        e_cm = 0

    b = 20.0
    h = 20.0

    # Activación del ajuste
    if e_cm > b:
        h = b + 0.1 * e_cm

    # Redondear a 0.5 cm para practicidad
    h = round(h, 1)

    return b, h, e_cm

def predimensionar_suncho(Pu_kN, Mu_kNm, fck_MPa, diam_min=20):
    Pu = abs(Pu_kN) * 1000  # N
    Mu = abs(Mu_kNm) * 1e6  # Nmm

    # Área mínima por axial
    Ag = Pu / (0.35 * fck_MPa * 10)  # cm²

    # Diámetro inicial
    d = math.sqrt(4 * Ag / math.pi)

    # Excentricidad
    e_cm = Mu / Pu / 10 if Pu > 0 else 0

    # Ajuste moderado si e > d
    if e_cm > d:
        d = d + 0.1 * e_cm

    # Redondeo a moldes estándar
    moldes = [20, 25, 30, 35, 40]
    d_final = max(diam_min, min(moldes, key=lambda m: abs(m - d)))

    return d_final, e_cm
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
# cargar diagrama JSON
# ==============================

def cargar_diagrama(nombre):
    try:
        with open(nombre, encoding="utf-8") as f:
            return json.load(f)
    except UnicodeDecodeError:
        with open(nombre, encoding="latin-1") as f:
            return json.load(f)

# Ejemplo de uso
diagrama = cargar_diagrama("datos/diagramas_interaccion/diagramaI4_fc20_gamma080.json")

print("Título:", diagrama["title"])
print("Resumen verificación:", diagrama["verificationSummary"])

# ==============================
# Función para convertir a MPa

def convertir_a_mpa(Pu_kN, Mu_kNm, b_cm, h_cm, phi):
    # pasar a N y N·mm
    Pu_N = Pu_kN * 1000
    Mu_Nmm = Mu_kNm * 1e6

    # dimensiones en mm
    b_mm = b_cm * 10
    h_mm = h_cm * 10

    # coordenadas normalizadas en MPa
    x_val = phi * Mu_Nmm / (b_mm * h_mm**2)
    y_val = phi * Pu_N / (b_mm * h_mm)

    return round(x_val, 3), round(y_val, 3)

# ==============================
# Función para cargar diagramas de interacción  

def calcular_gamma(h_cm, recubrimiento_cm):
    """Calcula gamma = (h - 2d') / h"""
    return (h_cm - 2*recubrimiento_cm) / h_cm

def interpola_y(points, x_obj):
    """Interpola el valor de Y en una curva dada para un X objetivo."""
    for i in range(len(points)-1):
        x1, y1 = points[i]["x"], points[i]["y"]
        x2, y2 = points[i+1]["x"], points[i+1]["y"]
        if (x1 <= x_obj <= x2) or (x2 <= x_obj <= x1):
            t = (x_obj - x1) / (x2 - x1)
            return y1 + t * (y2 - y1)
    return None

def buscar_rho_por_punto(x_obj, y_obj, curves):
    """
    Devuelve la cuantía rho_g mínima que alcanza Y >= y_obj
    para un X dado. Interpola entre curvas si es necesario.
    Devuelve con 3 decimales si coincide con un dato de la curva.
    """
    curves_sorted = sorted(curves, key=lambda c: c["parameterValue"])

    # candidatos: (rho, y_res)
    candidatos = []
    for c in curves_sorted:
        y_res = interpola_y(c["points"], x_obj)
        if y_res is not None:
            candidatos.append((c["parameterValue"], y_res))

    if not candidatos:
        return None

    # buscar curvas entre las que interpolar
    abajo = None
    arriba = None
    for rho, y in candidatos:
        if y < y_obj:
            abajo = (rho, y)
        elif y >= y_obj and arriba is None:
            arriba = (rho, y)
            break

    # fuera de rango
    if arriba is None:
        return round(candidatos[-1][0], 3)
    if abajo is None:
        return round(candidatos[0][0], 3)

    # si coincide exactamente con un valor de la curva, devolver ese
    if abs(arriba[1] - y_obj) < 1e-6:
        return round(arriba[0], 3)
    if abs(abajo[1] - y_obj) < 1e-6:
        return round(abajo[0], 3)

    # interpolación lineal entre abajo y arriba
    rho_calc = abajo[0] + (y_obj - abajo[1]) * (arriba[0]-abajo[0]) / (arriba[1]-abajo[1])
    return round(rho_calc, 3)

def interp_rho(curves, rho_target):
    """Interpola una curva completa entre dos cuantías rho_g."""
    curves_sorted = sorted(curves, key=lambda c: c["parameterValue"])
    # control de rango
    if rho_target < curves_sorted[0]["parameterValue"] or rho_target > curves_sorted[-1]["parameterValue"]:
        raise ValueError(f"rho_target={rho_target} fuera del rango de curvas")
    # buscar curvas inferior y superior
    c_inf = max([c for c in curves_sorted if c["parameterValue"] <= rho_target], key=lambda c: c["parameterValue"])
    c_sup = min([c for c in curves_sorted if c["parameterValue"] >= rho_target], key=lambda c: c["parameterValue"])
    if c_inf["parameterValue"] == c_sup["parameterValue"]:
        return c_inf["points"]  # coincide exacto
    t_rho = (rho_target - c_inf["parameterValue"]) / (c_sup["parameterValue"] - c_inf["parameterValue"])
    puntos_interp = []
    for p_inf, p_sup in zip(c_inf["points"], c_sup["points"]):
        x = p_inf["x"] + t_rho*(p_sup["x"] - p_inf["x"])
        y = p_inf["y"] + t_rho*(p_sup["y"] - p_inf["y"])
        puntos_interp.append({"x": x, "y": y})
    return puntos_interp

def cargar(seccion_tipo, fc, gamma):
    """Carga el diagrama JSON correspondiente a fc y gamma."""
    gamma_str = f"{int(gamma*100):03}"
    carpeta = "datos/diagramas_interaccion"
    for archivo in os.listdir(carpeta):
        if archivo.startswith(f"diagrama{seccion_tipo}") and f"_fc{fc}_" in archivo and f"gamma{gamma_str}" in archivo:
            with open(os.path.join(carpeta, archivo), encoding="utf-8") as f:
                return json.load(f)
    raise FileNotFoundError(f"No se encontró diagrama para fc={fc}, γ={gamma}")

def cargar_diagrama_interpolado(seccion_tipo, fck, gamma_calc, x_obj, y_obj):
    gamma_normativos = [0.50, 0.60, 0.70, 0.80, 0.90]

    g_inf = max([g for g in gamma_normativos if g <= gamma_calc], default=None)
    g_sup = min([g for g in gamma_normativos if g >= gamma_calc], default=None)

    # Validación de rango
    if g_inf is None or g_sup is None:
        raise ValueError(f"γ={gamma_calc:.2f} fuera del rango normativo {gamma_normativos}")

    # Coincidencia exacta
    if g_inf == g_sup:
        print(f"[DEBUG] γ={gamma_calc:.2f} coincide con valor normativo → se usa directamente el diagrama de {g_inf}")
        diag = cargar(seccion_tipo, fck, g_inf)
        rho_calc = buscar_rho_por_punto(x_obj, y_obj, diag["curves"])
        if rho_calc is None:
            raise ValueError(f"Punto (X={x_obj:.3f}, Y={y_obj:.3f}) fuera del rango del diagrama γ={g_inf}")
        puntos_final = interp_rho(diag["curves"], rho_calc)
        return {
            "title": f"Diagrama directo: γ={gamma_calc:.2f}, ρg={rho_calc:.3f}",
            "points": puntos_final,
            "rho_calc": rho_calc
        }

    # Interpolación
    t_gamma = (gamma_calc - g_inf) / (g_sup - g_inf)
    print(f"[DEBUG] γ_calc={gamma_calc:.3f}, γ_inf={g_inf:.2f}, γ_sup={g_sup:.2f}, t_gamma={t_gamma:.3f}")
    print(f"[DEBUG] Punto objetivo: X={x_obj:.3f}, Y={y_obj:.3f}")

    diag_inf = cargar(seccion_tipo, fck, g_inf)
    diag_sup = cargar(seccion_tipo, fck, g_sup)
    print(f"[DEBUG] Diagramas cargados: inf={g_inf}, sup={g_sup}")

    rho_inf = buscar_rho_por_punto(x_obj, y_obj, diag_inf["curves"])
    rho_sup = buscar_rho_por_punto(x_obj, y_obj, diag_sup["curves"])
    print(f"[DEBUG] rho_inf={rho_inf}, rho_sup={rho_sup}")

    if rho_inf is None or rho_sup is None:
        raise ValueError(f"Punto (X={x_obj:.3f}, Y={y_obj:.3f}) fuera del rango de los diagramas")

    rho_calc = rho_inf + t_gamma * (rho_sup - rho_inf)
    print(f"[DEBUG] rho_calc={rho_calc}")

    curva_inf = interp_rho(diag_inf["curves"], rho_calc)
    curva_sup = interp_rho(diag_sup["curves"], rho_calc)
    print(f"[DEBUG] Curvas interpoladas en rho_calc")

    puntos_final = []
    for p_inf, p_sup in zip(curva_inf, curva_sup):
        x = p_inf["x"] + t_gamma*(p_sup["x"] - p_inf["x"])
        y = p_inf["y"] + t_gamma*(p_sup["y"] - p_inf["y"])
        puntos_final.append({"x": x, "y": y})

    return {
        "title": f"Interpolación doble: γ={gamma_calc:.2f}, ρg={rho_calc:.3f}",
        "points": puntos_final,
        "rho_calc": rho_calc
    }


def calcular_x_normalizado(Pu, Mu, b, h):
    """
    Coordenada X normalizada en MPa.
    Pu: carga axial solicitante [kN]
    Mu: momento solicitante [kNm]
    b, h: dimensiones de la sección [cm]
    """
    Mu_cm = Mu * 100  # kNm → kN·cm
    return 10 * Mu_cm / (b * (h**2))

def calcular_y_normalizado(Pu, b, h):
    """
    Coordenada Y normalizada en MPa.
    Pu: carga axial solicitante [kN]
    b, h: dimensiones de la sección [cm]
    """
    return 10 * Pu / (b * h)
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
    Ag_calc = Pn_max / (alpha * (0.85*fck*(1-rho_g) + fy*rho_g)) / 100  # cm²
    Ast_min_calc = rho_g * Ag_calc

    # Ajuste normativo de armadura mínima
    if tipo == "estribos":
        Ast_min_adopt = max(Ast_min_calc, 4 * (math.pi*(12/10)**2/4))  # 4Ø12
    else:
        Ast_min_adopt = max(Ast_min_calc, 6 * (math.pi*(12/10)**2/4))  # 6Ø12
    # Resistencia axial reducida φPn (kN)
    Ag_mm2 = Ag_calc * 100
    As_mm2 = Ast_min_adopt * 100
    Pn = phi * (0.85*fck*(Ag_mm2 - As_mm2) + fy*As_mm2) / 1000  # kN

    return Ag_calc, Ast_min_calc, Ast_min_adopt, Pn

# ==============================
# 2. Dimensiones según tipo
# ==============================
def dimensiones(Ag, tipo):
    if tipo == "estribos":
        lado_calc = math.sqrt(Ag)  # cm
        lado_adopt = max(lado_calc, 20) # mínimo 20 cm
        return lado_calc, lado_adopt, f"Sección cuadrada ≈ {lado_calc:.1f} cm (calc), adoptada ≥ {lado_adopt:.1f} cm"
    elif tipo == "sunchos":
        diam_calc = math.sqrt(4*Ag/math.pi)  # cm
        diam_adopt = max(diam_calc, 20) # mínimo 20 cm
        return diam_calc, diam_adopt, f"Sección circular ≈ Ø{diam_calc:.1f} cm (calc), adoptada ≥ Ø{diam_adopt:.1f} cm"
# ==============================
# 3. Condiciones de armado
# ==============================
# ==============================
# Función para componer armadura y calcular d'
# ==============================
def composicion_armadura(seccion_tipo, n_barras, diam_long, b=None, h=None,
                         d_usado=None, diam_estribo=8, agregado=19):
    """
    seccion_tipo: 'I' (rectangular 2 lados), 'II' (rectangular 4 lados), 'III' (circular)
    n_barras: número de barras longitudinales
    diam_long: diámetro de barra longitudinal (mm)
    b, h: dimensiones de la sección (cm)
    d_usado: recubrimiento ingresado por el usuario (cm)
    diam_estribo: diámetro de estribo (mm)
    agregado: tamaño máximo de agregado (mm)
    """

    # d' normativo (20 mm + estribo + barra/2)
    d_norma = 2.0 + diam_estribo/10 + diam_long/20  # cm

    # d' real efectivo = el normativo, porque es el que gobierna la geometría
    d_real = d_norma

    texto_rec = (
        f"d' usado (nominal ingresado) = {d_usado:.2f} cm\n"
        f"d' mínimo normativo = {d_norma:.2f} cm (20 mm + Øestribo={diam_estribo/10:.2f} + Øbarra/2={diam_long/20:.2f})\n"
        f"d' real (efectivo para cálculo) = {d_real:.2f} cm"
    )

    # cálculo de separación con d_real
    if seccion_tipo == "I":
        s_ejes = b - 2*d_real
        sep_libre = s_ejes - diam_long/10
        texto = f"{n_barras}Ø{diam_long} → {n_barras//2} y {n_barras - n_barras//2} en las dos caras"

    elif seccion_tipo == "II":
        esquinas = 4
        resto = max(n_barras - esquinas, 0)
        texto = f"{n_barras}Ø{diam_long} → {esquinas} en esquinas + {resto} distribuidas en lados"

        # distribución en lados
        barras_por_lado = resto // 2 if resto > 0 else 0

        if barras_por_lado > 0:
            s_ejes = b - 2*d_real
            # separación entre ejes = longitud útil / (barras_por_lado+1)
            sep_ejes = s_ejes / (barras_por_lado+1)
            sep_libre = sep_ejes - diam_long/10
        else:
            # si solo hay esquinas, separación libre = lado útil - Øbarra
            s_ejes = b - 2*d_real
            sep_libre = s_ejes - diam_long/10


    elif seccion_tipo == "III":
        texto = f"{n_barras}Ø{diam_long} → {n_barras} barras equidistantes en perímetro circular"
        circ = math.pi * b
        sep_libre = circ/n_barras - diam_long/10

    # chequeo normativo de separación
    sep_min = max(diam_long/10, 2.5, (agregado+0.5)/10)
    cumple_sep = sep_libre >= sep_min

    texto_sep = f"Separación libre ≈ {sep_libre:.1f} cm (mínimo normativo ≈ {sep_min:.1f} cm → cumple={cumple_sep})"

    return texto_rec + "\n" + texto + "\n" + texto_sep
    
def armado_estribos(Ast_min, lado_cm, b=None, h=None, d_usado=None,
                    seccion_tipo="II", diametros=[12,16,20,25,32], agregado=19):
    opciones = []
    descartes = []
    Ast_min_cm2 = Ast_min

    # mínimo normativo: 4Ø12
    area_barra_min = math.pi*(1.2**2)/4   # Ø12 → 1.2 cm
    Ast_min_norma = 4 * area_barra_min
    texto_norma = f"Ast mínimo normativo ≈ {Ast_min_norma:.2f} cm² → 4Ø12"
    texto_calc = f"Ast cálculo teórico ≈ {Ast_min_cm2:.2f} cm²"

    if d_usado is None:
        d_usado = 3.0  # valor por defecto

    # filtro dinámico de diámetros
    diam_max = 20 if d_usado < 4.0 else 25
    diametros_filtrados = [d for d in diametros if 12 <= d <= diam_max]

    for d in diametros_filtrados:
        diam_long_cm = d/10   # mm → cm
        area_barra = math.pi*(diam_long_cm**2)/4

        if seccion_tipo == "I":
            Ast_por_lado = Ast_min_cm2 / 2
            n_por_lado = max(2, math.ceil(Ast_por_lado / area_barra))
            n_barras = 2 * n_por_lado
            s_ejes = b - 2*d_usado
            if n_por_lado == 2:
                sep_libre = s_ejes - diam_long_cm
            else:
                sep_libre = s_ejes/(n_por_lado-1) - diam_long_cm

        elif seccion_tipo == "II":
            Ast_por_lado = Ast_min_cm2 / 4
            n_por_lado = max(1, math.ceil(Ast_por_lado / area_barra))
            n_barras = 4 * n_por_lado
            s_ejes = b - 2*d_usado
            if n_por_lado == 1:
                sep_libre = s_ejes - diam_long_cm
            else:
                sep_libre = s_ejes/(n_por_lado-1) - diam_long_cm

        elif seccion_tipo == "III":
            n_barras = max(4, math.ceil(Ast_min_cm2 / area_barra))
            circ = math.pi*b
            sep_libre = circ/n_barras - diam_long_cm

        # chequeo normativo
        sep_min = max(diam_long_cm, 2.0, (agregado+0.5)/10)
        if sep_libre < sep_min:
            descartes.append((n_barras, d, sep_libre, sep_min))
            continue

        Ast_real = n_barras * area_barra
        # tolerancia: aceptar ≥90% del Ast requerido
        if Ast_real >= 0.9*Ast_min_cm2:
            opciones.append((n_barras, d, Ast_real))
        else:
            descartes.append((n_barras, d, sep_libre, sep_min))

    # salida si no hay opciones viables
    if not opciones:
        texto_desc = "\n".join(
            [f"{n}Ø{d} → sep_libre={sep:.1f} cm < mínimo={minimo:.1f} cm"
             for (n,d,sep,minimo) in descartes]
        )
        return {
            "texto": "⚠️ No hay combinación viable.\n" + texto_desc,
            "n_barras": 0,
            "diam_long": 0,
            "Ast_real": 0,
            "diam_estribo": 0,
            "sep_max": 0
        }

    # elegir la opción más ajustada y con menos barras
    mejor = min(opciones, key=lambda x: (abs(x[2]-Ast_min_cm2), x[0]))
    n_barras_calc, diam_long, Ast_real_calc = mejor
    texto_prop = f"Armadura propuesta: {n_barras_calc}Ø{diam_long} → Ast real ≈ {Ast_real_calc:.2f} cm²"

    texto_comp = composicion_armadura(
        seccion_tipo,
        n_barras_calc,
        diam_long,
        b=b,
        h=h,
        d_usado=d_usado
    )

    # diámetro mínimo de estribo
    if diam_long <= 16: diam_estribo = 6
    elif diam_long <= 25: diam_estribo = 8
    elif diam_long <= 32: diam_estribo = 10
    else: diam_estribo = 12

    # límites de separación (convertidos a cm)
    limite1 = (12*diam_long)/10
    limite2 = (48*diam_estribo)/10
    limite3 = lado_cm
    s_max = min(limite1, limite2, limite3)

    if s_max == limite1: gobernante = f"12·db = {limite1:.1f} cm"
    elif s_max == limite2: gobernante = f"48·Øestribo = {limite2:.1f} cm"
    else: gobernante = f"lado menor = {limite3:.1f} cm"

    return {
        "texto": (
            f"{texto_calc}\n"
            f"{texto_norma}\n"
            f"{texto_prop}\n"
            f"{texto_comp}\n"
            f"Estribos Ø{diam_estribo} mm (mínimo normativo Ø{diam_estribo}, cumple=True)\n"
            f"Separación máxima ≈ {s_max:.1f} cm (gobierna {gobernante})"
        ),
        "n_barras": n_barras_calc,
        "diam_long": diam_long,
        "Ast_real": Ast_real_calc,
        "diam_estribo": diam_estribo,
        "sep_max": s_max
    }

def armado_sunchos(Ast_min, diam_cm, fck, fy, diametros=[12, 16, 20], diam_espiral=10, recubrimiento=40):
    opciones = []
    Ast_min_cm2 = Ast_min

    # mínimo normativo: 6Ø12
    area_barra_min = math.pi*(12/10)**2/4
    Ast_min_norma = 6 * area_barra_min
    texto_norma = f"Ast mínimo normativo ≈ {Ast_min_norma:.2f} cm² → 6Ø12"
    texto_calc = f"Ast cálculo teórico ≈ {Ast_min_cm2:.2f} cm²"

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

    return {
        "texto": (
            f"{texto_calc}\n"
            f"{texto_norma}\n"
            f"{texto_prop}\n"
            f"Espiral Ø{diam_espiral} mm → paso ≈ {s:.1f} cm "
            f"(cumple rango 2.5–8.0 cm: {cumple})\n"
            f"Asp/s requerido ≈ {Asp_s:.2f} cm²/m"
        ),
        "n_barras": n_barras_calc,
        "diam_long": diam_long,
        "Ast_real": Ast_real_calc,
        "diam_espiral": diam_espiral,
        "paso": s,
        "cumple": cumple
    }
 
 # ==============================

# 4. Cálculo de pesos de armadura
# ==============================
def peso_barra(diam_mm, largo_m):
    # área en cm²
    diam_cm = diam_mm/10
    area_cm2 = math.pi*(diam_cm**2)/4
    # densidad acero ≈ 0.00617 kg/cm²·cm
    return area_cm2 * (largo_m*100) * 0.00617

def peso_estribos(diam_mm, b_cm, h_cm, altura_m, paso_cm):
    diam_cm = diam_mm/10
    area_cm2 = math.pi*(diam_cm**2)/4
    largo_estribo_cm = 2*(b_cm+h_cm) + 20  # +20 cm ganchos aprox
    cantidad = int(altura_m*100/paso_cm)+1
    return area_cm2 * largo_estribo_cm * cantidad * 0.00617
# ==============================
# Función para GUARDAR resultados en CSV
# ==============================

def guardar_resultados_csv(
    portico, id_columna, tipo, tipo_seccion, dimensiones_texto,
    altura_libre_m, d_recubrimiento_cm, Pu, Mu, fck, lambda_val,
    clasificacion, Ast, detalle_armadura, nota_diagrama,
    carpeta="salidas/columnas"
):
    import os, csv

    archivo = os.path.join(carpeta, "planilla_columnas.csv")
    campos = [
        "Pórtico","Columna","Tipo","Tipo_seccion","Dimensiones",
        "Altura_libre(m)","Recubrimiento(cm)",
        "Pu(kN)","Mu(kNm)","f'c(MPa)","λ","Clasificación",
        "Ast(cm²)","n_barras","diam_long(cm)","diam_estribo(cm)",
        "paso_estribo(cm)","cumple_estribo","Nota_diagrama"
    ]

    os.makedirs(carpeta, exist_ok=True)

    # Leer contenido existente
    filas = []
    if os.path.exists(archivo):
        with open(archivo, "r", newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter=";")
            filas = list(reader)

    # Datos de armadura
    n_barras = detalle_armadura.get('n_barras', "")
    diam_long = detalle_armadura.get('diam_long', "")
    diam_estribo = detalle_armadura.get('diam_estribo', detalle_armadura.get('diam_espiral', ""))
    paso_estribo = round(detalle_armadura.get('sep_max', detalle_armadura.get('paso', 0)), 1) if diam_estribo else ""
    cumple_estribo = detalle_armadura.get('cumple', "") if diam_estribo else ""

    nueva_fila = {
        "Pórtico": portico,
        "Columna": id_columna,
        "Tipo": tipo,
        "Tipo_seccion": tipo_seccion,
        "Dimensiones": dimensiones_texto,
        "Altura_libre(m)": f"{altura_libre_m:.2f}",
        "Recubrimiento(cm)": f"{d_recubrimiento_cm:.2f}",
        "Pu(kN)": Pu,
        "Mu(kNm)": Mu,
        "f'c(MPa)": fck,
        "λ": f"{lambda_val:.2f}",
        "Clasificación": clasificacion,
        "Ast(cm²)": f"{Ast:.2f}",
        "n_barras": n_barras,
        "diam_long(cm)": diam_long,
        "diam_estribo(cm)": diam_estribo,
        "paso_estribo(cm)": paso_estribo,
        "cumple_estribo": cumple_estribo,
        "Nota_diagrama": nota_diagrama
    }

    # Reemplazar si ya existe esa columna en ese pórtico
    reemplazado = False
    for i, fila in enumerate(filas):
        if fila["Pórtico"] == str(portico) and fila["Columna"] == id_columna:
            filas[i] = nueva_fila
            reemplazado = True
            break

    if not reemplazado:
        filas.append(nueva_fila)

    # Reescribir todo el archivo
    with open(archivo, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=campos, delimiter=";")
        writer.writeheader()
        writer.writerows(filas)
# ==============================
# Ejemplo de uso completo
# ==============================
# Determinar factor k según altura o nivel
def factor_k(nivel_columna):
    """
    Devuelve k según el nivel de la columna.
    Por ejemplo:
    - Planta baja: 1.0
    - Pisos superiores: 0.7
    """
    if nivel_columna.lower() in ["pb", "planta baja", "0"]:
        return 1.0
    else:
        return 0.7


from pathlib import Path
import json

# =========================================================
# PROGRAMA PRINCIPAL – análisis de COLUMNAS
# =========================================================

BASE = Path(__file__).parent

with open(BASE / "datos" / "estructura.json", encoding="utf-8") as f:
    estructura = json.load(f)

print("Pórticos encontrados:")
porticos = list(estructura.keys())
for i, nombre in enumerate(porticos, start=1):
    print(f"{i}: {nombre}")

seleccion = input("Ingrese el número o nombre del pórtico a calcular: ")

# 👉 permitir tanto número como nombre
if seleccion.isdigit():
    idx = int(seleccion) - 1
    if idx < 0 or idx >= len(porticos):
        raise ValueError("⚠️ Número de pórtico inválido")
    nombre_portico = porticos[idx]
else:
    if seleccion not in estructura:
        raise ValueError("⚠️ Pórtico no encontrado en el archivo JSON")
    nombre_portico = seleccion

memoria_txt = []
memoria_txt.append(f"\n=== MEMORIA DE CÁLCULO ===\n")
memoria_txt.append(f"PÓRTICO {nombre_portico}\n")


print(f"Procesando pórtico: {nombre_portico}")
columnas = estructura[nombre_portico].get("columnas", {})

if not columnas:
    print("⚠️  No hay columnas en este pórtico")
else:
    print(f"\n=== CÁLCULO DE COLUMNAS PARA PÓRTICO: {nombre_portico} ===\n")

    # -----------------------------
    # Loop sobre columnas
    # -----------------------------
    for col_id, col in columnas.items():
        print(f"\nProcesando columna {col_id} del pórtico {nombre_portico}")

        # 👉 preguntar tipo de columna para ESTA columna
        tipo_input = input(
            f"Tipo de columna {col_id} (1=estribos/rectangulares, 2=sunchos/circulares): "
        ).lower()
        TIPOS = {"1": "estribos", "2": "sunchos", "estribos": "estribos", "sunchos": "sunchos"}
        tipo = TIPOS.get(tipo_input, "estribos")

        # -----------------------------
        # Datos básicos de la columna
        # -----------------------------
        Pu = abs(col["P_kN"])            # carga mayorada
        altura = col["altura_m"]
        fck = definir_fck(col_id, Pu)    # criterio de resistencia
        fy = 420                         # MPa fijo
        Mu_inf = col.get("Mu_kNm_inf", 0)
        Mu_sup = col.get("Mu_kNm_sup", 0)
        Mu = max(abs(Mu_inf), abs(Mu_sup))

        # -----------------------------
        # Predimensionado inicial
        # -----------------------------
        b, h, e_cm = predimensionar_seccion(Pu, Mu, fck)

        # -----------------------------
        # Impresión de resultados
        # -----------------------------
        print(f"{col_id:6s} | Pu={Pu:7.0f} kN | h={altura:.2f} m | fck={fck:2d} | fy={fy} "
            f"| sec=({b}, {h}) | tipo={tipo} | Mu={Mu:.2f} kNm")

        if tipo == "estribos":
            print(f"Sección sugerida: b={b} cm, h≈{h:.1f} cm (e={e_cm:.1f} cm)")
        else:  # sunchos
            print(f"Sección sugerida: Ø{max(b,h)} cm (e={e_cm:.1f} cm)")

        # -----------------------------
        # 🚩 FLEXO-COMPRESIÓN
        # -----------------------------
        if tipo == "sunchos":
            seccion_tipo = "III"  # fijo, circular
            b = h = max(b, h)     # diámetro = lado mayor del predimensionado
            dprima_cm = 4.0       # fijo en circulares
        else:
            seccion_tipo_input = input(
                f"Tipo de sección para {col_id} (1=rectangular armadura a 2 lados, 2=rectangular armadura a 4 lados): "
            )
            TIPOS_SECCION = {"1": "I", "2": "II"}
            seccion_tipo = TIPOS_SECCION.get(seccion_tipo_input, "I")

            b = 20.0  # fijo reglamentario
            dprima_cm = 4.0 if e_cm >= 20 else 3.6

        h_sugerido = h
        b_sugerido = b

        gamma = calcular_gamma(h_sugerido, dprima_cm)

        while gamma < 0.70:
            print(f"⚠️ γ={gamma:.2f} < 0.70 → redimensionando sección...")
            h_sugerido += 5  # aumentar lado mayor o diámetro en 5 cm
            if seccion_tipo == "III":
                b_sugerido = h_sugerido  # circular
            gamma = calcular_gamma(h_sugerido, dprima_cm)

        print(f"👉 Sección final predimensionada: "
            f"{'Ø' if seccion_tipo=='III' else 'b=20, h='}{h_sugerido:.1f} cm "
            f"(γ={gamma:.2f})")
        h = h_sugerido
        b = b_sugerido
        Ag = b * h

        # altura libre para esbeltez (ya cargada en col)
        Lc = col["altura_m"]
        nivel = col.get("nivel", "alto")
        k = factor_k(nivel)
        lc_cm = k * Lc * 100


        # momento de inercia
        if seccion_tipo == "III":
            I = math.pi * (b**4) / 64
        else:
            I = b * h**3 / 12
        r = math.sqrt(I / Ag)
        lambda_val = lc_cm / r
        limite = 34

        # chequeo esbeltez y amplificación
        if lambda_val > limite:
            E = 4700 * math.sqrt(fck) * 1000
            Pcr = (math.pi**2 * E * I) / (lc_cm**2)
            Pu_N = Pu * 1000
            delta = 1 / (1 - Pu_N / Pcr)
            Mu_amp = Mu * delta
            print(f"⚠️ Columna esbelta → δ = {delta:.2f}, Mu amplificado = {Mu_amp:.2f} kNm")
            x_obj = calcular_x_normalizado(Pu, Mu_amp, b, h)
        else:
            x_obj = calcular_x_normalizado(Pu, Mu, b, h)

        y_obj = calcular_y_normalizado(Pu, b, h)
        diagrama = cargar_diagrama_interpolado(seccion_tipo, fck, gamma, x_obj, y_obj)
        rho_g = diagrama["rho_calc"]
        Ast_adoptado = rho_g * Ag

        if tipo == "estribos":
            detalle_armadura = armado_estribos(Ast_adoptado, min(b, h), b=b, h=h, d_usado=dprima_cm, seccion_tipo=seccion_tipo)
        else:
            detalle_armadura = armado_sunchos(Ast_adoptado, b, fck, fy)

        # -----------------------------
        # Salida y exportación
        # -----------------------------
        print("\n=== MEMORIA DE CÁLCULO DE COLUMNA ===")
        print(f"Tipo: {tipo}")
        print(f"Pu = {Pu} kN, Mu = {Mu} kNm")
        print(f"f'c = {fck} MPa, fy = {fy} MPa")
        print(f"Dimensiones: b={b} cm, h={h} cm")
        print(f"d' = {dprima_cm:.2f} cm → γ = {gamma:.2f}")
        print(f"Cuantía geométrica ρg = {rho_g:.3f}")
        print(f"Ast adoptado = {Ast_adoptado:.2f} cm²")
        print(detalle_armadura)
        print(f"Diagrama usado: {diagrama['title']}")
        print(f"Esbeltez λ = {lambda_val:.2f} → {'CORTA' if lambda_val <= limite else 'ESBELTA'}")

        nota_diagrama = f"γ={gamma:.2f}, X={x_obj:.3f}, Y={y_obj:.3f}, ρg={rho_g:.3f}"
        # Definir dimensiones_texto
        if seccion_tipo == "III":
            dimensiones_texto = f"Ø{b:.1f}"
        else:
            dimensiones_texto = f"b={b:.1f}, h={h:.1f}"

        # Definir clasificación
        clasificacion = "ESBELTA" if lambda_val > 35 else "CORTA"

        # Llamada corregida
        guardar_resultados_csv(
            portico=nombre_portico,
            id_columna=col_id,
            tipo=tipo,
            tipo_seccion=seccion_tipo,
            dimensiones_texto=dimensiones_texto,
            altura_libre_m=Lc,
            d_recubrimiento_cm=dprima_cm,
            Pu=Pu,
            Mu=Mu,
            fck=fck,
            lambda_val=lambda_val,
            clasificacion=clasificacion,
            Ast=Ast_adoptado,
            detalle_armadura=detalle_armadura,
            nota_diagrama=nota_diagrama
        )
        print("\n✅ Resultados guardados en salida/planilla_columnas.csv")

        # -----------------------------
        # Acumular en TXT
        # -----------------------------
        empalme_m = (48 * detalle_armadura['diam_long']) / 1000  # en metros
        largo_total = Lc + empalme_m
        if seccion_tipo == "I":
            descripcion_diagrama = (
                f"DIAGRAMA I - rectangulares con barras en las caras extremas "
                f"(f'c={fck} MPa, fy={fy} MPa, gamma={gamma:.2f})"
            )
        elif seccion_tipo == "II":
            descripcion_diagrama = (
                f"DIAGRAMA II - rectangulares con barras en las cuatro caras "
                f"(f'c={fck} MPa, fy={fy} MPa, gamma={gamma:.2f})"
            )
        else:  # III
            descripcion_diagrama = (
                f"DIAGRAMA III - circulares "
                f"(f'c={fck} MPa, fy={fy} MPa, gamma={gamma:.2f})"
            )
        
        texto_col = f"""
                Columna: {col_id}
                Sección: {dimensiones_texto}
                Altura libre: {Lc:.2f} m
                Recubrimiento: {dprima_cm:.2f} cm
                Pu = {Pu:.2f} kN, Mu = {Mu:.2f} kNm
                f'c = {fck} MPa
                Esbeltez λ = {lambda_val:.2f} → {clasificacion}
                Ast adoptado = {Ast_adoptado:.2f} cm²
                Armadura longitudinal: {detalle_armadura['n_barras']}Ø{detalle_armadura['diam_long']}      
                Descripción: {descripcion_diagrama}
                Nota diagrama: {nota_diagrama}
        """

        # Diferenciar según tipo
        if tipo == "estribos":
            texto_col += f"        Estribos: Ø{detalle_armadura['diam_estribo']} mm @ {detalle_armadura['sep_max']:.1f} cm\n"

            # Despiece completo (barras + estribos)
            texto_col += f"""
                Despiece:
                - Barras longitudinales: {detalle_armadura['n_barras']}Ø{detalle_armadura['diam_long']}
                Largo individual ≈ {Lc:.2f} m + empalme (48·db ≈ {empalme_m:.2f} m) → {largo_total:.2f} m
                Peso por barra ≈ {peso_barra(detalle_armadura['diam_long'], largo_total):.2f} kg
                Peso total ≈ {peso_barra(detalle_armadura['diam_long'], largo_total)*detalle_armadura['n_barras']:.2f} kg
                - Estribos: Ø{detalle_armadura['diam_estribo']} mm
                Largo individual ≈ {2*(b+h)/100 + 0.2:.2f} m (incluye ganchos)
                Cantidad ≈ {int(Lc*100/detalle_armadura['sep_max'])+1}
                Peso total ≈ {peso_estribos(detalle_armadura['diam_estribo'], b, h, Lc, detalle_armadura['sep_max']):.2f} kg
                -------------------------------------------------
            """

        elif tipo == "sunchos":
            texto_col += f"        Espiral Ø{detalle_armadura['diam_espiral']} mm → paso ≈ {detalle_armadura['paso']:.1f} cm (cumple: {detalle_armadura['cumple']})\n"
            for linea in detalle_armadura['texto'].splitlines():
                texto_col += f"                {linea}\n"

            # Calcular largo total del espiral
            Dc_m = (b/100) - 2*(dprima_cm/100)  # diámetro confinado en m
            vueltas = Lc / (detalle_armadura['paso']/100)  # cantidad de vueltas
            largo_espiral = vueltas * (math.pi * Dc_m)     # largo total en m
            peso_espiral_total = peso_barra(detalle_armadura['diam_espiral'], largo_espiral)

            texto_col += f"""
                Despiece:
                - Barras longitudinales: {detalle_armadura['n_barras']}Ø{detalle_armadura['diam_long']}
                Largo individual ≈ {Lc:.2f} m + empalme (48·db ≈ {empalme_m:.2f} m) → {largo_total:.2f} m
                Peso por barra ≈ {peso_barra(detalle_armadura['diam_long'], largo_total):.2f} kg
                Peso total ≈ {peso_barra(detalle_armadura['diam_long'], largo_total)*detalle_armadura['n_barras']:.2f} kg
                - Espiral: Ø{detalle_armadura['diam_espiral']} mm
                Largo total ≈ {largo_espiral:.2f} m
                Peso total ≈ {peso_espiral_total:.2f} kg
                -------------------------------------------------
            """


        # Volumen de hormigón
        if seccion_tipo in ["I","II"]:  # rectangulares
            volumen_hormigon = (b/100) * (h/100) * Lc  # b,h en cm → pasamos a m
        else:  # circulares
            diam_m = b/100  # diámetro en m
            volumen_hormigon = math.pi * (diam_m/2)**2 * Lc

        texto_col += f"    Volumen de hormigón ≈ {volumen_hormigon:.2f} m³\n"
        # ==============================
        # Verificación normativa φPn
        # ==============================
        # Área bruta Ag (m²)
        if seccion_tipo in ["I","II"]:
            Ag = (b/100) * (h/100)
        else:
            diam_m = b/100
            Ag = math.pi * (diam_m/2)**2

        # Área de acero longitudinal Ast (cm² → m²)
        Ast_m2 = (detalle_armadura['n_barras'] *
                (math.pi * (detalle_armadura['diam_long']/10)**2 / 4)) / 10000

        phi = PHI[tipo]
        alpha = ALPHA[tipo]

        Pn = alpha * fck * 1e6 * (Ag - Ast_m2) + fy * 1e6 * Ast_m2  # N
        Pn_red = phi * Pn / 1000  # kN

        texto_col += f"        Resistencia nominal reducida φPn ≈ {Pn_red:.2f} kN\n"
        # Comparación con la carga solicitada Pu
        if Pn_red < Pu:
            texto_col += f"        Verificación: NO CUMPLE (φPn = {Pn_red:.2f} kN < Pu = {Pu:.2f} kN)\n"
        else:
            texto_col += f"        Verificación: CUMPLE (φPn = {Pn_red:.2f} kN ≥ Pu = {Pu:.2f} kN)\n"
        # Excentricidad de cálculo
        excentricidad = Mu / Pu  # en m

        # Momento resistente nominal
        Mn = Pn * excentricidad  # N·m
        Mn_red = phi * Mn / 1000  # kNm

        if Mn_red < Mu:
            texto_col += f"        Verificación a flexo-compresión: NO CUMPLE (φMn = {Mn_red:.2f} kNm < Mu = {Mu:.2f} kNm)\n"
        else:
            texto_col += f"        Verificación a flexo-compresión: CUMPLE (φMn = {Mn_red:.2f} kNm ≥ Mu = {Mu:.2f} kNm)\n"

        
        memoria_txt.append(texto_col)

        # Guardar memoria en TXT
        SALIDA = BASE / "salidas" / "columnas"
        SALIDA.mkdir(parents=True, exist_ok=True)

        with open(SALIDA / f"memoria_{nombre_portico}.txt", "w", encoding="utf-8") as f:
            f.writelines(memoria_txt)

        print(f"✅ Memoria y despiece guardados en {SALIDA}/memoria_{nombre_portico}.txt")
import math
import json
import os
import csv, os

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
    gamma_normativos = [0.70, 0.80, 0.90]

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
    descartes = []   # aquí guardamos las combinaciones descartadas
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
        diam_long_cm = d/10   # pasar mm → cm
        area_barra = math.pi*(diam_long_cm**2)/4
        n_barras = max(4, math.ceil(Ast_min_cm2 / area_barra))
        Ast_real = n_barras * area_barra

        # separación libre en cm
        if seccion_tipo == "I":
            s_ejes = b - 2*d_usado
            sep_libre = s_ejes/(n_barras//2) - diam_long_cm if n_barras > 2 else s_ejes - diam_long_cm

        elif seccion_tipo == "II":
            esquinas = 4
            resto = max(n_barras - esquinas, 0)
            barras_por_lado = resto//2 if resto > 0 else 0
            s_ejes = b - 2*d_usado
            if barras_por_lado > 0:
                sep_ejes = s_ejes/(barras_por_lado+1)
                sep_libre = sep_ejes - diam_long_cm
            else:
                sep_libre = s_ejes - diam_long_cm

        elif seccion_tipo == "III":
            circ = math.pi*b
            sep_libre = circ/n_barras - diam_long_cm

        # chequeo normativo
        sep_min = max(diam_long_cm, 2.5, (agregado+0.5)/10)
        if sep_libre < sep_min:
            descartes.append((n_barras, d, sep_libre, sep_min))
            continue

        opciones.append((n_barras, d, Ast_real))

    # salida si no hay opciones viables
    if not opciones:
        texto_desc = "\n".join(
            [f"{n}Ø{d} → sep_libre={sep:.1f} cm < mínimo={minimo:.1f} cm"
             for (n,d,sep,minimo) in descartes]
        )
        return "⚠️ No hay combinación viable.\n" + texto_desc

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

# ==============================
# Función para GUARDAR resultados en CSV
# ==============================

def guardar_resultados_csv(
    id_columna, Pu, Mu, fck, b, h, lambda_val,
    clasificacion, Ast, detalle_armadura, nota_diagrama,
    carpeta="salidas"
):
    archivo = os.path.join(carpeta, "planilla_columnas.csv")
    campos = [
        "ID", "Pu(kN)", "Mu(kNm)", "f'c(MPa)",
        "b(cm)", "h(cm)", "λ", "Clasificación",
        "Ast(cm²)", "Armadura", "Estribos/Sunchos", "Nota_diagrama"
    ]

    os.makedirs(carpeta, exist_ok=True)

    # Crear encabezados si el archivo no existe
    if not os.path.exists(archivo):
        with open(archivo, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f, delimiter=";")
            writer.writerow(campos)

    # Construir campos compactos
    armadura = f"{detalle_armadura.get('n_barras','')}Ø{detalle_armadura.get('diam_long','')}"

    if "diam_estribo" in detalle_armadura:
        refuerzo = f"Ø{detalle_armadura.get('diam_estribo','')} c/{detalle_armadura.get('sep_max','')}cm"
    elif "diam_espiral" in detalle_armadura:
        cumple = detalle_armadura.get("cumple", "")
        refuerzo = f"Ø{detalle_armadura.get('diam_espiral','')} c/{detalle_armadura.get('paso','')}cm (cumple={cumple})"
    else:
        refuerzo = ""

    # Agregar fila
    with open(archivo, "a", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter=";")
        writer.writerow([
            id_columna, Pu, Mu, fck, b, h,
            f"{lambda_val:.2f}", clasificacion,
            f"{Ast:.2f}", armadura, refuerzo, nota_diagrama
        ])
# ==============================
# Ejemplo de uso completo
# ==============================

if __name__ == "__main__":
    id_columna = input("Ingrese ID de la columna: ")
    Pu = float(input("Ingrese Pu (kN): "))
    Mu = float(input("Ingrese Mu (kNm): "))
    fck = int(input("Ingrese f'c (20,25,30 MPa): "))
    fy = 420  # MPa
    tipo = input("Tipo de columna (1=estribos / 2=sunchos): ").lower()

    TIPOS = {"1": "estribos", "2": "sunchos", "estribos": "estribos", "sunchos": "sunchos"}
    tipo = TIPOS.get(tipo, "estribos")

    # 🚩 Bifurcación inmediata según momento
    if Mu == 0:
        # === AXIAL PURO ===
        Ag, Ast_min_calc, Ast_min_adopt, Pn = calcular_seccion(Pu, fck, fy=fy, tipo=tipo)
        lado_calc, lado_adopt, descripcion_geom = dimensiones(Ag, tipo)

        # 🚩 Pedir altura libre para esbeltez
        Lc = float(input("Ingrese altura libre de la columna (m): "))
        k = 1.0  # caso más común en tu obra
        lc_cm = k * Lc * 100  # pasar a cm

        # Radio de giro (rectangular o circular)
        if tipo == "sunchos":  # circular
            I = math.pi * (lado_adopt**4) / 64
        else:  # rectangular
            I = lado_adopt * (lado_adopt**3) / 12
        r = math.sqrt(I / Ag)
        lambda_val = lc_cm / r
        limite = 34  # valor típico ACI para columnas indesplazables

        # 🚩 Chequeo de esbeltez y amplificación de momentos
        if lambda_val > limite:
            E = 4700 * math.sqrt(fck) * 1000  # N/cm²
            Pcr = (math.pi**2 * E * I) / (lc_cm**2)
            Pu_N = Pu * 1000
            if Pu_N >= Pcr:
                print("⚠️ Columna esbelta → riesgo de pandeo, Pu ≥ Pcr → NO CUMPLE")
            else:
                print(f"⚠️ Columna esbelta → considerar reducción de resistencia por pandeo (Pcr = {Pcr/1000:.2f} kN)")


        print("\n=== MEMORIA DE CÁLCULO DE COLUMNA (AXIAL PURO) ===")
        print(f"Tipo: {tipo}")
        print(f"Pu = {Pu} kN")
        print(f"f'c = {fck} MPa, fy = {fy} MPa")
        print(f"Ag = {Ag:.2f} cm², As calc = {Ast_min_calc:.2f} cm², As adoptado = {Ast_min_adopt:.2f} cm²")
        print(descripcion_geom)
        print(f"Resistencia axial reducida φPn = {Pn:.2f} kN")
        print(f"Resultado: {'CUMPLE' if Pu <= Pn else 'NO CUMPLE'}")
        print(f"Esbeltez λ = {lambda_val:.2f} → {'CORTA' if lambda_val <= limite else 'ESBELTA'}")

        if tipo == "estribos":
            detalle_armadura = armado_estribos(
                Ast_min_adopt,
                lado_cm=lado_adopt,
                b=lado_adopt,
                h=lado_adopt,
                d_usado=3.8,
                seccion_tipo="II"
            )
            print(detalle_armadura)
        else:
            detalle_armadura = armado_sunchos(
                Ast_min_adopt,
                diam_cm=lado_adopt,
                fck=fck,
                fy=fy,
                diametros=[12, 16, 20],
                diam_espiral=10,
                recubrimiento=40
            )
            print(detalle_armadura)
        # Nota al pie vacía porque es axial puro
        nota_diagrama = ""

        guardar_resultados_csv(
            id_columna=id_columna,
            Pu=Pu,
            Mu=Mu,
            fck=fck,
            b=lado_adopt,
            h=lado_adopt,
            lambda_val=lambda_val,
            clasificacion="CORTA" if lambda_val <= limite else "ESBELTA",
            Ast=Ast_min_adopt,
            detalle_armadura=detalle_armadura,
            nota_diagrama=nota_diagrama
        )
        print("\n✅ Resultados guardados en salida/planilla_columnas.csv")

    else:
        # === FLEXO-COMPRESIÓN ===
        dprima_cm = float(input("Ingrese d' (cm): "))
        h = float(input("Altura h (cm): "))
        b = float(input("Base b (cm, si circular poner igual a h): "))
        seccion_tipo = input("Tipo de sección (I=rectangular 2 lados, II=rectangular 4 lados, III=circular): ").upper()

        gamma = calcular_gamma(h, dprima_cm)
        Ag = b * h  # área bruta en cm²

        # 🚩 Pedir altura libre para esbeltez
        Lc = float(input("Ingrese altura libre de la columna (m): "))
        k = 1.0
        lc_cm = k * Lc * 100

        if seccion_tipo == "III":  # circular
            I = math.pi * (b**4) / 64
        else:  # rectangular
            I = b * h**3 / 12
        r = math.sqrt(I / Ag)
        lambda_val = lc_cm / r
        limite = 34

        # 🚩 Chequeo de esbeltez y amplificación de momentos
        if lambda_val > limite:
            # módulo de elasticidad aproximado del hormigón (puedes ajustar según fck)
            E = 4700 * math.sqrt(fck)  # MPa
            E = E * 1000  # pasar a N/cm²

            # momento de inercia ya calculado (I)
            Pcr = (math.pi**2 * E * I) / (lc_cm**2)  # carga crítica de pandeo
            Pu_N = Pu * 1000  # pasar Pu a N
            delta = 1 / (1 - Pu_N / Pcr)
            Mu_amp = Mu * delta

            print(f"⚠️ Columna esbelta → δ = {delta:.2f}, Mu amplificado = {Mu_amp:.2f} kNm")

            # volver a calcular coordenadas con Mu amplificado
            x_obj = calcular_x_normalizado(Pu, Mu_amp, b, h)
            y_obj = calcular_y_normalizado(Pu, b, h)
        else:
            # usar Mu original
            x_obj = calcular_x_normalizado(Pu, Mu, b, h)
            y_obj = calcular_y_normalizado(Pu, b, h)


        # Coordenadas normalizadas del punto de solicitación
        x_obj = calcular_x_normalizado(Pu, Mu, b, h)
        y_obj = calcular_y_normalizado(Pu, b, h)

        diagrama = cargar_diagrama_interpolado(seccion_tipo, fck, gamma, x_obj, y_obj)
        rho_g = diagrama["rho_calc"]
        Ast_adoptado = rho_g * Ag

        if tipo == "estribos":
            detalle_armadura = armado_estribos(
                Ast_adoptado,
                min(b, h),
                b=b,
                h=h,
                d_usado=dprima_cm,
                seccion_tipo=seccion_tipo
            )
        else:
            detalle_armadura = armado_sunchos(Ast_adoptado, b, fck, fy)

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
        # Nota al pie con valores normativos
        nota_diagrama = f"γ={gamma:.2f}, X={x_obj:.3f}, Y={y_obj:.3f}, ρg={rho_g:.3f}"

        guardar_resultados_csv(
            id_columna=id_columna,
            Pu=Pu,
            Mu=Mu,
            fck=fck,
            b=b,
            h=h,
            lambda_val=lambda_val,
            clasificacion="CORTA" if lambda_val <= limite else "ESBELTA",
            Ast=Ast_adoptado,
            detalle_armadura=detalle_armadura,
            nota_diagrama=nota_diagrama
        )
        print("\n✅ Resultados guardados en salida/planilla_columnas.csv")
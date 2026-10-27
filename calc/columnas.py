import json
import os
import math
from pathlib import Path

from . import rutas
from . import columnas_biaxial

# ============================================================
# Bloque 1 (portado): predimensión y parámetros normativos
# Extraídos de Legacy/P04_Columnas_portico.py
# ============================================================

PHI = {
    "estribos": 0.65,   # columnas con estribos rectos
    "sunchos": 0.70     # columnas zunchadas (espiral)
}
ALPHA = {
    "estribos": 0.80,
    "sunchos": 0.85
}


def definir_fck(nombre_columna, Pu_kN):
    """Define el fck de la columna. C0 = planta baja, C1/C2 = pisos superiores."""
    if str(nombre_columna).startswith("C0"):
        if Pu_kN > 600:
            return 25
        return 20
    return 20


def predimensionar_seccion(Pu_kN, Mu_kNm, fck_MPa):
    """Predimensionado controlado para columnas rectangulares."""
    Pu = abs(Pu_kN) * 1000  # N
    Mu = abs(Mu_kNm) * 1e6  # Nmm

    if Pu > 0:
        e_cm = Mu / Pu / 10.0  # cm
    else:
        e_cm = 0

    if abs(Mu_kNm) > 200:
        b = 30.0
    elif abs(Mu_kNm) > 120:
        b = 25.0
    else:
        b = 20.0

    h = b + 0.15 * e_cm
    h = max(30.0, min(h, 55.0))
    h = math.ceil(h / 5.0) * 5.0

    return b, h, e_cm


def predimensionar_suncho(Pu_kN, Mu_kNm, fck_MPa, diam_min=20):
    """Predimensionado para columnas circulares (zunchadas)."""
    Pu = abs(Pu_kN) * 1000  # N
    Mu = abs(Mu_kNm) * 1e6  # Nmm

    Ag = Pu / (0.35 * fck_MPa * 10)  # cm²
    d = math.sqrt(4 * Ag / math.pi)

    e_cm = Mu / Pu / 10 if Pu > 0 else 0

    if e_cm > d:
        d = d + 0.1 * e_cm

    moldes = [20, 25, 30, 35, 40]
    d_final = max(diam_min, min(moldes, key=lambda m: abs(m - d)))

    return d_final, e_cm


def cargar_diagrama(nombre):
    try:
        with open(nombre, encoding="utf-8") as f:
            return json.load(f)
    except UnicodeDecodeError:
        with open(nombre, encoding="latin-1") as f:
            return json.load(f)


def calcular(datos, portico="Portico 4"):
    entrada = os.path.join("salidas", "solicitaciones", f"{portico}.json")
    salida = os.path.join("salidas", "columnas", f"resultados_{portico}_columnas.json")

    with open(entrada, "r") as f:
        solicitaciones = json.load(f)

    resultados = []

    for cid, valores in solicitaciones.get("envolvente", {}).get("columnas", {}).items():
        Mu = max(valores.get("M_inf", 0), valores.get("M_sup", 0))
        N = valores.get("N", 0)
        chequeo = {
            "id": cid,
            "Mu": Mu,
            "N": N,
            "estado": "OK" if Mu < 200 and N < 400 else "Revisar",
        }
        resultados.append(chequeo)

    os.makedirs(os.path.dirname(salida), exist_ok=True)
    with open(salida, "w") as f:
        json.dump(resultados, f, indent=4)

    return resultados


def calcular_global() -> list[dict]:
    acumulado = rutas.acumular_solicitaciones_columnas() or {}
    estructura = rutas.cargar_estructura()
    resultados = []
    for gid in sorted(acumulado):
        archivo = Path(acumulado[gid])
        datos = rutas.leer_json(archivo, {}) or {}
        env = datos.get("envolvente_global", {}) or {}
        catalogo = datos.get("catalogo", {}) or {}
        n_val = float(env.get("N", 0.0) or 0.0)
        deriva = float(env.get("deriva", 0.0) or 0.0)
        porticos = None
        for clave in ("pórticos", "porticos", "pÃ³rticos"):
            val = datos.get(clave)
            if isinstance(val, dict):
                porticos = val
                break
        if porticos is None:
            porticos = {}
        mx = 0.0
        my = 0.0
        for nombre_portico, pdata in (porticos or {}).items():
            env_p = (pdata or {}).get("envolvente_pórtico", {}) or {}
            m_inf = float(env_p.get("M_inf", 0.0) or 0.0)
            m_sup = float(env_p.get("M_sup", 0.0) or 0.0)
            m = max(abs(m_inf), abs(m_sup))
            try:
                p = rutas.pórtico(estructura, nombre_portico) or {}
                direccion = (p.get("referencia_planta", {}) or {}).get("direccion")
            except Exception:
                direccion = None
            if direccion == "x":
                my = max(my, m)
            else:
                mx = max(mx, m)
        Mu = max(mx, my)
        porticos_nombres = list(porticos.keys())
        # Bloque 2: predimensión de sección (unidireccional por Mu)
        fck = definir_fck(gid, n_val)
        b, h, e_cm = predimensionar_seccion(n_val, Mu, fck)
        biaxial = columnas_biaxial.combinar_momentos_globales(mx, my)
        resultados.append({
            "id": gid,
            "eje_id": catalogo.get("eje_id"),
            "eje_nombre": catalogo.get("eje_nombre"),
            "x": catalogo.get("x"),
            "nivel": catalogo.get("nivel"),
            "Mx": mx,
            "My": my,
            "Mu": Mu,
            "N": n_val,
            "deriva": deriva,
            "porticos": porticos_nombres,
            "n_porticos": len(porticos_nombres),
            "biaxial": biaxial,
            "fck_MPa": fck,
            "b_cm": b,
            "h_cm": h,
            "excentricidad_cm": e_cm,
            "estado": "OK" if Mu < 200 and n_val < 400 else "Revisar",
        })

    salida = rutas.SAL_COLUMNAS / "resultados_columnas_globales.json"
    rutas.SAL_COLUMNAS.mkdir(parents=True, exist_ok=True)
    rutas.guardar_json(salida, resultados)
    return resultados


if __name__ == "__main__":
    print(calcular({}))

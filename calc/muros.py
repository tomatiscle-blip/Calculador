"""calc.muros — Cálculo de muros ESTRUCTURALES (los de la estructura mixta).

No confundir con los muros del catálogo de Cargas: esos son elementos de carga
(mortero/ladrillo como carga lineal sobre un pórtico) y se calculan en
calc/cargas.py. Acá se calcula el muro portante que se crea en Inicio sobre un
pórtico (estructura mixta).

Estado del cálculo:
    * Peso propio por capas: SÍ (suma γ·e de cada capa con datos de materiales).
    * Peso total y carga lineal resultante: SÍ.
    * Resistencia térmica: SÍ cuando todas las capas tienen conductividad λ
      (``materiales.json``, campo ``lambda_w_mk``) o son cámaras de aire (su
      resistencia se tabula por espesor según ISO 6946). Si falta alguna λ, el
      peso propio SIGUE calculando y la térmica queda en ``None`` con aviso.
    * Diseño estructural y cimiento: PENDIENTES (requieren reglamentación de
      mampostería). No se inventan: quedan como claves con ``None`` y un aviso.
"""

from __future__ import annotations

from . import materiales

# ---------------------------------------------------------------------------
# Resistencias superficiales (ISO 6946) para un muro vertical, flujo horizontal
# del calor. R_si = cara interior (ambiente), R_se = cara exterior (intemperie).
# ---------------------------------------------------------------------------
R_SI = 0.13   # m²·K/W
R_SE = 0.04   # m²·K/W

# Resistencia térmica de cámaras de aire NO ventiladas (ISO 6946, flujo
# horizontal), en función del espesor. El aire no tiene λ propia: su resistencia
# no crece lineal con el espesor, por eso se tabula y se interpola.
# Pares (espesor_m, R en m²·K/W).
CAMARAS_AIRE = (
    (0.005, 0.11),
    (0.007, 0.13),
    (0.010, 0.15),
    (0.015, 0.17),
    (0.025, 0.18),
    (0.050, 0.18),
    (0.100, 0.18),
    (0.300, 0.18),
)


def _gamma_de(clave: str, grupo: str) -> float:
    """Peso específico (kN/m³) de una clave dentro de un grupo de materiales.

    ``grupo`` usa la misma clave que ``materiales.pesos_especificos()``:
    «Mamposteria», «Morteros» o «Aislantes».
    """
    pesos = materiales.pesos_especificos()
    grupo_tabla = pesos.get(grupo, {})
    item = grupo_tabla.get(clave)
    if item is None:
        raise ValueError(
            f"No existe la clave «{clave}» en el grupo «{grupo}» de materiales."
        )
    return float(item.get("gamma", item.get("q", 0.0)))


def _lambda_de(clave: str, grupo: str) -> float:
    """Conductividad térmica λ (W/m·K) de una clave dentro de un grupo.

    Usa ``materiales.conductividades()`` (lee ``lambda_w_mk`` de
    ``materiales.json``). Fala claro si el material no tiene λ cargada.
    """
    tabla = materiales.conductividades()
    grupo_tabla = tabla.get(grupo, {})
    item = grupo_tabla.get(clave)
    if item is None:
        raise ValueError(
            f"No hay conductividad térmica λ para «{clave}» en el grupo «{grupo}»."
        )
    return float(item["lambda_w_mk"])


def _es_camara_aire(clave: str) -> bool:
    return str(clave) == "camara_aire"


def _resistencia_camara_aire(espesor_m: float) -> float:
    """Resistencia térmica (m²·K/W) de una cámara de aire según su espesor.

    Interpola linealmente la tabla ISO 6946 ``CAMARAS_AIRE``. Por debajo del
    mínimo y por encima del máximo se satura en los valores extremos.
    """
    if espesor_m <= CAMARAS_AIRE[0][0]:
        return CAMARAS_AIRE[0][1]
    for (e0, r0), (e1, r1) in zip(CAMARAS_AIRE, CAMARAS_AIRE[1:]):
        if espesor_m <= e1:
            if e1 == e0:
                return r1
            t = (espesor_m - e0) / (e1 - e0)
            return r0 + t * (r1 - r0)
    return CAMARAS_AIRE[-1][1]


def grupos_disponibles() -> dict[str, dict]:
    """Grupos de materiales útiles para armar capas de un muro (por clave)."""
    pesos = materiales.pesos_especificos()
    return {
        "Mamposteria": pesos.get("Mamposteria", {}),
        "Morteros": pesos.get("Morteros", {}),
        "Aislantes": pesos.get("Aislantes", {}),
    }


def peso_propio_capas(capas: list[dict]) -> dict:
    """Peso propio de un muro compuesto a partir de sus capas.

    Cada capa es un dict con ``clave``, ``grupo`` y ``espesor_m`` (m).
    Devuelve el detalle por capa, el espesor total y el peso propio en kN/m²
    (suma de γ·e) y en kN/m (por metro lineal de muro y por metro de altura no
    aplica: es kN/m² de superficie vertical).
    """
    detalle: list[dict] = []
    espesor_total = 0.0
    peso_kN_m2 = 0.0
    for capa in capas:
        clave = str(capa.get("clave", ""))
        grupo = str(capa.get("grupo", ""))
        espesor = float(capa.get("espesor_m", 0.0))
        if espesor <= 0.0:
            continue
        gamma = _gamma_de(clave, grupo)
        aporte = gamma * espesor  # kN/m²
        espesor_total += espesor
        peso_kN_m2 += aporte
        detalle.append({
            "clave": clave,
            "grupo": grupo,
            "espesor_m": espesor,
            "gamma_kN_m3": gamma,
            "peso_kN_m2": aporte,
        })
    return {
        "capas": detalle,
        "espesor_total_m": espesor_total,
        "peso_propio_kN_m2": peso_kN_m2,
    }


def resistencia_termica_capas(capas: list[dict]) -> dict:
    """Resistencia térmica de un muro compuesto por capas (ISO 6946).

    Cada capa aporta R = e/λ (m²·K/W), salvo las cámaras de aire, cuya
    resistencia se tabula por espesor (no tienen λ). Se suman las resistencias
    superficiales interior (R_si) y exterior (R_se) de un muro vertical con
    flujo horizontal del calor. Devuelve:

        * ``capas``: detalle por capa (espesor, λ, R de cada una).
        * ``R_capas_m2K_W``: suma de las resistencias de las capas.
        * ``R_total_m2K_W``: R_si + ΣR_capas + R_se.
        * ``U_W_m2K``: transmitancia térmica = 1/R_total (menor es mejor aislación).

    Fala con ``ValueError`` si una capa no tiene λ cargada y no es cámara de
    aire: así ``calcular_muro_compuesto`` puede decidir dejar la térmica en
    ``None`` sin impedir el cálculo del peso propio.
    """
    detalle: list[dict] = []
    R_capas = 0.0
    for capa in capas:
        clave = str(capa.get("clave", ""))
        grupo = str(capa.get("grupo", ""))
        espesor = float(capa.get("espesor_m", 0.0))
        if espesor <= 0.0:
            continue
        if _es_camara_aire(clave):
            R = _resistencia_camara_aire(espesor)
            lam = None  # el aire no tiene λ: su R es tabulada
        else:
            lam = _lambda_de(clave, grupo)
            R = espesor / lam if lam > 0 else float("inf")
        R_capas += R
        detalle.append({
            "clave": clave,
            "grupo": grupo,
            "espesor_m": espesor,
            "lambda_w_mk": lam,
            "R_m2K_W": R,
        })
    R_total = R_SI + R_capas + R_SE
    U = 1.0 / R_total if R_total > 0 else float("inf")
    return {
        "capas": detalle,
        "R_si_m2K_W": R_SI,
        "R_se_m2K_W": R_SE,
        "R_capas_m2K_W": R_capas,
        "R_total_m2K_W": R_total,
        "U_W_m2K": U,
    }


def calcular_muro_compuesto(
    nombre: str,
    capas: list[dict],
    altura_m: float,
    longitud_m: float,
) -> dict:
    """Calcula un muro estructural compuesto por capas.

    Devuelve peso propio (kN/m² y kN/m lineales), peso total del paño
    (kN) con la geometría dada, y la resistencia térmica (ISO 6946) cuando
    todas las capas tienen conductividad λ o son cámaras de aire. Si falta
    alguna λ, el peso propio se calcula igual y la térmica queda en ``None``
    con aviso. El diseño estructural y el cimiento siguen pendientes.
    """
    base = peso_propio_capas(capas)
    if altura_m <= 0.0:
        raise ValueError("La altura del muro debe ser mayor que cero.")
    peso_lineal_kN_m = base["peso_propio_kN_m2"] * altura_m  # kN por metro lineal
    peso_total_kN = peso_lineal_kN_m * longitud_m

    # Resistencia térmica: se calcula solo si TODAS las capas tienen λ (o son
    # cámaras de aire). Si falta alguna, no se inventa: se deja en None y se
    # avisa, pero el peso propio de arriba sigue siendo válido.
    resistencia_termica: dict | None
    pendientes: list[str] = []
    try:
        resistencia_termica = resistencia_termica_capas(capas)
    except ValueError as exc:
        resistencia_termica = None
        pendientes.append(f"Resistencia térmica incompleta: {exc}")
    pendientes.append(
        "Diseño estructural y cimiento: pendientes (requieren reglamentación "
        "de mampostería)."
    )

    return {
        "nombre": nombre,
        "tipo": "muro",
        "capas": base["capas"],
        "espesor_total_m": base["espesor_total_m"],
        "altura_m": altura_m,
        "longitud_m": longitud_m,
        "peso_propio_kN_m2": base["peso_propio_kN_m2"],
        "peso_lineal_kN_m": peso_lineal_kN_m,
        "peso_total_kN": peso_total_kN,
        # Resistencia térmica (ISO 6946) o None si falta alguna λ:
        "resistencia_termica": resistencia_termica,
        # Pendientes de implementar (requieren reglamentación de mampostería):
        "diseno_estructural": None,
        "cimiento": None,
        "pendiente": " ".join(pendientes),
    }


def calcular_muro(elemento: dict) -> dict:
    """Composición simple (retrocompatibilidad con el cálculo previo)."""
    densidad = materiales.buscar(elemento["clave"])  # kg/m³
    espesor = float(elemento["espesor_m"])
    altura = float(elemento["altura_m"])
    q_kN_m = densidad * espesor * altura * 9.81 / 1000
    return {
        "id": elemento.get("id"),
        "tipo": "muro",
        "carga_kN_m": q_kN_m,
    }


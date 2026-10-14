"""calc.muros — Cálculo de muros ESTRUCTURALES (los de la estructura mixta).

No confundir con los muros del catálogo de Cargas: esos son elementos de carga
(mortero/ladrillo como carga lineal sobre un pórtico) y se calculan en
calc/cargas.py. Acá se calcula el muro portante que se crea en Inicio sobre un
pórtico (estructura mixta).

Estado del cálculo:
    * Peso propio por capas: SÍ (suma γ·e de cada capa con datos de materiales).
    * Peso total y carga lineal resultante: SÍ.
    * Resistencia térmica, resistencia estructural, dimensionado y cimiento:
      PENDIENTES de datos/fórmulas (conductividades λ y reglamentación de
      mampostería). No se inventan: quedan como claves con ``None`` y un aviso.
"""

from __future__ import annotations

from . import materiales


def _gamma_de(clave: str, grupo: str) -> float:
    """Peso específico (kN/m³) de una clave dentro de un grupo de materiales.

    ``grupo`` usa la misma clave que ``materiales.pesos_especificos()``:
    «Mamposteria» o «Morteros».
    """
    pesos = materiales.pesos_especificos()
    grupo_tabla = pesos.get(grupo, {})
    item = grupo_tabla.get(clave)
    if item is None:
        raise ValueError(
            f"No existe la clave «{clave}» en el grupo «{grupo}» de materiales."
        )
    return float(item.get("gamma", item.get("q", 0.0)))


def grupos_disponibles() -> dict[str, dict]:
    """Grupos de materiales útiles para armar capas de un muro (por clave)."""
    pesos = materiales.pesos_especificos()
    return {
        "Mamposteria": pesos.get("Mamposteria", {}),
        "Morteros": pesos.get("Morteros", {}),
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


def calcular_muro_compuesto(
    nombre: str,
    capas: list[dict],
    altura_m: float,
    longitud_m: float,
) -> dict:
    """Calcula un muro estructural compuesto por capas.

    Devuelve peso propio (kN/m² y kN/m lineales), peso total del paño
    (kN) con la geometría dada, y placeholders para los cálculos que aún
    requieren datos/fórmulas (termorresistencia, diseño estructural, cimiento).
    """
    base = peso_propio_capas(capas)
    if altura_m <= 0.0:
        raise ValueError("La altura del muro debe ser mayor que cero.")
    peso_lineal_kN_m = base["peso_propio_kN_m2"] * altura_m  # kN por metro lineal
    peso_total_kN = peso_lineal_kN_m * longitud_m
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
        # Pendientes de implementar (requieren λ y fórmulas de mampostería):
        "resistencia_termica": None,
        "diseno_estructural": None,
        "cimiento": None,
        "pendiente": (
            "Resistencia térmica, diseño estructural y cimiento requieren "
            "conductividades λ y reglamentación de mampostería (aún no cargadas)."
        ),
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


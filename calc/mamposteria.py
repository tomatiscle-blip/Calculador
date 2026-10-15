"""calc.mamposteria — Verificación estructural de muros de mampostería (CIRSOC 501).

Reglamento: CIRSOC 501 (Edición Septiembre 2007, INTI). Los datos (factores φ,
resistencias fm Tabla 6.3, fórmulas de Pn y Vn, cuantías de armadura) se leen de
``datos/mamposteria.json`` (biblioteca GLOBAL). No se inventan valores: si falta
un dato, se avisa.

Unidades del reglamento: fm y fy en MPa (N/mm²), áreas en mm², lo que da fuerzas
en N (se convierten a kN al presentar). El módulo acepta metros donde tiene
sentido (espesor, altura, longitud) y convierte internamente.

Estado del cálculo:
    * Compresión axial Pn (Art. 9.1.6): SÍ, dos fórmulas según esbeltez h/r.
    * Corte Vn = Vm + Vs: SÍ (con topes).
    * Recomendación de armadura (horizontal, vertical/columnas, dinteles): SÍ.
    * Fuera del plano (viento/sismo): Mcr = Sn·fr disponible.

Chequeo orientativo: usa el PESO PROPIO del muro (calc.muros) como carga axial.
"""

from __future__ import annotations

import math

from . import rutas


def cargar(biblioteca: dict | None = None) -> dict:
    """Lee datos/mamposteria.json. Avisa claro si no está."""
    if biblioteca is not None:
        return biblioteca
    datos = rutas.leer_json(rutas.MAMPOSTERIA, por_defecto=None)
    if not datos:
        raise FileNotFoundError(
            f"No se pudo leer la biblioteca de mampostería: {rutas.MAMPOSTERIA}"
        )
    return datos


def _mapeo_unidad(datos: dict) -> dict:
    return {
        k: v for k, v in (datos.get("mapeo_material_a_unidad", {}) or {}).items()
        if not str(k).startswith("_")
    }


def fm_de(clave_material: str, mortero: str = "mortero_normal_n",
          datos: dict | None = None) -> float:
    """Resistencia característica a compresión fm (MPa) de una hoja de mampostería.

    ``clave_material`` es la clave de materiales.json (ej. ``ladrillo_hueco_portante``).
    ``mortero`` es ``mortero_normal_n`` (por defecto), ``mortero_intermedio_i`` o
    ``mortero_elevado_e``.

    Fala con ``ValueError`` si el material no está mapeado a un tipo de unidad de
    la Tabla 6.3 (no se inventa fm para materiales no tabulados).
    """
    datos = cargar(datos)
    unidad = _mapeo_unidad(datos).get(str(clave_material))
    if unidad is None:
        raise ValueError(
            f"El material «{clave_material}» no está mapeado a un tipo de unidad "
            f"de la Tabla 6.3 del CIRSOC 501. Agregalo en datos/mamposteria.json "
            f"(mapeo_material_a_unidad)."
        )
    mortero_clave = mortero if mortero.startswith("mortero_") else f"mortero_{mortero}"
    tabla = datos["resistencias_caracteristicas_fm_indicativas_tabla_6_3"][unidad]
    if mortero_clave not in tabla:
        raise ValueError(f"No hay fm para el mortero «{mortero_clave}» en «{unidad}».")
    return float(tabla[mortero_clave])


def radio_de_giro_placa(espesor_mm: float) -> float:
    """Radio de giro r (mm) de la sección rectangular del muro (dirección del espesor)."""
    return espesor_mm / math.sqrt(12.0)


def esbeltez(altura_mm: float, espesor_mm: float) -> float:
    """Relación de esbeltez h/r del muro (adimensional)."""
    r = radio_de_giro_placa(espesor_mm)
    return altura_mm / r if r > 0 else float("inf")


def resistencia_axial_Pn(
    fm_mpa: float,
    an_mm2: float,
    fy_mpa: float,
    as_mm2: float,
    altura_mm: float,
    espesor_mm: float,
    datos: dict | None = None,
) -> dict:
    """Resistencia nominal a compresión Pn (N) del muro (Art. 9.1.6).

    Usa (1 - (h/140r)²) si h/r <= 99, o (70r/h)² si h/r > 99. Para mampostería
    simple, pasar ``as_mm2 = 0``.
    """
    datos = cargar(datos)
    c = datos["calculo_estructural"]["compresion_axial"]
    factor_global = float(c["factor_global"])       # 0.80
    factor_mamp = float(c["factor_mamposteria"])    # 0.80
    limite_hr = float(c["limite_esbeltez_hr"])      # 99

    hr = esbeltez(altura_mm, espesor_mm)
    if hr <= limite_hr:
        reduccion = 1.0 - (hr / 140.0) ** 2
        rama = "h/r <= 99"
    else:
        r = radio_de_giro_placa(espesor_mm)
        reduccion = (70.0 * r / altura_mm) ** 2 if altura_mm > 0 else 0.0
        rama = "h/r > 99"

    pn = factor_global * (factor_mamp * fm_mpa * (an_mm2 - as_mm2) + fy_mpa * as_mm2) * reduccion
    return {"Pn_N": pn, "Pn_kN": pn / 1000.0, "esbeltez_hr": hr, "rama": rama,
            "factor_reduccion": reduccion}


def phi(datos: dict | None = None, armada: bool = False) -> float:
    """Factor de reducción φ según sea mampostería armada o simple (flexión/axial)."""
    datos = cargar(datos)
    phi_tabla = datos["factores_reduccion_resistencia_phi"]
    clave = ("flexion_axial_mamposteria_armada" if armada
             else "flexion_axial_mamposteria_simple")
    return float(phi_tabla[clave])


def resistencia_corte_Vn(
    fm_mpa: float,
    an_mm2: float,
    p_axial_N: float,
    m_kNm: float,
    v_kN: float,
    dv_mm: float,
    fy_mpa: float,
    av_mm2: float = 0.0,
    s_mm: float = 0.0,
    datos: dict | None = None,
) -> dict:
    """Resistencia nominal al corte Vn = Vm + Vs (N), con topes (Art. 9.1.6).

    ``m_kNm``/``v_kN`` son el momento y corte concomitantes; con M=V=0 queda el
    término M/(V·dv) = 0 (caso solo axial).
    """
    datos = cargar(datos)
    c = datos["calculo_estructural"]["corte"]
    coef_vm = float(c["coef_vm"])                 # 0.083
    aporte_p = float(c["aporte_carga_axial_Vm"])  # 0.25
    coef_vs = float(c["coef_vs"])                 # 0.5
    mv_inf = float(c["mv_limite_inferior"])       # 0.25
    mv_sup = float(c["mv_limite_superior"])       # 1.0
    lim_bajo = float(c["limite_Vn_mv_bajo_0_25"])  # 0.5
    lim_alto = float(c["limite_Vn_mv_sobre_1_0"])  # 0.3

    v_N = v_kN * 1000.0
    m_Nmm = m_kNm * 1e6
    mv = (m_Nmm / (v_N * dv_mm)) if (v_N > 0 and dv_mm > 0) else 0.0
    termino = max(0.0, 4.0 - 1.75 * mv)
    vm = coef_vm * termino * an_mm2 * math.sqrt(fm_mpa) + aporte_p * p_axial_N
    vs = coef_vs * (av_mm2 / s_mm) * fy_mpa * dv_mm if s_mm > 0 else 0.0
    vn = vm + vs

    sqrt_fm = math.sqrt(fm_mpa)
    if mv <= mv_inf:
        tope = lim_bajo * an_mm2 * sqrt_fm
    elif mv >= mv_sup:
        tope = lim_alto * an_mm2 * sqrt_fm
    else:
        t = (mv - mv_inf) / (mv_sup - mv_inf)
        coef = lim_bajo + t * (lim_alto - lim_bajo)
        tope = coef * an_mm2 * sqrt_fm

    return {"Vm_N": vm, "Vs_N": vs, "Vn_N": vn, "Vn_kN": vn / 1000.0,
            "tope_Vn_N": tope, "mv": mv, "supera_tope": vn > tope}


def _hoja_portante(capas: list[dict]) -> dict | None:
    """Hoja resistente = la capa de Mampostería de mayor espesor.

    Los revoques (Morteros) y aislantes no se cuentan como portantes: de la hoja
    de mampostería salen fm y el espesor t del muro.
    """
    mamp = [c for c in capas
            if str(c.get("grupo", "")).lower() == "mamposteria"
            and float(c.get("espesor_m", 0.0) or 0.0) > 0.0]
    if not mamp:
        return None
    return max(mamp, key=lambda c: float(c.get("espesor_m", 0.0)))


def verificar_muro(
    capas: list[dict],
    altura_m: float,
    mortero: str = "mortero_normal_n",
    armada: bool = False,
    carga_axial_kN_m: float | None = None,
    longitud_referencia_m: float = 1.0,
    fy_mpa: float | None = None,
    datos: dict | None = None,
) -> dict:
    """Verificación estructural orientativa de un muro (CIRSOC 501).

    Toma la hoja de Mampostería de mayor espesor como portante (de ahí fm y t),
    calcula esbeltez h/r, capacidad a compresión φ·Pn y a corte φ·Vn, y arma
    recomendaciones de armadura (horizontal, vertical/columnas, dintel).

    ``carga_axial_kN_m``: carga axial de compresión por metro de muro en la
    sección crítica. Si es ``None``, no se sentencia pasa/no pasa (solo capacidad).

    Devuelve un dict con el detalle y, si no se puede calcular (falta fm mapeado
    u otro dato), deja el campo en ``None`` y explica en ``pendiente``.
    """
    datos = cargar(datos)
    ce = datos["calculo_estructural"]
    fy = float(fy_mpa) if fy_mpa is not None else float(ce["compresion_axial"]["fy_default_mpa"])

    hoja = _hoja_portante(capas)
    if hoja is None or altura_m <= 0:
        return {"ok": False,
                "pendiente": "El muro no tiene una hoja de Mampostería o falta la altura."}

    clave = str(hoja["clave"])
    espesor_m = float(hoja["espesor_m"])
    t_mm = espesor_m * 1000.0
    h_mm = altura_m * 1000.0

    try:
        fm = fm_de(clave, mortero, datos)
    except ValueError as exc:
        return {"ok": False, "pendiente": str(exc)}

    hr = esbeltez(h_mm, t_mm)
    es_esbelta = hr > float(datos["limitaciones_geometricas_y_esbeltez"]["relacion_esbeltez_maxima_hr_99"])

    # Franja de referencia (por defecto 1 m de longitud de muro).
    L_ref_mm = longitud_referencia_m * 1000.0
    An = t_mm * L_ref_mm  # mm²

    # Mampostería simple: sin armadura (As = 0).
    pn = resistencia_axial_Pn(fm, An, fy, 0.0, h_mm, t_mm, datos)
    factor = phi(datos, armada=armada)
    phi_Pn_kN = factor * pn["Pn_kN"]

    # Corte basal de la franja: solo aporta la mampostería (sin armadura de corte).
    p_axial_N = (carga_axial_kN_m * 1000.0) if carga_axial_kN_m is not None else 0.0
    vn = resistencia_corte_Vn(fm, An, p_axial_N, 0.0, 0.0, 0.8 * h_mm, fy, 0.0, 0.0, datos)
    phi_corte = float(datos["factores_reduccion_resistencia_phi"]["corte"])
    phi_Vn_kN = phi_corte * vn["Vn_kN"]

    compresion_ok = None
    if carga_axial_kN_m is not None:
        compresion_ok = carga_axial_kN_m <= phi_Pn_kN

    return {
        "ok": True,
        "hoja_portante": {
            "clave": clave,
            "fm_mpa": fm,
            "espesor_m": espesor_m,
            "mortero": mortero,
        },
        "esbeltez_hr": hr,
        "es_esbelta": es_esbelta,
        "area_franja_m2": An / 1e6,
        "compresion": {
            "Pn_kN": pn["Pn_kN"],
            "phi": factor,
            "phi_Pn_kN": phi_Pn_kN,
            "rama": pn["rama"],
            "factor_reduccion": pn["factor_reduccion"],
            "carga_axial_kN_m": carga_axial_kN_m,
            "pasa": compresion_ok,
        },
        "corte": {
            "Vn_kN": vn["Vn_kN"],
            "phi": phi_corte,
            "phi_Vn_kN": phi_Vn_kN,
        },
        "armada": armada,
        "fy_mpa": fy,
        "recomendaciones": _recomendar_armadura(datos, armada, es_esbelta, compresion_ok, An, fy),
        "pendiente": None,
    }


def _recomendar_armadura(datos: dict, armada: bool, es_esbelta: bool,
                         compresion_ok: bool | None, an_mm2: float, fy_mpa: float) -> list[str]:
    """Recomendaciones de armadura con los datos del reglamento. No inventa:
    lo que requiere un dato no cargado queda marcado como pendiente."""
    arma = datos["calculo_estructural"]["armadura"]
    horiz = arma["refuerzo_horizontal"]
    vert = arma["armadura_vertical_columnas"]
    recos: list[str] = []

    # Armadura horizontal en juntas.
    recos.append(
        f"Refuerzo horizontal en juntas: Ø máximo {horiz.get('diametro_maximo_juntas_mm', '?')} mm "
        f"(no más de la mitad del espesor de junta). ρ/espaciamiento mínimo: pendiente dato."
    )
    # Columnas de refuerzo (vertical).
    rho_min = float(vert["cuantia_minima"])
    rho_max = float(vert["cuantia_maxima"])
    as_min_cm2 = rho_min * an_mm2 / 100.0
    recos.append(
        f"Columnas de refuerzo (extremos y encuentros): cuantía {rho_min}–{rho_max} de An "
        f"→ As mín ≳ {as_min_cm2:.2f} cm² por franja; mínimo {vert.get('minimo_barras', 4)} barras; "
        f"estribos Ø ≥ {vert.get('estribos_diametro_min_mm', 6)} mm."
    )
    # Dintel sobre vanos.
    recos.append(
        "Dintel sobre vanos: " + str(arma["dinteles"]["criterio"])
        + " Luz libre máxima que dispara dintel: pendiente dato."
    )
    # Sentencias según el chequeo.
    if es_esbelta:
        recos.append("⚠ Muro ESBELTO (h/r > 99): reducir altura, aumentar espesor o armar/confinar.")
    if compresion_ok is False:
        recos.append("⚠ NO pasa compresión con la carga axial dada: aumentar espesor o usar mampostería armada.")
    if not armada:
        recos.append("Nota: calculado como mampostería SIMPLE (sin aporte de armadura a Pn). "
                     "En zona 0 el armado no es obligatorio por sismo, pero verificar flexión fuera del plano.")
    return recos



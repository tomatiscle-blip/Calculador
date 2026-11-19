"""Análisis elástico de paños macizos con elementos de placa de PyNite."""

from __future__ import annotations

import math

import numpy as np
from Pynite import FEModel3D

from . import cargas, configuracion_obra


def analizar_paño_macizo(elemento: dict, valores: dict, espesor_m: float,
                         malla_m: float | None = None) -> dict:
    """Malla y resuelve una placa rectangular cargada con D/L.

    El caso unidireccional apoya dos bordes opuestos. Las placas cruzadas
    aceptan tres o cuatro bordes apoyados. Una casetonada puede modelarse en
    una dirección cuando así se indica en ``panel_ejes``.
    """
    luz = float(elemento.get("luz_transversal_m", 0.0) or 0.0)
    ancho = float(elemento.get("ancho_losa_m", 0.0) or 0.0)
    if luz <= 0 or ancho <= 0 or espesor_m <= 0:
        raise ValueError("La losa necesita luz, ancho y espesor estructural positivos.")
    tipologia = elemento.get("tipologia")
    if tipologia not in ("maciza", "maciza_cruzada", "casetonada"):
        raise ValueError("El análisis de placa requiere una losa maciza o casetonada.")
    apoyos = elemento.get("apoya_en", {}) or {}
    apoyos_seleccionados = [
        lado for lado in ("izq", "der", "inf", "sup") if apoyos.get(lado)
    ]
    bordes_apoyados = (
        ["izq", "der"]
        if tipologia == "maciza"
        else apoyos_seleccionados
    )
    if not all(apoyos.get(lado) for lado in ("izq", "der")):
        raise ValueError("La losa necesita apoyo en los bordes opuestos izq y der.")
    if tipologia == "maciza" and len(apoyos_seleccionados) > 2:
        raise ValueError(
            "Una losa maciza con tres o cuatro bordes apoyados debe usar tipología cruzada."
        )
    modelo_placa = (elemento.get("panel_ejes") or {}).get("modelo_placa")
    casetonada_unidireccional = (
        tipologia == "casetonada" and modelo_placa == "casetonada_unidireccional"
    )
    if tipologia == "maciza_cruzada" and len(bordes_apoyados) < 3:
        raise ValueError("La losa cruzada necesita al menos tres bordes apoyados.")
    if tipologia == "casetonada":
        if casetonada_unidireccional and len(bordes_apoyados) != 2:
            raise ValueError(
                "La casetonada unidireccional necesita exactamente dos bordes opuestos."
            )
        if not casetonada_unidireccional and len(bordes_apoyados) < 3:
            raise ValueError(
                "La casetonada cruzada necesita al menos tres bordes apoyados."
            )

    seccion_casetonada = elemento.get("casetonada", {}) or {}
    kx_mod = float(seccion_casetonada.get("kx_mod", 1.0)) if tipologia == "casetonada" else 1.0
    ky_mod = float(seccion_casetonada.get("ky_mod", 1.0)) if tipologia == "casetonada" else 1.0
    if not all(math.isfinite(factor) and 0 < factor <= 1 for factor in (kx_mod, ky_mod)):
        raise ValueError("Los factores de rigidez kx y ky deben ser mayores que 0 y no mayores que 1.")

    ec = configuracion_obra.modulo_elasticidad_hormigon()
    E = float(ec["Ec_MPa"]) * 1000.0  # kN/m²
    nu = 0.20
    modelo = FEModel3D()
    modelo.add_material("Hormigon_losa", E, E / (2 * (1 + nu)), nu, 0.0)
    if malla_m is None:
        malla_m = max(0.25, min(0.75, min(luz, ancho) / 8.0))
    modelo.add_rectangle_mesh(
        "Paño", float(malla_m), luz, ancho, espesor_m, "Hormigon_losa",
        kx_mod=kx_mod, ky_mod=ky_mod,
        origin=(0.0, 0.0, 0.0), plane="XY", element_type="Quad",
    )
    modelo.meshes["Paño"].generate()

    tol = 1e-7
    nodos = list(modelo.nodes.values())
    bordes = {
        "izq": lambda n: abs(n.X) <= tol,
        "der": lambda n: abs(n.X - luz) <= tol,
        "inf": lambda n: abs(n.Y) <= tol,
        "sup": lambda n: abs(n.Y - ancho) <= tol,
    }
    # Evitar modos rígidos en el plano sin empotrar la flexión de placa.
    esquina = min(nodos, key=lambda n: n.X**2 + n.Y**2)
    segundo = min((n for n in nodos if n.name != esquina.name),
                  key=lambda n: abs(n.Y) + abs(n.X - min(luz, ancho)))
    for nodo in nodos:
        vertical = any(bordes[lado](nodo) for lado in bordes_apoyados)
        modelo.def_support(
            nodo.name,
            support_DX=nodo.name in {esquina.name, segundo.name},
            support_DY=nodo.name == esquina.name,
            support_DZ=vertical,
        )

    D = float(valores["D_kNm2"])
    L = float(valores["L_kNm2"])
    for quad in modelo.quads.values():
        modelo.add_quad_surface_pressure(quad.name, -D, case="D")
        if L:
            modelo.add_quad_surface_pressure(quad.name, -L, case="L")
    combinaciones = {
        nombre: fD * D + fL * L
        for nombre, (fD, fL, fW) in cargas.FACTORES_COMBINACIONES.items()
        if fW == 0
    }
    modelo.add_load_combo("D_servicio", {"D": 1.0})
    modelo.add_load_combo("L_servicio", {"L": 1.0})
    for nombre, factores in cargas.FACTORES_COMBINACIONES.items():
        fD, fL, fW = factores
        if not fW:
            modelo.add_load_combo(nombre, {"D": fD, "L": fL})
    modelo.analyze_linear(check_stability=True, check_statics=False)

    resultados = {}
    momentos_apoyos = {}
    coordenadas = (-1 / math.sqrt(3), 0.0, 1 / math.sqrt(3))
    for nombre_combo in combinaciones:
        if nombre_combo not in modelo.load_combos:
            continue
        componentes = {k: [] for k in ("Mx", "My", "Mxy", "Qx", "Qy")}
        for quad in modelo.quads.values():
            for xi in coordenadas:
                for eta in coordenadas:
                    m = np.asarray(quad.moment(xi, eta, local=True, combo_name=nombre_combo)).reshape(-1)
                    q = np.asarray(quad.shear(xi, eta, local=True, combo_name=nombre_combo)).reshape(-1)
                    for clave, valor in zip(("Mx", "My", "Mxy"), m[:3]):
                        componentes[clave].append(float(valor))
                    for clave, valor in zip(("Qx", "Qy"), q[:2]):
                        componentes[clave].append(float(valor))
        resultados[nombre_combo] = {
            k: {"max": max(v, default=0.0), "min": min(v, default=0.0),
                "max_abs": max((abs(x) for x in v), default=0.0)}
            for k, v in componentes.items()
        }
        momentos_apoyos[nombre_combo] = {
            lado: {
                "componente": "Mx" if lado in ("izq", "der") else "My",
                "max": 0.0,
                "min": 0.0,
                "max_abs": 0.0,
                "metodo": "apoyo simple idealizado",
            }
            for lado in bordes_apoyados
        }

    reacciones = []
    for lado in bordes_apoyados:
        nodos_borde = [n for n in nodos if bordes[lado](n)]
        total_por_caso = {}
        for caso, combo in (("D", "D_servicio"), ("L", "L_servicio")):
            total = 0.0
            for n in nodos_borde:
                # Las esquinas compartidas se distribuyen entre los bordes que
                # realmente las sostienen para evitar contarlas dos veces.
                cantidad_bordes = sum(bordes[b](n) for b in bordes_apoyados)
                rxns = n.RxnFZ or {}
                total += float(rxns.get(combo, 0.0)) / max(1, cantidad_bordes)
            total_por_caso[caso] = abs(total)
        longitud_borde = ancho if lado in ("izq", "der") else luz
        reacciones.append({
            "apoyo": lado,
            "apoya_en": apoyos.get(lado),
            "reaccion_total_D_kN": total_por_caso["D"],
            "reaccion_total_L_kN": total_por_caso["L"],
            "cargas_kN_m": {
                "D": total_por_caso["D"] / longitud_borde,
                "L": total_por_caso["L"] / longitud_borde,
            },
        })
    nombres_modelo = {
        "maciza": "maciza unidireccional (dos bordes apoyados)",
        "maciza_cruzada": f"maciza cruzada ({len(bordes_apoyados)} bordes apoyados)",
        "casetonada": (
            "casetonada equivalente ortótropa unidireccional (dos bordes apoyados)"
            if casetonada_unidireccional
            else f"casetonada equivalente ortótropa ({len(bordes_apoyados)} bordes apoyados)"
        ),
    }
    return {
        "metodo": "PyNite FEModel3D · placa Quad mallada · análisis lineal",
        "modelo": nombres_modelo[tipologia],
        "Ec_MPa": float(ec["Ec_MPa"]),
        "Ec_referencia": ec["articulo"],
        "kx_mod": kx_mod,
        "ky_mod": ky_mod,
        "malla_aproximada_m": float(malla_m),
        "cantidad_elementos": len(modelo.quads),
        "cargas": {"D_kNm2": D, "L_kNm2": L},
        "combinaciones_kNm2": combinaciones,
        "solicitaciones": resultados,
        "momentos_apoyos": momentos_apoyos,
        "reacciones": reacciones,
        "alcance": "Placa maciza elástica; apoyos simples verticales en los bordes indicados. No dimensiona armaduras, no modela continuidad ni verifica flecha, corte o fisuración.",
    }

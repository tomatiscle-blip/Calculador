"""Vinculacion y analisis conjunto de paños macizos continuos."""

from __future__ import annotations

import math

import numpy as np
from Pynite import FEModel3D

from . import cargas, configuracion_obra

_TOL = 1e-7
_TIPOLOGIAS_CONTINUAS = {"maciza", "maciza_cruzada"}


def _espesor(elemento: dict) -> float | None:
    espesores = []
    for componente in elemento.get("componentes", []) or []:
        if (componente.get("grupo"), componente.get("clave")) in {
            ("Forjados", "Losa_maciza"),
            ("Hormigon", "armado"),
        }:
            try:
                espesor = float(componente.get("espesor_m", 0.0) or 0.0)
            except (TypeError, ValueError):
                return None
            if not math.isfinite(espesor) or espesor < 0:
                return None
            espesores.append(espesor)
    total = sum(espesores)
    return total if total > 0 else None


def _geometria(elemento: dict) -> dict | None:
    panel = elemento.get("panel_ejes") or {}
    ejes = panel.get("ejes") or {}
    try:
        x0 = float(ejes["x"]["inicio"]["coordenada_m"])
        x1 = float(ejes["x"]["fin"]["coordenada_m"])
        y0 = float(ejes["y"]["inicio"]["coordenada_m"])
        y1 = float(ejes["y"]["fin"]["coordenada_m"])
        libre0 = float(panel["coordenada_libre_inicio_m"])
        libre1 = float(panel["coordenada_libre_fin_m"])
        cota = float(panel["cota_apoyo_m"])
    except (KeyError, TypeError, ValueError):
        return None
    direccion = panel.get("direccion_luz")
    if direccion == "x":
        x_min, x_max = sorted((x0, x1))
        y_min, y_max = sorted((libre0, libre1))
        lados = {
            "izq": ("x", x_min),
            "der": ("x", x_max),
            "inf": ("y", y_min),
            "sup": ("y", y_max),
        }
    elif direccion == "y":
        x_min, x_max = sorted((libre0, libre1))
        y_min, y_max = sorted((y0, y1))
        lados = {
            "izq": ("y", y_min),
            "der": ("y", y_max),
            "inf": ("x", x_min),
            "sup": ("x", x_max),
        }
    else:
        return None
    if x_max - x_min <= _TOL or y_max - y_min <= _TOL or not math.isfinite(cota):
        return None
    return {
        "x_min": x_min,
        "x_max": x_max,
        "y_min": y_min,
        "y_max": y_max,
        "cota": cota,
        "lados": lados,
    }


def _apoyo_en_borde(elemento: dict, geometria: dict, eje: str, coordenada: float) -> dict | None:
    apoyos = elemento.get("apoya_en", {}) or {}
    for lado, (eje_lado, coord_lado) in geometria["lados"].items():
        if eje_lado == eje and abs(coord_lado - coordenada) <= _TOL:
            apoyo = apoyos.get(lado)
            return apoyo if isinstance(apoyo, dict) else None
    return None


def _clave_apoyo(apoyo: dict | None) -> tuple[str, ...] | None:
    if not apoyo:
        return None
    tipo = str(apoyo.get("tipo", ""))
    if tipo == "portico":
        return tipo, str(apoyo.get("portico", "")), str(apoyo.get("viga", ""))
    if tipo == "muro":
        return tipo, str(apoyo.get("nombre", ""))
    return None


def borde_compartido(elemento_a: dict, elemento_b: dict) -> dict | None:
    """Devuelve el borde comun si ambos paños son compatibles para continuidad."""
    if elemento_a.get("tipologia") not in _TIPOLOGIAS_CONTINUAS:
        return None
    if elemento_b.get("tipologia") not in _TIPOLOGIAS_CONTINUAS:
        return None
    geometria_a = _geometria(elemento_a)
    geometria_b = _geometria(elemento_b)
    espesor_a = _espesor(elemento_a)
    espesor_b = _espesor(elemento_b)
    if not geometria_a or not geometria_b or espesor_a is None or espesor_b is None:
        return None
    if abs(geometria_a["cota"] - geometria_b["cota"]) > _TOL:
        return None
    if abs(espesor_a - espesor_b) > _TOL:
        return None

    bordes = (
        ("x", geometria_a["x_max"], geometria_b["x_min"], "y"),
        ("x", geometria_a["x_min"], geometria_b["x_max"], "y"),
        ("y", geometria_a["y_max"], geometria_b["y_min"], "x"),
        ("y", geometria_a["y_min"], geometria_b["y_max"], "x"),
    )
    for eje, coord_a, coord_b, eje_longitudinal in bordes:
        if abs(coord_a - coord_b) > _TOL:
            continue
        minimo = f"{eje_longitudinal}_min"
        maximo = f"{eje_longitudinal}_max"
        if (
            abs(geometria_a[minimo] - geometria_b[minimo]) > _TOL
            or abs(geometria_a[maximo] - geometria_b[maximo]) > _TOL
        ):
            continue
        apoyo_a = _apoyo_en_borde(elemento_a, geometria_a, eje, coord_a)
        apoyo_b = _apoyo_en_borde(elemento_b, geometria_b, eje, coord_b)
        clave_a = _clave_apoyo(apoyo_a)
        if clave_a is None or clave_a != _clave_apoyo(apoyo_b):
            continue
        return {
            "eje": eje,
            "coordenada_m": coord_a,
            "eje_longitudinal": eje_longitudinal,
            "inicio_m": geometria_a[minimo],
            "fin_m": geometria_a[maximo],
            "apoyo": apoyo_a,
        }
    return None


def candidatos(elemento: dict, elementos: dict[str, dict], excluir: str | None = None) -> list[str]:
    """Lista paños del catálogo que pueden compartir continuidad con el dado."""
    salida = []
    for nombre, otro in elementos.items():
        if nombre == excluir or otro.get("tipo") != "losa":
            continue
        if borde_compartido(elemento, otro):
            salida.append(nombre)
    return salida


def actualizar_vinculos(
    elementos: dict[str, dict], nombre: str, nombres_continuos: set[str]
) -> None:
    """Guarda vinculos simetricos por ID para que sobrevivan a renombres."""
    elemento = elementos[nombre]
    identificador = str(elemento.get("id", ""))
    if not identificador:
        raise ValueError(f"La losa {nombre} no tiene un ID para guardar su continuidad.")
    ids_aceptados = set()
    for otro_nombre in nombres_continuos:
        otro = elementos.get(otro_nombre)
        if not otro or not otro.get("id"):
            raise ValueError(f"No se puede vincular la losa {otro_nombre}: falta su ID.")
        if not borde_compartido(elemento, otro):
            raise ValueError(f"Los paños {nombre} y {otro_nombre} ya no son compatibles.")
        ids_aceptados.add(str(otro["id"]))
    panel = elemento.get("panel_ejes")
    if not isinstance(panel, dict) and ids_aceptados:
        raise ValueError(f"La losa {nombre} no fue definida desde ejes.")

    for otro in elementos.values():
        panel_otro = otro.get("panel_ejes")
        if not isinstance(panel_otro, dict):
            continue
        enlaces = {
            str(vinculo)
            for vinculo in panel_otro.get("continuidad_con", [])
            if str(vinculo) != identificador
        }
        otro_id = str(otro.get("id", ""))
        if otro is not elemento and otro_id in ids_aceptados:
            enlaces.add(identificador)
        if enlaces:
            panel_otro["continuidad_con"] = sorted(enlaces)
        else:
            panel_otro.pop("continuidad_con", None)

    if not isinstance(panel, dict):
        return
    if ids_aceptados:
        panel["continuidad_con"] = sorted(ids_aceptados)
    else:
        panel.pop("continuidad_con", None)


def quitar_vinculos(elementos: dict[str, dict], identificador: str) -> None:
    """Elimina referencias a una losa borrada de los paños restantes."""
    for elemento in elementos.values():
        panel = elemento.get("panel_ejes")
        if not isinstance(panel, dict):
            continue
        enlaces = [
            str(vinculo)
            for vinculo in panel.get("continuidad_con", [])
            if str(vinculo) != str(identificador)
        ]
        if enlaces:
            panel["continuidad_con"] = enlaces
        else:
            panel.pop("continuidad_con", None)


def grupo_continuo(nombre: str, elementos: dict[str, dict]) -> dict[str, dict]:
    """Resuelve el grupo de paños vinculados a una losa y valida sus enlaces."""
    por_id: dict[str, str] = {}
    for otro_nombre, elemento in elementos.items():
        if elemento.get("tipo") != "losa":
            continue
        identificador = str(elemento.get("id", ""))
        if identificador:
            if identificador in por_id:
                raise ValueError(f"Hay IDs duplicados de losa: {identificador}.")
            por_id[identificador] = otro_nombre
    if nombre not in elementos:
        raise KeyError(f"No existe la losa {nombre}.")
    grupo = {nombre}
    pendientes = [nombre]
    while pendientes:
        actual = pendientes.pop()
        elemento = elementos[actual]
        panel = elemento.get("panel_ejes") or {}
        for id_vecino in panel.get("continuidad_con", []) or []:
            vecino = por_id.get(str(id_vecino))
            if vecino is None:
                raise ValueError(
                    f"La continuidad de {actual} refiere a una losa inexistente; "
                    "editá el vínculo de continuidad."
                )
            elemento_vecino = elementos[vecino]
            enlaces_vecino = {
                str(v) for v in (elemento_vecino.get("panel_ejes") or {}).get(
                    "continuidad_con", []
                ) or []
            }
            if str(elemento.get("id", "")) not in enlaces_vecino:
                raise ValueError(
                    f"La continuidad entre {actual} y {vecino} no es simétrica; "
                    "editá el vínculo de continuidad."
                )
            if not borde_compartido(elemento, elemento_vecino):
                raise ValueError(
                    f"La continuidad entre {actual} y {vecino} ya no es compatible "
                    "por geometría, nivel, espesor o apoyo compartido."
                )
            if vecino not in grupo:
                grupo.add(vecino)
                pendientes.append(vecino)
    return {elemento_nombre: elementos[elemento_nombre] for elemento_nombre in sorted(grupo)}


def analizar_grupo(
    paños: list[dict], malla_m: float | None = None
) -> dict[str, dict]:
    """Resuelve paños macizos coplanares con malla y nodos compartidos."""
    if len(paños) < 2:
        raise ValueError("El análisis conjunto requiere al menos dos paños.")
    geometria = {}
    for paño in paños:
        nombre = paño["nombre"]
        elemento = paño["elemento"]
        apoyos = elemento.get("apoya_en", {}) or {}
        if not apoyos.get("izq") or not apoyos.get("der"):
            raise ValueError(f"La losa {nombre} necesita sus dos apoyos principales.")
        if elemento.get("tipologia") == "maciza":
            if apoyos.get("inf") or apoyos.get("sup"):
                raise ValueError(
                    f"La losa maciza unidireccional {nombre} no puede tener apoyos "
                    "perpendiculares adicionales."
                )
        elif sum(bool(apoyos.get(lado)) for lado in ("izq", "der", "inf", "sup")) < 3:
            raise ValueError(f"La losa cruzada {nombre} necesita al menos tres apoyos.")
        geom = _geometria(elemento)
        if geom is None:
            raise ValueError(f"La losa {nombre} no tiene geometría compatible desde ejes.")
        geometria[nombre] = geom
    espesores = [float(paño["espesor_m"]) for paño in paños]
    if max(espesores) - min(espesores) > _TOL:
        raise ValueError("El análisis continuo requiere el mismo espesor estructural en todos los paños.")
    for i, paño_a in enumerate(paños):
        ga = geometria[paño_a["nombre"]]
        for paño_b in paños[i + 1:]:
            gb = geometria[paño_b["nombre"]]
            overlap_x = min(ga["x_max"], gb["x_max"]) - max(ga["x_min"], gb["x_min"])
            overlap_y = min(ga["y_max"], gb["y_max"]) - max(ga["y_min"], gb["y_min"])
            if overlap_x > _TOL and overlap_y > _TOL:
                raise ValueError("Los paños continuos no pueden superponerse en superficie.")
            se_tocan = (
                abs(ga["x_max"] - gb["x_min"]) <= _TOL
                or abs(ga["x_min"] - gb["x_max"]) <= _TOL
                or abs(ga["y_max"] - gb["y_min"]) <= _TOL
                or abs(ga["y_min"] - gb["y_max"]) <= _TOL
            )
            if se_tocan:
                borde = borde_compartido(paño_a["elemento"], paño_b["elemento"])
                if borde is None:
                    raise ValueError(
                        f"Los paños {paño_a['nombre']} y {paño_b['nombre']} se tocan, "
                        "pero no comparten un borde de apoyo compatible."
                    )
                id_a = str(paño_a["elemento"].get("id", ""))
                id_b = str(paño_b["elemento"].get("id", ""))
                enlaces_a = {
                    str(v)
                    for v in (paño_a["elemento"].get("panel_ejes") or {}).get(
                        "continuidad_con", []
                    ) or []
                }
                enlaces_b = {
                    str(v)
                    for v in (paño_b["elemento"].get("panel_ejes") or {}).get(
                        "continuidad_con", []
                    ) or []
                }
                if id_b not in enlaces_a or id_a not in enlaces_b:
                    raise ValueError(
                        f"Los paños {paño_a['nombre']} y {paño_b['nombre']} comparten "
                        "geometría, pero su continuidad no fue confirmada."
                    )

    if malla_m is None:
        dimension_minima = min(
            min(geom["x_max"] - geom["x_min"], geom["y_max"] - geom["y_min"])
            for geom in geometria.values()
        )
        malla_m = max(0.25, min(0.75, dimension_minima / 8.0))
    if not math.isfinite(malla_m) or malla_m <= 0:
        raise ValueError("El paso de malla debe ser positivo y finito.")

    xs = set()
    ys = set()
    for geom in geometria.values():
        for eje, minimo, maximo in (
            ("x", geom["x_min"], geom["x_max"]),
            ("y", geom["y_min"], geom["y_max"]),
        ):
            n = max(1, math.ceil((maximo - minimo) / malla_m))
            coordenadas = np.linspace(minimo, maximo, n + 1)
            (xs if eje == "x" else ys).update(round(float(c), 9) for c in coordenadas)
    xs_ordenados = sorted(xs)
    ys_ordenados = sorted(ys)
    indice_x = {x: i for i, x in enumerate(xs_ordenados)}
    indice_y = {y: i for i, y in enumerate(ys_ordenados)}

    ec = configuracion_obra.modulo_elasticidad_hormigon()
    E = float(ec["Ec_MPa"]) * 1000.0
    nu = 0.20
    modelo = FEModel3D()
    modelo.add_material("Hormigon_losa", E, E / (2 * (1 + nu)), nu, 0.0)

    def nombre_nodo(x: float, y: float) -> str:
        return f"N{indice_x[x]}_{indice_y[y]}"

    for x in xs_ordenados:
        for y in ys_ordenados:
            if any(
                geom["x_min"] - _TOL <= x <= geom["x_max"] + _TOL
                and geom["y_min"] - _TOL <= y <= geom["y_max"] + _TOL
                for geom in geometria.values()
            ):
                modelo.add_node(nombre_nodo(x, y), x, y, 0.0)

    quads_por_pano: dict[str, list[tuple[str, float, float, float, float]]] = {
        p["nombre"]: [] for p in paños
    }
    apoyos_por_pano: dict[str, list[dict]] = {p["nombre"]: [] for p in paños}
    lados_continuos: dict[str, set[str]] = {
        p["nombre"]: set() for p in paños
    }
    for i, paño_a in enumerate(paños):
        for paño_b in paños[i + 1:]:
            borde = borde_compartido(paño_a["elemento"], paño_b["elemento"])
            if borde is None:
                continue
            for paño in (paño_a, paño_b):
                geometria_pano = geometria[paño["nombre"]]
                for lado, (eje, coordenada) in geometria_pano["lados"].items():
                    if (
                        eje == borde["eje"]
                        and abs(coordenada - borde["coordenada_m"]) <= _TOL
                    ):
                        lados_continuos[paño["nombre"]].add(lado)
                        break
    nodos_apoyados = set()
    quad_numero = 0
    for paño in paños:
        nombre = paño["nombre"]
        geom = geometria[nombre]
        elemento = paño["elemento"]
        xs_local = [x for x in xs_ordenados if geom["x_min"] - _TOL <= x <= geom["x_max"] + _TOL]
        ys_local = [y for y in ys_ordenados if geom["y_min"] - _TOL <= y <= geom["y_max"] + _TOL]
        for ix in range(len(xs_local) - 1):
            for iy in range(len(ys_local) - 1):
                x0, x1 = xs_local[ix], xs_local[ix + 1]
                y0, y1 = ys_local[iy], ys_local[iy + 1]
                quad_name = f"Q{quad_numero}"
                quad_numero += 1
                modelo.add_quad(
                    quad_name,
                    nombre_nodo(x0, y0),
                    nombre_nodo(x1, y0),
                    nombre_nodo(x1, y1),
                    nombre_nodo(x0, y1),
                    float(paño["espesor_m"]),
                    "Hormigon_losa",
                )
                quads_por_pano[nombre].append((quad_name, x0, x1, y0, y1))
                modelo.add_quad_surface_pressure(quad_name, -float(paño["D_kNm2"]), case="D")
                if float(paño["L_kNm2"]):
                    modelo.add_quad_surface_pressure(
                        quad_name, -float(paño["L_kNm2"]), case="L"
                    )

        apoyos = elemento.get("apoya_en", {}) or {}
        for lado, (eje, coordenada) in geom["lados"].items():
            apoyo = apoyos.get(lado)
            if not apoyo:
                continue
            inicio = geom[f"{'y' if eje == 'x' else 'x'}_min"]
            fin = geom[f"{'y' if eje == 'x' else 'x'}_max"]
            borde = {
                "lado": lado,
                "eje": eje,
                "coordenada_m": coordenada,
                "inicio_m": inicio,
                "fin_m": fin,
                "apoya_en": apoyo,
            }
            apoyos_por_pano[nombre].append(borde)
            for x in xs_local:
                for y in ys_local:
                    coord = x if eje == "x" else y
                    longitudinal = y if eje == "x" else x
                    if abs(coord - coordenada) <= _TOL and inicio - _TOL <= longitudinal <= fin + _TOL:
                        nodos_apoyados.add(nombre_nodo(x, y))

    nodos_modelo = list(modelo.nodes.values())
    nodo_pin = min(nodos_modelo, key=lambda n: (n.X, n.Y))
    nodo_segundo = min(
        (
            n for n in nodos_modelo
            if abs(n.X - nodo_pin.X) <= _TOL and n.Y > nodo_pin.Y + _TOL
        ),
        key=lambda n: n.Y - nodo_pin.Y,
        default=None,
    )
    if nodo_segundo is not None:
        nodo_lateral = nodo_segundo.name
    else:
        nodo_segundo = min(
            (
                n for n in nodos_modelo
                if abs(n.Y - nodo_pin.Y) <= _TOL and n.X > nodo_pin.X + _TOL
            ),
            key=lambda n: n.X - nodo_pin.X,
        )
        nodo_lateral = nodo_segundo.name
    for nodo in nodos_modelo:
        modelo.def_support(
            nodo.name,
            support_DX=nodo.name in {nodo_pin.name, nodo_lateral},
            support_DY=nodo.name == nodo_pin.name,
            support_DZ=nodo.name in nodos_apoyados,
        )

    for nombre, (fD, fL, fW) in cargas.FACTORES_COMBINACIONES.items():
        if fW == 0:
            modelo.add_load_combo(nombre, {"D": fD, "L": fL})
    modelo.add_load_combo("D_servicio", {"D": 1.0})
    modelo.add_load_combo("L_servicio", {"L": 1.0})
    modelo.analyze_linear(check_stability=True, check_statics=False)

    resultados = {}
    gauss = (-1 / math.sqrt(3), 0.0, 1 / math.sqrt(3))
    quad_por_nombre = modelo.quads
    for paño in paños:
        nombre = paño["nombre"]
        valores_combo = {}
        momentos_apoyos_combo = {}
        for nombre_combo, (fD, fL, fW) in cargas.FACTORES_COMBINACIONES.items():
            if fW != 0:
                continue
            componentes = {k: [] for k in ("Mx", "My", "Mxy", "Qx", "Qy")}
            for quad_name, _, _, _, _ in quads_por_pano[nombre]:
                quad = quad_por_nombre[quad_name]
                for xi in gauss:
                    for eta in gauss:
                        momentos = np.asarray(
                            quad.moment(xi, eta, local=True, combo_name=nombre_combo)
                        ).reshape(-1)
                        cortes = np.asarray(
                            quad.shear(xi, eta, local=True, combo_name=nombre_combo)
                        ).reshape(-1)
                        for clave, valor in zip(("Mx", "My", "Mxy"), momentos[:3]):
                            componentes[clave].append(float(valor))
                        for clave, valor in zip(("Qx", "Qy"), cortes[:2]):
                            componentes[clave].append(float(valor))
            valores_combo[nombre_combo] = {
                clave: {
                    "max": max(valores, default=0.0),
                    "min": min(valores, default=0.0),
                    "max_abs": max((abs(v) for v in valores), default=0.0),
                }
                for clave, valores in componentes.items()
            }
            momentos_apoyos = {}
            for borde in apoyos_por_pano[nombre]:
                if borde["lado"] not in lados_continuos[nombre]:
                    continue
                valores_borde = []
                componente = "Mx" if borde["eje"] == "x" else "My"
                for quad_name, x0, x1, y0, y1 in quads_por_pano[nombre]:
                    if borde["eje"] == "x":
                        if abs(x0 - borde["coordenada_m"]) <= _TOL:
                            xi = -1 / math.sqrt(3)
                        elif abs(x1 - borde["coordenada_m"]) <= _TOL:
                            xi = 1 / math.sqrt(3)
                        else:
                            continue
                        puntos_borde = ((xi, eta) for eta in gauss)
                    else:
                        if abs(y0 - borde["coordenada_m"]) <= _TOL:
                            eta = -1 / math.sqrt(3)
                        elif abs(y1 - borde["coordenada_m"]) <= _TOL:
                            eta = 1 / math.sqrt(3)
                        else:
                            continue
                        puntos_borde = ((xi, eta) for xi in gauss)
                    quad = quad_por_nombre[quad_name]
                    for xi, eta in puntos_borde:
                        momentos = np.asarray(
                            quad.moment(
                                xi, eta, local=True, combo_name=nombre_combo
                            )
                        ).reshape(-1)
                        valores_borde.append(float(momentos[0 if componente == "Mx" else 1]))
                momentos_apoyos[borde["lado"]] = {
                    "componente": componente,
                    "max": max(valores_borde, default=0.0),
                    "min": min(valores_borde, default=0.0),
                    "max_abs": max((abs(v) for v in valores_borde), default=0.0),
                    "metodo": "puntos de integración adyacentes al borde",
                }
            momentos_apoyos_combo[nombre_combo] = momentos_apoyos
        reacciones = []
        for borde in apoyos_por_pano[nombre]:
            nodos_borde = []
            for nodo in modelo.nodes.values():
                coord = nodo.X if borde["eje"] == "x" else nodo.Y
                longitudinal = nodo.Y if borde["eje"] == "x" else nodo.X
                if (
                    abs(coord - borde["coordenada_m"]) <= _TOL
                    and borde["inicio_m"] - _TOL <= longitudinal <= borde["fin_m"] + _TOL
                ):
                    nodos_borde.append(nodo)
            carga_total = {}
            for caso in ("D", "L"):
                total = 0.0
                for nodo in nodos_borde:
                    rxns = nodo.RxnFZ or {}
                    bordes_nodo = [
                        otro_borde
                        for bordes_otro in apoyos_por_pano.values()
                        for otro_borde in bordes_otro
                        if (
                            abs(
                                (nodo.X if otro_borde["eje"] == "x" else nodo.Y)
                                - otro_borde["coordenada_m"]
                            )
                            <= _TOL
                            and otro_borde["inicio_m"] - _TOL
                            <= (nodo.Y if otro_borde["eje"] == "x" else nodo.X)
                            <= otro_borde["fin_m"] + _TOL
                        )
                    ]
                    lineas_nodo = {
                        (
                            otro_borde["eje"],
                            round(otro_borde["coordenada_m"], 7),
                            _clave_apoyo(otro_borde["apoya_en"]),
                        )
                        for otro_borde in bordes_nodo
                    }
                    linea_actual = (
                        borde["eje"],
                        round(borde["coordenada_m"], 7),
                        _clave_apoyo(borde["apoya_en"]),
                    )
                    segmentos_linea = sum(
                        1
                        for otro_borde in bordes_nodo
                        if (
                            otro_borde["eje"],
                            round(otro_borde["coordenada_m"], 7),
                            _clave_apoyo(otro_borde["apoya_en"]),
                        ) == linea_actual
                    )
                    total += float(rxns.get(f"{caso}_servicio", 0.0)) / (
                        max(1, len(lineas_nodo)) * max(1, segmentos_linea)
                    )
                carga_total[caso] = abs(total)
            longitud = borde["fin_m"] - borde["inicio_m"]
            reacciones.append({
                "apoyo": borde["lado"],
                "apoya_en": borde["apoya_en"],
                "reaccion_total_D_kN": carga_total["D"],
                "reaccion_total_L_kN": carga_total["L"],
                "cargas_kN_m": {
                    "D": carga_total["D"] / longitud,
                    "L": carga_total["L"] / longitud,
                },
            })

        D = float(paño["D_kNm2"])
        L = float(paño["L_kNm2"])
        resultados[nombre] = {
            "metodo": "PyNite FEModel3D · malla Quad conjunta con nodos compartidos",
            "modelo": f"maciza continua ({len(paños)} paños)",
            "Ec_MPa": float(ec["Ec_MPa"]),
            "Ec_referencia": ec["articulo"],
            "kx_mod": 1.0,
            "ky_mod": 1.0,
            "malla_aproximada_m": float(malla_m),
            "cantidad_elementos": len(modelo.quads),
            "cargas": {"D_kNm2": D, "L_kNm2": L},
            "combinaciones_kNm2": {
                nombre_combo: fD * D + fL * L
                for nombre_combo, (fD, fL, fW) in cargas.FACTORES_COMBINACIONES.items()
                if fW == 0
            },
            "solicitaciones": valores_combo,
            "momentos_apoyos": momentos_apoyos_combo,
            "reacciones": reacciones,
            "alcance": (
                "Placas macizas coplanares modeladas en conjunto con nodos compartidos "
                "en los bordes declarados continuos; los apoyos conservan desplazamiento "
                "vertical restringido y giro compatible entre paños. No modela la rigidez "
                "de vigas/muros ni dimensiona armaduras, flecha, corte o fisuración."
            ),
        }
    return resultados

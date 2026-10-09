"""
calc.portico — El MOTOR de cálculo: resuelve el pórtico 2D y devuelve SOLO mecánica.

Qué es (y qué no es)
--------------------
Es la capa 2 de la arquitectura (ver ARQUITECTURA.md): el MOTOR. No sabe de
hormigón ni de acero, no dimensiona nada y no pide datos por teclado. Recibe el
pórtico (geometría + cargas) y devuelve las solicitaciones (M, V, N), las
reacciones en las bases y los desplazamientos, POR BARRA Y POR COMBINACIÓN. El
dimensionado (vigas, columnas, bases) se hace después leyendo esto.

Reemplaza a anaStruct (el viejo P01). Se validó que Pynite reproduce a anaStruct
en <0,1 % con los pórticos ya calculados (tools/comparar_motores.py), y esta
versión POR COMBINACIÓN dio idéntica a lo guardado en datos/estructura.json.

Para que los números coincidan se usan las MISMAS rigideces relativas que usaba
anaStruct (EA=15000, EI=5000): las fuerzas internas de un pórtico elástico no
dependen de E, solo del cociente entre la rigidez axial y la de flexión.

Convención de signos (la de Pynite): M en kN·m, V y N en kN. En las columnas
M_inf es el momento en la base y M_sup en la punta; en las vigas M_izq/M_der son
los extremos y M_campo el momento mínimo de la luz.

Uso:
    py -m calc.portico                        qué pórticos hay y qué está resuelto
    py -m calc.portico "Portico 1"            resuelve ese pórtico e imprime el informe
    py -m calc.portico "Portico 1" --guardar  además guarda salidas/solicitaciones
"""

from __future__ import annotations

from copy import deepcopy
import hashlib
import json
import sys
from datetime import datetime
from pathlib import Path

from Pynite import FEModel3D

from . import cargas, materiales, rutas

# ---------------------------------------------------------------------------
# Rigidez (los MISMOS valores con que se validó Pynite contra anaStruct)
# ---------------------------------------------------------------------------
E_MOD = 1.0
G_MOD = 1.0
EA = 15000.0
EI = 5000.0
A_SEC = EA / E_MOD
I_SEC = EI / E_MOD
J_SEC = EI / E_MOD
POISSON = 0.3
PESO_ESPECIFICO_HORMIGON_KN_M3 = 25.0
VERSION_MOTOR = 2


def combinaciones() -> dict[str, tuple[float, float, float, float]]:
    """
    Las combinaciones del pórtico: (fD, fL, fW, fP).

    Las 4 de CIRSOC salen de `cargas.FACTORES_COMBINACIONES` (una sola cuenta).
    'servicio' (D + L, sin viento) es el estado límite de servicio. Las cargas
    puntuales ('P') se aplican enteras en todas las combinaciones, igual que en
    el cálculo viejo (P01).
    """
    tabla: dict[str, tuple[float, float, float, float]] = {"servicio": (1.0, 1.0, 0.0, 1.0)}
    for nombre, (fD, fL, fW) in cargas.FACTORES_COMBINACIONES.items():
        tabla[nombre] = (fD, fL, fW, 1.0)
    return tabla


# ---------------------------------------------------------------------------
# Nudos y nombres
# ---------------------------------------------------------------------------
def _k(x) -> float:
    return round(float(x), 4)


def _nodo(x: float, y: float) -> str:
    return f"N{_k(x)}_{_k(y)}"


def _piso_de_viga(viga_id: str) -> int | None:
    """De 'V0-1' saca 0; de 'V12-2' saca 12. None si el nombre no sigue el patrón."""
    cabeza = viga_id.split("-")[0]
    if cabeza[:1].upper() == "V" and cabeza[1:].isdigit():
        return int(cabeza[1:])
    return None


def _peso_propio_viga_kN_m(viga: dict) -> tuple[float, float, bool]:
    """Devuelve q propio, peralte usado y si el peralte se predimensionó."""
    try:
        ancho_cm = float(viga["b_cm"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("falta un ancho b_cm válido para calcular su peso propio") from exc
    if ancho_cm <= 0:
        raise ValueError("el ancho b_cm debe ser mayor que cero para calcular su peso propio")

    peralte = viga.get("h_cm")
    peralte_predimensionado = peralte is None
    if peralte_predimensionado:
        longitudes = [
            float(tramo.get("longitud_m", 0.0))
            for tramo in viga.get("tramos", [])
        ]
        luz_maxima = max(longitudes, default=0.0)
        if luz_maxima <= 0:
            raise ValueError("no hay una luz válida para predimensionar el peralte")
        peralte = round(luz_maxima * 100 / 12.5)
    try:
        peralte_cm = float(peralte)
    except (TypeError, ValueError) as exc:
        raise ValueError("el peralte h_cm debe ser un número válido") from exc
    if peralte_cm <= 0:
        raise ValueError("el peralte h_cm debe ser mayor que cero")

    area_m2 = ancho_cm * peralte_cm / 10000
    return area_m2 * PESO_ESPECIFICO_HORMIGON_KN_M3, peralte_cm, peralte_predimensionado


def _aplicaciones_losas_portico(
    nombre_portico: str, elementos: dict, estructura_completa: dict
) -> tuple[list[dict], list[str]]:
    """Convierte reacciones de paños por ejes en cargas lineales sobre sus vigas."""
    estructura = estructura_completa.get(nombre_portico, {})
    referencia = estructura.get("referencia_planta", {}) or {}
    direccion_portico = str(referencia.get("direccion", "x")).lower()
    losas_asociadas = []
    for nombre_losa, elemento in elementos.items():
        panel = elemento.get("panel_ejes")
        apoyos = elemento.get("apoya_en")
        if (
            elemento.get("tipo") != "losa"
            or not elemento.get("activo", True)
            or not isinstance(panel, dict)
            or not isinstance(apoyos, dict)
            or direccion_portico
            != ("y" if str(panel.get("direccion_luz", "x")) == "x" else "x")
            or not any(
                (apoyos.get(lado) or {}).get("portico") == nombre_portico
                for lado in ("izq", "der")
            )
        ):
            continue
        losas_asociadas.append((nombre_losa, elemento, panel, apoyos))
    if not losas_asociadas:
        return [], []

    eje_longitudinal_id = str(referencia.get("eje_longitudinal_id", ""))
    ejes = rutas.leer_json(rutas.DATOS / "ejes.json", {}) or {}
    eje_longitudinal = next(
        (
            eje for eje in ejes.get(direccion_portico, [])
            if str(eje.get("id", "")) == eje_longitudinal_id
        ),
        None,
    )
    origen = referencia.get("origen_longitudinal_m")
    if origen is None and eje_longitudinal is not None:
        origen = float(eje_longitudinal["coordenada_m"])
    if origen is None:
        return [], [
            f"{nombre_portico}: no se pudo ubicar el origen longitudinal para cargas de losa."
        ]
    origen = float(origen)
    vigas = {
        str(viga_id): viga
        for viga_id, viga in (estructura.get("vigas", {}) or {}).items()
    }
    aplicaciones: list[dict] = []
    avisos: list[str] = []
    for nombre_losa, elemento, panel, apoyos in losas_asociadas:
        direccion_luz = str(panel.get("direccion_luz", "x"))
        ejes_apoyo = (panel.get("ejes", {}) or {}).get(direccion_luz, {})
        inicio_libre = float(panel.get("coordenada_libre_inicio_m", 0.0))
        fin_libre = float(panel.get("coordenada_libre_fin_m", 0.0))
        luz = float(elemento.get("luz_transversal_m", 0.0) or 0.0)
        if luz <= 0 or fin_libre <= inicio_libre:
            avisos.append(f"{nombre_losa}: geometría de paño inválida; no se aplicó.")
            continue

        for lado, borde in zip(("izq", "der"), ("inicio", "fin")):
            destino = apoyos.get(lado) or {}
            if destino.get("portico") != nombre_portico:
                continue
            viga_id = str(destino.get("viga", ""))
            viga = vigas.get(viga_id)
            if viga is None:
                avisos.append(
                    f"{nombre_losa}: la viga {viga_id or '?'} de apoyo ya no existe en "
                    f"{nombre_portico}; no se aplicó."
                )
                continue
            cota_apoyo = panel.get("cota_apoyo_m")
            if cota_apoyo is not None and (
                viga.get("cota_m") is None
                or abs(float(viga["cota_m"]) - float(cota_apoyo)) > 1e-4
            ):
                avisos.append(
                    f"{nombre_losa}: la viga {viga_id} de {nombre_portico} ya no está "
                    "en el nivel elegido; no se aplicó."
                )
                continue

            eje_guardado = ejes_apoyo.get(borde, {}) or {}
            eje_actual_id = str(referencia.get("eje_id", ""))
            familia_transversal = "y" if direccion_portico == "x" else "x"
            if (
                str(eje_guardado.get("id", "")) != eje_actual_id
                or abs(float(referencia.get("desfase_m", 0.0) or 0.0)) > 1e-8
            ):
                avisos.append(
                    f"{nombre_losa}: el pórtico {nombre_portico} ya no coincide con el "
                    "eje de apoyo guardado; no se aplicó."
                )
                continue
            eje_actual = next(
                (
                    eje for eje in ejes.get(familia_transversal, [])
                    if str(eje.get("id", "")) == eje_actual_id
                ),
                None,
            )
            if eje_actual is None or abs(
                float(eje_actual["coordenada_m"])
                - float(eje_guardado.get("coordenada_m", 0.0))
            ) > 1e-4:
                avisos.append(
                    f"{nombre_losa}: cambió la posición del eje de apoyo de "
                    f"{nombre_portico}; revisá la geometría del paño."
                )
                continue

            intervalos: list[tuple[float, float]] = []
            aplicaciones_lado = []
            for tramo in viga.get("tramos", []):
                tramo_id = str(tramo.get("id", ""))
                longitud = float(tramo.get("longitud_m", 0.0))
                local_inicio = float(tramo.get("x_inicio", 0.0))
                local_fin = float(tramo.get("x_fin", local_inicio + longitud))
                global_inicio = origen + min(local_inicio, local_fin)
                global_fin = origen + max(local_inicio, local_fin)
                inicio = max(inicio_libre, global_inicio)
                fin = min(fin_libre, global_fin)
                if fin <= inicio + 1e-8:
                    continue
                if tramo.get("es_voladizo") and "der" in tramo_id.lower():
                    x_inicio = longitud - (fin - global_inicio)
                    x_fin = longitud - (inicio - global_inicio)
                else:
                    x_inicio = inicio - global_inicio
                    x_fin = fin - global_inicio
                intervalos.append((inicio, fin))
                aplicaciones_lado.append({
                    "id": f"losa:{elemento.get('id', nombre_losa)}:{lado}:{tramo_id}",
                    "carga_id": elemento.get("id"),
                    "portico": nombre_portico,
                    "tramo_id": tramo_id,
                    "x_inicio_m": x_inicio,
                    "x_fin_m": x_fin,
                    "ancho_modo": "manual",
                    "ancho_tributario_m": luz / 2.0,
                    "modo_cargas_previas": "sumar",
                    "activa": True,
                })
            intervalos.sort()
            cubierto_hasta = inicio_libre
            for inicio, fin in intervalos:
                if inicio > cubierto_hasta + 1e-6:
                    break
                cubierto_hasta = max(cubierto_hasta, fin)
            if cubierto_hasta < fin_libre - 1e-6:
                avisos.append(
                    f"{nombre_losa}: la viga {viga_id} de {nombre_portico} no cubre "
                    "todo el borde libre; no se aplicó su reacción."
                )
                continue
            aplicaciones.extend(aplicaciones_lado)
    return aplicaciones, avisos


def _cargas_asignadas(nombre_portico: str) -> tuple[dict[str, list[dict]], list[str]]:
    """Calcula las aplicaciones activas de `datos/cargas.json` para un pórtico."""
    if not nombre_portico:
        return {}, []
    datos = cargas.datos_cargas()
    estructura_completa = rutas.cargar_estructura()
    estructura = estructura_completa.get(nombre_portico, {})
    tramos: dict[str, dict] = {}
    for viga in estructura.get("vigas", {}).values():
        for tramo in viga.get("tramos", []):
            if tramo.get("id"):
                tramos[tramo["id"]] = tramo

    por_id = {
        elemento.get("id"): (nombre, elemento)
        for nombre, elemento in datos.get("elementos", {}).items()
        if elemento.get("id")
    }
    biblioteca = None
    resultado: dict[str, list[dict]] = {}
    aplicaciones_losas, avisos = _aplicaciones_losas_portico(
        nombre_portico, datos.get("elementos", {}), estructura_completa
    )
    aplicaciones_existentes = list(datos.get("aplicaciones", []))
    ids_aplicados_manualmente = {
        aplicacion.get("carga_id")
        for aplicacion in aplicaciones_existentes
        if aplicacion.get("activa", True)
        and aplicacion.get("portico") == nombre_portico
        and aplicacion.get("carga_id")
    }
    ids_losa_manual = {
        aplicacion.get("carga_id")
        for aplicacion in aplicaciones_losas
        if aplicacion.get("carga_id") in ids_aplicados_manualmente
    }
    if ids_losa_manual:
        nombres_por_id = {
            elemento.get("id"): nombre
            for nombre, elemento in datos.get("elementos", {}).items()
        }
        for carga_id in sorted(ids_losa_manual, key=str):
            avisos.append(
                f"{nombres_por_id.get(carga_id, carga_id)}: ya tiene una aplicación "
                f"manual en {nombre_portico}; no se sumó además su reacción automática."
            )
        aplicaciones_losas = [
            aplicacion for aplicacion in aplicaciones_losas
            if aplicacion.get("carga_id") not in ids_losa_manual
        ]
    aplicaciones = aplicaciones_existentes + aplicaciones_losas
    for aplicacion in aplicaciones:
        if not aplicacion.get("activa", True):
            continue
        reacciones = aplicacion.get("reacciones", [])
        if reacciones:
            apoyo = next(
                (r for r in reacciones if r.get("portico") == nombre_portico),
                None,
            )
            if apoyo is None:
                continue
            tramo_id = apoyo.get("tramo_id", "")
            x_puntual = float(apoyo.get("x_m", 0.0))
        else:
            if aplicacion.get("portico") != nombre_portico:
                continue
            apoyo = None
            tramo_id = aplicacion.get("tramo_id", "")

        encontrado = por_id.get(aplicacion.get("carga_id"))
        if not encontrado:
            avisos.append(f"Aplicación {aplicacion.get('id', '?')}: no existe su carga; se ignora")
            continue
        nombre, fuente = encontrado
        if not fuente.get("activo", True):
            continue
        tramo = tramos.get(tramo_id)
        if not tramo:
            avisos.append(f"{aplicacion.get('id', '?')}: no existe el tramo {tramo_id}; se ignora")
            continue
        longitud = float(tramo.get("longitud_m", 0.0))
        if apoyo is not None:
            if longitud <= 0 or x_puntual < 0 or x_puntual > longitud + 1e-8:
                avisos.append(f"{aplicacion.get('id', '?')}: punto fuera del tramo {tramo_id}; se ignora")
                continue
        else:
            x0 = float(aplicacion.get("x_inicio_m", 0.0))
            x1 = float(aplicacion.get("x_fin_m", longitud))
            if longitud <= 0 or x0 < 0 or x1 <= x0 or x1 > longitud + 1e-8:
                avisos.append(f"{aplicacion.get('id', '?')}: intervalo fuera del tramo {tramo_id}; se ignora")
                continue

        elemento = deepcopy(fuente)
        ancho = None
        if elemento.get("tipo") in ("losa", "cubierta"):
            modo = aplicacion.get("ancho_modo", "manual")
            luz = float(elemento.get("luz_transversal_m", 0.0))
            if modo == "media_luz":
                ancho = luz / 2
            elif modo == "luz_completa":
                ancho = luz
            elif modo == "entre_porticos":
                ancho = cargas.ancho_tributario_entre_porticos(
                    estructura_completa, nombre_portico
                )
            else:
                ancho = aplicacion.get("ancho_tributario_m")
            if ancho is None or float(ancho) <= 0:
                avisos.append(f"{aplicacion.get('id', '?')}: falta ancho tributario; se ignora")
                continue
            elemento["ancho_tributario_m"] = float(ancho)
        try:
            if biblioteca is None:
                biblioteca = materiales.cargar()
            items = cargas.items_de_elemento(f"{fuente.get('id')} · {nombre}", elemento, biblioteca)
        except (KeyError, ValueError) as exc:
            avisos.append(f"{aplicacion.get('id', '?')}: no se pudo calcular {nombre}: {exc}")
            continue

        if apoyo is not None:
            influencia = float(apoyo.get("influencia_m", 0.0))
            if influencia <= 0 or elemento.get("tipo") != "muro":
                avisos.append(
                    f"{aplicacion.get('id', '?')}: la reacción transversal requiere muro y ancho positivo; se ignora"
                )
                continue
            x_local = x_puntual
            if tramo.get("es_voladizo") and "der" in str(tramo_id).lower():
                x_local = longitud - x_puntual
            for item in items:
                if item["tipo"] != "D":
                    continue
                resultado.setdefault(tramo_id, []).append({
                    "aplicacion_id": aplicacion.get("id", ""),
                    "carga_id": fuente.get("id", ""),
                    "descripcion": nombre,
                    "tipo": "D",
                    "tipo_aplicacion": "puntual",
                    "valor_kN": float(item["valor"]) * influencia,
                    "x_m": x_local,
                    "influencia_m": influencia,
                    "separacion_m": float(aplicacion.get("separacion_m", influencia * 2)),
                    "modo_cargas_previas": "sumar",
                    "signo": -1.0,
                })
            continue

        # Pynite toma las cargas parciales desde el eje local del miembro. En el
        # voladizo derecho ese eje apunta desde la punta hacia la columna.
        a, b = x0, x1
        if tramo.get("es_voladizo") and "der" in str(tramo_id).lower():
            a, b = longitud - x1, longitud - x0
        for item in items:
            tipo = item["tipo"]
            if tipo not in ("D", "L", "W"):
                continue
            valor = float(item["valor"])
            # D y L actúan hacia abajo. La succión de cubierta actúa hacia arriba.
            signo = 1.0 if tipo == "W" and "succión" in item["descripcion"].lower() else -1.0
            resultado.setdefault(tramo_id, []).append({
                "aplicacion_id": aplicacion.get("id", ""),
                "carga_id": fuente.get("id", ""),
                "descripcion": nombre,
                "tipo": tipo,
                "valor_kN_m": valor,
                "x_inicio_m": a,
                "x_fin_m": b,
                "ancho_tributario_m": float(ancho) if ancho is not None else None,
                "modo_cargas_previas": aplicacion.get("modo_cargas_previas", "reemplazar"),
                "signo": signo,
            })
    return resultado, avisos


# ---------------------------------------------------------------------------
# Armar el modelo (pórtico plano)
# ---------------------------------------------------------------------------
def construir(portico: dict, nombre: str = "") -> tuple[FEModel3D, dict, list[str]]:
    """
    Arma el modelo Pynite del pórtico como PÓRTICO PLANO (libre en el plano,
    fijo fuera de él: DZ, RX y RY bloqueados en todos los nudos). Aplica las
    cargas como casos base D, L, W y P, incluido el peso propio de las vigas,
    y deja cargadas las combinaciones.

    Devuelve (modelo, mapa_de_barras, avisos). Columnas conservan un miembro FE
    por ID; vigas y voladizos exponen listas de segmentos FE bajo el ID lógico
    original, para que las transferencias y cargas interiores no cambien la API
    de solicitaciones.
    """
    aplicaciones, avisos = _cargas_asignadas(nombre)
    m = FEModel3D()
    m.add_material("mat", E_MOD, G_MOD, POISSON, 0.0)
    m.add_section("sec", A_SEC, I_SEC, I_SEC, J_SEC)

    mapa: dict[str, dict] = {
        "columnas": {}, "vigas": {}, "voladizos": {}, "segmentos": {},
        "nudos_columna_inferior": set(), "nudos_viga": set(), "apoyos": [],
    }
    # Solo se informan como aplicadas las cargas de miembros que se pudieron
    # crear en el modelo. Una asignación a un tramo sin apoyos no debe aparecer
    # en el resultado como si hubiera entrado al cálculo.
    mapa["cargas_aplicadas"] = []
    columnas = portico.get("columnas", {})
    datos_carga_proyecto = cargas.datos_cargas() if nombre else {}
    viento_config = datos_carga_proyecto.get("viento", {})
    w_general = 0.0
    if viento_config.get("activo"):
        items_viento = cargas.items_viento_general(datos_carga_proyecto, materiales.cargar())
        w_general = sum(float(item["valor"]) for item in items_viento)
    mapa["viento_general_kN_m"] = w_general

    # --- Columnas (de abajo hacia arriba: i-end = base, j-end = punta) -----
    for cid, col in columnas.items():
        x = _k(col["x"])
        y0 = _k(col.get("nivel", 0.0))
        h = _k(col["altura_m"])
        nb, nt = _nodo(x, y0), _nodo(x, y0 + h)
        if nb not in m.nodes:
            m.add_node(nb, x, y0, 0.0)
        if nt not in m.nodes:
            m.add_node(nt, x, y0 + h, 0.0)
        m.def_support(nb, False, False, True, True, True, False)
        m.def_support(nt, False, False, True, True, True, False)
        m.add_member(f"C-{cid}", nb, nt, "mat", "sec")
        mapa["columnas"][cid] = f"C-{cid}"
        if y0 > 1e-6:
            mapa["nudos_columna_inferior"].add(nb)

    # --- Bases (apoyos) ----------------------------------------------------
    for base_id, base in portico.get("bases", {}).items():
        nodo = _nodo(base["x"], 0.0)
        if nodo not in m.nodes:
            avisos.append(f"{base_id}: no hay columna en x={base['x']}, se ignora")
            continue
        if base.get("tipo") == "articulado":
            m.def_support(nodo, True, True, True, False, False, False)
        else:  # empotramiento (por defecto)
            m.def_support(nodo, True, True, True, True, True, True)
        mapa["apoyos"].append({
            "x": float(base["x"]),
            "y": 0.0,
            "tipo": "articulado" if base.get("tipo") == "articulado" else "empotramiento",
        })

    # --- Vigas (por piso), voladizos y sus cargas --------------------------
    for viga_id, viga in portico.get("vigas", {}).items():
        try:
            peso_propio_kN_m, peralte_cm, peralte_predimensionado = (
                _peso_propio_viga_kN_m(viga)
            )
        except ValueError as exc:
            peso_propio_kN_m = None
            avisos.append(f"{viga_id}: {exc}; no se incorpora el peso propio.")
        piso = _piso_de_viga(viga_id)
        cols = sorted(
            (c for cid, c in columnas.items() if piso is not None and cid.startswith(f"C{piso}-")),
            key=lambda c: c["x"],
        )
        if not cols:
            avisos.append(f"{viga_id}: no encuentro sus columnas (piso {piso}); se ignora")
            continue
        y_viga = _k(cols[0].get("nivel", 0.0) + cols[0]["altura_m"])

        for i, tramo in enumerate(viga.get("tramos", [])):
            cargas_tramo = tramo.get("cargas", {})
            tramo_id = str(tramo["id"])
            miembro_ids = []
            if tramo.get("es_voladizo"):
                if "izq" in tramo["id"]:
                    camino_inicio = _k(cols[0]["x"] - tramo["longitud_m"])
                    camino_fin = _k(cols[0]["x"])
                else:
                    camino_inicio = _k(cols[-1]["x"] + tramo["longitud_m"])
                    camino_fin = _k(cols[-1]["x"])
            else:
                if i >= len(cols) - 1:
                    avisos.append(f"{tramo['id']}: faltan columnas para el tramo; se ignora")
                    continue
                camino_inicio = _k(tramo.get("x_inicio", cols[i]["x"]))
                camino_fin = _k(tramo.get("x_fin", cols[i + 1]["x"]))
            if abs(camino_fin - camino_inicio) <= 1e-8:
                avisos.append(f"{tramo_id}: longitud nula; se ignora")
                continue
            signo_x = 1.0 if camino_fin > camino_inicio else -1.0
            longitud_camino = abs(camino_fin - camino_inicio)
            puntos_nativos = tramo.get("cargas_puntuales", [])
            puntos_aplicados = [
                carga for carga in aplicaciones.get(tramo_id, [])
                if carga.get("tipo_aplicacion") == "puntual"
            ]
            estaciones_carga = []
            for cp in puntos_nativos:
                coordenadas = cp.get("coordenadas", {})
                x_cp = _k(coordenadas.get("x", camino_inicio))
                y_cp = _k(coordenadas.get("y", y_viga))
                estacion = (x_cp - camino_inicio) * signo_x
                if estacion < -1e-6 or estacion > longitud_camino + 1e-6 or abs(y_cp - y_viga) > 1e-4:
                    avisos.append(f"{tramo_id}: carga puntual fuera de la viga; se ignora")
                    continue
                estaciones_carga.append((
                    max(0.0, min(longitud_camino, estacion)), x_cp, cp, "nativa"
                ))
            for carga in puntos_aplicados:
                estacion = float(carga["x_m"])
                if estacion < -1e-6 or estacion > longitud_camino + 1e-6:
                    avisos.append(f"{tramo_id}: carga puntual aplicada fuera de la viga; se ignora")
                    continue
                estacion = max(0.0, min(longitud_camino, estacion))
                x_cp = _k(camino_inicio + signo_x * estacion)
                estaciones_carga.append((estacion, x_cp, carga, "aplicacion"))

            estaciones_nodos = [0.0, longitud_camino]
            estaciones_nodos.extend(
                (float(c["x"]) - camino_inicio) * signo_x
                for c in columnas.values()
                if abs(_k(c.get("nivel", 0.0)) - y_viga) <= 1e-4
                and -1e-6 < (float(c["x"]) - camino_inicio) * signo_x < longitud_camino - 1e-6
            )
            estaciones_nodos.extend(
                carga_nodal[0] for carga_nodal in estaciones_carga
                if 1e-6 < carga_nodal[0] < longitud_camino - 1e-6
            )
            estaciones_nodos = sorted({
                round(max(0.0, min(longitud_camino, estacion)), 6)
                for estacion in estaciones_nodos
            })
            nodos_camino = []
            for estacion in estaciones_nodos:
                x_nodo = _k(camino_inicio + signo_x * estacion)
                n_nodo = _nodo(x_nodo, y_viga)
                if n_nodo not in m.nodes:
                    m.add_node(n_nodo, x_nodo, y_viga, 0.0)
                    m.def_support(n_nodo, False, False, True, True, True, False)
                nodos_camino.append((estacion, n_nodo))
                mapa["nudos_viga"].add(n_nodo)

            # Cada apoyo superior y cada carga puntual crea un nudo compartido;
            # se divide la barra FE, pero se conserva el ID lógico del tramo.
            for indice_segmento, (
                (estacion_i, nodo_i), (estacion_j, nodo_j)
            ) in enumerate(zip(nodos_camino, nodos_camino[1:]), 1):
                miembro_id = f"V-{tramo_id}__s{indice_segmento}"
                m.add_member(miembro_id, nodo_i, nodo_j, "mat", "sec")
                miembro_ids.append(miembro_id)
                mapa["segmentos"].setdefault(tramo_id, []).append({
                    "miembro_id": miembro_id,
                    "estacion_inicio_m": estacion_i,
                    "estacion_fin_m": estacion_j,
                })
            if not miembro_ids:
                avisos.append(f"{tramo_id}: no se pudo generar ningún segmento FE; se ignora")
                continue
            if tramo.get("es_voladizo"):
                mapa["voladizos"][tramo_id] = miembro_ids
            else:
                mapa["vigas"][tramo_id] = miembro_ids

            cargas_nuevas = aplicaciones.get(tramo_id, [])
            mapa["cargas_aplicadas"].extend(
                dict(carga, tramo_id=tramo_id) for carga in cargas_nuevas
            )
            reemplaza_anteriores = any(
                c.get("modo_cargas_previas") == "reemplazar" for c in cargas_nuevas
            )
            # Las aplicaciones nuevas reemplazan las cargas distribuidas de P00
            # cuando así se indicó; las cargas puntuales se conservan.
            wD = float(cargas_tramo.get("D_total", 0.0))
            wL = float(cargas_tramo.get("L_total", 0.0))
            for indice_segmento, segmento in enumerate(mapa["segmentos"][tramo_id]):
                miembro_id = segmento["miembro_id"]
                estacion_i = segmento["estacion_inicio_m"]
                estacion_j = segmento["estacion_fin_m"]
                if not reemplaza_anteriores and abs(wD) > 1e-9:
                    m.add_member_dist_load(miembro_id, "FY", -wD, -wD, case="D")
                if not reemplaza_anteriores and abs(wL) > 1e-9:
                    m.add_member_dist_load(miembro_id, "FY", -wL, -wL, case="L")
                if peso_propio_kN_m is not None:
                    m.add_member_dist_load(
                        miembro_id, "FY", -peso_propio_kN_m, -peso_propio_kN_m,
                        case="D",
                    )
                for carga_aplicada in cargas_nuevas:
                    if carga_aplicada.get("tipo_aplicacion") == "puntual":
                        continue
                    inicio = max(estacion_i, float(carga_aplicada["x_inicio_m"]))
                    fin = min(estacion_j, float(carga_aplicada["x_fin_m"]))
                    if fin - inicio <= 1e-8:
                        continue
                    x_local_inicio = inicio - estacion_i
                    x_local_fin = fin - estacion_i
                    m.add_member_dist_load(
                        miembro_id, "FY",
                        carga_aplicada["signo"] * carga_aplicada["valor_kN_m"],
                        carga_aplicada["signo"] * carga_aplicada["valor_kN_m"],
                        x1=x_local_inicio, x2=x_local_fin,
                        case=carga_aplicada["tipo"],
                    )
            if peso_propio_kN_m is not None:
                mapa["cargas_aplicadas"].append({
                    "aplicacion_id": f"peso-propio:{tramo_id}",
                    "carga_id": f"peso-propio:{viga_id}",
                    "descripcion": (
                        f"Peso propio {viga_id} · "
                        f"{float(viga['b_cm']):g}×{peralte_cm:g} cm"
                        + (" · h predimensionado L/12,5" if peralte_predimensionado else "")
                    ),
                    "tipo": "D",
                    "valor_kN_m": peso_propio_kN_m,
                    "x_inicio_m": 0.0,
                    "x_fin_m": longitud_camino,
                    "ancho_tributario_m": None,
                    "modo_cargas_previas": "sumar",
                    "signo": -1.0,
                    "tramo_id": tramo_id,
                })
            for carga_nodal in estaciones_carga:
                x_cp, carga_puntual, origen = carga_nodal[1], carga_nodal[2], carga_nodal[3]
                n_cp = _nodo(x_cp, y_viga)
                if origen == "nativa":
                    valor = -float(carga_puntual["valor_kN"])
                    caso = "P"
                else:
                    valor = (
                        float(carga_puntual["signo"])
                        * float(carga_puntual["valor_kN"])
                    )
                    caso = str(carga_puntual["tipo"])
                m.add_node_load(n_cp, "FY", valor, case=caso)

    nudos_viga = mapa["nudos_viga"]
    nudos_sin_apoyo = mapa["nudos_columna_inferior"] - nudos_viga
    if nudos_sin_apoyo:
        ids = ", ".join(sorted(nudos_sin_apoyo))
        raise ValueError(
            "Hay columnas de niveles superiores sin viga o voladizo inferior en su base: "
            f"{ids}. No se resuelve un modelo con transferencias desconectadas."
        )

    # --- Viento (caso W): fuerza horizontal en la punta de cada columna ----
    for viga_id, viga in portico.get("vigas", {}).items():
        piso = _piso_de_viga(viga_id)
        if piso is None or not viga.get("tramos"):
            continue
        cols = [c for cid, c in columnas.items() if cid.startswith(f"C{piso}-")]
        if not cols:
            continue
        carga_legacy = viga["tramos"][0].get("cargas", {})
        if "viento" in datos_carga_proyecto:
            # La configuración del proyecto es la única fuente del viento.
            # Evita sumar el viento general nuevo con el total copiado por P00.
            w_total = w_general
        else:
            w_total = float(carga_legacy.get("W_total", 0.0))
        total = w_total * float(cols[0]["altura_m"])
        if abs(total) < 1e-9:
            continue
        por_col = total / len(cols)  # reparto igual entre las columnas del piso
        for cid, col in columnas.items():
            if not cid.startswith(f"C{piso}-"):
                continue
            n_punta = _nodo(col["x"], col.get("nivel", 0.0) + col["altura_m"])
            if n_punta in m.nodes:
                m.add_node_load(n_punta, "FX", por_col, case="W")

    # --- Combinaciones -----------------------------------------------------
    for nombre, (fD, fL, fW, fP) in combinaciones().items():
        factores: dict[str, float] = {}
        if fD:
            factores["D"] = fD
        if fL:
            factores["L"] = fL
        if fW:
            factores["W"] = fW
        if fP:
            factores["P"] = fP
        m.add_load_combo(nombre, factores)

    return m, mapa, avisos


# ---------------------------------------------------------------------------
# Resolver y extraer las solicitaciones
# ---------------------------------------------------------------------------
def _envolvente_vacia(mapa: dict) -> dict:
    return {
        "columnas": {cid: {} for cid in mapa["columnas"]},
        "vigas": {t: {} for t in mapa["vigas"]},
        "voladizos": {t: {} for t in mapa["voladizos"]},
        "bases": {},
    }


def _acumular(bolsa: dict, clave: str, valores: dict) -> None:
    """Guarda, por barra, el máximo EN MÓDULO de cada magnitud (la envolvente)."""
    actual = bolsa.setdefault(clave, {})
    for campo, valor in valores.items():
        v = abs(float(valor))
        if campo not in actual or v > actual[campo]:
            actual[campo] = v


def _resultados_segmentos(
    modelo: FEModel3D, miembro_ids: list[str], combinacion: str,
) -> dict:
    miembros = [modelo.members[miembro_id] for miembro_id in miembro_ids]
    primero, ultimo = miembros[0], miembros[-1]
    return {
        "M_izq": float(primero.moment("Mz", 0.0, combinacion)),
        "M_der": float(ultimo.moment("Mz", ultimo.L(), combinacion)),
        "M_campo": min(float(miembro.min_moment("Mz", combinacion)) for miembro in miembros),
        "V_izq": float(primero.shear("Fy", 0.0, combinacion)),
        "V_der": float(ultimo.shear("Fy", ultimo.L(), combinacion)),
    }


def _diagrama_momento(modelo: FEModel3D, mapa: dict, combinacion: str) -> dict:
    """Muestrea momento y corte de todas las barras para una combinación."""
    barras = []
    grupos = (
        ("columna", mapa["columnas"].items()),
        ("viga", mapa["vigas"].items()),
        ("voladizo", mapa["voladizos"].items()),
    )
    maximo = 0.0
    maximo_corte = 0.0
    for tipo, elementos in grupos:
        for identificador, valor in elementos:
            miembro_ids = [valor] if tipo == "columna" else valor
            for indice, miembro_id in enumerate(miembro_ids, start=1):
                miembro = modelo.members[miembro_id]
                estaciones, momentos = miembro.moment_array("Mz", 31, combinacion)
                estaciones_corte, cortes = miembro.shear_array("Fy", 31, combinacion)
                estaciones_m = [float(x) for x in estaciones]
                momentos_kNm = [float(m) for m in momentos]
                estaciones_corte_m = [float(x) for x in estaciones_corte]
                cortes_kN = [float(v) for v in cortes]
                if (
                    len(estaciones_m) != len(momentos_kNm)
                    or not estaciones_m
                    or len(estaciones_corte_m) != len(estaciones_m)
                    or len(cortes_kN) != len(estaciones_m)
                ):
                    raise ValueError(
                        f"Pynite devolvió diagramas inválidos para {identificador}."
                    )
                maximo = max(maximo, *(abs(m) for m in momentos_kNm))
                maximo_corte = max(maximo_corte, *(abs(v) for v in cortes_kN))
                etiqueta = str(identificador)
                if len(miembro_ids) > 1:
                    etiqueta = f"{etiqueta} · segmento {indice}"
                barras.append({
                    "id": etiqueta,
                    "tipo": tipo,
                    "inicio": {
                        "x": float(miembro.i_node.X),
                        "y": float(miembro.i_node.Y),
                    },
                    "fin": {
                        "x": float(miembro.j_node.X),
                        "y": float(miembro.j_node.Y),
                    },
                    "x_m": estaciones_m,
                    "M_kNm": momentos_kNm,
                    "V_kN": cortes_kN,
                })
    return {
        "combinacion": combinacion,
        "max_abs_kNm": maximo,
        "max_abs_kN": maximo_corte,
        "apoyos": mapa.get("apoyos", []),
        "barras": barras,
    }


def calcular(portico: dict, nombre: str = "") -> dict:
    """
    Resuelve el pórtico y devuelve las solicitaciones, reacciones y
    desplazamientos por barra y por combinación, más una ENVOLVENTE con los
    máximos en módulo (la que usarán los dimensionadores). No guarda nada: eso
    lo hace `guardar`.
    """
    m, mapa, avisos = construir(portico, nombre)
    m.analyze()

    combos = list(combinaciones())
    solicitaciones: dict[str, dict] = {}
    envolvente = _envolvente_vacia(mapa)
    diagrama_gobernante = None
    maximo_momento_global = -1.0

    for combo in combos:
        bloque: dict[str, dict] = {
            "columnas": {}, "vigas": {}, "voladizos": {}, "bases": {}, "desplazamientos": {},
        }

        for cid, miembro_id in mapa["columnas"].items():
            mem = m.members[miembro_id]
            col = portico["columnas"][cid]
            h = float(col.get("altura_m", 1.0)) or 1.0
            n_punta = _nodo(col["x"], col.get("nivel", 0.0) + col["altura_m"])
            datos = {
                "M_inf": float(mem.moment("Mz", 0.0, combo)),
                "M_sup": float(mem.moment("Mz", mem.L(), combo)),
                "N": float(mem.axial(0.0, combo)),
                "deriva": float(m.nodes[n_punta].DX[combo]) / h,
            }
            bloque["columnas"][cid] = datos
            _acumular(envolvente["columnas"], cid, datos)

        for grupo in ("vigas", "voladizos"):
            for clave, miembro_ids in mapa[grupo].items():
                datos = _resultados_segmentos(m, miembro_ids, combo)
                bloque[grupo][clave] = datos
                _acumular(envolvente[grupo], clave, datos)

        for base_id, base in portico.get("bases", {}).items():
            nodo = _nodo(base["x"], 0.0)
            if nodo not in m.nodes:
                continue
            n = m.nodes[nodo]
            datos = {
                "Fx": float(n.RxnFX.get(combo, 0.0)),
                "Fy": float(n.RxnFY.get(combo, 0.0)),
                "Mz": float(n.RxnMZ.get(combo, 0.0)),
            }
            bloque["bases"][base_id] = datos
            _acumular(envolvente["bases"], base_id, datos)

        for nodo_id, nodo in m.nodes.items():
            bloque["desplazamientos"][nodo_id] = {
                "DX": float(nodo.DX[combo]),
                "DY": float(nodo.DY[combo]),
                "RZ": float(nodo.RZ[combo]),
            }

        solicitaciones[combo] = bloque
        diagrama = _diagrama_momento(m, mapa, combo)
        if diagrama["max_abs_kNm"] > maximo_momento_global:
            diagrama_gobernante = diagrama
            maximo_momento_global = diagrama["max_abs_kNm"]

    return {
        "portico": nombre,
        "generado": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "rigidez": {"EA": EA, "EI": EI},
        "combinaciones": combos,
        "avisos": avisos,
        "cargas_aplicadas": mapa.get("cargas_aplicadas", []),
        "viento_general_kN_m": mapa.get("viento_general_kN_m", 0.0),
        "solicitaciones": solicitaciones,
        "envolvente": envolvente,
        "diagrama_momentos": diagrama_gobernante,
    }


# ---------------------------------------------------------------------------
# Informe y guardado
# ---------------------------------------------------------------------------
def informe(resultado: dict) -> str:
    lineas = [
        f"SOLICITACIONES DEL PÓRTICO: {resultado['portico']}",
        "=" * 74,
        f"Motor Pynite (pórtico plano) · rigidez EA={EA:g} EI={EI:g} · {resultado['generado']}",
        "Unidades: M [kN·m], N y V [kN], deriva [m/m]. Signos según Pynite.",
    ]
    for aviso in resultado["avisos"]:
        lineas.append(f"[AVISO] {aviso}")

    aplicaciones = resultado.get("cargas_aplicadas", [])
    if aplicaciones:
        lineas.extend(["", "CARGAS APLICADAS DESDE EL PROYECTO:"])
        for carga in aplicaciones:
            linea = (
                f"{carga['carga_id']} · {carga['descripcion']} → {carga['tramo_id']} "
                f"x={carga['x_inicio_m']:.2f}–{carga['x_fin_m']:.2f} m: "
                f"{carga['valor_kN_m']:.2f} kN/m ({carga['tipo']}); "
                f"P00: {carga['modo_cargas_previas']}"
            )
            if carga.get("ancho_tributario_m") is not None:
                linea += f"; b trib.={carga['ancho_tributario_m']:.2f} m"
            lineas.append(linea)
        lineas.append("Las cargas puntuales guardadas previamente en P00 se conservaron.")

    tabla = resultado["solicitaciones"]
    for combo in resultado["combinaciones"]:
        bloque = tabla[combo]
        lineas.append(f"\n=== {combo} ===")

        if bloque["columnas"]:
            lineas.append("  Columnas:   barra        M_base     M_punta          N      deriva")
            for cid, d in bloque["columnas"].items():
                lineas.append(
                    f"    {cid:<10} {d['M_inf']:>10.2f} {d['M_sup']:>11.2f} "
                    f"{d['N']:>10.2f} {d['deriva']:>11.5f}"
                )

        for titulo, grupo in (("Vigas", "vigas"), ("Voladizos", "voladizos")):
            if not bloque[grupo]:
                continue
            lineas.append(f"  {titulo}:     tramo        M_izq      M_der    M_campo      V_izq      V_der")
            for clave, d in bloque[grupo].items():
                lineas.append(
                    f"    {clave:<12} {d['M_izq']:>8.2f} {d['M_der']:>10.2f} {d['M_campo']:>10.2f} "
                    f"{d['V_izq']:>10.2f} {d['V_der']:>10.2f}"
                )

        if bloque["bases"]:
            lineas.append("  Bases:      base             Fx         Fy         Mz")
            for base_id, d in bloque["bases"].items():
                lineas.append(f"    {base_id:<10} {d['Fx']:>10.2f} {d['Fy']:>10.2f} {d['Mz']:>10.2f}")

    lineas.append("\nENVOLVENTE (máximos en módulo, para el dimensionado)")
    lineas.append("-" * 74)
    for cid, d in resultado["envolvente"]["columnas"].items():
        if d:
            lineas.append(f"    {cid:<10} |M| máx = {max(d.get('M_inf', 0.0), d.get('M_sup', 0.0)):8.2f} kNm"
                          f"   |N| máx = {d.get('N', 0.0):8.2f} kN")
    for grupo in ("vigas", "voladizos"):
        for clave, d in resultado["envolvente"][grupo].items():
            if d:
                lineas.append(f"    {clave:<12} |M| máx = {d.get('M_campo', 0.0):8.2f} kNm"
                              f"   |V| máx = {d.get('V_izq', 0.0):8.2f} kN")
    return "\n".join(lineas)


def guardar(resultado: dict) -> Path:
    """Guarda las solicitaciones en salidas/solicitaciones/<portico>.json."""
    rutas.SAL_SOLICITACIONES.mkdir(parents=True, exist_ok=True)
    nombre = rutas.nombre_seguro(resultado["portico"]) or "portico"
    archivo = rutas.SAL_SOLICITACIONES / f"{nombre}.json"
    resultado["_firma_entradas_motor"] = huella_entrada(resultado["portico"])
    return rutas.guardar_json(archivo, resultado)


def huella_entrada(nombre_portico: str) -> str:
    """Firma las entradas mecánicas y la versión del modelo."""
    estructura = rutas.cargar_estructura()
    if nombre_portico not in estructura:
        raise ValueError(f"No existe {nombre_portico} en datos/estructura.json")

    return _huella_entradas(
        estructura[nombre_portico],
        rutas.leer_json(rutas.CARGAS, {}) or {},
        rutas.leer_json(rutas.MATERIALES, {}) or {},
    )


def _huella_entradas(datos_portico: dict, cargas: dict, materiales: dict) -> str:
    datos_portico = deepcopy(datos_portico)
    for viga in datos_portico.get("vigas", {}).values():
        # b y h determinan el peso propio; los demás datos son de dimensionado.
        for clave in ("recubrimiento_cm", "fc_MPa", "fy_MPa"):
            viga.pop(clave, None)

    entradas = {
        "version_motor": VERSION_MOTOR,
        "portico": datos_portico,
        "cargas": cargas,
        "materiales": materiales,
    }
    contenido = json.dumps(
        entradas, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")
    return hashlib.sha256(contenido).hexdigest()


# ---------------------------------------------------------------------------
# Listado y consola
# ---------------------------------------------------------------------------
def estado_texto() -> str:
    estructura = rutas.cargar_estructura()
    lineas = ["PÓRTICOS DISPONIBLES (datos/estructura.json)", "=" * 74]
    if not estructura:
        lineas.append("(no hay pórticos cargados)")
    for nombre, p in estructura.items():
        n_col = len(p.get("columnas", {}))
        n_tra = sum(len(v.get("tramos", [])) for v in p.get("vigas", {}).values())
        lineas.append(f" - {nombre:<26} {n_col} columnas · {n_tra} tramos")

    archivos = rutas.listar(rutas.SAL_SOLICITACIONES, "*.json")
    lineas.append("")
    if archivos:
        lineas.append(f"Solicitaciones ya guardadas ({len(archivos)}):")
        for a in archivos:
            lineas.append(f"   {a.name}")
    else:
        lineas.append("(todavía no hay solicitaciones guardadas)")
    lineas.append("")
    lineas.append('Resolver uno:   py -m calc.portico "<nombre o número>" [--guardar]')
    return "\n".join(lineas)


def _main(argv: list[str]) -> int:
    opciones = [a for a in argv if a.startswith("-")]
    nombres = [a for a in argv if not a.startswith("-")]

    if not nombres:
        print(estado_texto())
        return 0

    estructura = rutas.cargar_estructura()
    codigo = 0
    for referencia in nombres:
        try:
            nombre = rutas.resolver_portico(referencia, estructura)
        except KeyError as exc:
            print(f"{exc} Ver: py -m calc.portico")
            codigo = 1
            continue
        resultado = calcular(estructura[nombre], nombre)
        print(informe(resultado))
        if "--guardar" in opciones:
            print(f"\nGuardado en: {guardar(resultado)}")
    return codigo


if __name__ == "__main__":
    raise SystemExit(_main(sys.argv[1:]))

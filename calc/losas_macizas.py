"""Análisis elástico inicial de losas macizas y casetonadas.

Representa paños rectangulares cargados uniformemente, con dos o cuatro bordes
simplemente apoyados. Calcula solicitaciones y reacciones; no dimensiona la
armadura ni representa continuidad entre paños.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path

from . import cargas, materiales, rutas
from .formato import numero_local
from .losas_pynite import analizar_paño_macizo
from .losas_continuidad import analizar_grupo, grupo_continuo
from .losas_seccion import (
    componentes_con_hormigon_casetonada,
    espesor_equivalente_hormigon_casetonada,
)


@dataclass
class ResultadoLosaMaciza:
    nombre: str
    carga_id: str
    metodo: str
    criterio_espesor: str
    luz_m: float
    ancho_m: float
    espesor_m: float
    peso_propio_kNm2: float
    componentes: list[dict]
    D2_kNm2: float
    uso_sobrecarga: str
    D_kNm2: float
    L_kNm2: float
    q_servicio_kNm2: float
    combinaciones_kNm2: dict[str, float]
    combo_critica: str
    momento_ultimo_kNm_m: float
    corte_ultimo_kN_m: float
    reacciones: list[dict]
    tipologia: str = "maciza"
    seccion_casetonada: dict | None = None
    analisis_pynite: dict | None = None
    memoria: str = ""


def validar_seccion_casetonada(datos: dict) -> dict:
    """Valida la geometría y los factores de una placa casetonada ortótropa."""
    if not isinstance(datos, dict):
        raise ValueError("Definí la geometría y los factores de rigidez de la casetonada.")
    claves = (
        "altura_total_m",
        "espesor_capa_m",
        "ancho_nervio_x_m",
        "ancho_nervio_y_m",
        "separacion_libre_x_m",
        "separacion_libre_y_m",
        "kx_mod",
        "ky_mod",
    )
    valores = {}
    for clave in claves:
        try:
            valor = float(datos.get(clave, 0.0))
        except (TypeError, ValueError):
            raise ValueError(f"El dato de casetonada '{clave}' debe ser numérico.") from None
        if not math.isfinite(valor):
            raise ValueError(f"El dato de casetonada '{clave}' debe ser finito.")
        valores[clave] = valor

    h = valores["altura_total_m"]
    capa = valores["espesor_capa_m"]
    bx = valores["ancho_nervio_x_m"]
    by = valores["ancho_nervio_y_m"]
    sx = valores["separacion_libre_x_m"]
    sy = valores["separacion_libre_y_m"]
    kx = valores["kx_mod"]
    ky = valores["ky_mod"]
    if h <= 0 or capa <= 0 or capa >= h:
        raise ValueError("La casetonada necesita canto total y capa positivos, con capa menor que el canto.")
    if min(bx, by, sx, sy) <= 0:
        raise ValueError("Los anchos y las separaciones libres de los nervios deben ser positivos.")
    if not all(0 < factor <= 1 for factor in (kx, ky)):
        raise ValueError("Los factores de rigidez kx y ky deben ser mayores que 0 y no mayores que 1.")

    aligeramiento_estructural = bool(datos.get("aligeramiento_estructural", False))
    espesor_minimo_capa = max(
        sx / 12,
        sy / 12,
        0.04 if aligeramiento_estructural else 0.05,
    )
    incumplimientos = []
    if min(bx, by) < 0.10:
        incumplimientos.append("ancho mínimo de nervio de 100 mm")
    if h > 3.5 * min(bx, by):
        incumplimientos.append("canto total no mayor que 3,5 veces el ancho mínimo del nervio")
    if max(sx, sy) > 0.75:
        incumplimientos.append("separación libre entre nervios no mayor que 750 mm")
    if capa < espesor_minimo_capa:
        incumplimientos.append(
            f"capa superior mínima de {numero_local(espesor_minimo_capa * 100, 1)} cm"
        )
    if incumplimientos:
        articulos = "8.8.1.2–8.8.1.4 y 8.8.2.1.1" if aligeramiento_estructural else "8.8.1.2–8.8.1.4 y 8.8.3.1"
        raise ValueError(
            "La sección no cumple los límites geométricos del Proyecto CIRSOC 201-24 "
            f"(arts. {articulos}): " + "; ".join(incumplimientos) + ". "
            "El art. 8.8.1.8 requiere tratarla como losa y vigas."
        )

    paso_x = sx + bx
    paso_y = sy + by
    fraccion_nervios = bx / paso_x + by / paso_y - (bx * by) / (paso_x * paso_y)
    articulo_aligeramiento = "8.8.2.1.1" if aligeramiento_estructural else "8.8.3.1"
    return {
        **valores,
        "aligeramiento_estructural": aligeramiento_estructural,
        "espesor_minimo_capa_m": espesor_minimo_capa,
        "paso_x_m": paso_x,
        "paso_y_m": paso_y,
        "fraccion_nervios": fraccion_nervios,
        "espesor_equivalente_hormigon_m": espesor_equivalente_hormigon_casetonada(datos),
        "referencia_cirsoc": (
            "Proyecto CIRSOC 201-24, arts. "
            f"8.8.1.1–8.8.1.4 y {articulo_aligeramiento}"
        ),
    }


def _espesor_hormigon(elemento: dict) -> float:
    """Toma el espesor de la capa de hormigón armado de la composición."""
    if elemento.get("tipologia") == "casetonada":
        seccion = elemento.get("casetonada", {}) or {}
        return float(seccion.get("altura_total_m", 0.0) or 0.0)
    espesores = [
        float(c.get("espesor_m", 0.0))
        for c in elemento.get("componentes", [])
        if (c.get("grupo"), c.get("clave")) in {
            ("Forjados", "Losa_maciza"), ("Hormigon", "armado")
        }
    ]
    espesor = sum(espesores)
    if espesor <= 0:
        raise ValueError(
            "Para analizar la losa maciza, definí el forjado base con su espesor "
            "estructural en la composición."
        )
    return espesor


def calcular_desde_elemento_carga(
    nombre: str,
    elemento: dict,
    biblioteca: dict | None = None,
    analisis_override: dict | None = None,
) -> ResultadoLosaMaciza:
    """Calcula solicitaciones y reacciones de losas macizas o casetonadas."""
    tipologia = elemento.get("tipologia")
    if elemento.get("tipo") != "losa" or tipologia not in (
        "maciza", "maciza_cruzada", "casetonada"
    ):
        raise ValueError("Este análisis corresponde a una losa maciza o casetonada.")

    luz = float(elemento.get("luz_transversal_m", 0.0) or 0.0)
    ancho = float(elemento.get("ancho_losa_m", 0.0) or 0.0)
    if not math.isfinite(luz) or not math.isfinite(ancho) or luz <= 0 or ancho <= 0:
        raise ValueError("Completá la luz entre apoyos y el ancho del paño.")

    seccion_casetonada = (
        validar_seccion_casetonada(elemento.get("casetonada", {}) or {})
        if tipologia == "casetonada"
        else None
    )
    espesor = _espesor_hormigon(elemento)
    componentes_elemento = list(elemento.get("componentes", []))
    componente_estructural = None
    if seccion_casetonada:
        componentes_elemento = componentes_con_hormigon_casetonada(elemento)
        componente_estructural = componentes_elemento[-1]
    else:
        componente_estructural = next((
            c for c in componentes_elemento
            if (c.get("grupo"), c.get("clave")) in {
                ("Forjados", "Losa_maciza"), ("Hormigon", "armado")
            }
        ), None)
    elemento_calculo = {**elemento, "componentes": componentes_elemento}
    bib = biblioteca if biblioteca is not None else materiales.cargar()
    valores = cargas.valores_superficiales(elemento_calculo, bib)
    D = float(valores["D_kNm2"])
    L = float(valores["L_kNm2"])
    if componente_estructural is None:
        raise ValueError("Definí el espesor estructural de hormigón de la losa.")
    peso_propio = cargas.valor_componente_superficial(componente_estructural, bib)[1]
    componentes = []
    D2 = 0.0
    for componente, (descripcion, valor) in zip(componentes_elemento, valores["componentes"]):
        estructural = (componente.get("grupo"), componente.get("clave")) in {
            ("Forjados", "Losa_maciza"), ("Forjados", "Losa_casetonada"),
            ("Hormigon", "armado")
        }
        componentes.append({"descripcion": descripcion, "valor_kNm2": valor, "estructural": estructural})
        if not estructural:
            D2 += valor
    uso_sobrecarga = next((
        str(item.get("uso", elemento.get("sobrecarga", "")))
        for item in bib.get("Sobrecargas", [])
        if item.get("clave") == elemento.get("sobrecarga")
    ), str(elemento.get("sobrecarga", "")))
    analisis = analisis_override or analizar_paño_macizo(
        elemento_calculo, valores, espesor
    )
    combinaciones = analisis["combinaciones_kNm2"]
    combo = max(
        analisis["solicitaciones"],
        key=lambda nombre_combo: max(
            analisis["solicitaciones"][nombre_combo]["Mx"]["max_abs"],
            analisis["solicitaciones"][nombre_combo]["My"]["max_abs"],
        ),
    )
    gobernante = analisis["solicitaciones"][combo]
    momento = max(gobernante["Mx"]["max_abs"], gobernante["My"]["max_abs"])
    corte = max(gobernante["Qx"]["max_abs"], gobernante["Qy"]["max_abs"])
    reacciones = analisis["reacciones"]

    resultado = ResultadoLosaMaciza(
        nombre=nombre,
        carga_id=str(elemento.get("id", "")),
        metodo=analisis["metodo"] + " · " + analisis["modelo"],
        criterio_espesor=(
            "CIRSOC, Tabla 7.3.1.1: L/20, redondeado hacia arriba al cm "
            "(losa maciza unidireccional simplemente apoyada; hormigón normal, fy=420 MPa)"
            if tipologia == "maciza" and elemento.get("espesor_automatico", False)
            else (
                "Espesor ingresado manualmente; verificar límites de espesor y flecha "
                "para losa en dos direcciones."
                if tipologia == "maciza_cruzada"
                else (
                    "Canto ingresado; sección geométrica contrastada con "
                    f"{seccion_casetonada['referencia_cirsoc']}."
                )
                if tipologia == "casetonada"
                else "Espesor ingresado manualmente"
            )
        ),
        luz_m=luz,
        ancho_m=ancho,
        espesor_m=espesor,
        peso_propio_kNm2=peso_propio,
        componentes=componentes,
        D2_kNm2=D2,
        uso_sobrecarga=uso_sobrecarga,
        D_kNm2=D,
        L_kNm2=L,
        q_servicio_kNm2=D + L,
        combinaciones_kNm2=combinaciones,
        combo_critica=combo,
        momento_ultimo_kNm_m=momento,
        corte_ultimo_kN_m=corte,
        reacciones=reacciones,
        tipologia=str(tipologia),
        seccion_casetonada=seccion_casetonada,
        analisis_pynite=analisis,
    )
    resultado.memoria = generar_memoria(resultado)
    return resultado


def calcular_grupo_desde_elementos(
    elementos: dict[str, dict], biblioteca: dict | None = None
) -> dict[str, ResultadoLosaMaciza]:
    """Calcula juntas losas macizas enlazadas como continuas."""
    if len(elementos) < 2:
        raise ValueError("El grupo continuo debe contener al menos dos paños.")
    bib = biblioteca if biblioteca is not None else materiales.cargar()
    panos = []
    for nombre, elemento in elementos.items():
        if elemento.get("tipo") != "losa" or elemento.get("tipologia") not in (
            "maciza", "maciza_cruzada"
        ):
            raise ValueError(
                f"El grupo continuo solo admite losas macizas; revisar {nombre}."
            )
        espesor = _espesor_hormigon(elemento)
        valores = cargas.valores_superficiales(elemento, bib)
        panos.append({
            "nombre": nombre,
            "elemento": elemento,
            "espesor_m": espesor,
            "D_kNm2": float(valores["D_kNm2"]),
            "L_kNm2": float(valores["L_kNm2"]),
        })
    analisis_por_nombre = analizar_grupo(panos)
    return {
        nombre: calcular_desde_elemento_carga(
            nombre,
            elemento,
            bib,
            analisis_override=analisis_por_nombre[nombre],
        )
        for nombre, elemento in elementos.items()
    }


def reacciones_y_anchos_por_borde(
    elemento: dict,
    biblioteca: dict | None = None,
    elementos: dict[str, dict] | None = None,
) -> dict[str, dict]:
    """Devuelve reacciones lineales y anchos equivalentes para cada borde apoyado."""
    nombre = str(elemento.get("id", "losa"))
    if elementos is not None:
        nombre = next(
            (
                nombre_elemento
                for nombre_elemento, otro in elementos.items()
                if str(otro.get("id", "")) == str(elemento.get("id", ""))
                and elemento.get("id")
            ),
            nombre,
        )
        grupo = grupo_continuo(nombre, elementos) if nombre in elementos else {}
        grupo_activo = {
            nombre_grupo: losa
            for nombre_grupo, losa in grupo.items()
            if losa.get("activo", True)
        }
        if nombre in grupo_activo and len(grupo_activo) > 1:
            resultado = calcular_grupo_desde_elementos(
                grupo_activo, biblioteca
            )[nombre]
        else:
            resultado = calcular_desde_elemento_carga(nombre, elemento, biblioteca)
    else:
        resultado = calcular_desde_elemento_carga(nombre, elemento, biblioteca)
    carga_referencia = (
        resultado.D_kNm2 if resultado.D_kNm2 > 0
        else resultado.L_kNm2
    )
    if carga_referencia <= 0:
        raise ValueError("No hay carga superficial positiva para transferir las reacciones.")
    clave_referencia = "D" if resultado.D_kNm2 > 0 else "L"
    return {
        apoyo["apoyo"]: {
            "cargas_kN_m": dict(apoyo["cargas_kN_m"]),
            "ancho_tributario_m": (
                apoyo["cargas_kN_m"][clave_referencia] / carga_referencia
            ),
        }
        for apoyo in resultado.reacciones
    }


def generar_memoria(r: ResultadoLosaMaciza) -> str:
    casetonada_unidireccional = (
        r.tipologia == "casetonada"
        and "unidireccional" in str(r.analisis_pynite.get("modelo", "")).lower()
    )
    if r.tipologia == "casetonada":
        titulo = (
            "LOSA CASETONADA UNIDIRECCIONAL"
            if casetonada_unidireccional
            else "LOSA CASETONADA BIDIRECCIONAL"
        )
    elif r.tipologia == "maciza_cruzada":
        titulo = "LOSA MACIZA CRUZADA"
    else:
        titulo = "LOSA MACIZA UNIDIRECCIONAL"
    lineas = [
        f"MEMORIA DE CÁLCULO - {titulo}",
        "=" * (23 + len(titulo)),
        "",
        f"Losa: {r.nombre}",
        f"Carga del proyecto: {r.carga_id or '(sin ID)' }",
        f"Luz entre apoyos: {numero_local(r.luz_m)} m",
        f"Ancho del paño: {numero_local(r.ancho_m)} m",
        f"Espesor adoptado: {numero_local(r.espesor_m * 100, 1)} cm",
        f"Criterio de espesor: {r.criterio_espesor}",
        "",
        "Composición y cargas:",
    ]
    for item in r.componentes:
        tipo = "hormigón estructural" if item["estructural"] else "carga accesoria"
        lineas.append(f" - {item['descripcion']} ({tipo}) → {numero_local(item['valor_kNm2'], 3)} kN/m²")
    if r.seccion_casetonada:
        s = r.seccion_casetonada
        lineas.extend([
            "",
            "Geometría de la sección casetonada:",
            f" - Canto total: {numero_local(s['altura_total_m'] * 100, 1)} cm",
            f" - Capa superior: {numero_local(s['espesor_capa_m'] * 100, 1)} cm "
            f"(mínimo aplicado: {numero_local(s['espesor_minimo_capa_m'] * 100, 1)} cm)",
            f" - Nervios X/Y: {numero_local(s['ancho_nervio_x_m'] * 100, 1)} / "
            f"{numero_local(s['ancho_nervio_y_m'] * 100, 1)} cm",
            f" - Separación libre X/Y: {numero_local(s['separacion_libre_x_m'] * 100, 1)} / "
            f"{numero_local(s['separacion_libre_y_m'] * 100, 1)} cm",
            f" - Rigidez equivalente de hormigón: {numero_local(s['espesor_equivalente_hormigon_m'] * 100)} cm "
            "(para peso propio)",
            f" - Aligeramientos estructurales según 8.8.2.1: "
            f"{'sí' if s['aligeramiento_estructural'] else 'no'}",
            f" - Factores ortótropos ingresados por el usuario: "
            f"kx={numero_local(s['kx_mod'], 3)}; ky={numero_local(s['ky_mod'], 3)}",
            f" - Hipótesis: nervios y capa superior monolíticos; "
            f"{'trabajo unidireccional equivalente' if casetonada_unidireccional else 'trabajo en dos direcciones (8.8.1.1)'}",
            f" - Verificación geométrica: {s['referencia_cirsoc']}",
        ])
    lineas.extend([
        f"Subtotal cargas accesorias (D2): {numero_local(r.D2_kNm2, 3)} kN/m²",
        f"Total cargas permanentes (D1+D2): {numero_local(r.D_kNm2, 3)} kN/m²",
        "",
        "Cargas consideradas:",
        f" - D1 (peso propio): {numero_local(r.peso_propio_kNm2, 3)} kN/m²",
        f" - D2 (accesoria): {numero_local(r.D2_kNm2, 3)} kN/m²",
        f" - D total: {numero_local(r.D_kNm2, 3)} kN/m²",
        f" - Sobrecarga seleccionada: {r.uso_sobrecarga} → {numero_local(r.L_kNm2, 3)} kN/m²",
        f" - Carga de servicio (D+L): {numero_local(r.q_servicio_kNm2, 3)} kN/m²",
        "",
        f"Modelo de análisis: {r.metodo}",
        "Combinaciones consideradas (kN/m²):",
    ])
    lineas.extend(f" - {nombre}: {numero_local(valor, 3)}" for nombre, valor in r.combinaciones_kNm2.items())
    lineas.extend([
        "",
        f"Combinación crítica: {r.combo_critica}",
        f"Momento último requerido: {numero_local(r.momento_ultimo_kNm_m, 3)} kN·m/m",
        f"Corte último: {numero_local(r.corte_ultimo_kN_m, 3)} kN/m",
    ])
    continuidad_modelada = "conjunta" in str(r.analisis_pynite.get("metodo", ""))
    momentos_apoyos = (
        r.analisis_pynite.get("momentos_apoyos", {})
        if continuidad_modelada else {}
    )
    if momentos_apoyos:
        lineas.extend([
            "",
            f"Momentos próximos a los apoyos ({r.combo_critica}; kN·m/m):",
        ])
        for lado, momento in momentos_apoyos.get(r.combo_critica, {}).items():
            lineas.append(
                f" - {lado} ({momento['componente']}): "
                f"máx. {numero_local(momento['max'], 3)}; mín. {numero_local(momento['min'], 3)}; "
                f"máx. abs. {numero_local(momento['max_abs'], 3)} "
                f"[{momento.get('metodo', 'resultado de placa')}]"
            )
    lineas.extend([
        "",
        "Análisis de placa en PyNite:",
        f" - Ec = {numero_local(r.analisis_pynite['Ec_MPa'], 0)} MPa ({r.analisis_pynite['Ec_referencia']})",
        f" - Malla: {r.analisis_pynite['cantidad_elementos']} elementos Quad; paso aproximado "
        f"{numero_local(r.analisis_pynite['malla_aproximada_m'])} m",
        f" - Alcance: {r.analisis_pynite['alcance']}",
        "",
        "Reacciones de apoyo (sin combinar; cargas lineales):",
    ])
    for apoyo in r.reacciones:
        c = apoyo["cargas_kN_m"]
        nombre_apoyo = apoyo.get("apoya_en")
        referencia = ""
        if isinstance(nombre_apoyo, dict):
            referencia = (
                f" → muro {nombre_apoyo['nombre']}"
                if nombre_apoyo.get("tipo") == "muro"
                else (
                    f" → pórtico {nombre_apoyo['portico']}, viga {nombre_apoyo['viga']}"
                    if nombre_apoyo.get("tipo") == "portico"
                    else ""
                )
            )
        longitud = r.ancho_m if apoyo["apoyo"] in ("izq", "der") else r.luz_m
        lineas.append(
            f" - {apoyo['apoyo']}{referencia}: D={numero_local(c['D'], 3)}; L={numero_local(c['L'], 3)} kN/m "
            f"(reacción total D={numero_local(apoyo['reaccion_total_D_kN'], 3)}, "
            f"L={numero_local(apoyo['reaccion_total_L_kN'], 3)} kN en {numero_local(longitud)} m)"
        )
    lineas.extend([
        "",
        "=================================",
        (
            "Alcance: solicitaciones elásticas de un grupo de paños conectados."
            if continuidad_modelada
            else "Alcance: solicitaciones elásticas de un paño aislado."
        ),
    ])
    if r.tipologia == "casetonada":
        lineas.append(
            "Los factores kx/ky de la casetonada son parámetros ingresados por el usuario; "
            "CIRSOC no prescribe esos factores."
        )
    lineas.extend([
        (
            "No dimensiona armadura ni verifica flecha, corte o fisuración; "
            "la continuidad declarada se modela solo mediante placas conectadas."
            if continuidad_modelada
            else "No dimensiona armadura ni verifica flecha, corte o continuidad entre paños."
        ),
        "Fin de la memoria de cálculo",
    ])
    return "\n".join(lineas) + "\n"

def guardar_resultado(resultado: ResultadoLosaMaciza) -> tuple[Path, Path]:
    tipo_archivo = (
        "casetonada" if resultado.tipologia == "casetonada"
        else "maciza_cruzada" if resultado.tipologia == "maciza_cruzada"
        else "maciza"
    )
    carpeta = rutas.SAL_LOSAS / ("casetonadas" if resultado.tipologia == "casetonada" else "macizas")
    carpeta.mkdir(parents=True, exist_ok=True)
    rutas.SAL_LOSAS.mkdir(parents=True, exist_ok=True)
    base = rutas.nombre_seguro(resultado.nombre)
    json_path = carpeta / f"{base}.json"
    txt_path = rutas.SAL_LOSAS / f"memoria_losa_{base}.txt"
    datos = asdict(resultado)
    datos.pop("memoria", None)
    rutas.guardar_json(json_path, datos)
    rutas.guardar_texto(txt_path, resultado.memoria)
    return txt_path, json_path

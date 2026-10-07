"""Análisis elástico inicial de losas macizas unidireccionales.

Este primer modelo representa un paño rectangular, simplemente apoyado en dos
bordes opuestos y cargado uniformemente. Calcula solicitaciones y reacciones;
no dimensiona la armadura ni representa continuidad entre paños.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path

from . import cargas, materiales, rutas


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
    memoria: str = ""


def _espesor_hormigon(elemento: dict) -> float:
    """Toma el espesor de la capa de hormigón armado de la composición."""
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
    nombre: str, elemento: dict, biblioteca: dict | None = None
) -> ResultadoLosaMaciza:
    """Calcula M, V y reacciones de un paño macizo unidireccional simplemente apoyado."""
    if elemento.get("tipo") != "losa" or elemento.get("tipologia") != "maciza":
        raise ValueError("Este análisis corresponde a una losa maciza.")

    luz = float(elemento.get("luz_transversal_m", 0.0) or 0.0)
    ancho = float(elemento.get("ancho_losa_m", 0.0) or 0.0)
    if luz <= 0 or ancho <= 0:
        raise ValueError("Completá la luz entre apoyos y el ancho del paño.")

    espesor = _espesor_hormigon(elemento)
    bib = biblioteca if biblioteca is not None else materiales.cargar()
    valores = cargas.valores_superficiales(elemento, bib)
    D = float(valores["D_kNm2"])
    L = float(valores["L_kNm2"])
    componente_estructural = next((
        c for c in elemento.get("componentes", [])
        if (c.get("grupo"), c.get("clave")) in {
            ("Forjados", "Losa_maciza"), ("Hormigon", "armado")
        }
    ), None)
    peso_propio = cargas.valor_componente_superficial(componente_estructural, bib)[1]
    componentes = []
    D2 = 0.0
    for componente, (descripcion, valor) in zip(elemento.get("componentes", []), valores["componentes"]):
        estructural = (componente.get("grupo"), componente.get("clave")) in {
            ("Forjados", "Losa_maciza"), ("Hormigon", "armado")
        }
        componentes.append({"descripcion": descripcion, "valor_kNm2": valor, "estructural": estructural})
        if not estructural:
            D2 += valor
    uso_sobrecarga = next((
        str(item.get("uso", elemento.get("sobrecarga", "")))
        for item in bib.get("Sobrecargas", [])
        if item.get("clave") == elemento.get("sobrecarga")
    ), str(elemento.get("sobrecarga", "")))
    combinaciones = cargas.combinaciones_de_componentes(D, L, 0.0)
    combo, q_u = max(combinaciones.items(), key=lambda par: par[1])

    # Fórmulas elásticas para una franja de un metro, simplemente apoyada.
    momento = q_u * luz**2 / 8.0
    corte = q_u * luz / 2.0
    reacciones = []
    apoya = elemento.get("apoya_en", {})
    for lado in ("izq", "der"):
        reacciones.append({
            "apoyo": lado,
            "apoya_en": apoya.get(lado) if isinstance(apoya, dict) else None,
            "cargas_kN_m": {"D": D * luz / 2.0, "L": L * luz / 2.0},
            "combinaciones_kN_m": {
                nombre_combo: q * luz / 2.0
                for nombre_combo, q in combinaciones.items()
            },
            "reaccion_total_D_kN": D * luz * ancho / 2.0,
            "reaccion_total_L_kN": L * luz * ancho / 2.0,
        })

    resultado = ResultadoLosaMaciza(
        nombre=nombre,
        carga_id=str(elemento.get("id", "")),
        metodo="Franja de 1 m; unidireccional, simplemente apoyada",
        criterio_espesor=(
            "CIRSOC 201-25, Tabla 7.3.1.1: L/20, redondeado hacia arriba al cm "
            "(hormigón normal, fy=420 MPa)"
            if elemento.get("espesor_automatico", False)
            else "Espesor ingresado manualmente"
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
    )
    resultado.memoria = generar_memoria(resultado)
    return resultado


def generar_memoria(r: ResultadoLosaMaciza) -> str:
    lineas = [
        "MEMORIA DE CÁLCULO - LOSA MACIZA",
        "=================================",
        "",
        f"Losa: {r.nombre}",
        f"Carga del proyecto: {r.carga_id or '(sin ID)' }",
        f"Luz entre apoyos: {r.luz_m:.2f} m",
        f"Ancho del paño: {r.ancho_m:.2f} m",
        f"Espesor adoptado: {r.espesor_m * 100:.1f} cm",
        f"Criterio de espesor: {r.criterio_espesor}",
        "",
        "Materiales seleccionados:",
    ]
    for item in r.componentes:
        tipo = "peso propio estructural" if item["estructural"] else "carga accesoria"
        lineas.append(f" - {item['descripcion']} ({tipo}) → {item['valor_kNm2']:.3f} kN/m²")
    lineas.extend([
        f"Subtotal cargas accesorias (D2): {r.D2_kNm2:.3f} kN/m²",
        f"Total cargas permanentes (D1+D2): {r.D_kNm2:.3f} kN/m²",
        "",
        "Cargas consideradas:",
        f" - D1 (peso propio): {r.peso_propio_kNm2:.3f} kN/m²",
        f" - D2 (accesoria): {r.D2_kNm2:.3f} kN/m²",
        f" - D total: {r.D_kNm2:.3f} kN/m²",
        f" - Sobrecarga seleccionada: {r.uso_sobrecarga} → {r.L_kNm2:.3f} kN/m²",
        f" - Carga de servicio (D+L): {r.q_servicio_kNm2:.3f} kN/m²",
        "",
        f"Modelo de análisis: {r.metodo}",
        "Combinaciones consideradas (kN/m²):",
    ])
    lineas.extend(f" - {nombre}: {valor:.3f}" for nombre, valor in r.combinaciones_kNm2.items())
    lineas.extend([
        "",
        f"Combinación crítica: {r.combo_critica}",
        f"Momento último requerido: {r.momento_ultimo_kNm_m:.3f} kN·m/m",
        f"Corte último: {r.corte_ultimo_kN_m:.3f} kN/m",
        "",
        "Reacciones por borde largo:",
    ])
    for apoyo in r.reacciones:
        c = apoyo["cargas_kN_m"]
        lineas.append(f" - {apoyo['apoyo']}: D={c['D']:.3f}, L={c['L']:.3f} kN/m")
    lineas.extend([
        "",
        "=================================",
        "Alcance: solicitaciones elásticas de un paño aislado.",
        "No dimensiona armadura ni verifica flecha, corte o continuidad entre paños.",
        "Fin de la memoria de cálculo",
    ])
    return "\n".join(lineas) + "\n"

def guardar_resultado(resultado: ResultadoLosaMaciza) -> tuple[Path, Path]:
    carpeta = rutas.SAL_LOSAS / "macizas"
    carpeta.mkdir(parents=True, exist_ok=True)
    rutas.SAL_LOSAS.mkdir(parents=True, exist_ok=True)
    base = rutas.nombre_seguro(resultado.nombre)
    json_path = carpeta / f"{base}.json"
    txt_path = rutas.SAL_LOSAS / f"memoria_losa_maciza_{base}_{datetime.now():%Y%m%d_%H%M%S}.txt"
    datos = asdict(resultado)
    datos.pop("memoria", None)
    rutas.guardar_json(json_path, datos)
    rutas.guardar_texto(txt_path, resultado.memoria)
    return txt_path, json_path

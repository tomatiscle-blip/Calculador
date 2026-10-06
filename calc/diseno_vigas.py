"""Adaptador entre el motor de pórticos actual y el dimensionador legacy P02."""

from __future__ import annotations

import json
from pathlib import Path

from . import cargas, rutas


def _mayor_absoluto(bloques: list[dict], grupo: str, tramo_id: str, campo: str) -> float:
    valores = [
        float(bloque.get(grupo, {}).get(tramo_id, {}).get(campo, 0.0))
        for bloque in bloques
        if bloque.get(grupo, {}).get(tramo_id, {}).get(campo) is not None
    ]
    return max(valores, key=abs) if valores else 0.0


def preparar_datos(portico: str, secciones: dict | None = None) -> tuple[dict, list[str]]:
    """Construye la entrada que P02 esperaba usando geometría y motor actuales."""
    estructura = rutas.cargar_estructura()
    if portico not in estructura:
        raise ValueError(f"No existe {portico} en datos/estructura.json")
    archivo_motor = rutas.SAL_SOLICITACIONES / f"{rutas.nombre_seguro(portico)}.json"
    motor = rutas.leer_json(archivo_motor)
    if not motor:
        raise ValueError(f"Primero resolvé las solicitaciones de {portico}.")
    if motor.get("avisos"):
        raise ValueError("El motor tiene avisos; corregilos antes de dimensionar vigas: " + "; ".join(motor["avisos"]))

    bloques = motor.get("solicitaciones", {})
    servicio = bloques.get("servicio")
    if not servicio:
        raise ValueError("El resultado del motor no contiene la combinación de servicio.")
    ultimos = [bloque for nombre, bloque in bloques.items() if nombre != "servicio"]
    if not ultimos:
        raise ValueError("El resultado del motor no contiene combinaciones de dimensionado.")

    datos_portico = estructura[portico]
    aplicadas = [a for a in motor.get("cargas_aplicadas", []) if a.get("tipo") in ("D", "L")]
    advertencias: list[str] = []
    datos_vigas: dict = {}
    secciones = secciones or {}
    for viga_id, viga in datos_portico.get("vigas", {}).items():
        entrada = {k: v for k, v in viga.items() if k in ("b_cm", "h_cm", "recubrimiento_cm", "fc_MPa", "fy_MPa")}
        entrada.update(secciones.get(viga_id, {}))
        if entrada.get("b_cm") is None or entrada.get("fc_MPa") is None:
            raise ValueError(f"Faltan ancho b y/o resistencia fc para {viga_id}; completalos en la ventana de vigas.")
        entrada["tramos"] = []
        for tramo in viga.get("tramos", []):
            tramo_id = tramo.get("id", "")
            grupo = "voladizos" if tramo.get("es_voladizo") else "vigas"
            bloque_serv = servicio.get(grupo, {}).get(tramo_id)
            if not bloque_serv:
                raise ValueError(f"El motor no tiene resultados de servicio para {tramo_id}.")

            aplicadas_tramo = [a for a in aplicadas if a.get("tramo_id") == tramo_id]
            longitud = float(tramo.get("longitud_m", 0.0))
            if longitud <= 0:
                raise ValueError(f"{tramo_id}: la longitud debe ser mayor que cero.")
            q = {"D": 0.0, "L": 0.0}
            for carga in aplicadas_tramo:
                ancho_aplicado = max(0.0, min(longitud, float(carga.get("x_fin_m", longitud))) - max(0.0, float(carga.get("x_inicio_m", 0.0))))
                q[carga["tipo"]] += float(carga["valor_kN_m"]) * ancho_aplicado / longitud
                if ancho_aplicado < longitud - 1e-8:
                    advertencias.append(f"{tramo_id}: P02 estima flecha con carga uniforme equivalente; hay cargas parciales.")
            carga_puntual = tramo.get("cargas_puntuales", [])
            if carga_puntual:
                advertencias.append(f"{tramo_id}: revisar cargas puntuales; P02 no las incluye en su estimación de flecha.")
            combinaciones = cargas.combinaciones_de_componentes(q["D"], q["L"], 0.0)
            datos = {
                "id": tramo_id,
                "longitud_m": longitud,
                "es_voladizo": bool(tramo.get("es_voladizo", False)),
                "cargas": {"D_total": q["D"], "L_total": q["L"], "W_total": 0.0, "combinaciones": combinaciones},
                "Mu_kNm_apoyo_izq": abs(_mayor_absoluto(ultimos, grupo, tramo_id, "M_izq")),
                "Mu_kNm_apoyo_der": abs(_mayor_absoluto(ultimos, grupo, tramo_id, "M_der")),
                "Mu_kNm_campo": abs(_mayor_absoluto(ultimos, grupo, tramo_id, "M_campo")),
                "Vu_kN_izq": abs(_mayor_absoluto(ultimos, grupo, tramo_id, "V_izq")),
                "Vu_kN_der": abs(_mayor_absoluto(ultimos, grupo, tramo_id, "V_der")),
                "Mu_servicio_izq": abs(float(bloque_serv.get("M_izq", 0.0))),
                "Mu_servicio_der": abs(float(bloque_serv.get("M_der", 0.0))),
                "Mu_servicio_campo": abs(float(bloque_serv.get("M_campo", 0.0))),
                "Vu_servicio_izq": abs(float(bloque_serv.get("V_izq", 0.0))),
                "Vu_servicio_der": abs(float(bloque_serv.get("V_der", 0.0))),
            }
            if tramo.get("es_voladizo"):
                datos["Mu_kNm"] = datos["Mu_kNm_apoyo_izq"]
                datos["Vu_kN_emp"] = datos["Vu_kN_izq"]
            entrada["tramos"].append(datos)
        datos_vigas[viga_id] = entrada
    return datos_vigas, sorted(set(advertencias))


def dimensionar(portico: str, secciones: dict | None = None) -> tuple[Path, list[str]]:
    """Ejecuta los cálculos conservados en P02 y genera sus salidas habituales."""
    import P02_Viga_portico as legacy

    datos_vigas, advertencias = preparar_datos(portico, secciones)
    coef_kd = rutas.leer_json(rutas.COEFICIENTES_KD)
    if not coef_kd:
        raise FileNotFoundError(f"Falta la tabla {rutas.COEFICIENTES_KD}")
    legacy.RESULTADOS = legacy._resultados_vacios()
    legacy.procesar_vigas(datos_vigas, coef_kd, rutas.SAL_VIGAS, portico)
    salida = rutas.SAL_VIGAS / f"resultados_{rutas.nombre_seguro(portico)}_vigas.json"
    legacy.guardar_resultados_json(salida)
    return salida, advertencias

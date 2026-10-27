"""Cálculo consolidado de vigas globales (M/V por tramo de pórtico)."""

from __future__ import annotations

import os
from pathlib import Path

try:
    from . import rutas
except ImportError:  # ejecución directa como script
    from calc import rutas


def _nivel_de_viga(viga_id: str) -> str:
    texto = str(viga_id or "")
    if texto.startswith("V") and len(texto) > 1 and texto[1].isdigit():
        return f"V{texto[1]}"
    return ""


def _mapa_tramos_portico(portico: str, estructura: dict) -> dict[str, str]:
    """Devuelve {tramo_id: viga_id} para un pórtico."""
    resultado: dict[str, str] = {}
    datos = (estructura or {}).get(portico, {}) or {}
    for viga_id, viga in (datos.get("vigas", {}) or {}).items():
        for tramo in (viga.get("tramos", []) or []):
            tid = tramo.get("id")
            if tid:
                resultado[str(tid)] = str(viga_id)
    return resultado


def calcular_global() -> list[dict]:
    estructura = rutas.cargar_estructura()
    resultados: list[dict] = []
    for archivo in sorted(rutas.listar(rutas.SAL_SOLICITACIONES, "*.json")):
        portico = archivo.stem
        datos = rutas.leer_json(archivo, {}) or {}
        if datos.get("avisos"):
            continue
        mapa_tramos = _mapa_tramos_portico(portico, estructura)
        envolvente = datos.get("envolvente", {}) or {}
        for grupo in ("vigas", "voladizos"):
            for tramo_id, valores in (envolvente.get(grupo, {}) or {}).items():
                valores = valores or {}
                m_izq = float(valores.get("M_izq", 0.0) or 0.0)
                m_der = float(valores.get("M_der", 0.0) or 0.0)
                m_campo = float(valores.get("M_campo", 0.0) or 0.0)
                m = max(abs(m_izq), abs(m_der), abs(m_campo))
                v_izq = float(valores.get("V_izq", 0.0) or 0.0)
                v_der = float(valores.get("V_der", 0.0) or 0.0)
                v = max(abs(v_izq), abs(v_der))
                viga_id = mapa_tramos.get(str(tramo_id), "")
                resultados.append({
                    "id": tramo_id,
                    "portico": portico,
                    "nivel": _nivel_de_viga(viga_id),
                    "viga": viga_id,
                    "tipo": grupo,
                    "M": m,
                    "V": v,
                    "M_izq": m_izq,
                    "M_der": m_der,
                    "M_campo": m_campo,
                    "V_izq": v_izq,
                    "V_der": v_der,
                    "estado": "OK" if m < 200 and v < 400 else "Revisar",
                })

    salida = rutas.SAL_VIGAS / "resultados_vigas_globales.json"
    rutas.SAL_VIGAS.mkdir(parents=True, exist_ok=True)
    rutas.guardar_json(salida, resultados)
    return resultados


if __name__ == "__main__":
    print(calcular_global())
"""Calculo de vigas de fundacion (arriostramiento entre bases globales)."""
from __future__ import annotations
import math
import sys
from pathlib import Path
_RAIZ = Path(__file__).resolve().parent.parent
if str(_RAIZ) not in sys.path:
    sys.path.insert(0, str(_RAIZ))

try:
    from . import rutas, bases
except ImportError:
    from calc import rutas, bases


def _coord(bid, estructura):
    info = (estructura.get("bases", {}) or {}).get(bid, {}) or {}
    return float(info.get("X", 0.0) or 0.0), float(info.get("Y", 0.0) or 0.0)


def calcular() -> list[dict]:
    datos = rutas.leer_json(rutas.SAL_BASES / "resultados_bases_globales.json", []) or []
    estructura = rutas.cargar_estructura()
    resultados = []
    vistos = set()
    for portico in sorted({p for b in datos for p in (b.get("porticos") or [])}):
        bases_portico = [b for b in datos if portico in (b.get("porticos") or [])]
        # ordenar por X luego Y
        bases_portico.sort(key=lambda b: _coord(b.get("id"), estructura))
        for i in range(len(bases_portico) - 1):
            b1 = bases_portico[i]
            b2 = bases_portico[i + 1]
            id1 = b1.get("id")
            id2 = b2.get("id")
            clave = tuple(sorted((id1, id2)))
            if clave in vistos:
                continue
            vistos.add(clave)
            x1, y1 = _coord(id1, estructura)
            x2, y2 = _coord(id2, estructura)
            L = math.hypot(x2 - x1, y2 - y1)
            if L < 0.1:
                continue
            h = max(0.30, round(L / 10.0, 2))
            b = max(0.25, round(h * 0.5, 2))
            armadura = "4φ16" if L <= 4.0 else "6φ16"
            estado = "OK" if L <= 6.0 else "Revisar"
            resultados.append({
                "id": f"VF-{len(resultados)+1:02d}",
                "bases": [id1, id2],
                "portico": portico,
                "longitud_m": round(L, 2),
                "ancho_m": b,
                "altura_m": h,
                "armadura": armadura,
                "estado": estado,
            })
    salida = rutas.SAL_BASES / "resultados_vigas_fundacion.json"
    rutas.SAL_BASES.mkdir(parents=True, exist_ok=True)
    rutas.guardar_json(salida, resultados)
    return resultados


if __name__ == "__main__":
    print(calcular())
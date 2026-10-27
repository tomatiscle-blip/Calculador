import json
import os
from pathlib import Path

from . import rutas


def calcular(datos, portico="Portico 4"):
    entrada = os.path.join("salidas", "solicitaciones", f"{portico}.json")
    salida = os.path.join("salidas", "bases", f"resultados_{portico}_bases.json")

    with open(entrada, "r") as f:
        solicitaciones = json.load(f)

    resultados = []
    for bid, valores in solicitaciones.get("envolvente", {}).get("bases", {}).items():
        resultados.append({"id": bid, **valores})

    os.makedirs(os.path.dirname(salida), exist_ok=True)
    with open(salida, "w") as f:
        json.dump(resultados, f, indent=4)

    return resultados


def calcular_global() -> list[dict]:
    consolidado = rutas.leer_json(rutas.SOLICITACIONES_BASES, {}) or {}
    estructura = rutas.cargar_estructura()
    resultados = []
    for bid in sorted(consolidado):
        entradas = consolidado[bid] or []
        combinado: dict = {}
        # Cargas verticales/horizontales: se toma el máximo en módulo.
        for clave in ("Fx", "Fy"):
            valores = [float(e.get(clave, 0.0) or 0.0) for e in entradas]
            if valores:
                combinado[clave] = max(valores, key=abs)
        # Momentos por eje: Mz del pórtico se asigna a Mx o My según dirección.
        mx = 0.0
        my = 0.0
        for entrada in entradas:
            mz = float(entrada.get("Mz", 0.0) or 0.0)
            try:
                portico = rutas.pórtico(estructura, entrada.get("portico")) or {}
                direccion = (portico.get("referencia_planta", {}) or {}).get("direccion")
            except Exception:
                direccion = None
            if direccion == "x":
                if abs(mz) > abs(my):
                    my = mz
            else:
                if abs(mz) > abs(mx):
                    mx = mz
        combinado["Mx"] = mx
        combinado["My"] = my
        combinado["M"] = max(abs(mx), abs(my))
        resultados.append({
            "id": bid,
            "porticos": [e.get("portico") for e in entradas],
            "n_porticos": len(entradas),
            "momentos_otro_sentido": any(bool(e.get("momentos_otro_sentido")) for e in entradas),
            **combinado,
            "estado": "OK" if combinado.get("M", 0) < 200 and abs(combinado.get("Fy", 0)) < 400 else "Revisar",
        })

    salida = rutas.SAL_BASES / "resultados_bases_globales.json"
    rutas.SAL_BASES.mkdir(parents=True, exist_ok=True)
    rutas.guardar_json(salida, resultados)
    return resultados


if __name__ == "__main__":
    print(calcular({}))

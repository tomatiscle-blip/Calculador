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
    resultados = []
    for bid in sorted(consolidado):
        entradas = consolidado[bid] or []
        combinado: dict = {}
        for entrada in entradas:
            for clave, valor in entrada.items():
                if clave == "portico":
                    continue
                if isinstance(valor, (int, float)):
                    if clave not in combinado or abs(float(valor)) > abs(float(combinado[clave])):
                        combinado[clave] = valor
        resultados.append({
            "id": bid,
            "porticos": [e.get("portico") for e in entradas],
            "n_porticos": len(entradas),
            "momentos_otro_sentido": any(bool(e.get("momentos_otro_sentido")) for e in entradas),
            **combinado,
        })

    salida = rutas.SAL_BASES / "resultados_bases_globales.json"
    rutas.SAL_BASES.mkdir(parents=True, exist_ok=True)
    rutas.guardar_json(salida, resultados)
    return resultados


if __name__ == "__main__":
    print(calcular({}))

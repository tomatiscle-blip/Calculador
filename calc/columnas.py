import json
import os
from . import rutas

def calcular(datos, portico="Portico 4"):
    entrada = os.path.join("salidas", "solicitaciones", f"{portico}.json")
    salida = os.path.join("salidas", "columnas", f"resultados_{portico}_columnas.json")

    with open(entrada, "r") as f:
        solicitaciones = json.load(f)

    resultados = []

    # Leer la envolvente de columnas
    for cid, valores in solicitaciones.get("envolvente", {}).get("columnas", {}).items():
        Mu = max(valores.get("M_inf", 0), valores.get("M_sup", 0))
        N = valores.get("N", 0)
        chequeo = {
            "id": cid,
            "Mu": Mu,
            "N": N,
            "estado": "OK" if Mu < 200 and N < 400 else "Revisar"
        }
        resultados.append(chequeo)

    os.makedirs(os.path.dirname(salida), exist_ok=True)
    with open(salida, "w") as f:
        json.dump(resultados, f, indent=4)

    return resultados


def calcular_global() -> list[dict]:
    """Dimensiona columnas usando el consolidado por ID global."""
    consolidado = rutas.leer_json(rutas.SOLICITACIONES_COLUMNAS, {}) or {}
    resultados = []
    for cid, lista in consolidado.items():
        M_inf = sum(float(e.get("M_inf", 0) or 0) for e in lista)
        M_sup = sum(float(e.get("M_sup", 0) or 0) for e in lista)
        N = sum(float(e.get("N", 0) or 0) for e in lista)
        Mu = max(abs(M_inf), abs(M_sup))
        resultados.append({
            "id": cid,
            "Mu": Mu,
            "N": N,
            "porticos": [e.get("portico") for e in lista],
            "estado": "OK" if Mu < 200 and N < 400 else "Revisar",
        })
    salida = rutas.SAL_COLUMNAS / "resultados_columnas_globales.json"
    rutas.SAL_COLUMNAS.mkdir(parents=True, exist_ok=True)
    rutas.guardar_json(salida, resultados)
    return resultados


if __name__ == "__main__":
    print(calcular({}))



import json
import os

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

if __name__ == "__main__":
    print(calcular({}))


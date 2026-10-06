import json
import os

def calcular(datos, portico="Portico 4"):
    entrada = os.path.join("salidas", "solicitaciones", f"{portico}.json")
    salida = os.path.join("salidas", "vigas", f"resultados_{portico}_vigas.json")

    with open(entrada, "r") as f:
        solicitaciones = json.load(f)

    resultados = []

    # Leer la envolvente de vigas
    for tramo, valores in solicitaciones.get("envolvente", {}).get("vigas", {}).items():
        Mu = valores.get("M_campo", 0)
        V = valores.get("V_izq", 0)
        chequeo = {
            "id": tramo,
            "Mu": Mu,
            "V": V,
            "estado": "OK" if Mu < 200 else "Revisar"
        }
        resultados.append(chequeo)

    os.makedirs(os.path.dirname(salida), exist_ok=True)
    with open(salida, "w") as f:
        json.dump(resultados, f, indent=4)

    return resultados

if __name__ == "__main__":
    print(calcular({}))

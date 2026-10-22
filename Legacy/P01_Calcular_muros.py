from calc import cargas

def calcular_muros():
    # Leer muros definidos en datos/cargas.json
    muros = cargas.leer_elemento("Muros")
    
    resultados = []
    for muro in muros:
        carga = cargas.calcular(muro)  # usa la cuenta única
        resultados.append({"id": muro["id"], "carga": carga})
    
    return resultados

if __name__ == "__main__":
    for r in calcular_muros():
        print(f"Muro {r['id']}: {r['carga']} kN/m")

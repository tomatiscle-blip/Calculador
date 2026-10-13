from calc import materiales

def calcular_muro(elemento: dict) -> dict:
    densidad = materiales.buscar(elemento["clave"])  # kg/m³
    espesor = float(elemento["espesor_m"])
    altura = float(elemento["altura_m"])
    q_kN_m = densidad * espesor * altura * 9.81 / 1000
    return {
        "id": elemento.get("id"),
        "tipo": "muro",
        "carga_kN_m": q_kN_m
    }

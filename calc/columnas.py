import json
import os
from pathlib import Path

from . import rutas


def calcular(datos, portico="Portico 4"):
    entrada = os.path.join("salidas", "solicitaciones", f"{portico}.json")
    salida = os.path.join("salidas", "columnas", f"resultados_{portico}_columnas.json")

    with open(entrada, "r") as f:
        solicitaciones = json.load(f)

    resultados = []

    for cid, valores in solicitaciones.get("envolvente", {}).get("columnas", {}).items():
        Mu = max(valores.get("M_inf", 0), valores.get("M_sup", 0))
        N = valores.get("N", 0)
        chequeo = {
            "id": cid,
            "Mu": Mu,
            "N": N,
            "estado": "OK" if Mu < 200 and N < 400 else "Revisar",
        }
        resultados.append(chequeo)

    os.makedirs(os.path.dirname(salida), exist_ok=True)
    with open(salida, "w") as f:
        json.dump(resultados, f, indent=4)

    return resultados


def calcular_global() -> list[dict]:
    acumulado = rutas.acumular_solicitaciones_columnas() or {}
    estructura = rutas.cargar_estructura()
    resultados = []
    for gid in sorted(acumulado):
        archivo = Path(acumulado[gid])
        datos = rutas.leer_json(archivo, {}) or {}
        env = datos.get("envolvente_global", {}) or {}
        catalogo = datos.get("catalogo", {}) or {}
        n_val = float(env.get("N", 0.0) or 0.0)
        deriva = float(env.get("deriva", 0.0) or 0.0)
        porticos = None
        for clave in ("pórticos", "porticos", "pÃ³rticos"):
            val = datos.get(clave)
            if isinstance(val, dict):
                porticos = val
                break
        if porticos is None:
            porticos = {}
        mx = 0.0
        my = 0.0
        for nombre_portico, pdata in (porticos or {}).items():
            env_p = (pdata or {}).get("envolvente_pórtico", {}) or {}
            m_inf = float(env_p.get("M_inf", 0.0) or 0.0)
            m_sup = float(env_p.get("M_sup", 0.0) or 0.0)
            m = max(abs(m_inf), abs(m_sup))
            try:
                p = rutas.pórtico(estructura, nombre_portico) or {}
                direccion = (p.get("referencia_planta", {}) or {}).get("direccion")
            except Exception:
                direccion = None
            if direccion == "x":
                my = max(my, m)
            else:
                mx = max(mx, m)
        Mu = max(mx, my)
        porticos_nombres = list(porticos.keys())
        resultados.append({
            "id": gid,
            "eje_id": catalogo.get("eje_id"),
            "eje_nombre": catalogo.get("eje_nombre"),
            "x": catalogo.get("x"),
            "nivel": catalogo.get("nivel"),
            "Mx": mx,
            "My": my,
            "Mu": Mu,
            "N": n_val,
            "deriva": deriva,
            "porticos": porticos_nombres,
            "n_porticos": len(porticos_nombres),
            "estado": "OK" if Mu < 200 and n_val < 400 else "Revisar",
        })

    salida = rutas.SAL_COLUMNAS / "resultados_columnas_globales.json"
    rutas.SAL_COLUMNAS.mkdir(parents=True, exist_ok=True)
    rutas.guardar_json(salida, resultados)
    return resultados


if __name__ == "__main__":
    print(calcular({}))

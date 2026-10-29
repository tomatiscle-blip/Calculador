p = 'calc/cargas.py'
s = open(p, encoding='utf-8').read()
old = '''def items_viento_general(datos: dict, biblioteca: dict) -> list[dict]:
    """
    El viento que empuja al portico (carga horizontal). No es de un elemento:
    es del conjunto, por eso se suma aparte.
    """
    configuracion = datos.get("viento", {})
    if not configuracion.get("activo"):
        return []
    viento = materiales.viento(biblioteca)
    ancho = float(configuracion.get("ancho_tributario_m", 0.0))
    w_viento = carga_lineal_viento(viento["V"], ancho, viento["Cd"], viento["rho"])
    return [{
        "descripcion": f"Viento diseno {viento['ciudad']} (V={viento['V']:.0f} m/s) – general",
        "tipo": "W",
        "valor": w_viento,
        "unidad": "kN/m",
        "composicion": f"b = {ancho:.2f} m · Cd = {viento['Cd']}",
    }]'''
new = '''def items_viento_general(datos: dict, biblioteca: dict, ancho_override: float | None = None) -> list[dict]:
    """
    El viento que empuja al portico (carga horizontal). No es de un elemento:
    es del conjunto, por eso se suma aparte.
    """
    configuracion = datos.get("viento", {})
    if not configuracion.get("activo"):
        return []
    viento = materiales.viento(biblioteca)
    if ancho_override is not None and ancho_override > 0:
        ancho = float(ancho_override)
    else:
        ancho = float(configuracion.get("ancho_tributario_m", 0.0))
    w_viento = carga_lineal_viento(viento["V"], ancho, viento["Cd"], viento["rho"])
    return [{
        "descripcion": f"Viento diseno {viento['ciudad']} (V={viento['V']:.0f} m/s) – general",
        "tipo": "W",
        "valor": w_viento,
        "unidad": "kN/m",
        "composicion": f"b = {ancho:.2f} m · Cd = {viento['Cd']}",
    }]'''
assert old in s, "old not found"
s = s.replace(old, new, 1)
open(p, "w", encoding="utf-8").write(s)
print("patched")

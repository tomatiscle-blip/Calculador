p = "script.py"
s = open(p, encoding="utf-8").read()
old = """    porticos_activos = viento_config.get("porticos_activos")
    anchos_por_portico = viento_config.get("anchos_por_portico", {}) or {}
    if not porticos_activos:
        if viento_config.get("activo"):
            direccion = str(viento_config.get("direccion_viento", "x")).strip().lower()
            if direccion == direccion_portico(portico):
                items_viento = cargas.items_viento_general(datos_carga_proyecto, materiales.cargar(), ancho_override=anchos_por_portico.get(nombre))
                w_general = sum(float(item["valor"]) for item in items_viento)
    else:
        if nombre in (porticos_activos or []):
            items_viento = cargas.items_viento_general(datos_carga_proyecto, materiales.cargar(), ancho_override=anchos_por_portico.get(nombre))
            w_general = sum(float(item["valor"]) for item in items_viento)
"""
new = """    # Viento global: se aplica a porticos paralelos a la direccion del viento.
    if viento_config.get("activo"):
        direccion = str(viento_config.get("direccion_viento", "x")).strip().lower()
        if direccion == direccion_portico(portico):
            anchos_por_portico = viento_config.get("anchos_por_portico", {}) or {}
            ancho = anchos_por_portico.get(nombre, viento_config.get("ancho_tributario_m", 0.0))
            items_viento = cargas.items_viento_general(datos_carga_proyecto, materiales.cargar(), ancho_override=ancho)
            w_general = sum(float(item["valor"]) for item in items_viento)
"""
assert old in s, "old not found in script.py"
s = s.replace(old, new, 1)
open(p, "w", encoding="utf-8").write(s)
print("done")

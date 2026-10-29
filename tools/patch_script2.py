p = "script.py"
s = open(p, encoding="utf-8").read()
old = """    # Viento global: se aplica a porticos paralelos a la direccion del viento.
    if viento_config.get("activo"):
        direccion = str(viento_config.get("direccion_viento", "x")).strip().lower()
        if direccion == direccion_portico(portico):
            anchos_por_portico = viento_config.get("anchos_por_portico", {}) or {}
            ancho = anchos_por_portico.get(nombre, viento_config.get("ancho_tributario_m", 0.0))
            items_viento = cargas.items_viento_general(datos_carga_proyecto, materiales.cargar(), ancho_override=ancho)
            w_general = sum(float(item["valor"]) for item in items_viento)
"""
new = """    # Viento: se aplica a todos los porticos si esta activo (sin filtrar por direccion).
    if viento_config.get("activo"):
        anchos_por_portico = viento_config.get("anchos_por_portico", {}) or {}
        ancho = anchos_por_portico.get(nombre, viento_config.get("ancho_tributario_m", 0.0))
        if not ancho or ancho <= 0:
            ancho = viento_config.get("ancho_tributario_m", 0.0)
        items_viento = cargas.items_viento_general(datos_carga_proyecto, materiales.cargar(), ancho_override=ancho)
        w_general = sum(float(item["valor"]) for item in items_viento)
"""
assert old in s, "old not found in script.py"
s = s.replace(old, new, 1)
open(p, "w", encoding="utf-8").write(s)
print("done")

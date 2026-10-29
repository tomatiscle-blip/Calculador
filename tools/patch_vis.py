p = "app/visualizacion.py"
s = open(p, encoding="utf-8").read()
old = '''        porticos_activos = configuracion_viento.get("porticos_activos")
        activo = bool(configuracion_viento.get("activo"))
        if porticos_activos is None:
            # modo global: solo pórticos paralelos a la dirección del viento
            from calc import rutas
            estructura = rutas.cargar_estructura() or {}
            dir_portico = str((estructura.get(self.nombre, {}).get("referencia_planta") or {}).get("direccion", "x")).lower()
            activo = activo and (dir_portico == self.direccion_viento)
        else:
            activo = bool(porticos_activos) and (self.nombre in (porticos_activos or []))
        if activo:'''
new = '''        # El visor muestra viento si está activo globalmente (el motor decide si aplica).
        activo = bool(configuracion_viento.get("activo"))
        if activo:'''
assert old in s, "old not found"
s = s.replace(old, new, 1)
open(p, "w", encoding="utf-8").write(s)
print("done")
'''
open("tools/patch_vis.py","w",encoding="utf-8").write(s_patch) if False else None
'''
open("tools/patch_vis.py", "w", encoding="utf-8").write("""p = 'app/visualizacion.py'
s = open(p, encoding='utf-8').read()
old = '''        porticos_activos = configuracion_viento.get("porticos_activos")
        activo = bool(configuracion_viento.get("activo"))
        if porticos_activos is None:
            # modo global: solo pórticos paralelos a la dirección del viento
            from calc import rutas
            estructura = rutas.cargar_estructura() or {}
            dir_portico = str((estructura.get(self.nombre, {}).get("referencia_planta") or {}).get("direccion", "x")).lower()
            activo = activo and (dir_portico == self.direccion_viento)
        else:
            activo = bool(porticos_activos) and (self.nombre in (porticos_activos or []))
        if activo:'''
new = '''        # El visor muestra viento si esta activo globalmente (el motor decide si aplica).
        activo = bool(configuracion_viento.get("activo"))
        if activo:'''
assert old in s, "old not found"
s = s.replace(old, new, 1)
open(p, 'w', encoding='utf-8').write(s)
print('done')
""")
print("written")

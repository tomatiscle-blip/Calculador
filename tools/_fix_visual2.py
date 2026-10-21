# -*- coding: utf-8 -*-
"""Remove remaining debug paint_log.txt block from visualizacion.py"""
import io

path = r"d:/03 Ingenieria/Calculador/app/visualizacion.py"

with io.open(path, "r", encoding="utf-8") as f:
    content = f.read()

# Fix: Remove debug paint_log.txt block before legend
old = '''        # --- leyenda visual: muestra de color + símbolo por tipo de elemento ---
        try:
            with open("paint_log.txt", "a", encoding="utf-8") as _f:
                _f.write("before legend\\\\n")
        except Exception:
            pass
        tipos_presentes = sorted({t for t in escena["tipos_losa"].values()})'''

new = '''        # --- leyenda visual: muestra de color + símbolo por tipo de elemento ---
        tipos_presentes = sorted({t for t in escena["tipos_losa"].values()})'''

if old in content:
    content = content.replace(old, new, 1)
    print("Removed debug paint_log block before legend")
else:
    print("OLD TEXT NOT FOUND")
    # Find what's actually there
    idx = content.find("paint_log")
    if idx >= 0:
        start = max(0, idx - 100)
        end = min(len(content), idx + 200)
        print(repr(content[start:end]))

with io.open(path, "w", encoding="utf-8") as f:
    f.write(content)

print("Done")
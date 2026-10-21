# -*- coding: utf-8 -*-
"""Fix incomplete implementation in visualizacion.py"""
import io

path = r"d:/03 Ingenieria/Calculador/app/visualizacion.py"

with io.open(path, "r", encoding="utf-8") as f:
    content = f.read()

# Fix 1: Remove debug paint_log.txt blocks and early return in paintEvent
old1 = '''        # --- rótulos orientados según el eje, sin superposición ---
        try:
            with open("paint_log.txt", "a", encoding="utf-8") as _f:
                _f.write("before rotulos\\\\n")
        except Exception:
            pass
        # self._pintar_rotulos(painter, escena["rotulos"], px, py)
        self._pintar_rotulos(painter, escena["rotulos"], px, py)
        try:
            with open("paint_log.txt", "a", encoding="utf-8") as _f:
                _f.write("after rotulos call\\\\n")
        except Exception:
            pass
        return

        # --- textos: nivel, leyenda y conteo de elementos del nivel ---'''

new1 = '''        # --- rótulos orientados según el eje, sin superposición ---
        self._pintar_rotulos(painter, escena["rotulos"], px, py)

        # --- textos: nivel, leyenda y conteo de elementos del nivel ---'''

if old1 in content:
    content = content.replace(old1, new1, 1)
    print("Fix 1: Removed debug paint_log blocks and early return")
else:
    print("Fix 1: OLD TEXT NOT FOUND")
    idx = content.find("rótulos orientados")
    if idx >= 0:
        print(repr(content[idx:idx+600]))

# Fix 2: Remove debug paint_log.txt block before legend
old2 = '''        # --- leyenda visual: muestra de color + sÃmbolo por tipo de elemento ---
        try:
            with open("paint_log.txt", "a", encoding="utf-8") as _f:
                _f.write("before legend\\\\n")
        except Exception:
            pass
        tipos_presentes = sorted({t for t in escena["tipos_losa"].values()})'''

new2 = '''        # --- leyenda visual: muestra de color + sÃmbolo por tipo de elemento ---
        tipos_presentes = sorted({t for t in escena["tipos_losa"].values()})'''

if old2 in content:
    content = content.replace(old2, new2, 1)
    print("Fix 2: Removed debug paint_log block before legend")
else:
    print("Fix 2: OLD TEXT NOT FOUND")
    idx = content.find("leyenda visual")
    if idx >= 0:
        print(repr(content[idx:idx+400]))

with io.open(path, "w", encoding="utf-8") as f:
    f.write(content)

print("Done")
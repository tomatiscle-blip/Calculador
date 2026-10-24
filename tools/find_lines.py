# Find the exact line number of 'def partes' and 'def proximo_nombre_portico'
lines = open("calc/rutas.py", encoding="utf-8").read().splitlines()
for i, ln in enumerate(lines):
    if ln.startswith("def partes") or ln.startswith("def proximo_nombre_portico") or ln.startswith("    return f"):
        print(i + 1, "|", ln[:90])
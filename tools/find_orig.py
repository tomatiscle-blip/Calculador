lines = open(r"d:\03 Ingenieria\Calculador\tools\rutas_original.txt", encoding="utf-8").read().splitlines()
for i, ln in enumerate(lines):
    if "proximo_nombre_portico" in ln or "def partes" in ln or "Pórtico" in ln or "Usa" in ln:
        print(i+1, "|", ln[:100])
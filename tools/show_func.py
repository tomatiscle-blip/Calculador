lines = open(r"d:\03 Ingenieria\Calculador\calc\rutas.py", encoding="utf-8").read().splitlines()
for i in range(626, 640):
    print(i + 1, "|", repr(lines[i]))

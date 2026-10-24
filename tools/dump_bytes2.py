lines = open("calc/rutas.py", encoding="utf-8").read().splitlines(keepends=True)
for i in range(632, 642):
    print(f"{i+1}| {lines[i]!r}")

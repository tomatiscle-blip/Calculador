import ast

for f in [
    "calc/rutas.py",
    "calc/portico.py",
    "calc/pipeline.py",
]:
    src = open(f, encoding="utf-8").read()
    try:
        ast.parse(src)
        print(f"{f}: SYNTAX OK")
    except SyntaxError as e:
        print(f"{f}: SYNTAX ERROR {e}")

import ast

for f in [
    "calc/rutas.py",
]:
    src = open(f, encoding="utf-8").read()
    try:
        ast.parse(src)
        print(f"{f}: SYNTAX OK")
    except SyntaxError as e:
        print(f"{f}: SYNTAX ERROR at line {e.lineno}: {e.msg}")
        print("CONTEXT:")
        lines = src.splitlines()
        for i in range(max(0, e.lineno - 4), min(len(lines), e.lineno + 2)):
            print(f"{i+1:4d}| {lines[i]}")

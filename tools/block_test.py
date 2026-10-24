import ast

src = open("calc/rutas.py", encoding="utf-8").read().splitlines()
block = "\n".join(src[627:645])
print("--- block ---")
print(block)
try:
    ast.parse(block)
    print("BLOCK OK")
except SyntaxError as e:
    print("BLOCK ERROR:", e.msg, "at", e.lineno, e.offset)


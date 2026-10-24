import ast

src = open("calc/rutas.py", encoding="utf-8").read()
try:
    ast.parse(src)
    print("OK")
except SyntaxError as e:
    print("msg:", e.msg)
    print("lineno:", e.lineno)
    print("offset:", e.offset)
    print("text:", repr(e.text))
    print("end_offset:", e.end_lineno)
    # print a wider region
    lines = src.splitlines()
    lo = max(0, (e.lineno or 0) - 6)
    hi = min(len(lines), (e.end_lineno or 0) + 3)
    for i in range(lo, hi):
        print(f"{i+1}| {lines[i]}")

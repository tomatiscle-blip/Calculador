import difflib

original = open(r"d:\03 Ingenieria\Calculador\tools\rutas_original.txt", encoding="utf-8").read().splitlines(keepends=True)
current = open(r"d:\03 Ingenieria\Calculador\calc\rutas.py", encoding="utf-8").read().splitlines(keepends=True)

diff = list(difflib.unified_diff(original, current, lineterm="", fromfile="git/rutas.py", tofile="calc/rutas.py"))
print("diff lines:", len(diff))
for ln in diff[:120]:
    print(ln)
if len(diff) > 120:
    print("... (truncated)")

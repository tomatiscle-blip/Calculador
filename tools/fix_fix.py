# This script fixes the syntax error in fix_rutas.py
path = "fix_rutas.py"
lines = open(path, encoding="utf-8").read().splitlines(keepends=True)

for i, ln in enumerate(lines):
    if ln.strip().startswith("if cur_lines[i].startswith"):
        lines[i] = '    if cur_lines[i].startswith("    return f\\"P" + "\\u00f3" + "rtico"):\n'
        print("replaced line", i + 1)
        break

open(path, "w", encoding="utf-8").writelines(lines)
print("done")

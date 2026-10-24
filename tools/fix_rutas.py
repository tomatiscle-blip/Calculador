# -*- coding: utf-8 -*-
"""Fix calc/rutas.py: restore 'proximo_nombre_portico' and remove mangled duplicate."""

backup_path = r"d:\03 Ingenieria\Calculador\tools\rutas_original.txt"
current_path = r"d:\03 Ingenieria\Calculador\calc\rutas.py"

with open(backup_path, encoding="utf-8") as f:
    backup_lines = f.read().splitlines(keepends=True)

with open(current_path, encoding="utf-8") as f:
    cur_lines = f.read().splitlines(keepends=True)

# ---------- Step 1: remove the mangled duplicate at module level ----------
s_dup = None
for i, ln in enumerate(cur_lines):
    if ln.startswith('    Usa "Pórtico N" con tilde:'):
        s_dup = i
        break

if s_dup is None:
    raise SystemExit("duplicate start not found")

# find the 'return f"Pórtico {numero}"' right after the duplicate start
e_dup = None
for i in range(s_dup, len(cur_lines)):
    if cur_lines[i].startswith("    return f\"P" + "\u00f3" + "rtico"):
        e_dup = i
        break

if e_dup is None:
    raise SystemExit("duplicate end not found")

print("removing module-level duplicate lines", s_dup + 1, "to", e_dup + 1)
del cur_lines[s_dup:e_dup + 1]

# ---------- Step 2: restore 'proximo_nombre_portico' from the git backup ----------
start = None
for i, ln in enumerate(backup_lines):
    if ln.startswith("def proximo_nombre_portico()"):
        start = i
        break

if start is None:
    raise SystemExit("original function not found in backup")

orig_func = "".join(backup_lines[start:])
print("original_func lines:", orig_func.count(chr(10)))

s = None
e = None
for i, ln in enumerate(cur_lines):
    if ln.startswith("def proximo_nombre_portico()"):
        s = i
    if s is not None and ln.startswith("def partes("):
        e = i
        break

if s is None or e is None:
    raise SystemExit("broken region not found")

print("replacing lines", s + 1, "to", e, "with restored function")
cur_lines[s:e] = [orig_func, "\n"]

# ---------- Write back ----------
with open(current_path, "w", encoding="utf-8") as f:
    f.writelines(cur_lines)

print("written; total lines:", len(cur_lines))

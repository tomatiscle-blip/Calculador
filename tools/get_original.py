import subprocess, os
os.chdir(r"d:\03 Ingenieria\Calculador")
out = subprocess.run(
    ["git", "show", "main:calc/rutas.py"],
    capture_output=True, text=True
)
if out.returncode != 0:
    print("ERR:", out.stderr)
else:
    with open(r"d:\03 Ingenieria\Calculador\tools\rutas_original.txt", "w", encoding="utf-8") as f:
        f.write(out.stdout)
    print("saved", len(out.stdout), "bytes")


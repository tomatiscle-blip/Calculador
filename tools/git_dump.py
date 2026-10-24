import subprocess
p = r"d:\03 Ingenieria\Calculador"
with open(p + r"\tools\git_orig.txt", "wb") as f:
    subprocess.run(["git", "-C", p, "show", "main:calc/rutas.py"], stdout=f)
print("done")
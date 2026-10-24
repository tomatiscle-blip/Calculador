# Fix the duplicated encoding in calc/rutas.py for proximo_nombre_portico()
# Also discover and correct any other double-encoded strings in the file.

path = r"d:\03 Ingenieria\Calculador\calc\rutas.py"
src = open(path, encoding="utf-8").read()

# Markers in the current (double-encoded) file
start_marker = 'def proximo_nombre_portico() -> str:'
end_marker = 'def partes(estructura: dict) -> tuple[dict, dict, dict, dict]:'

i = src.find(start_marker)
j = src.find(end_marker)
if i == -1 or j == -1 or j < i:
    raise SystemExit("markers not found")

# Extract the broken function (inclusive of its docstring + body up to before 'def partes')
broken_func = src[i:j]

# Build a corrected version of the function
corrected_func = '''def proximo_nombre_portico() -> str:
    """
    Próximo nombre libre de pórtico. Respeta el criterio viejo de P00
    (cantidad de pórticos + 1) sin pisar nombres ya usados.

    Usa "Pórtico N" con tilde: es la convención del resto del proyecto
    (la app crea "Pórtico N" y los datos existentes usan esa forma).
    """
    existentes = set(listar_porticos())
    numero = len(existentes) + 1
    while f"Pórtico {numero}" in existentes or f"Portico {numero}" in existentes:
        numero += 1
    return f"Pórtico {numero}"

'''

new_src = src[:i] + corrected_func + src[j:]

# Optional: scan for other double-encoded words and fix them (in-place)
# Patterns like 'PÃ³rtico' (encoding of 'Pórtico') etc.
import re

pairs = [
    ("Ã³", "ó"),
    ("Ã©", "é"),
    ("Ã±", "ñ"),
    ("Ã¢", "â"),
    ("Ã¨", "è"),
    ("Ã¯", "ï"),
    ("Ã²", "ò"),
    ("Ãº", "ú"),
]

fixed = new_src
for bad, good in pairs:
    fixed = fixed.replace(bad, good)

# Write back
with open(path, "w", encoding="utf-8") as f:
    f.write(fixed)

print("done; new length:", len(fixed))

import glob

KEYS = ("portico", "construir", "calcular", "cargar_estructura", "columnas", "estructura")

for f in sorted(glob.glob("tests/test_*.py")):
    with open(f, encoding="utf-8") as fh:
        lines = fh.readlines()
    print("===", f, f"({len(lines)} lines) ===")
    for i, linea in enumerate(lines, 1):
        if any(k in linea for k in KEYS):
            print(f"{i}: {linea.rstrip()}")
    print()

import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent
if str(RAIZ) not in sys.path:
    sys.path.insert(0, str(RAIZ))
from calc import rutas  # noqa: E402
estructura = rutas.cargar_estructura() or {}

print("=== Pórtico 7 (completo) ===")
p7 = estructura.get('Pórtico 7', {})
for k, v in p7.items():
    print(f"  {k} = {v}")

print("\n=== Muro 7 (completo) ===")
m7 = estructura.get('Muro 7', {})
for k, v in m7.items():
    print(f"  {k} = {v}")

print("\n=== Todos los ejes ===")
ejes = estructura.get('ejes', {})
if isinstance(ejes, dict):
    for k, v in ejes.items():
        print(f"  {k} = {v}")
else:
    print("  No hay clave 'ejes' en estructura")

print("\n=== Claves de estructura ===")
print(list(estructura.keys()))
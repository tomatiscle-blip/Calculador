"""Calculo consolidado de columnas globales (Mx/My por ID global)."""
from __future__ import annotations
import sys
from pathlib import Path
_RAIZ = Path(__file__).resolve().parent.parent
if str(_RAIZ) not in sys.path:
    sys.path.insert(0, str(_RAIZ))

try:
    from . import columnas
except ImportError:
    from calc import columnas


def main() -> None:
    resultados = columnas.calcular_global()
    print(f"Columnas globales calculadas: {len(resultados)}")


if __name__ == "__main__":
    main()
"""Calculo consolidado de bases globales (Mx/My por ID global)."""
from __future__ import annotations
import sys
from pathlib import Path
_RAIZ = Path(__file__).resolve().parent.parent
if str(_RAIZ) not in sys.path:
    sys.path.insert(0, str(_RAIZ))

try:
    from . import bases
except ImportError:
    from calc import bases


def main() -> None:
    resultados = bases.calcular_global()
    print(f"Bases globales calculadas: {len(resultados)}")


if __name__ == "__main__":
    main()
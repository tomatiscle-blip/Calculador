"""Cálculo consolidado de columnas globales (Mx/My por ID global)."""
from __future__ import annotations

try:
    from . import columnas
except ImportError:  # ejecución directa como script
    from calc import columnas


def main() -> None:
    resultados = columnas.calcular_global()
    print(f"Columnas globales calculadas: {len(resultados)}")


if __name__ == "__main__":
    main()

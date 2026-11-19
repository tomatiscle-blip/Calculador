"""Formato numérico para informes y exportaciones orientadas a Excel."""

from __future__ import annotations

from numbers import Real


def numero_local(valor: Real, decimales: int = 2) -> str:
    """Devuelve un número con punto de miles y coma decimal."""
    return f"{valor:,.{decimales}f}".replace(",", "_").replace(".", ",").replace("_", ".")


def valor_csv(valor: object, decimales: int = 2) -> object:
    """Localiza números para CSV sin alterar texto ni valores vacíos."""
    if isinstance(valor, float):
        return numero_local(valor, decimales)
    return valor


def parse_numero(valor: object, defecto: float = 0.0) -> float:
    """Lee números antiguos con punto decimal y los nuevos con coma decimal."""
    texto = "" if valor is None else str(valor).strip()
    if not texto:
        return defecto
    if "," in texto:
        texto = texto.replace(".", "").replace(",", ".")
    return float(texto)

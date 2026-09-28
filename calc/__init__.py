"""
calc — Núcleo de cálculo del proyecto Calculador.

Idea central (regla de oro del proyecto):

    El cálculo NO sabe que existe una pantalla, y la pantalla NO sabe
    ninguna fórmula.

Cada módulo de este paquete:
    * no llama a input() ni imprime carteles,
    * no depende del directorio desde el que se ejecuta el programa,
    * recibe datos y devuelve datos (diccionarios / listas),
    * deja los archivos en su lugar usando calc.rutas.

Así el mismo código sirve para:
    - la app de escritorio (PySide6),
    - los scripts de consola viejos (P00..P06, L00, C00, V0x),
    - pruebas automáticas.

Módulos:
    rutas       -> dónde está cada archivo (§ única fuente de verdad de rutas)
    pipeline    -> orden de las etapas, dependencias entre ellas y estado (semáforo)
"""

__all__ = ["rutas", "pipeline"]

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
    materiales  -> la biblioteca del proyecto (datos/materiales.json)
    cargas      -> de los elementos que reciben carga a las cargas D/L/W y las combinaciones
    losas       -> el cálculo de una losa alivianada (por tipología)
    losas_macizas -> análisis elástico inicial de paños macizos unidireccionales
    portico     -> el MOTOR: resuelve el pórtico 2D y devuelve las solicitaciones (M, V, N)
"""

__all__ = ["rutas", "pipeline", "materiales", "cargas", "losas", "losas_macizas", "portico"]

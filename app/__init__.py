"""
app — Interfaz gráfica del Calculador (ventana de Windows, PySide6).

Regla de oro del proyecto: la pantalla NO calcula nada.
Solo:
  * lee el estado del proyecto (calc.pipeline),
  * muestra los resultados que ya existen en `salidas/`,
  * abre los archivos con el programa que corresponda,
  * lanza las etapas que ya se pueden ejecutar sin escribir datos.

Para abrirla: doble clic en `calculador.bat` (o `py -m app`).
"""

__all__ = ["principal"]

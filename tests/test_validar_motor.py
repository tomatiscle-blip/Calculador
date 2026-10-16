"""
ARCHIVADO (16/10/2026): control histórico, ya no se corre.

Qué era
-------
El primer control de la migración del motor: se hizo ANTES de decidir si
Pynite reemplazaba a anaStruct (el P01). Leía un pórtico de
datos/estructura.json y lo resolvía por dos caminos independientes:

  1. anaStruct FRESCO, cargando a mano D/L/W sobre los tramos (igual que P01);
  2. calc.portico (Pynite), con el mismo pórtico.

Comparaba por combinación reacciones de base, M/V de la viga y M/N de la
columna (tolerancia 0,05 con malla 50). Concluyó que los dos motores daban
lo mismo (la única diferencia, ~0,04 kN·m, era el momento de campo por la
malla) y con eso se decidió migrar a Pynite.

Por qué no corre más
--------------------
Comparaba contra el pórtico "Portico 1" con las cargas D_total/L_total/
W_total guardadas DENTRO de cada tramo (el esquema de datos de P00/P01).
Las obras actuales ya no guardan así las cargas: viven en datos/cargas.json
como aplicaciones sobre tramos, así que el control quedó sin datos con los
que trabajar. Replicarlo al modelo actual no aporta nada nuevo.

Qué lo reemplaza
----------------
* tools/comparar_motores.py — comparación permanente entre los resultados
  anaStruct guardados en estructura.json y los que calcula Pynite.
* tools/regresion.py + tests/golden/ — red de seguridad de resultados.
"""

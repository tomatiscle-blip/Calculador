# README – Script columnas.py

## Descripción
El script `columnas.py` permite calcular y registrar resultados de columnas de hormigón armado (axiales, flexocompresión, cortas o esbeltas) en una planilla CSV compacta y auditable.  
Incluye la propuesta de armadura longitudinal y de confinamiento (estribos o sunchos), junto con parámetros normativos y notas técnicas.

## Funcionalidad principal
- Entrada de datos:
  - ID de la columna
  - Pu (carga axial en kN)
  - Mu (momento en kNm)
  - f’c (resistencia del hormigón en MPa)
  - Dimensiones b, h (cm)
  - d’ (recubrimiento en cm)
  - Tipo de columna (estribos o sunchos)
  - Tipo de sección (rectangular o circular)
  - Altura libre (m)

- Cálculo automático:
  - Esbeltez λ y clasificación (CORTA / ESBELTA)
  - Amplificación de momento si corresponde
  - Cuantía geométrica ρg
  - Ast teórico, mínimo normativo y adoptado
  - Propuesta de armadura longitudinal (número y diámetro de barras)
  - Verificación de estribos o sunchos (diámetro, separación/paso, cumplimiento normativo)

- Salida compacta en CSV:
  Cada fila contiene:
  ID;Pu(kN);Mu(kNm);f'c(MPa);b(cm);h(cm);λ;Clasificación;Ast(cm²);Armadura;Estribos/Sunchos;Nota_diagrama

## Ejemplo de salida
ID;Pu(kN);Mu(kNm);f'c(MPa);b(cm);h(cm);λ;Clasificación;Ast(cm²);Armadura;Estribos/Sunchos;Nota_diagrama  
C0-1;2113.0;211.3;25;35.0;50.0;32.91;CORTA;38.50;8Ø25;Ø8 c/30.0cm;γ=0.84, X=2.415, Y=12.074, ρg=0.022  
C0-4;481.0;31.05;20;30.0;30.0;43.63;ESBELTA;9.00;8Ø12;Ø10 c/8.2cm (cumple=False);γ=0.73, X=1.144, Y=5.344, ρg=0.010  

## Archivos generados
- `salidas/planilla_columnas.csv` → planilla compacta en formato CSV UTF‑8.
- Opcional: puede exportarse a ODS/XLSX desde LibreOffice Calc para visualización avanzada.

## Uso
1. Ejecutar el script `columnas.py`.
2. Ingresar los datos solicitados por consola.
3. El programa calcula y guarda automáticamente la fila en `planilla_columnas.csv`.
4. Abrir el archivo en LibreOffice Calc, WPS o Excel para revisar.

## Notas
- El CSV se guarda en UTF‑8 para soportar símbolos técnicos (Ø, λ, γ, ρg).
- El campo *Estribos/Sunchos* se adapta según el tipo de columna.
- La planilla es compacta y auditable, sin bloques de texto largos.
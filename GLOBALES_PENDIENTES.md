# Mejoras pendientes: entidades globales y cálculo consolidado

## Contexto
Las columnas y bases son entidades globales por coordenadas (`x`, nivel). Cada pórtico guarda listas de IDs globales (`columnas`, `bases`). Esto evita duplicación de datos geométricos.

## Problema actual
El dimensionado de columnas se ejecuta por pórtico y asume que cada pórtico puede resolver sus columnas aisladamente. Si una columna global es compartida por varios pórticos, cada uno puede aportar solicitaciones distintas. El cálculo actual no consolida esas solicitaciones antes de dimensionar.

Bases pueden ser independientes, pero si comparten coordenadas deberían acumular cargas de varios pórticos.

## Objetivo
Permitir que columnas/bases globales se dimensionen con solicitaciones consolidadas de todos los pórticos que las usan.

## Etapas propuestas

### Etapa 1: Solicitaciones por columna global
- Modificar el motor de pórticos para guardar solicitaciones por ID global de columna/base, además del resultado por pórtico.
- Archivo sugerido: `calc/portico.py` y `pipeline` de etapas.

### Etapa 2: Consolidación antes de dimensionar columnas
- En la etapa de columnas, recorrer `estructura["columnas"]`.
- Para cada columna global, reunir solicitaciones de todos los pórticos que la referencian.
- Combinar (sumar/escalar según criterio) las solicitaciones consolidadas.

### Etapa 3: Dimensionado global de columnas
- Dimensionar cada columna global una sola vez con solicitaciones consolidadas.
- Guardar resultado global en `salidas/columnas/` y referenciarlo desde cada pórtico.

### Etapa 4: Bases
- Aplicar el mismo criterio de consolidación para bases globales compartidas.

### Etapa 5: UI y pipeline
- Ajustar `app/principal.py` para que la etapa de columnas no se ejecute “plana” por pórtico cuando haya compartición.
- Mostrar en la UI qué columnas/bases son compartidas y su estado.

## Criterio de aceptación
- Una columna compartida por dos pórticos se dimensiona una sola vez.
- Las solicitaciones de ambos pórticos se reflejan en el resultado.
- No se duplican entidades globales.

## Nombres legibles para humanos
- El ID global (`COL-002`) es clave interna.
- Se puede agregar un campo `nombre` legible basado en eje + nivel (ej. `C0-1`, `Eje A - Nivel 0`).
- La UI y resultados deberían mostrar `nombre` + coordenadas, no solo el ID.

## Plan para mañana — portado de columnas/bases a Python
- Bloque 1: extraer predimensión y parámetros normativos del legacy `P04`.
- Bloque 2: portar bucle principal de dimensionado a `calc/columnas.py`.
- Bloque 3: adaptar escritura de resultados/CSV y memoria.
- Bloque 4: cablear UI (`app/principal.py`) a `calcular_global()`.
- Repetir para bases (`calc/bases.py` inexistente) con datos de terreno únicos.

## Viento y ancho activo en planta (futuro)
- Definir dirección de viento por proyecto y porticos afectados (según orientación).
- Calcular ancho activo por eje a partir de la planta y ejes.
- Aplicar cargas de viento solo a pórticos con el sentido correspondiente.
- Si una columna recibe excesivo momento por viento, evaluar agregar pórtico transversal para compensar.
- Integrar con columnas globales para consolidar solicitaciones de viento.
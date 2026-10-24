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

## Actualización: UI y planilla local (Octubre 2026)
- Pestaña "Columnas globales" habilitada cuando todos los porticos que comparten la columna estan resueltos.
- En la planilla de columnas por portico, si una columna es global se calcula igual pero se avisa que debe revisarse en Columnas globales; el dato local no vale, vale el global.
- Preferido: el listado de columnas por portico se acumula/actualiza en columnas globales, quedando siempre activo y una planilla de todas las columnas de la obra.
- Bases siempre globales: acumular por BASE-XXX y flag momentos_otro_sentido para activar/desactivar contribucion transversal.

## Actualización: identidad global por posición y momentos por eje

### Implementado
- Columnas globales identificadas por coordenada global `(X, Y, nivel)`, no por `(x, nivel)` local del pórtico.
- Numeración por piso: `C0-x` para PB, `C1-x` para N1, etc.
- Catálogo global (`estructura["columnas"]`) guarda `X`, `Y`, `x`, `nivel`, `eje_id`, `eje_nombre`, `nivel_inicio`, `piso`, `momentos_otro_sentido`.
- `obtener_o_crear_columna_global(..., piso=...)` reutiliza columna global si ya existe en esa posición y piso.
- Momentos por eje consolidados: las columnas globales acumulan `Mx`, `My`, `Mx_otro_sentido`, `My_otro_sentido`.
- Opción en UI para eliminar la obra actual y crear una obra nueva (limpia `estructura.json`, `cargas.json` y salidas).

### Pendiente
- El visor de pórtico (`app/visualizacion.py::VistaPortico`) usa `columna["x"]` (coordenada local del pórtico que creó la columna). En pórticos transversales que reutilizan columnas globales, esa coordenada no corresponde y el dibujo colapsa (línea vertical, vigas como puntos).
- Corregir `_geometria_tramos()` para calcular la coordenada local de cada columna a partir de su `X/Y` global y la `referencia_planta` + `posicion_planta_m` del pórtico activo, en lugar de `columna["x"]`.
- Podría ser también el armado de estructura del pórtico (asociación de columnas) lo que deba normalizar la coordenada local por pórtico.
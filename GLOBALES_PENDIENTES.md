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
- Bloque 1: extraer predimensión y parámetros normativos del legacy `P04`.  ✅ IMPLEMENTADO en `calc/columnas.py` (`definir_fck`, `predimensionar_seccion`, `predimensionar_suncho`, `PHI`, `ALPHA`, `cargar_diagrama`).
- Bloque 2: portar bucle principal de dimensionado a `calc/columnas.py`.  (en progreso: predimensión + fck + b/h integrados; falta armadura/diagramas y memoria)
- Bloque 3: adaptar escritura de resultados/CSV y memoria.  (pendiente)
- Bloque 4: cablear UI (`app/principal.py`) a `calcular_global()`.  (ya cableado; falta dimensionado real)
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
\n## Actualización: vigas globales (Octubre 2026)\n\n### Decisión de acumulado\n- Las vigas globales se acumulan por tramo (cada viga física de cada pórtico), no se deduplican por elemento estructural compartido.\n- El ID global de tramo conserva el nombre del pórtico + tramo (ej. `Portico 2.V-01.T1`), por lo que no se duplican nombres: el pórtico y el tramo identifican unívocamente cada fila.\n- El nivel/piso se deriva del identificador de viga (`V{piso}-...`) cuando está disponible.\n- El cálculo es manual por botón (igual que columnas/bases globales); no se recalcula automáticamente al cambiar cargas.\n\n### Plan implementado\n- `calc/vigas_globales.py` (nuevo): recorre `salidas/solicitaciones/*.json`, acumula envolvente de vigas y voladizos por pórtico/tramo, calcula M/V básicos, exporta `salidas/vigas/resultados_vigas_globales.json`.\n- `calc/pipeline.py`: nueva etapa `vigas_globales` dependiente de `vigas` y `columnas_globales`.\n- `app/principal.py`: pestaña \"Vigas globales\" con tabla (ID, Pórtico, Nivel, M, V, Estado) y botón de cálculo.\n"
## Actualización: bases globales con Mx/My (Octubre 2026)

### Implementado
- `calc/bases.py::calcular_global()` ahora acumula solicitaciones por `BASE-XXX` y descompone el momento en plano `Mz` de cada pórtico en `Mx`/`My` según `referencia_planta.direccion` (análogo a columnas globales).
- Se generan `Mx`, `My`, `M` (máximo en módulo) y se conservan `Fx`, `Fy`, porticos y `momentos_otro_sentido`.
- Nuevo wrapper `calc/bases_globales.py` para ejecución directa.
- Pipeline/UI ya exponen la etapa "Bases globales" (`resultados_bases_globales.json`).

### Pendiente
- Integrar el visor de pórtico para bases transversales (coordenada local por pórtico activo).
- Validar dimensionado de zapatas con momentos bidireccionales consolidados.

## Actualización: reorganización UI y dimensionado global con Bresler (Octubre 2026)

### Decisión de diseño
- Las pestañas "Columnas globales", "Vigas globales" y "Bases globales" **deben salir de la subpestaña de pórticos** y ubicarse al final del proceso (sección de consolidación), porque dependen de que todos los pórticos estén resueltos.
- Orden natural: geometría → cargas → pórticos → consolidación global → dimensionado definitivo.
- Las columnas no compartidas (un solo pórtico) ya tienen su dimensionado definitivo en la etapa "7. Columnas" (P04, unidireccional). Solo las compartidas deben re-dimensionarse con datos consolidados.

### Pendiente de implementación
1. **Reorganización UI** (`app/principal.py`):
   - Mover las pestañas globales fuera de `_pagina_portico` a una sección final (nueva pestaña "Consolidación" o al final de "Estado y etapas").
   - Habilitarlas solo cuando todos los pórticos estén resueltos.
2. **Dimensionado global de columnas con Bresler**:
   - Migrar el bucle principal de `Legacy/P04_Columnas_portico.py` a `calc/columnas.py` (Bloque 2 pendiente).
   - Para cada columna global compartida con `Mx` y `My` ambos no-cero, usar `calc/columnas_biaxial.py` (`verificar_biaxial` / `bresler_pn`) en lugar de `Mu = max(Mx, My)`.
   - Columnas con un solo momento no-cero → dimensionado unidireccional (caso actual).
   - Agregar botón "Dimensionar columnas globales" en la pestaña global que ejecute este cálculo.
3. **Bases globales**: aplicar el mismo criterio bidireccional consolidado al dimensionado de zapatas (cuando se porte a `calc/bases.py`).
4. **Vigas globales**: mantener acumulado por tramo; no requiere cambio de método, solo reubicación UI.

### Criterio de aceptación
- Las pestañas globales se calculan al final, con todos los pórticos resueltos.
- Una columna compartida con Mx y My no-cero se dimensiona con Bresler.
- Las columnas no compartidas quedan definitivas con el dato del pórtico.
- No se duplican entidades globales.
## Actualización: propuesta de pestaña "Estructura global" (Octubre 2026)

### Decisión de diseño
- Crear una pestaña de nivel superior **"Estructura global"** con subpestañas internas:
  - Columnas globales
  - Vigas globales
  - Bases globales
- Sacar las pestañas globales de la subpestaña de pórticos y moverlas ahí.
- Administrar toda la estructura global desde un solo lugar, separado del flujo por pórtico.
- Habilitar/deshabilitarlas según el estado de resolución de los pórticos.
## Actualización: pestaña "Estructura global" (UI) — Octubre 2026

### Implementado
- Nueva pestaña de nivel superior "Estructura global" con subpestañas:
  - Columnas globales
  - Vigas globales
  - Bases globales (nueva, con botón de cálculo)
- Las pestañas globales se sacaron de la subpestaña de pórticos.
- `calc/bases.calcular_global()` ya disponible desde la UI.

### Pendiente
- Habilitar/deshabilitar "Estructura global" según resolución de todos los pórticos.
- Dimensionado global de columnas con Bresler (pendiente Bloque 2).
- Bases globales bidireccional definitivo.



## Bug reportado (Octubre 2026): diagrama de momentos en pórticos con apoyos intermedios
- En `calc/portico.py`, las columnas se agregan como nudos de la viga solo si `abs(_k(c.get("nivel")) - y_viga) <= 1e-4` y su x cae dentro del tramo.
- Las columnas centrales cuyo nivel no coincide exactamente con la cota de la viga (redondeo) o cuyo x cae fuera del rango por tolerancia quedan fuera, por lo que el diagrama no muestra el quiebre en ese apoyo.
- Pendiente: revisar tolerancia y la fuente de `columnas` del pórtico (IDs globales vs. datos locales) para confirmar que las columnas intermedias estén incluidas.

### Validación
- El ancho tributario automático como mitad de separaciones entre pórticos vecinos es válido para pórticos con el mismo plano de aplicación.
- La dirección del viento se define a nivel de obra (X o Y); solo reciben carga los pórticos cuya normal coincide con esa dirección.
- Pórticos paralelos a la dirección de viento no reciben carga directa en ese sentido.

### Implementado (base)
- `calc/cargas.py`: `porticos_afectados_por_viento(estructura, direccion_viento)` y `anchos_tributarios_viento(estructura, direccion_viento)`.

### Pendiente
- Guardar `direccion_viento` en `datos/cargas.json`.
- Aplicar automáticamente los anchos a los pórticos afectados desde la UI de cargas/viento.
- Considerar bordes de obra y pórticos irregulares (distancia real, no módulo uniforme).
## Pendiente (Octubre 2026): viento clásico por pórtico
- Activación manual: checkbox "Aplicar viento" por pórtico en la pestaña de cargas.
- Ancho tributario: siempre manual (no calcular automáticamente desde ejes), para evitar cálculos excesivos.
- El viento se aplica solo a los pórticos marcados, independientemente de la dirección automática.
- Quitar o ignorar la regla automática por dirección (porticos_afectados_por_viento) cuando se use el modo manual.

## Bug reportado (Octubre 2026): diagrama de momentos en pórticos con apoyos intermedios — RESUELTO
- Causa: en `calc/portico.py` ~línea 710, al agregar columnas como nudos de la viga se comparaba `nivel` (cota inferior) con `y_viga` (cota superior = nivel+altura), por lo que ninguna columna interior coincidía y quedaba fuera del diagrama.
- Fix aplicado: comparar `nivel + altura_m` con `y_viga`, consistente con el filtro de `cols` (~línea 650).
- Solución general (no manual): aplica a todos los pórticos/obras automáticamente.


## Actualizacion: Fase 1 - Persistencia de dimensionado por portico (Octubre 2026)

### Implementado
- Nuevo modulo `calc/columnas_persistencia.py` con `generar_dimensionados_json()`.
- Genera, por cada columna global, un archivo `salidas/columnas/dimensionado_{GID}.json`.
- Cada archivo acumula, por portico que usa la columna, los datos de dimensionado uniaxial: Pu, Mu, b, h, fck, Ast, armadura longitudinal, estribos/sunchos, gamma, X, Y, rho, clasificacion, lambda, altura libre, recubrimiento.
- No modifica el calculo por portico ni la UI; solo agrega persistencia estructurada.
- Ejecucion: `python calc/columnas_persistencia.py` (post-dimensionado de columnas).

### Pendiente (Fases 3-4)
- Unificacion de seccion/armadura por columna global.
- Planilla global unica (una fila por columna fisica).
- Verificacion biaxial (Bresler) para columnas con Mx y My no-cero.

## Actualizacion: Fase 2 - Unificacion por columna global (Octubre 2026)

### Implementado
- Nuevo modulo `calc/columnas_unificacion.py` con `unificar()`.
- Genera `salidas/columnas/unificado_{GID}.json` con seccion/armadura candidata unificada por GID.
- Genera `salidas/columnas/planilla_columnas_globales.csv` (una fila por columna fisica).
- Criterio de unificacion: seccion del portico con mayor Ast (desempate mayor area b*h); fck maximo; armadura del portico con mayor Ast; estribos mayor diametro / menor paso; clasificacion ESBELTA si algun portico lo es.
- Incluye solicitaciones globales N, Mx, My, Mu_max y flag biaxial (sin verificacion Bresler).
- Ejecucion: `python calc/columnas_unificacion.py` (requiere Fase 1 ejecutada).

### Pendiente (Fase 4)
- Verificacion biaxial (Bresler) para columnas con Mx y My no-cero.

## Actualizacion: Fase 4 - Verificacion biaxial global (Octubre 2026)

### Implementado
- Nuevo modulo `calc/columnas_biaxial_global.py` con `ejecutar()`.
- Para cada columna global unificada:
  - Si es uniaxial: se marca como RESUELTA (uniaxial).
  - Si es biaxial: calcula Pnx y Pny (capacidad uniaxial por eje via diagrama de interaccion), P0 (carga axial pura) y ejecuta `verificar_biaxial()` (Bresler).
- Genera `salidas/columnas/resultado_final_{GID}.json` con el resultado final de cada columna global.
- Ejecucion: `python calc/columnas_biaxial_global.py` (requiere Fase 2 ejecutada).

## Actualizacion: Planilla global final (Octubre 2026)

### Implementado
- Nuevo modulo `calc/columnas_planilla_final.py` con `generar()`.
- Genera `salidas/columnas/planilla_columnas_globales_final.csv` con una fila por columna global, incluyendo resultados de verificacion biaxial (Pnx, Pny, P0, Pn_biaxial, Cumple) para las biaxiales.
- Ejecucion: `python calc/columnas_planilla_final.py` (requiere Fases 2 y 4 ejecutadas).

## Actualizacion final: integracion UI + pipeline (Octubre 2026)

### Implementado
- Boton "Calcular columnas globales (completo)" en la pestana Estructura global -> Columnas globales (app/principal.py).
- Tabla de planilla global final en la misma pestana (20 columnas, incluye Pnx/Pny/P0/Cumple).
- Etapa "columnas_globales_completo" agregada a calc/pipeline.py (orquesta Fases 1->2->4->planilla).
- Orquestador calc/columnas_globales_completo.py para ejecucion manual.
- Tests: tests/test_columnas_globales_fases.py (5 tests OK).

### Estado
- Columnas globales: Fases 1, 2, 4 + planilla final + UI + pipeline COMPLETAS.
- Bases globales: pendiente (siguiente etapa).

## Actualizacion: Viento manual por portico (Octubre 2026)

### Implementado
- Modo manual de viento por portico en calc/portico.py y script.py.
- UI en app/cargas_proyecto.py: tabla de porticos con checkbox 'Aplicar viento' y ancho tributario manual.
- Visor de planta (app/visualizacion.py) muestra viento con ancho por portico o fallback.
- Fix: presion_viento faltante en calc/cargas.py.
- Fix: viento aplica FX o FY segun direccion del portico.
- Fix: visor muestra viento si esta activo globalmente.

### Notas
- Viento se aplica a todos los porticos si esta activo (sin filtrar por direccion en el calculo).
- Visor muestra flecha en la direccion configurada.

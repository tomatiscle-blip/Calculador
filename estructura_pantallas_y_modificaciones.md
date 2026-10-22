# Estructura de ventanas y pestañas — CALCULADOR (app principal PySide6)

Archivo de referencia para dar indicaciones. La app se abre con `calculador.bat` (`pyw -m app`). Ventana principal: **"Calculador — Estructuras"** (QMainWindow). Tamaño inicial 1340×820, mínimo 1050×680.

## Banda superior (siempre visible)

- Título de la aplicación (CALCULADOR).
- Combo **"Obra"** (selector de obra activa) + botón **"Nueva obra"** + **"Actualizar"**.
- Etiqueta resumen de etapas del pórtico activo.
- (Los datos de carátula se editan en la pestaña Inicio.)

## Ventana principal: pestañas (QTabWidget, documentMode)

### 1) INICIO
Pantalla carátula + navegación rápida. Visión global de toda la obra.

**Bloques:**
- **"Carátula de la obra"**: nombre, ID, lugar, comitente, responsable, notas. Botón **"Guardar"**.
- **"Pórticos y muros de esta obra"**: botones **"Crear pórtico"**, **"Eliminar pórtico / muro"**.
- **"Muros de esta obra"**: botón **"Crear muro"**.
- **"Accesos rápidos"**: botones para ir a Cargas, Asignar cargas a barras, Losa desde ejes.
- **"Cargas del proyecto"**: etiqueta resumen + botón **"Ver cargas"**.
- **"Solicitaciones del pórtico"**: etiqueta + botones **"Resolver motor"**, **"Ver JSON"**, **"Diagrama de momentos"**.
- **"Planta de la obra"**: `grupo_planta_inicio` con selector **"Nivel visible"** (Todos los niveles / PB / N1 / …) en fila compacta ARRIBA del visor, y debajo el visor de planta por nivel (`VistaPlantaPorticos` dentro de `scroll_vista_portico_p`), que ocupa todo el ancho restante. Incluye la tabla de resumen del pórtico activo (`tabla_inicio`). No hay solapas ni pestañas internas dentro de Inicio.
- **Panel lateral derecho**: datos de ejes/niveles y referencias (`resumen_ejes`, `resumen_niveles`), botones **"Definir ejes"**, **"Definir niveles"**, **"Asociar pórticos"**, referencia de cargas.
- También muestra referencias de cargas y botón para ver el último análisis.

### 2) CARGAS (`PaginaCargas`)
Pestaña con 2 solapas (QTabWidget `"pestanias_cargas"`):

- **a) "Definir elementos"** (scroll): catálogo de elementos de carga de la obra.
  - Tabla de elementos: ID, Descripción, Aplicación, Transferencia a apoyos, Categoría, D, L, Incluir, Aplicaciones manuales.
  - Botones: Nueva carga, Nueva losa desde ejes, Guardar cargas, Editar composición, Eliminar elemento, Aplicar a tramos, Editar composición del elemento, Quitar asignación, Diseñar losa, Editar paño de losa, Generar informe, Abrir informe.
  - Grupo **"Aplicar y revisar cargas en barras"**: selector de elemento, tabla de aplicaciones (Pórtico(s), Carga, Tipo, Barra/tramo, Desde/Hasta, Ancho).
- **b) "Elementos y asignaciones"** (scroll con catálogo+aplicaciones).

**Diálogos asociados:**
- `EditorElemento`: **"Componer carga — {nombre}"**.
- `EditorAplicaciones`: **"Aplicaciones de {carga_id}"**.

### 3) PÓRTICOS
- Combo **"Pórtico activo"** (pórtico seleccionado para toda la app).
- **"Esquema del pórtico seleccionado"**: `VistaPortico`.
- Tabla de resumen del pórtico (`tabla_portico_resumen`).
- Selector de nivel visible (`combo_nivel_planta`) para la planta del pórtico.

### 4) 1 · ESTADO Y ETAPAS
- Tabla de etapas (semáforo: ok / desactualizada / pendiente / sin_datos / a_desarrollar).
- Detalle de la etapa seleccionada + botones **"Ejecutar"** y **"Abrir salida"**.

### 5) 2 · VIGAS
- Etiqueta resumen.
- Botón **"Dimensionar vigas"**.
- Tabla de vigas.
- Panel de planillas: lista de planillas TXT + visor de texto (splitter).

### 6) 3 · COLUMNAS
- Etiqueta resumen.
- Botón **"Dimensionar columnas"**.
- Tabla de columnas.
- Panel de planillas: lista + visor (splitter).

### 7) 4 · BASES
- Etiqueta resumen.
- Botón **"Dimensionar bases"**.
- Tabla de bases.
- Panel de planillas: lista + visor (splitter).

### 8) 5 · LOSAS
- Etiqueta resumen.
- Botón **"Dimensionar losas"**.
- Tabla de losas.

### 9) 6 · MUROS
- Tabla de muros de la obra.
- Botones para crear/editar muros.
- Visor de muro (`VistaPortico`).

### 10) 7 · ARCHIVOS
- Gestión de archivos del proyecto.
- Importar/exportar datos.

---
> Si existe una implementación previa en `tests/`, `tools/`, archivos temporales, `_tmp_diff_app.txt`, `_test_*.py`: analizar primero si contiene funcionalidad reutilizable.
> No reimplementar una funcionalidad que ya exista.
> Si se encuentra una implementación funcional de `VistaPlantaPorticos`, `VistaPortico`, pestaña Pórticos o reorganización de Inicio: integrarla en la aplicación principal antes de crear una nueva versión.

---

## Cambio de interfaz — completado (implementación real)

### Objetivo
Reorganizar la interfaz principal para separar claramente:
- Vista global de la obra.
- Vista particular del pórtico seleccionado.

### Importante
Cambios aplicados sobre la aplicación real ejecutada mediante `calculador.bat`. Se integraron los componentes existentes (no se crearon tests, demos, ventanas temporales ni prototipos). Los cambios quedan visibles al ejecutar `calculador.bat`.

### Estado actual (ya implementado)

La pestaña Inicio contiene ahora:
- `VistaPlantaPorticos` integrada en el grupo **"Planta de la obra"** con `scroll_vista_portico_p` y `combo_nivel_planta`.
- `tabla_inicio` para el resumen del pórtico activo.
- Referencias de ejes y niveles actualizadas.

La pestaña **"Pórticos"** fue creada con:
- Combo **"Pórtico activo"**.
- `VistaPortico` (esquema del pórtico seleccionado).
- Controles asociados al pórtico seleccionado.

### Archivos modificados
- `app/principal.py`: integración de `VistaPlantaPorticos`, `combo_nivel_planta`, `tabla_inicio`, pestaña "Pórticos", métodos `_cargar_inicio`, `_cambiar_nivel_planta`, `_pagina_porticos`, `_llenar`.

### Clases modificadas
- `VentanaPrincipal` (en `app/principal.py`): reorganización de la pestaña Inicio y creación de la pestaña Pórticos.

### Métodos modificados / agregados
- `__init__`: instancia de `VistaPlantaPorticos`, `scroll_vista_portico_p`, `combo_nivel_planta`, `tabla_inicio`; `addTab(..., "Pórticos")`.
- `_pagina_inicio`: `grupo_planta_inicio` con selector de nivel en fila compacta arriba del visor; visor ocupa todo el ancho restante. Sin solapas internas.
- `_cargar_inicio`: repuebla combo de niveles conservando nivel visible; la planta se dibuja con TODOS los ejes y niveles de la obra.
- `_cambiar_nivel_planta`: actualización de nivel visible en `VistaPlantaPorticos`.
- `_pagina_porticos`: construcción de la pestaña "Pórticos" con combo, `VistaPortico` y tabla de resumen.
- `_llenar`: método genérico para llenar `QTableWidget` con filas (celdas, `esta_mal`).
- `refrescar`: relee la estructura desde disco (`rutas.cargar_estructura`) y repinta la planta.

### Verificación obligatoria
Antes de finalizar informar:
- archivos modificados;
- clases modificadas;
- métodos modificados.

`calculador.bat` ejecutado con código 0 (sin tracebacks). `VentanaPrincipal()` instanciada sin errores.

- `VistaPlantaPorticos._encuadre`: encuadre de la planta separado del `paintEvent` (márgenes 78/62 px para las burbujas de eje).
- `VistaPlantaPorticos._colocar_rotulos`: cálculo de posición de rótulos con evasión de solapamiento en espiral y ajuste dentro del widget.
- `VistaPlantaPorticos._pintar_rotulos`: dibuja solo los rótulos colocados.
- `VistaPlantaPorticos._escena`: los nombres de columnas y bases pasan al sistema de rótulos.
- `VistaPlantaPorticos.paintEvent`: burbujas de eje con radio según el texto y escalonado cuando dos ejes están juntos.


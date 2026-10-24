# # INTERFAZ.md

Este documento describe la implementación real de la interfaz gráfica del calculador.

Contiene:
- Pestañas y navegación.
- Widgets principales.
- Ventanas y diálogos.
- Clases y métodos relevantes de la GUI.
- Restricciones de implementación de la interfaz.

No describe cálculos ni arquitectura general del sistema: esos temas se documentan en `ARQUITECTURA.md`.

> Si el cambio afecta navegación, pestañas, widgets, diálogos o visualización, consultar este documento antes de modificar código.

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
- Subpestañas integradas: **"Pórtico"**, **"Vigas"**, **"Columnas"** y **"Bases"**.
- **"Pórtico"**: `VistaPortico` (esquema), resumen y envolvente del pórtico seleccionado.
- **"Vigas" / "Columnas" / "Bases"**: tablas, botones de dimensionado y planillas correspondientes.

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
## Reglas para modificar la interfaz

1. Buscar primero componentes existentes antes de crear nuevos.
2. Integrar vistas existentes antes de duplicarlas.
3. No crear versiones paralelas de `VistaPortico` o `VistaPlantaPorticos`.
4. Mantener la navegación principal centralizada en `VentanaPrincipal`.
5. Los cálculos no deben implementarse en la GUI.
6. La interfaz muestra datos y lanza procesos; la lógica vive en `calc/`.
7. Si existe una implementación funcional en `tests/`, `tools/` o archivos temporales recientes, evaluar reutilizarla antes de crear una nueva.
8. Antes de modificar una pantalla existente, verificar si ya existe documentación o implementación funcional dentro de la aplicación principal.
9. Las modificaciones deben integrarse sobre la aplicación real ejecutada por `calculador.bat`; evitar demos, prototipos o ventanas paralelas.

## Organización actual de la interfaz

### Criterio de diseño
Reorganizar la interfaz principal para separar claramente:
- Vista global de la obra.
- Vista particular del pórtico seleccionado.

### Importante
Cambios aplicados sobre la aplicación real ejecutada mediante `calculador.bat`. Se integraron los componentes existentes (no se crearon tests, demos, ventanas temporales ni prototipos). Los cambios quedan visibles al ejecutar `calculador.bat`.

### Estado actual 

La pestaña Inicio contiene ahora:
- `VistaPlantaPorticos` integrada en el grupo **"Planta de la obra"** con `scroll_vista_portico_p` y `combo_nivel_planta`.
- `tabla_inicio` para el resumen del pórtico activo.
- Referencias de ejes y niveles actualizadas.

La pestaña **"Pórticos"** fue creada con:
- Combo **"Pórtico activo"**.
- `VistaPortico` (esquema del pórtico seleccionado).
- Controles asociados al pórtico seleccionado.

## Problemas conocidos

### Migración a columnas y bases globales

La arquitectura actual utiliza columnas y bases globales identificadas mediante:

- `COL-XXX`
- `BASE-XXX`

Los pórticos almacenan referencias a estos elementos en lugar de mantener definiciones locales propias.

Ejemplo:

```json
"columnas": [
  "COL-001",
  "COL-002"
]


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

Archivo de referencia para dar indicaciones. La app se abre con `calculador.bat`, que inicia `-m app` usando `.venv\Scripts\pythonw.exe` (sin consola). También puede iniciarse desde una terminal con `py -m app`. Ventana principal: **"Calculador — Estructuras"** (QMainWindow). Tamaño inicial 1340×820, mínimo 1050×680.

## Ciclo de arranque y actualización

1. `calculador.bat` cambia el directorio a la raíz del proyecto y lanza el Python gráfico del entorno `.venv`. La consola del archivo por lotes se cierra; los errores de inicio pueden consultarse ejecutando `py -m app` desde una terminal.
2. `app.__main__` llama a `app.principal.main()`. El arranque inicializa la selección/migración de obras, recarga el pipeline y asegura las carpetas necesarias.
3. Se crea `QApplication`, `VentanaPrincipal` y sus pestañas. La obra activa se lee desde la configuración de rutas; los datos de obra se guardan bajo `Obras/<obra>/datos/` y sus resultados bajo `Obras/<obra>/salidas/`.
4. La actualización de la ventana lee datos y resultados guardados para armar los resúmenes, tablas y visores. No debería resolver el modelo de losas ni el motor estructural completo por el solo hecho de abrir la aplicación.
5. Los cálculos se inician desde sus acciones específicas (por ejemplo, resolver solicitaciones o calcular una losa). Las reacciones automáticas de losas en el visor de pórtico se difieren durante la actualización general; el visor de muros solicita reacciones al seleccionar un muro.

El botón **Actualizar** y el cambio de obra refrescan los datos mostrados. Los resultados guardados no se recalculan automáticamente por editar o refrescar; se ejecuta el cálculo correspondiente desde la interfaz.

El orden de las nueve pestañas principales es: **Inicio**, **1 · Estado y etapas**,
**Cargas**, **Configuración de obra**, **Pórticos**, **Estructura global**,
**5 · Losas**, **6 · Muros** y **7 · Archivos**. Vigas, columnas y bases son
subpestañas de **Pórticos**, no pestañas principales.

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

### 2) 1 · ESTADO Y ETAPAS
- Tabla de etapas (semáforo: ok / desactualizada / pendiente / sin_datos / a_desarrollar).
- Detalle de la etapa seleccionada + botones **"Ejecutar"** y **"Abrir salida"**.

### 3) CARGAS (`PaginaCargas`)
Pestaña con 2 solapas (QTabWidget `"pestanias_cargas"`):

- **a) "Definir elementos"** (scroll): catálogo de elementos de carga de la obra.
  - Tabla de elementos: ID, Descripción, Aplicación, Transferencia a apoyos, Categoría, D, L, Incluir, Aplicaciones manuales.
  - Botones: Nueva carga, Nueva losa desde ejes, Guardar cargas, Editar composición, Eliminar elemento, Aplicar a tramos, Editar composición del elemento, Quitar asignación, Editar paño de losa, Generar informe, Abrir informe.
  - Grupo **"Aplicar y revisar cargas en barras"**: selector de elemento, tabla de aplicaciones (Pórtico(s), Carga, Tipo, Barra/tramo, Desde/Hasta, Ancho).
  - El cálculo de losas se inicia en **5 · Losas**, no desde este catálogo.
- **b) "Elementos y asignaciones"** (scroll con catálogo+aplicaciones).

**Diálogos asociados:**
- `EditorElemento`: **"Componer carga — {nombre}"**.
- `EditorAplicaciones`: **"Aplicaciones de {carga_id}"**.

### 4) CONFIGURACIÓN DE OBRA
- Configuración del material/modelo de análisis y datos de terreno/fundaciones.
- Los datos se guardan en la obra activa; no se comparten entre obras salvo la biblioteca global.

### 5) PÓRTICOS
- Combo **"Pórtico activo"** para elegir la estructura cuyos resultados se muestran y calculan.
- Subpestañas integradas: **"Pórtico"**, **"Vigas"**, **"Columnas"** y **"Bases"**.
- **"Pórtico"**: `VistaPortico` (esquema), resumen y envolvente del pórtico seleccionado.
- **"Vigas" / "Columnas" / "Bases"**: tablas, acciones de dimensionado y planillas correspondientes. Columnas y bases dependen de que existan resultados previos vigentes.

### 6) ESTRUCTURA GLOBAL
- Herramientas y resultados globales para columnas, vigas, bases y vigas de fundación que comparten elementos entre pórticos.

### 7) 5 · LOSAS
- Selector de paño activo y botón **"Calcular losa seleccionada"**.
- Lista de memorias vigentes (una por losa) y visor de la memoria seleccionada.
- Recalcular reemplaza la memoria vigente de esa losa, no crea otra versión en cada ejecución.
- Para una losa maciza sin continuidad declarada, la memoria no presenta momentos de apoyo como resultados calculados; si hay paños continuos, informa los momentos próximos a los bordes compartidos.

### Formato numérico de informes y planillas
- Los informes TXT usan coma decimal y punto para separar miles.
- Los CSV activos para Excel usan `;` como separador de campos y coma decimal. Los cómputos de losas existentes con el formato anterior se convierten al formato nuevo al guardar.
- Los archivos JSON mantienen números en formato estándar. Los scripts de `Legacy/` no forman parte de esta conversión.

### 8) 6 · MUROS
- Tabla de muros de la obra.
- Botones para crear, editar y eliminar muros.
- Visor de muro (`VistaPortico`).
- La lista y el resumen guardado se actualizan sin calcular las reacciones estructurales de las losas. Al elegir un muro en el selector se actualiza su visor y se consultan las reacciones asociadas.
- Al componer un muro, el informe detalla el peso por capa, las reacciones D/L de las losas y cubiertas, sus combinaciones gravitatorias, la resistencia/transmitancia térmica y la verificación de mampostería.
- El predimensionado de la zapata es solo por presión admisible; no calcula flexión ni armadura. La térmica es informativa si no se configuró un límite U.
- El informe se conserva al guardar la composición y queda visible en el resumen del muro.

### 9) 7 · ARCHIVOS
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



## Actualizacion: Columnas globales (UI) — Octubre 2026

### Implementado
- En la pestaña "Estructura global → Columnas globales" se agregó un segundo botón: "Calcular columnas globales (completo)".
- Al ejecutarlo, se corre el pipeline completo (Fases 1→2→4→planilla final) y se carga la planilla global final en una tabla de 20 columnas.
- La tabla muestra: Columna, Eje, Nivel, N, Mx, My, Mu_max, Biaxial, b, h, f'c, Ast, Armadura, Clasificacion, Pnx, Pny, P0, Pn_biaxial, Cumple, Porticos.

### Notas
- El botón "Calcular columnas globales" original sigue disponible para solo consolidar solicitaciones.
- Los cálculos siguen viviendo en calc/; la UI solo lanza y muestra resultados.

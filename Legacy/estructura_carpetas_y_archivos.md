# Estructura de carpetas y archivos — CALCULADOR

Guía general del programa para no perder tiempo. Indica qué conviene leer, qué ignorar y qué mover.

## Carpetas principales

| Carpeta | Qué es | Leer / ignorar |
|---|---|---|
| `app/` | Ventana PySide6 (principal). Código de la interfaz. | **Leer** — núcleo de la app. |
| `calc/` | Motores de cálculo (pipeline, portico, cargas, vigas, etc.). | **Leer** — núcleo de cálculo. |
| `datos/` | JSON de datos globales (ejes, niveles, materiales, cargas). | **Leer** — datos de entrada. |
| `Obras/` | Una carpeta por obra con sus datos y salidas. | **Leer** la obra activa; ignorar las demás. |
| `salidas/` | Resultados generados (JSON/TXT/planos). | Ignorar para entender lógica; útil para verificar. |
| `tests/` | Pruebas de humo. | **Leer** solo si se toca esa área. |
| `tools/` | Scripts auxiliares (regresión, comparación). | **Leer** solo si se necesita. |
| `Legacy/` | Código viejo (P00, P01, etc.). | **Ignorar** — no tocar; referencia histórica. |
| `respaldo_limpieza/` | Backups de limpiezas anteriores. | **Ignorar / mover** a un lugar externo. |
| `Bibliografia/` | Documentos de referencia. | **Ignorar** para código. |
| `imagenes/` | Recursos gráficos. | **Ignorar** salvo UI. |
| `.venv/` | Entorno virtual Python. | **Ignorar** totalmente. |
| `.vs/` | Configuración de Visual Studio. | **Ignorar**. |
| `__pycache__/` | Caché de Python. | **Ignorar**. |

## Archivos raíz importantes

- `calculador.bat` — abre la app.
- `estado.bat` — semáforo de etapas.
- `ver_cambios.bat` — cambios sin commitear.
- `ARQUITECTURA.md` — arquitectura y flujo.
- `README.md` — tablero de trabajo.
- `estructura_pantallas_y_modificaciones.md` — estructura de pantallas.
- `estructura_carpetas_y_archivos.md` — este archivo.
- `requirements.txt` — dependencias.

## Convenciones de movimiento

- **No mover** `app/`, `calc/`, `datos/`, `Obras/`, `salidas/`: la app depende de rutas relativas.
- **Se puede mover** `Legacy/`, `respaldo_limpieza/`, `Bibliografia/` fuera del repo para aligerar (no afectan ejecución).
- **No borrar** `.venv/`: sin él no corre la app.
- Los `.bat` deben quedar en la raíz.

---

## 🔧 ESTADO ACTUAL — PRINCIPIO 4 (biaxial / columnas globales)

### Qué se hizo (ya implementado y verificado)

| Archivo | Cambio |
|---|---|
| `calc/rutas.py` | Agregadas `obtener_o_crear_columna_global(x, nivel, altura_m, eje_id, eje_nombre, nivel_inicio, estructura)` y `obtener_o_crear_base_global(x, nivel, eje_id, eje_nombre, estructura)`. Dedupe por coordenadas `(x, nivel)` → `COL-001`, `COL-002`, … / `BASE-001`, `BASE-002`, … |
| `calc/rutas.py` | Agregados `siguiente_id_columna_global()` y `siguiente_id_base_global()`. |
| `app/principal.py` | Creación de pórticos ahora usa IDs globales en lugar de nombres locales acumulativos `C{piso}-N`. `datos["columnas"]` y `datos["bases"]` son listas de IDs globales (ej: `["COL-001", "COL-002"]`). Los datos reales viven en `estructura["columnas"][COL-XXX]` (catálogo global). |
| `app/visualizacion.py` | `actualizar()` ahora soporta ambos formatos: dict antiguo `{cid: datos}` y lista nueva de IDs globales. `_geometria_tramos()` fallback si no matchea prefijo `C{piso}-`. |
| `calc/rutas.py` | `acumular_solicitaciones_columnas()` ya implementado y verificado: lee `salidas/solicitaciones/*.json`, mapea cids a globales, escribe `salidas/columnas/<COL-00X>.json` con per-portico data y envelope global. Verificado `COL-002` acumulando entre Pórtico 2 y Pórtico 5. |
| Importación | `py -m calc.portico` y `py -m app` verificadas OK. |

### Qué FALTA para cerrar el principio 4

| Tarea | Archivo | Notas |
|---|---|---|
| Dimensionador de columnas lee `salidas/columnas/<COL-00X>.json` | `app/diseno_columnas_bases.py` | Todavía lee por pórtico individual. Para biaxial/Brescel debe leer el JSON global por columna que acumula todos los pórticos. |
| Dimensionador de bases lee `salidas/bases/<BASE-00X>.json` | `app/diseno_columnas_bases.py` | Idem para bases. |
| Migrar pórticos existentes | `datos/estructura.json` | Pórticos creados antes del cambio tienen `datos["columnas"]` como dict `{C0-15: {...}}` en lugar de lista `["COL-001", ...]`. El visor soporta ambos formatos (retrocompat), pero ideal migrar. |
| Smoke test completo | `py -m app` + `py -m calc.portico` | Verificar que crear un pórtico nuevo desde la GUI genere COL-XXX y que se reutilice en otro pórtico en misma coordenada. |

### Problema conocido: visor de inicio no grafica

Síntoma: al abrir `py -m app`, la pestaña Inicio / visor de pórtico no dibuja las columnas/vigas.

Causa probable: en `app/visualizacion.py`, `actualizar()` recibe `estructura.get("columnas", {})`. Si el pórtico activo tiene formato antiguo (dict con claves `C0-15`) o si `estructura_completa["columnas"]` no contiene los gid referenciados, `self.columnas` queda vacío y `_geometria_tramos()` no dibuja nada.

Revisar:
- `app/visualizacion.py` líneas ~398–420: construcción de `self.columnas` desde lista de IDs.
- `app/visualizacion.py` líneas ~504–528: `_geometria_tramos()` filtra por `C{piso}-`; fallback a todas las columnas si no matchea, pero solo si `self.columnas` no está vacío.
- Verificar que `estructura_completa.get("columnas", {})` tenga los datos de cada COL-XXX.

---

## Notas

- El análisis sigue siendo 2D por pórtico; la app muta hacia una planta de obra por nivel.
- Para tareas de UI, leer primero `estructura_pantallas_y_modificaciones.md`.
- Para tareas de cálculo, leer `calc/` y `ARQUITECTURA.md`.

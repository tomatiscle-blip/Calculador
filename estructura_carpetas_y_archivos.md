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
- `requirements.txt` — dependencias.

## Convenciones de movimiento

- **No mover** `app/`, `calc/`, `datos/`, `Obras/`, `salidas/`: la app depende de rutas relativas.
- **Se puede mover** `Legacy/`, `respaldo_limpieza/`, `Bibliografia/` fuera del repo para aligerar (no afectan ejecución).
- **No borrar** `.venv/`: sin él no corre la app.
- Los `.bat` deben quedar en la raíz.

## Notas

- El análisis sigue siendo 2D por pórtico; la app muta hacia una planta de obra por nivel.
- Para tareas de UI, leer primero `estructura_pantallas_y_modificaciones.md`.
- Para tareas de cálculo, leer `calc/` y `ARQUITECTURA.md`.

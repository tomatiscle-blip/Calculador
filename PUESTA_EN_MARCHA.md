# Puesta en marcha — Calculador

> **Empezar por `README.md`**: ahí está el tablero de trabajo (qué está hecho, qué falta
> y en qué orden). Este documento es el detalle de las decisiones, cómo reinstalar todo
> y la guía de git.

Documento de trabajo. Explica **qué se agregó**, **por qué**, **cómo se usa** y
**qué sigue**. El cálculo viejo (P00…P06, L00, C00, V0x, 00, 10) **no se tocó**:
sigue funcionando igual que siempre.

---

## 1. Decisiones tomadas

| Tema | Decisión | Motivo |
|---|---|---|
| Motor de cálculo | **Pynite (PyNiteFEA) 3.0.0** | 3D, P-Δ, placas para losas, resortes para el terreno, combinaciones de carga nativas, licencia **MIT** (podés empaquetar un `.exe` sin obligaciones). Ya está instalado. |
| anaStruct | Se retira (queda solo como control cruzado) | La comparación ya se hizo (02/10, `tools/comparar_motores.py`): en pórticos planos dan igual (<0,1 %). Se mantiene hasta terminar de migrar `P01` a Pynite. |
| Interfaz | **PySide6** (ventana de Windows, libre, ya instalada) | Programa de verdad: ícono, doble clic, sin navegador, empaquetable a `.exe`. |
| Scripts viejos | No se borra ninguno | Quedan como "envoltorios" de consola hasta que su cálculo esté migrado a `calc/`. |

---

## 2. Archivos nuevos y para qué sirve cada uno

| Archivo | Qué hace |
|---|---|
| `requirements.txt` | Lista de librerías del proyecto. **Esto es lo que evita que se pierdan otra vez.** |
| `.gitignore` | Le dice a git que ignore `.venv/`, `__pycache__/`, `.vs/` y las copias de prueba. |
| `calc/rutas.py` | **Todas** las rutas del proyecto en un solo lugar. Resuelve las rutas a partir de la ubicación real del archivo, así funciona la app, el `.exe` o un acceso directo (antes varias rutas eran relativas y se rompían si el programa arrancaba desde otra carpeta). |
| `calc/pipeline.py` | Las 10 etapas del cálculo, en orden, con sus entradas y salidas, más el **semáforo** de estado y la función para ejecutar una etapa. |
| `tools/regresion.py` | Red de seguridad: congela los resultados actuales y avisa si cambian. |
| `tools/comparar_motores.py` | Compara el mismo pórtico en **anaStruct y Pynite** (reacciones, momentos, cortantes, axiales) y chequea el equilibrio. Valida el cambio de motor. |
| `estado.bat` | **Doble clic** para ver el semáforo del proyecto. |
| `tests/golden/` | Las copias congeladas de referencia (no se versionan en git). |

---

## 3. Si se actualiza Windows o Python otra vez

Abrir PowerShell **en la carpeta del proyecto** y pegar:

```powershell
py -m venv .venv
.\.venv\Scripts\Activate.ps1
py -m pip install -r requirements.txt
py -m pip freeze > requirements.lock.txt
```

Con el entorno creado, los programas se lanzan con el `py` de la carpeta `.venv`
(o simplemente con `py`, como hasta ahora). Si algo falla, el comando
`py -m pip install -r requirements.txt` reinstala lo que falte.

---

## 4. Uso diario (lo que ya funciona hoy)

| Para… | Hacer |
|---|---|
| Ver el estado de todo | Doble clic en **`estado.bat`** |
| Ver el estado de un pórtico | `py -m calc.pipeline "Portico 1"` |
| Guardar la foto de resultados actuales | `py tools\regresion.py congelar` |
| Ver si algo cambió | `py tools\regresion.py comparar` |
| Ver las diferencias | `py tools\regresion.py comparar --diff` |
| Ver un resultado de referencia | `py tools\regresion.py ver salidas\vigas\resultados_Portico 3_vigas.json` |

### Cómo se lee el semáforo

| Ícono | Significa |
|---|---|
| `[ OK ]` | Calculado y al día. |
| `[OJO!]` | Se calculó, pero hay datos **más nuevos**: conviene recalcular. |
| `[  -  ]` | Falta calcular (la salida no existe todavía). |
| `[FALTA]` | Faltan los datos de entrada de esa etapa. |
| `[TODO ]` | Etapa planificada, todavía sin programar (ej. terreno). |

---

## 5. Qué sigue, en orden

1. ~~**Comparar motores (anaStruct vs Pynite)**:~~ **HECHO (02/10)** con
   `tools/comparar_motores.py`. Pórticos 1 y 3 coinciden **<0,1 %** → Pynite validado.
   Falta únicamente el Pórtico 2 (sus datos están incompletos: falta `nivel` en
   `C1-a`/`C1-b`). Con eso, `P01` pasa a Pynite y quedan habilitados el 3D, el P-Δ,
   las placas y los resortes.
2. **Datos que faltan**: `datos/terreno.json` (capas, nivel freático, `q_adm`,
   módulo de balasto) y `datos/tipos_losa.json` (vigueta / maciza / casetonada),
   más los diagramas de interacción como dato canónico.
3. **App PySide6**: pestañas por etapa, semáforo en color, tablas para cargar
   datos y botones para calcular y exportar.
4. **Migrar los scripts a `calc/`**: convertir `P00`, `P02`, `P04`, `P05`, `P06`
   y `L00` en funciones `calcular(datos)` — así se eliminan las 63 preguntas por
   teclado (`input()`) y las rutas relativas que hoy atan todo a una consola.

---

## 6. Detalles anotados (para no olvidarlos)

* `P06_Portico_dxf.py` tiene el pórtico fijo dentro del script: `PORTICO = "Portico 3"`.
* `P04_Columnas_portico.py` **reescribe entero** `salidas/columnas/planilla_columnas.csv` cada vez.
* `V03_guardar_vigas_excel.py` usa `pandas.to_excel` sin indicar motor: necesita
  `openpyxl` (falta instalarlo). `P03` ya usa `engine="xlsxwriter"`, por eso ese sí funciona.
* Los nombres de archivo con espacios y paréntesis (ej. `Portico 3(mercedes) `) son
  un problema en Windows: usar siempre `calc.rutas.nombre_seguro(...)`.
* `00_Analisis_cargas.py` guarda su configuración dentro del propio script; a futuro
  pasa a `datos/cargas.json` (la etapa "1. Análisis de cargas" del semáforo).
* **Limitación conocida del semáforo**: mientras los tres pórticos vivan en un mismo
  `datos/estructura.json`, cualquier cambio en uno de ellos marca a los otros como
  "a recalcular" (`[OJO!]`), porque la comparación es por fecha de archivo. Se
  resuelve el día que cada obra/pórtico tenga su propio archivo de datos.

---

## 7. Git sin ser programador (para no asustarse)

**Qué es**: git guarda "fotos" del proyecto. Cada foto es un *commit* y después
viaja a GitHub. Nada más que eso.

**Guardar una foto desde VS Code (lo más simple)**

1. Panel **Control de código fuente** → `Ctrl + Shift + G`.
2. Escribir el mensaje en la caja de arriba (ej. `se agrega bases portico 4`).
   **Sin mensaje, git NO guarda**: por eso parece que "quedó cargando".
3. `Ctrl + Enter` (o el botón ✓). Foto hecha.
4. Para subirla a GitHub: botón **Sync Changes** (o `git push`). La primera vez
   puede pedir el usuario y la contraseña de GitHub en una ventana aparte: hay
   que dejarla abierta y completarla, **esa ventana es la que parece "colgada"**.

**Ver qué falta guardar, sin saber git**: doble clic en `ver_cambios.bat`.

**Los dos comandos para mirar cómo viene todo**

```powershell
git status -sb        # qué falta guardar y si hay fotos sin subir (ahead N)
git log --oneline -5  # las últimas 5 fotos, con su mensaje
```

**Si te arrepentís de una foto que todavía NO subiste a GitHub**

```powershell
git reset --soft HEAD~1   # desarma la última foto y deja todo listo para rehacerla
                          # (no borra ningún archivo)
```

**Señales y qué significan**

| Lo que ves | Qué es | Qué hacer |
|---|---|---|
| `ahead 1` | Hay 1 foto local sin subir | Sync Changes / `git push` |
| Caja de mensaje vacía | Falta el título de la foto | Escribirlo y `Ctrl + Enter` |
| `index.lock` | Un programa quedó a mitad de camino | Cerrar VS Code y borrar `D:\03 Ingenieria\Calculador\.git\index.lock` |
| Pide usuario/contraseña | Es GitHub pidiendo permiso para subir | Completar esa ventana (no es un error) |



# Arquitectura — Calculador

Este documento dice **cómo está pensado el programa**: qué se calcula solo, qué se
combina con qué, dónde vive cada dato y en qué orden se construye lo que falta.
Sale de cómo se trabaja de verdad: primero se definen los elementos que reciben carga
(un muro, la losa de un piso, un techo), después el pórtico que los recibe, después las
solicitaciones y recién ahí el dimensionado de cada pieza según CIRSOC.

---

## 1. El requisito que manda: versatilidad

El programa tiene que servir para tres cosas muy distintas, **con el mismo motor y las
mismas cuentas**:

| Querés… | Definís | Y el programa te da… |
|---|---|---|
| Una **losa alivianada sola** | el elemento losa y nada más | su carga, su dimensionado y su cómputo (no pide pórtico) |
| Un **muro** que recibe el techo y la losa de un piso | el muro + los elementos que le llegan | la carga del muro y **su cimiento** |
| Un **pórtico** de una estructura complicada | los elementos + el pórtico | solicitaciones y dimensionado de cada barra |

La regla que hace posible esa versatilidad:

> **Todo es un elemento que se puede calcular solo, y los elementos se pueden combinar
> en un contenedor (el pórtico).**

El modelo permite calcular elementos sueltos o combinarlos dentro de un pórtico. En el
flujo actual de la app ya se definen y aplican cargas a tramos existentes; la creación
gráfica de la geometría del pórtico todavía está pendiente.

---

## 2. Las 4 capas

| Capa | Qué hace | Qué NO sabe |
|---|---|---|
| **1. Modelo** | Lo que definís: obra, materiales, elementos que reciben carga, pórticos, combinaciones | fórmulas |
| **2. Motor** | Solo mecánica: el pórtico → solicitaciones (M, V, N, flechas) | hormigón ni acero |
| **3. Dimensionadores** | Uno por familia: viga, columna, base, losa (por tipología), acero | de dónde vinieron las cargas |
| **4. Salidas** | Planillas, memorias, planos, cómputo: **vistas** del mismo resultado | nada de cálculo |

Cada capa se prueba sola. El motor no cambia cuando aparece el pórtico metálico: cambia
el dimensionador. Y la losa maciza o casetonada es la misma "losa" con otro `tipologia`.

---

## 3. Los 7 pasos del flujo de trabajo objetivo

1. **Obra** — nombre, ubicación (de ahí sale el viento de CIRSOC 102), reglamento.
2. **Materiales** — la biblioteca (`datos/materiales.json`): pesos específicos, cargas
   superficiales, sobrecargas, viento.
3. **Elementos que reciben carga** — muro, losa (alivianada / maciza / casetonada),
   techo/cubierta (+ encadenados). Cada uno declara sus componentes, su uso y su ancho
   tributario. De acá salen las cargas **D, L y W** y las combinaciones de CIRSOC.
4. **Pórticos** — geometría (naves, pisos, alturas, vinculaciones), **qué elementos
   cargan encima**, cargas puntuales y el viento como carga horizontal.
5. **Solicitaciones** — el motor resuelve el pórtico 2D y devuelve M, V, N y flechas por
   barra y por combinación. Sale a un archivo propio (`solicitaciones.json`).
6. **Dimensionamiento** — cada familia lee las solicitaciones y verifica: viga, columna,
   base, losa, acero (CIRSOC 201 / 301).
7. **Salidas** — planilla, memoria de cálculo, planos (DXF) y cómputo, armados desde los
   resultados (no recalculados).

**Estado de la interfaz:** este es el orden conceptual, pero la app aún no permite crear
la geometría en el paso 4. `datos/crear_estructura.py` la solicita por consola; la lista
superior de la app solo selecciona pórticos ya guardados. El editor visual de geometría
es el siguiente paso de producto.

---

## 4. Dónde vive cada dato

| Archivo | Qué contiene | Quién lo lee |
|---|---|---|
| `datos/materiales.json` | La biblioteca: pesos específicos, cargas superficiales, sobrecargas, viento | `calc/materiales.py` → todo el programa |
| `datos/cargas.json` | **Los elementos** (muro, losa, techo, encadenado) y el viento general | `calc/cargas.py` |
| `datos/estructura.json` | Los pórticos: columnas, vigas, tramos, bases y, tras dimensionar, datos de sección `b`/`fc` | `datos/crear_estructura.py`, `calc/portico.py`, `calc/diseno_vigas.py`, P02 |
| `datos/cargas.json` → `aplicaciones` | Qué carga va a qué pórtico y tramo, intervalo `x`, ancho tributario y modo de reemplazo | `app/cargas_proyecto.py`, `calc/cargas.py`, `calc/portico.py` |
| `datos/terreno.json` | Capas, nivel freático, q_adm, módulo de balasto | bases (falta crear) |
| `datos/tipos_losa.json` | Vigueta / maciza / casetonada con sus datos | losas (falta crear) |
| `salidas/…` | Todo lo que produce el cálculo | la ventana y los planos |

**Nada de esto se escribe dos veces.** Si un material aparece en la losa y en el pórtico,
es el mismo renglón del mismo archivo: lo que hace que el mismo contrapiso no pueda valer
0,80 en un lado y 0,85 en el otro.

---

## 5. Estado: qué está hecho y qué sigue

| Paso | Estado | Cómo se ve hoy |
|---|---|---|
| 1. Obra | Parcial | La app guarda identidad de proyecto y referencia de viento en `datos/cargas.json`; aún falta el manejo de varias obras/archivos aislados |
| 2. Materiales | **Hecho (30/09)** | `datos/materiales.json` es la única biblioteca y `calc/materiales.py` la reparte |
| 3. Elementos que reciben carga | **Hecho (30/09)** | `datos/cargas.json` + `calc/cargas.py` (calcula el conjunto **o un elemento solo**); 01/10: la **losa entrega sus reacciones** (`calc/losas.py`) |
| 4. Pórticos y asignación de cargas | Parcial | `datos/estructura.json` guarda geometría y `datos/cargas.json` ya guarda aplicaciones a pórtico/tramo/intervalo/ancho; falta crear y editar geometría desde la app |
| 5. Solicitaciones | **Hecho (03/10)** | `calc/portico.py` (Pynite) escribe `salidas/solicitaciones/<pórtico>.json`; Inicio puede ejecutar el motor y visualizar geometría/cargas |
| 6. Dimensionamiento | En marcha | Losas en `calc/losas.py`; vigas conectadas a la app mediante `calc/diseno_vigas.py`, que adapta solicitaciones al P02 legacy. Falta migrar la lógica P02 a `calc/vigas.py`; P04 columnas y P05 bases siguen pendientes |
| 7. Salidas | A migrar | `P03` Excel, `P06` DXF; las memorias están dentro de cada script |

Lo siguiente para retomar es **crear/editar pórticos desde la app** y visualizar su
geometría. Después, continuar la migración de dimensionadores (columnas y bases) y salidas.
La aplicación de cargas a tramos ya está modelada; revisar especialmente combinaciones,
cargas parciales y que no se dupliquen D/L.

---

## 6. Reglas de trabajo

1. **Un dato, un lugar.** Ningún número nuevo se escribe dentro de un script: va a `datos/`.
2. **Una cuenta, un lugar.** Las combinaciones de carga y las conversiones (m2 → m,
   γ × espesor) viven en `calc/cargas.py`; nadie las reescribe.
3. **Nada de `input()` en `calc/`.** Los datos entran por parámetro (los scripts viejos
   que piden por teclado se van reemplazando de a uno).
4. **Cada paso se valida con números**: antes de cambiar un cálculo, se corre
   `py tools\regresion.py comparar`. Si un número cambia, se explica por qué.
5. **El script viejo queda como envoltorio** hasta que su reemplazo dé los mismos
   números; recién ahí se archiva en `legacy\` (ver la sección 10 del README).

---

## 7. Alcance y archivos: cómo no se pierde nada (01/10)

Dos reglas nuevas, que son las que sostienen la versatilidad (desglosar el cálculo sin
que la info quede suelta ni se pise).

### 7.1. Todo cálculo tiene un alcance

| Alcance | Qué es | Dónde vive |
|---|---|---|
| **Obra** (con nombre) | Un proyecto: sus elementos, sus pórticos, sus salidas | `obras\<Obra>\` |
| **`_sueltos`** | Un elemento individual, para consultar rápido | `obras\_sueltos\` |

El `id` (`L0-1`) debería ser único **dentro** de su alcance: la identidad completa es
**`alcance/id`**. La separación real por obra y el selector de alcance todavía no están
implementados. Hoy se trabaja con los JSON compartidos de `datos/` y el inicio de la app
muestra la identidad del proyecto actual.

### 7.2. Un resultado, un archivo. Los agregados son vistas

- Cada elemento escribe **su** archivo: `salidas/<alcance>/<etapa>/<id>.json`. **Nunca**
  un JSON único donde todos escriben.
- Los **agregados se generan** leyendo los individuales y se pueden rehacer: el
  `_conjunto` de reacciones, los `.csv`, las planillas. **Nunca son la fuente.**

Implementado (01/10, losas): `salidas/losas/<id>.json` es la **fuente** (una por losa) y
`salidas/reacciones/_conjunto.json` es la **vista** (se arma leyendo esos archivos). La
memoria de la L00 quedó idéntica.

**Por qué** (el caso que lo motivó): `datos/estructura.json` **junta todos los pórticos**
con una clave por nombre (`"Portico 1"`, `"Portico 2"`…). El creador actual
`datos/crear_estructura.py` usa `calc/rutas.py::proximo_nombre_portico()` para proponer
un nombre libre, pero todavía trabaja por consola. La app no tiene aún editor para crear,
renombrar o modificar geometría; su selector superior solo elige los pórticos existentes.
La separación de datos por obra sigue pendiente.

### 7.3. Lo que es común va junto

La **biblioteca** (`materiales.json`, `viguetas.json`, `perfiles_metalicos.json`,
`coeficientes_kd.json`, `diagramas_interaccion\`) es **global**: es igual para todas las
obras. Por eso se queda arriba, compartida. Solo son **por obra** los datos del proyecto
(`estructura`, `cargas`, `losas`, `terreno`) y las **salidas**.

### 7.4. Estructura objetivo

```
Calculador\
├─ datos\                  biblioteca GLOBAL (común a todo)
│   ├─ materiales.json  viguetas.json  perfiles_metalicos.json
│   └─ coeficientes_kd.json  diagramas_interaccion\
├─ obras\
│   ├─ Casa Mercedes\
│   │   ├─ obra.json       nombre, ubicación (de ahí el viento), reglamento
│   │   ├─ datos\          estructura.json  cargas.json  losas.json  terreno.json
│   │   └─ salidas\        analisis_cargas\  vigas\  columnas\  bases\  losas\  reacciones\  dxf\
│   └─ _sueltos\
│       ├─ datos\          losas.json
│       └─ salidas\        losas\  reacciones\
```

### 7.5. Migración sin romper

`calc/rutas.py` pasa a ser "consciente del alcance" (`obra_actual()`, `datos()`,
`salidas()`), con el **alcance por defecto apuntando a lo de hoy** (`datos\` y `salidas\`
de la raíz). Así nada se rompe y se migra de a poco a `obras\`.

## 8. Estado de la app y cómo retomar — 05/10/2026

### Lo que ya se puede hacer desde `calculador.bat`

- **Inicio:** elegir el pórtico existente en la lista superior; ver su estado y la
  visualización de la geometría y las cargas. Desde ahí se pueden resolver las
  solicitaciones del motor Pynite.
- **Cargas:** crear/editar composiciones y asignarlas a un pórtico y tramo, con intervalo
  `x_inicio`–`x_fin`, ancho tributario y opción de reemplazar cargas previas. La carga
  permanente `D` y la sobrecarga `L` se conservan separadas para combinarlas una sola vez.
- **Vigas:** tras resolver el motor, el botón **Dimensionar vigas** pide `b` y `fc` por
  viga. `app/diseno_vigas.py` guarda esos datos en `estructura.json`; `calc/diseno_vigas.py`
  construye la entrada para P02 a partir de `estructura.json`, las aplicaciones de carga
  y `salidas/solicitaciones/<pórtico>.json`. P02 genera el JSON de resultados y la planilla
  TXT. `P02_Viga_portico.py` también conserva el uso por consola y ahora evita iniciar el
  diálogo al importarse.
- El diálogo de secciones fija el contraste de sus campos para evitar texto blanco sobre
  fondo blanco con el tema de Windows.

### Límites actuales que hay que tener presentes

- P02 sigue siendo el cálculo legacy; el adaptador lo conecta a los datos actuales, pero
  no reemplaza ni valida por sí mismo sus criterios de diseño.
- P02 aproxima la flecha con carga uniforme equivalente. Avisa cuando hay cargas parciales;
  las cargas puntuales no se incluyen en esa comprobación de flecha.
- Los valores de ubicación de viento CIRSOC 102-25 se guardan como referencia en el
  proyecto; la configuración todavía indica que el motor usa el cálculo de viento previo.
- La creación de geometría no está en la GUI. `datos/crear_estructura.py` sigue siendo el
  ingreso por consola; no se debe confundir el selector de pórtico con un editor.

### Próximo paso acordado

Crear una pantalla inicial para crear y editar uno o varios pórticos: nombre, cantidad de
pisos, luces de tramos, alturas, columnas/apoyos y voladizos. Debe guardar en el formato
actual de `datos/estructura.json`, evitar pisar nombres existentes y mostrar una vista
previa clara. Luego revisar con el usuario el flujo y la geometría antes de conectar más
dimensionadores. Mantener los scripts legacy disponibles durante esa transición.


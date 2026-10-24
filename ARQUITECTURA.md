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

El modelo permite calcular elementos sueltos o combinarlos dentro de un pórtico. En la
app, la ficha de `datos/cargas.json` representa la obra activa; `datos/estructura.json`
guarda varios pórticos bajo esa misma obra. Las cargas y sus elementos se comparten a
nivel de obra, y cada aplicación identifica su destino. Los resultados calculados son
salidas regenerables, no la fuente de las cargas.

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

**Estado de la interfaz:** Inicio permite crear geometría básica de un pórtico; se
seleccionan las columnas en cada tramo entre niveles y se puede indicar hasta qué nivel
llega el pórtico, aunque la obra tenga niveles superiores definidos. Todavía no incluye
edición posterior ni cargas puntuales. `datos/crear_estructura.py` sigue disponible por
consola; el selector de pórticos quedó dentro del bloque correspondiente en Inicio.

---

## 4. Dónde vive cada dato
### 4.1. Carpetas principales del proyecto

| Carpeta | Función |
|---|---|
| `app/` | Aplicación PySide6. Contiene ventanas, pestañas, diálogos y navegación. |
| `calc/` | Motores de cálculo y lógica estructural (cargas, pórticos, vigas, columnas, bases, losas, etc.). |
| `datos/` | Datos compartidos del sistema: materiales, coeficientes, tipologías y configuraciones generales. |
| `Obras/` | Una carpeta por obra. Cada obra contiene sus propios datos y resultados. |
| `salidas/` | Resultados generados por los cálculos. Son regenerables y no constituyen la fuente de datos. |
| `tests/` | Pruebas automáticas y verificaciones de funcionamiento. |
| `tools/` | Utilidades auxiliares (regresión, comparación de motores, diagnósticos y mantenimiento). |
| `legacy/` | Scripts históricos conservados como referencia durante la migración a `calc/`. |
| `Bibliografia/` | Reglamentos, documentación técnica y criterios de cálculo de consulta. |
| `imagenes/` | Recursos gráficos utilizados por la aplicación. |
| `.venv/` | Entorno virtual de Python. No forma parte del código del proyecto. |
| `.vs/` | Configuración local de Visual Studio. |
| `__pycache__/` | Caché generada automáticamente por Python. |
---
### 4.2. Archivos principales de la raíz

| Archivo | Función |
|---|---|
| `calculador.bat` | Abre la aplicación principal. |
| `estado.bat` | Muestra el estado de las etapas y el semáforo del proyecto. |
| `ver_cambios.bat` | Muestra los cambios pendientes de guardar en git. |
| `README.md` | Estado actual, hoja de ruta y trabajo en curso. |
| `ARQUITECTURA.md` | Diseño general del sistema y organización de datos. |
| `INTERFAZ.md` | Implementación real de la interfaz gráfica, pestañas y navegación. |
| `PUESTA_EN_MARCHA.md` | Instalación, recuperación del entorno y guía de git. |
| `requirements.txt` | Dependencias necesarias para ejecutar el proyecto. |
| `.gitignore` | Archivos y carpetas que no deben versionarse en git. |

## 5. Carátula y navegación del programa

La pantalla **Inicio** funciona como carátula de la obra activa y punto de entrada,
no como una pantalla técnica más. Su alcance mínimo:

- Identidad: nombre de obra, ID generado y ubicación (la ubicación alimenta la
  configuración de viento). Notas son opcionales; cliente/propietario puede agregarse
  cuando exista un uso concreto.
- Resumen: cantidad de pórticos, vigas aisladas y losas, más cálculos disponibles o
  pendientes.
- Acciones principales: crear obra, abrir obra y continuar con **Pórticos**, **Viga
  aislada** o **Losas**.
- El selector de pórtico aparece dentro de Pórticos, donde corresponde elegir qué
  estructura resolver. No ocupa la navegación global.

La navegación inicial agrupa el selector bajo Pórticos y ofrece accesos a Cargas y Losas.
La ficha y las rutas de cálculo usan la obra seleccionada en `Obras/`; la migración inicial
copia los datos y resultados existentes y conserva los originales. La biblioteca de
materiales permanece común.

Para controlar el alcance, no se agregan campos a la carátula hasta que un flujo los
necesite. El análisis de cargas se comparte por obra y puede asignarse a varios pórticos;
las salidas se regeneran desde los datos fuente. Los cálculos existentes de losas y vigas
aisladas se conectan desde la navegación, sin duplicar sus motores.

## 6. Estado: qué está hecho y qué sigue

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
5. **El script viejo queda como envoltorio** ver `legacy\` en caso de ser necesario se vuelce a el para implemetar soluciones pentientes. (ver la sección 10 del README).

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
**`alcance/id`**. La app ahora selecciona una obra en `Obras/`; cada obra mantiene sus
datos y resultados bajo su propia carpeta. `_sueltos` queda como idea futura para cálculos
individuales que no pertenecen a una obra.

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

### 7.5. Selección y migración de obras

`calculador.bat` abre la última obra activa. La primera apertura copia los datos fuente a
`Obras/<nombre>/datos`; los resultados previos se conservan bajo
`Obras/<nombre>/archivo_migracion/salidas` y los originales siguen en la raíz. El botón
**Nueva obra** crea una carpeta independiente y vacía; el selector cambia entre obras.
`calc/rutas.py` dirige cargas, estructura, losas, terreno, `moments_input.json` (cuando se genera) y resultados
a la obra seleccionada. Materiales, viguetas, perfiles, coeficientes y viento quedan en
`datos/` como bibliotecas compartidas. Las obras nuevas arrancan vacías: `cargas.json`,
`estructura.json` y `moments_input.json` se crean en la carpeta de la obra y ya no se copian
desde `datos/` de la raíz.

## 8. Evolución espacial: Columnas compartidas y pórticos cruzados (Octubre 2026)

Para superar la limitación de los pórticos planos independientes cuando se intersectan en planta (por ejemplo, una casa con un pórtico transversal en sentido Y y tres pórticos en sentido X que comparten columnas), el programa adopta los siguientes principios de modelado:

1. **Identidad Global de Columnas y Nudos:** 
   Las columnas y apoyos no son propiedad exclusiva de un único pórtico. Se definen a nivel de obra por su posición en la grilla y coordenadas $(X, Y)$, garantizando un ID único global.
2. **Pórticos como Planos de Análisis (Trazas):** 
   Un pórtico actúa como un plano virtual de corte para alimentar el motor PyNite. Cuando dos pórticos se cruzan en una esquina, **reutilizan el nudo y la columna existente** en lugar de duplicarlos o sobrescribirlos.
3. **Integridad al Eliminar:** 
   Eliminar un pórtico borra únicamente sus vigas y su definición geométrica como plano de análisis, preservando los elementos globales de apoyo (columnas y bases) si otros elementos de la obra los siguen referenciando.
4. **Acumulación de Solicitaciones (Biaxialidad):** 
   Las solicitaciones obtenidas al resolver de manera independiente los pórticos longitudinales y transversales se vinculan al ID global de la columna para que los dimensionadores posteriores evalúen el efecto bidireccional (flexión esviada).

### Estado actual de implementación

Los principios anteriores ya comenzaron a implementarse en la estructura de datos de la aplicación.

Actualmente:

- Las columnas utilizan identificadores globales (`COL-XXX`) a nivel de obra.
- Las bases utilizan identificadores globales (`BASE-XXX`) a nivel de obra.
- Los pórticos almacenan referencias a esos elementos globales en lugar de mantener copias independientes.
- La acumulación de solicitaciones por columna ya está prevista para que una misma columna física pueda recibir esfuerzos provenientes de varios pórticos.
- Se mantiene compatibilidad temporal con obras y pórticos creados bajo el esquema anterior mientras se completa la migración.
- Los dimensionadores legacy de columnas y bases (`P04`/`P05`) ya se ejecutan desde la app usando IDs globales y datos de la obra activa.

Este esquema permite avanzar hacia el dimensionamiento de columnas y bases considerando la interacción entre pórticos longitudinales y transversales, manteniendo el análisis estructural organizado por pórticos planos.
  

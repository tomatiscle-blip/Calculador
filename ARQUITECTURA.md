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

**Estado de la interfaz:** Inicio permite definir ejes y niveles, crear pórticos y muros,
y consultar la planta de la obra. El selector del pórtico activo y el esquema de elevación
están en **Pórticos**. La app también contiene flujos de cargas, losas, muros, dimensionado,
estructura global, estado de etapas y archivos. Los alcances específicos de cada pantalla
y sus acciones se describen en `INTERFAZ.md`.

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

### 5.1. Arranque y separación entre interfaz y cálculo

`calculador.bat` cambia a la carpeta del proyecto y abre `app` con el ejecutable
`.venv\Scripts\pythonw.exe` (sin consola). Para ver errores de arranque se puede usar
`py -m app` desde una terminal. El flujo de `app.principal.main()` inicializa las obras,
prepara directorios, crea `QApplication` y construye `VentanaPrincipal`.

La construcción y el refresco inicial leen datos y salidas para poblar la interfaz; no
deben lanzar el análisis del modelo de losas o de pórticos como trabajo de fondo implícito.
El visor de pórtico difiere la obtención de cargas automáticas durante el refresco general,
y el visor de muro calcula/consulta reacciones de losas al seleccionar el muro. Los cálculos
estructurales se disparan desde una acción explícita de la pantalla correspondiente. Esta
separación permite que abrir o refrescar la app no ejecute un análisis pesado sin que el
usuario lo haya pedido.

El refresco no invalida ni reemplaza resultados guardados. La validez se representa en el
semáforo del pipeline; cuando una entrada cambió, el usuario decide recalcular la etapa.

## 6. Funcionamiento implementado y evolución

| Área | Funcionamiento actual |
|---|---|
| Obras | El selector permite trabajar con obras bajo `Obras/`; la obra activa determina las rutas de sus datos y resultados. La biblioteca de materiales sigue siendo global. |
| Cargas | Los elementos y sus asignaciones se mantienen por obra. Se pueden consultar, editar y transferir cargas a pórticos, vigas o muros según el modelo. |
| Geometría | Desde la app se gestionan ejes, niveles, pórticos y muros. Los datos globales estructurales permiten compartir elementos entre pórticos. |
| Solicitaciones | `calc/portico.py` resuelve los pórticos con Pynite y genera resultados por pórtico que se usan en las vistas y dimensionadores. |
| Dimensionamiento | La interfaz ofrece flujos para losas, muros, vigas, columnas y bases, además de herramientas globales. Parte de la lógica sigue adaptando o ejecutando cálculos legacy; la migración se realiza gradualmente. |
| Resultados | La app muestra estados, memorias, planillas y otros archivos generados. Los resultados son salidas regenerables; las entradas permanecen en los datos de obra. |

El detalle de pestañas y acciones está en `INTERFAZ.md`; las limitaciones vigentes
de cada cálculo deben documentarse junto con su flujo y no inferirse de notas históricas
fechadas de esta sección.

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
  

## 10. Fuente unica de verdad

- datos/cargas.json
- datos/materiales.json
- salidas/solicitaciones/*.json

## 11. Componentes heredados / Legacy

- Legacy/P01_*
- Legacy/P02_Viga_portico.py
- Legacy/P04_* / P05_*
- script.py
- tools/patch_*.py
- calc/columnas_persistencia.py
- viento_activo en cubiertas
- porticos_afectados_por_viento

## 12. Arquitectura anterior / Historial de migraciones

(Contenido historico conservado.)
## 9. Arquitectura actual del flujo de cargas, viento y solicitaciones

Flujo real:

cargas.json -> calc.cargas -> calc.portico -> PyNite -> solicitaciones JSON -> envolventes -> dimensionado

### Aplicacion de cargas
- D: member loads distribuidos sobre tramos (q_D kN/m) en calc/portico.py construir.
- L: member loads distribuidos sobre tramos (q_L kN/m) en calc/portico.py construir.
- W: node loads FX en punta de columnas; w_total = items_viento_general(..., ancho_override=ancho). Exclusivo de calc/portico.py.

### Combinaciones
- calc/cargas.py FACTORES_COMBINACIONES + calc/portico.py combinaciones() (agrega servicio).

### Envolventes
- calc/portico.py calcular() via _acumular() (maximos en modulo por barra).

### Consumidores
- app/principal.py (tabla envolvente), calc/diseno_vigas.preparar_datos, Legacy/P02, app/visualizacion.

### Viento en informe de analisis de cargas
- cargas_del_conjunto no suma W a totales D/L ni combinaciones.
- informe_completo incluye seccion descriptiva VIENTO GLOBAL DEL PROYECTO.

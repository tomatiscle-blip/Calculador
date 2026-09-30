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

Por eso el programa **no arranca por el pórtico**: arranca por los elementos. El pórtico
es solo el que junta lo que ya se calculó por separado.

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

## 3. Los 7 pasos (en el orden del trabajo, y de la pantalla)

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

---

## 4. Dónde vive cada dato

| Archivo | Qué contiene | Quién lo lee |
|---|---|---|
| `datos/materiales.json` | La biblioteca: pesos específicos, cargas superficiales, sobrecargas, viento | `calc/materiales.py` → todo el programa |
| `datos/cargas.json` | **Los elementos** (muro, losa, techo, encadenado) y el viento general | `calc/cargas.py` |
| `datos/estructura.json` | Los pórticos: columnas, vigas, tramos, bases | `calc/motor.py` (hoy `P00`/`P01`) |
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
| 1. Obra | Falta | El nombre y el viento están en `datos/cargas.json` |
| 2. Materiales | **Hecho (30/09)** | `datos/materiales.json` es la única biblioteca y `calc/materiales.py` la reparte |
| 3. Elementos que reciben carga | **Hecho (30/09)** | `datos/cargas.json` + `calc/cargas.py` (calcula el conjunto **o un elemento solo**) |
| 4. Pórticos | A medias | La geometría está en `datos/estructura.json`, pero la carga **no** dice a qué pórtico va |
| 5. Solicitaciones | Falta | Hoy los M, V y N quedan mezclados dentro de `estructura.json` (`P01`) |
| 6. Dimensionamiento | A migrar | `P02` vigas, `P04` columnas, `P05` bases, `L00` losas; falta el acero |
| 7. Salidas | A migrar | `P03` Excel, `P06` DXF; las memorias están dentro de cada script |

Lo que sigue, en orden: **reparto de cargas al pórtico (paso 4)** → **solicitaciones en
archivo propio (paso 5)** → **dimensionadores uno por uno (paso 6)** → **salidas (7)**.

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

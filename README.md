# Calculador — Estructuras

Programa de cálculo de estructuras de hormigón para obras chicas (uno o dos pisos):
pórticos, vigas, columnas, bases, losas y planos. Es un programa de escritorio: se abre
con doble clic, no necesita internet ni servidor.

**Este archivo es el tablero de trabajo del proyecto.** Dice qué está terminado, qué
falta, en qué orden y cómo se retoma. Cada vez que volvemos al proyecto, se empieza acá.

- Última actualización: **07/10/2026**
- Trabajo actual: **ejes X/Y, niveles Z y referencias de pórticos por obra**. Desde Inicio se
  pueden definir líneas de ejes no uniformes en `datos/ejes.json`, cotas nombradas en
  `datos/niveles.json`. Al crear un pórtico nuevo se elige su dirección y el eje transversal
  exacto, y se seleccionan los cruces de ejes donde nace una columna para cada tramo entre niveles.
  Las alturas se obtienen de las cotas Z reales; los cruces omitidos no agregan columnas. Se
  pueden editar las coordenadas acumuladas individualmente para columnas de niveles superiores.
  Los pórticos existentes se ven en gris y ocupan su eje transversal. Se pregunta por voladizos
  y cargas puntuales en cada nivel de viga. Una columna superior que cae entre columnas inferiores
  divide la viga en nudos FE conectados; si cae fuera de la viga/voladizo, se impide crear la geometría.
  El análisis sigue siendo 2D por pórtico: esto no convierte el modelo en una estructura espacial.
  No se admiten desfases entre un pórtico y su eje de aplicación.
- Último trabajo: **columnas y bases accesibles desde la app** (07/10). Las dos pestañas
  tienen botones de dimensionado y formularios de entrada. Columnas toma N y M de la
  envolvente de Pynite; bases toma reacciones de servicio y los parámetros de terreno
  ingresados se guardan en `datos/terreno.json`. Se conservan los dimensionadores P04/P05;
  sus criterios y limitaciones siguen siendo los de esos scripts.
- Trabajo anterior: primera reorganización de **Inicio como carátula y navegación** (06/10):
  el selector quedó en el bloque de Pórticos, y hay accesos a Cargas y Losas. La ficha
  muestra la ubicación del proyecto, derivada de la configuración de viento.
  - **Obras:** la app ahora guarda cada obra en `Obras/<nombre>/datos` y
    `Obras/<nombre>/salidas`. Al crear o seleccionar una obra se preparan las subcarpetas
    de resultados (`vigas`, `columnas`, `bases` y las demás etapas), aunque todavía no
    haya cálculos. La obra activa se recuerda al volver a abrir el calculador; en el primer
    arranque se copian los datos fuente y las salidas anteriores se guardan aparte en
    `archivo_migracion/salidas`.
  - **Alcance:** cargas, estructura, losas, terreno y resultados son propios de cada obra.
    `moments_input.json` solo se copia al migrar datos existentes; una obra nueva no recibe
    los datos de ejemplo. Materiales y tablas de diseño siguen siendo compartidos. La viga
    aislada aún no tiene recorrido integrado en la app.
  - **Losas desde ejes:** desde Inicio o Cargas se delimita el paño seleccionando dos ejes X
    y dos Y; luego se define la dirección de luz, el nivel y las dos vigas de apoyo existentes
    sobre los ejes correspondientes. La composición y las cargas se editan después. Las losas
    unidireccionales transfieren automáticamente D/L como cargas lineales a las vigas de apoyo,
    con ancho tributario igual a la mitad de la luz. No hay que asignarlas manualmente además.
    En Cargas, la columna «Transferencia a apoyos» indica si la losa activa está aplicada
    automáticamente a ambas vigas; Inicio también informa cuántos paños llegan al pórtico elegido.
    Los apoyos deben coincidir exactamente con pórticos/vigas existentes; los bordes libres sí
    pueden desplazarse dentro del alcance común de ambas vigas. Mientras una losa use sus ejes,
    no se pueden mover ni borrar esos ejes: primero hay que editar o quitar el paño.
  - **Composición y cálculo:** cada losa conserva su composición y muestra D/L en kN/m².
    La composición alimenta el cálculo existente de viguetas para losas alivianadas y las
    solicitaciones de losas macizas. La tipología casetonada está disponible para definir
    cargas, pero todavía no tiene motor de cálculo; su peso debe incluirse en la composición.
- Trabajo anterior: conexión inicial de **P02 con el motor y la app** (05/10). El motor de
  `calc/portico.py` entrega los esfuerzos por combinación; `calc/diseno_vigas.py` adapta
  esos resultados al dimensionador conservado en `P02_Viga_portico.py`. Desde la pestaña
  **2 · Vigas** se ingresan `b` y `fc`, se generan el JSON y la planilla TXT. El script
  sigue funcionando por consola y detecta el formato actual de geometría.
  - **Importante:** P02 todavía no se migró como cálculo a `calc/vigas.py`; se conserva
    su lógica de diseño. La flecha usa carga uniforme equivalente cuando hay cargas
    parciales y no incluye cargas puntuales en esa estimación.
  - `calculador.bat` abre la app. Inicio muestra la obra y ubica el selector en el bloque
    **Pórticos**; también da accesos a Cargas y Losas.
  - La ficha de cada obra queda en `Obras/<nombre>/obra.json`; cargas y pórticos, en su
    subcarpeta `datos/`. El análisis de cargas pertenece a la obra y sus aplicaciones
    indican a qué pórtico/tramo llegan. La app recuerda la obra activa. Editar cargas no
    borra los informes anteriores: las salidas quedan dentro de la misma obra.
- Trabajo anterior: **el MOTOR de cálculo quedó armado en `calc/portico.py`** (03/10).
  Resuelve el pórtico 2D con **Pynite** (casos base D, L, W y P + las 5 combinaciones, con
  viento) y devuelve **M, V, N, reacciones y desplazamientos por barra y por combinación**,
  más una **envolvente** para el dimensionado. Es el reemplazo de `P01` (anaStruct).
  - Se usó con `py -m calc.portico "Portico 1" --guardar` → `salidas/solicitaciones/Portico 1.json`.
  - **Validado contra anaStruct fresco por combinación**: coinciden las 5 combinaciones
    (reacciones, M, V y N; el momento de campo difiere ≤ 0,05 kN·m por la malla de 50).
    Antes ya se había validado Pynite contra anaStruct en <0,1 % (`tools/comparar_motores.py`).
   - Vigas, columnas y bases ya leen los esfuerzos de ese resultado desde la app.
- Último commit: `cee9b34` — *"se implementa de apoco pynite"*.

---

## 1. Cómo se usa hoy

| Para… | Hacer |
|---|---|
| Abrir el programa (ventana) | doble clic en **`calculador.bat`** |
| Ver qué está calculado y qué falta | doble clic en **`estado.bat`** |
| Ver qué cambié y todavía no guardé en git | doble clic en **`ver_cambios.bat`** |
| Estado de un pórtico suelto | `py -m calc.pipeline "Portico 4"` |
| Calcular una **losa** sola | `py -m calc.losas "L00"` (ver `py -m calc.losas --lista`) |
| Reacciones de todas las losas | `py -m calc.losas --conjunto` |
| Guardar una foto de los resultados actuales | `py tools\regresion.py congelar` |
| Ver si algo cambió respecto de esa foto | `py tools\regresion.py comparar` |
| Comparar los dos motores (anaStruct vs Pynite) | `py tools\comparar_motores.py` |

La ventana tiene 8 pestañas: **Inicio** (estado del proyecto y solicitaciones),
**Cargas**, **1 · Estado y etapas** (el semáforo en colores), **2 · Vigas**,
**3 · Columnas**, **4 · Bases**, **5 · Losas** (lee las memorias) y **6 · Archivos**
(doble clic abre el archivo con Excel, el visor de DXF, etc.). Columnas requiere
solicitaciones vigentes; bases requiere además columnas calculadas. Los parámetros
geotécnicos se ingresan en su formulario y quedan guardados en `datos/terreno.json`.

Para crear un pórtico nuevo, definí primero ejes X/Y y al menos dos niveles Z desde **Inicio**.
El asistente permite elegir el eje transversal (sin desfase), los ejes de columna en cada
entrepiso, ver la ubicación en planta y distinguir pórticos existentes en gris. En niveles
superiores, cada columna seleccionada tiene su propia coordenada acumulada editable. También
permite definir voladizos en cada nivel de viga e ingresar cargas puntuales con su coordenada local X.
Si una columna de planta alta apoya entre columnas inferiores, la viga se subdivide en el
modelo FE para compartir ese nudo; no se acepta una base de columna que quede fuera de una
viga o voladizo inferior. Desde Inicio se puede eliminar un pórtico: se quitan sus aplicaciones
de carga, se conserva el catálogo y los resultados previos quedan como históricos. Si una losa
guardada lo referencia, primero hay que reasignar su destino.

En **Cargas**, las losas y cubiertas muestran sus cargas superficiales D/L sin exigir que
estén aplicadas a una viga. Al seleccionar una losa alivianada con luz, ancho y la capa
`Losa_alivianada`, el botón **Calcular viguetas de la losa** genera la memoria y el cómputo
en las salidas de la obra. Aplicarla a tramos sigue siendo opcional.

La pantalla **Inicio** tiene accesos separados a **Definir cargas del proyecto** (crear o
editar la composición de muros, losas, cubiertas y otros elementos) y **Aplicar / revisar
cargas en barras** (asignar cada elemento a uno o más tramos y consultar los intervalos y
anchos tributarios guardados para el pórtico seleccionado). Desde esta última vista se
puede editar la composición del elemento elegido, quitar una aplicación concreta sin
borrar el elemento del catálogo, o eliminar el elemento completo junto con sus aplicaciones.
Las aplicaciones nuevas se suman por defecto a las cargas distribuidas previas de P00.
Los muros pueden aplicarse como carga lineal sobre un tramo o perpendicularmente entre
dos pórticos adyacentes: se ingresan sus dos barras receptoras y posiciones, y el motor
modela el paño simplemente apoyado y aplica en cada una la reacción puntual `q·L/2`.
Para este modo y para el ancho tributario automático de losas/cubiertas, asociá los pórticos
a sus ejes transversales desde **Asociar pórticos a ejes…**; no se pueden asociar dos pórticos
paralelos al mismo plano. Sus coordenadas acumuladas definen las separaciones. El ancho tributario
automático suma medias separaciones hacia los pórticos vecinos: supone carga continua sobre el
pórtico y no identifica un módulo/paño individual. Las losas y cubiertas mantienen también las
opciones de media luz, luz completa y ancho tributario personalizado.

---

## 2. Las tres reglas del programa

Estas tres reglas son las que le dan "unidad" al programa. Todo lo nuevo se hace
respetándolas.

1. **La pantalla coordina el flujo y muestra resultados.** Las cuentas viven en `calc/` o,
   durante la migración, en scripts legacy como P02. La app lanza el motor, pasa datos al
   dimensionador y lee sus archivos de salida.
2. **Cada dato vive en un solo lugar** (`datos/*.json`). Si un peso específico, un
   espesor o una sobrecarga está escrito en dos archivos distintos, tarde o temprano
   dan números distintos.
3. **Cada cuenta vive en un solo lugar** (`calc/`), y se llama desde la ventana, desde
   la consola o desde otro cálculo. Nada de fórmulas repetidas en dos scripts.

Lo que hoy **no** cumple la regla 2 es el análisis de cargas: está escrito *adentro* de
un script y la losa lo repite por su cuenta (ver el punto 4).

## 3. Dónde quedamos

### Última sesión — 05/10/2026
- **P02 conectado a la app y al motor actual.** `P02_Viga_portico.py` ahora se puede
  importar sin iniciar preguntas de consola; ejecutado directamente, conserva el modo
  interactivo y deriva al adaptador cuando `estructura.json` tiene el formato actual.
  `calc/diseno_vigas.py` arma la entrada legacy desde geometría, cargas aplicadas y
  `salidas/solicitaciones/<pórtico>.json`, y conserva las combinaciones de esfuerzos del
  motor. La prueba manual del Pórtico 1 generó la planilla TXT y
  `resultados_Portico 1_vigas.json`.
- **Flujo en la app:** arriba se selecciona el pórtico; Inicio resuelve solicitaciones;
  la pestaña **2 · Vigas** pide ancho `b` y hormigón `fc` para cada viga y llama a P02.
  Esos dos datos quedan guardados en `datos/estructura.json`; la altura se predimensiona
  con el criterio existente de P02, y fy/recubrimiento mantienen sus valores por defecto.
- **Límite de la flecha:** P02 calcula esa comprobación con carga uniforme equivalente.
  La app informa si el tramo tiene cargas parciales; las cargas puntuales tampoco entran
  en esa estimación. Los esfuerzos de diseño se toman de las combinaciones del motor.
- **Contraste de campos corregido** en el diálogo de diseño de vigas para que el texto sea
  legible con el tema de Windows.
- **Pendiente de esta etapa:** edición posterior de la geometría multinivel desde la app,
  cargas de losas/muros aplicadas por paño y análisis espacial 3D. El asistente actual ya
  crea los pórticos por nivel, conserva las cargas puntuales en posiciones arbitrarias y
  conecta columnas desalineadas sobre vigas/voladizos inferiores en el modelo 2D.

### Dirección de interfaz acordada — carátula y navegación

La carátula representa **una obra abierta**. Debe mostrar nombre, identificador, ubicación
(necesaria para viento), notas opcionales y un resumen de los elementos/cálculos. Desde
allí se podrá crear o abrir una obra y entrar a tres recorridos: **Pórticos**, **Viga
aislada** y **Losas**. La selección del pórtico se hará dentro de Pórticos, no como
selector global. La ficha no será una pantalla de datos técnicos ni duplicará cargas.

El selector, la creación de carpetas y la migración inicial de la obra actual ya están
implementados. La carátula y el manejo de obras quedan conectados; los recorridos de viga
aislada y losa todavía deben integrarse sin rehacer sus motores. Evitar sumar módulos o
campos hasta que un flujo los necesite.

### Última sesión — 03/10/2026
- **El MOTOR quedó armado: `calc/portico.py`.** Resuelve el pórtico 2D con Pynite con los
  casos base **D, L, W y P** y las **5 combinaciones** (servicio + las 4 de CIRSOC, con
  viento), y devuelve **M, V, N, reacciones y desplazamientos** por barra y por combinación,
  más una **envolvente** (máximos en módulo) que es la que van a consumir los dimensionadores.
  - Se usa: `py -m calc.portico "Portico 1"` (informe) y `... --guardar` para escribir
    `salidas/solicitaciones/Portico 1.json`. `py -m calc.portico` lista qué hay.
  - Las 4 combinaciones de CIRSOC siguen viviendo en **un solo lugar**
    (`cargas.FACTORES_COMBINACIONES`): el motor las reusa, no las reescribe.
  - **Límite conocido (heredado de P01):** las **cargas puntuales** se aplican **enteras** en
    todas las combinaciones (no se escalan por fD/fL). Es a propósito, para reproducir P01.
- **Validado contra anaStruct fresco por combinación** (control `_tmp_validar_motor.py`):
  reacciones, M, V y N coinciden en las 5 combinaciones; la única diferencia (~0,04 kN·m) es
  el momento de campo, por la malla de 50 del control.
- **Pantalla de inicio en la app (`app/principal.py`)**: la ventana ahora abre en la pestaña
  **«Inicio»**, que muestra las dos piezas resueltas —**cargas** (etapa 1) y el **motor de
  solicitaciones**— con su estado, un botón para **resolver el pórtico** (lanza
  `py -m calc.portico "<pórtico>" --guardar`) y la **envolvente** del pórtico elegido.
  La pantalla no calcula: lanza el motor y lee `salidas/solicitaciones/<pórtico>.json`.

### Última sesión — 02/10/2026
- **Comparación de motores hecha y con resultado**: `py tools\comparar_motores.py`
  arma el pórtico en **Pynite** (como pórtico plano, bloqueando lo de fuera del plano)
  y lo compara contra los resultados **anaStruct ya guardados** en `datos/estructura.json`
  (reacciones, momentos en extremos de vigas y columnas, cortantes y axiales), con una
  tolerancia de ±2 %. También chequea el **equilibrio** (carga aplicada vs. suma de
  reacciones), para detectar referencias viejas.
  - **Pórticos 1 y 3**: coinciden en todo, **<0,1 %**. Pynite **reproduce** a anaStruct.
    Queda **validado** el cambio de motor (lo que pedía la Fase 4).
  - **Pórtico 2**: no se puede comparar todavía. Sus columnas **`C1-a`/`C1-b` no tienen
    el campo `nivel`**, así que el pórtico de 2 pisos queda degenerado (columnas
    superpuestas) y **una carga puntual de 25,3 kN se pierde**. Su resultado guardado
    **ni cierra el equilibrio**: Σ reacciones = 510,86 kN ≠ 536,16 kN aplicados.
    Hay que **completar los datos del Pórtico 2** antes de validarlo (va con la Fase 2).
- **Cómo se lee el informe**: cada fila tiene el valor `anaStruct` (guardado), el `Pynite`
  (calculado) y la diferencia; se comparan **módulos** porque las convenciones de signo
  difieren (anaStruct da la reacción `Fy` "hacia abajo" y el corte al revés que Pynite).

### Última sesión — 01/10/2026
- **La losa alivianada salió del script**: la cuenta clásica (cargas → momento →
  tabla de viguetas → cómputo) ahora vive en **`calc/losas.py`** como `calcular(datos)`,
  y `L00_…` quedó como **envoltorio** que pide por teclado y llama a esa función (igual
  que `00_…` con `calc/cargas.py`). El informe de la **L00 salió idéntico** a la memoria
  guardada (`memoria_losa_L00_20260928_202740`).
- **Se calcula UNA losa sola**, sin pórtico ni teclado: los datos de cada losa viven en
  **`datos/losas.json`** y se corre `py -m calc.losas "L00"` (o `--lista` para ver cuáles
  hay, `--guardar` para dejarla en `salidas/losas`). Es el primer elemento que queda
  **suelto y llamable** desde la ventana.
- Lo que sigue de la losa (Fase 1, abajo): sacarle el `1,881` fijo y las combinaciones
  propias, y que use la cuenta única de `calc/cargas.py`.
- **Reacciones de la losa (Esquema A)**: la losa entrega las **2 reacciones** (una por
  apoyo) en **D / L / W sin combinar** — porque es **unidireccional**. La **combinación se
  aplica UNA sola vez**: la cuenta vive en `calc/cargas.py` (`combinaciones_de_componentes`)
  y la usan el pórtico, la losa y las reacciones, así no pueden dar distinto. *(Esquema A:
  el pórtico combina; la losa combina solo para diseñarse.)*
- **Conjunto de losas**: `py -m calc.losas --conjunto` calcula todas las losas y guarda
  **un archivo por losa** (`salidas/losas/<id>.json`, la **fuente**) + la **vista**
  `salidas/reacciones/_conjunto.json` (armada leyendo esos archivos, se puede rehacer).
  Si una losa declara `apoya_en` (`portico`/`viga`), la vista además **suma por viga** —
  ese es el reparto al pórtico. La L00 da `D=9,72 · L=5,10 kN/m` y su archivo propio.
- **Diseño acordado (alcance y archivos)**: todo cálculo vive en un **alcance** — una
  **obra** o **`_sueltos`** (elemento individual) — y vale **un resultado, un archivo**
  (los agregados, como el `_conjunto`, son **vistas** que se regeneran). Quedó escrito en
  `ARQUITECTURA.md` (sección 7). Motivo: `estructura.json` junta todos los pórticos y
  **pisa si repetís el nombre** (`P00`: `len()+1`).

### Última sesión — 30/09/2026
- **Arreglada la ventana**: faltaba un import (`QListWidgetItem`) y se caía al abrirse.
- **La etapa 1 ahora se puede correr desde la ventana**: antes se caía al capturar su
  salida (imprime γ y Windows usaba cp1252). Se arregló en `calc/pipeline.py`.
- **Una sola biblioteca de materiales** (primer ítem de la Fase 1): las tablas que
  estaban adentro de `00_Analisis_cargas.py` se mudaron a `datos/materiales.json`, con
  una `clave` por material. Verificado: el informe del análisis da **los mismos
  números** (D=11,02 · L=3,00 · W=6,02 y las 4 combinaciones de CIRSOC); solo cambian
  2 nombres de material, que ahora son iguales en los dos módulos.
- **Quedan 3 valores repetidos para decidir** (anotados en el propio archivo, en la
  sección `_revisar`): contrapiso de cascotes (17,0 vs 16,0), cubierta de chapa
  (0,15 vs 0,07) y teja (0,65 vs 0,90), y peso propio de la losa (1,81 vs 1,881).
- **La cuenta de cargas salió del script**: ahora vive en `calc/cargas.py` y los
  elementos (losa, muro, techo, encadenado) en `datos/cargas.json`. `00_…` quedó como
  envoltorio y el informe guardado es **idéntico byte a byte** al anterior (mismo hash).
- **Versatilidad (lo que pediste)**: ya se calcula **un elemento solo**, sin armar un
  pórtico: `py -m calc.cargas "Losa Alivianada L0-1"`. `py -m calc.cargas --lista`
  muestra todos los elementos y cuáles están activos.
- **`ARQUITECTURA.md`**: las 4 capas, los 7 pasos del trabajo y dónde vive cada dato.

### Terminado y funcionando
- `calc/rutas.py` — todas las rutas del proyecto en un solo lugar (antes eran relativas
  y se rompían si el programa arrancaba desde otra carpeta).
- `calc/pipeline.py` — las 10 etapas en orden, con entradas, salidas y **semáforo**.
- `tools/regresion.py` — la red de seguridad (congela resultados y avisa si cambian).
- `app/` — la ventana (PySide6, 8 pestañas) + `calculador.bat`; `estado.bat`,
  `ver_cambios.bat`; `requirements.txt` y `.gitignore`.

### Sin guardar en git (quedó a medias el 28/09)
| Archivo | Qué es |
|---|---|
| `app/principal.py` | el arreglo de la ventana (el `import` que faltaba) |
| `datos/estructura.json` | ➕ **Pórtico 4** cargado (2 vigas de 6,00 m; columnas de 3 y 4 m) |
| `salidas/losas/computo_losas.csv` | ➕ fila de la losa **L00** |
| `salidas/analisis_cargas/…28-09-2026_2026.txt` | análisis de cargas nuevo |
| `salidas/losas/memoria_losa_L00_…txt` | memoria de la losa L00 |

### Pendiente de cálculo: el Pórtico 4
Está cargado pero **no calculado**. El semáforo lo dice así:

```
[  -  ] 4. Cálculo del pórtico     0/4 columnas calculadas
[  -  ] 5. Vigas de hormigón       falta resultados_Portico 4_vigas.json
[  -  ] 7. Columnas                el pórtico Portico 4 todavía no tiene columnas
[  -  ] 8. Bases (zapatas)         falta salidas/bases/bases_Portico 4.json
[FALTA] 9. Plano lateral (DXF)     falta vigas + planilla + bases del P4
```

> Nota: los Pórticos 1 a 3 figuran en amarillo ("recalcular") solo porque se tocó
> `estructura.json` al agregar el 4; sus números no cambiaron. Se arregla el día que
> cada obra tenga su propio archivo de datos.

---

## 4. El problema de fondo: hoy la misma carga se calcula dos veces

Es exactamente lo que se veía: **el análisis de cargas sirve para el pórtico, pero la
losa lo vuelve a hacer por su cuenta**. Y peor: son **dos bibliotecas de materiales
distintas**.

| Qué | Análisis de cargas (`00_Analisis_cargas.py`) | Losas (`L00_Losas_alivianadas.py`) |
|---|---|---|
| De dónde saca los materiales | `datos/materiales.json` ✅ (desde el 30/09) | `datos/materiales.json` |
| Peso propio de la losa alivianada | 1,81 kN/m² (`SISTEMAS["Forjados"]`) | **1,881 fijo en el código** (`D1 = 1.881`) |
| Combinaciones | las 4 de CIRSOC, **con viento** | solo 1,4D y 1,2D+1,6L, **sin viento** |
| Salida | un `.txt` | un `.txt` + una fila del CSV |
| Cómo llega al pórtico | `P00` **lee ese `.txt` con expresiones regulares** y lo copia a `estructura.json` | no llega: la losa vive aparte |

**Ejemplo concreto del daño** — contrapiso de cascotes y cal de 5 cm:

- `00_Analisis_cargas.py`: γ = **17,0** kN/m³ → 0,85 kN/m²
- `datos/materiales.json`: 0,80 kN/m² → γ = **16,0** kN/m³

Mismo material, 6 % de diferencia según quién lo mire. No es un error de nadie: es lo
que pasa cuando el mismo dato está en dos lugares.

### Cómo se está arreglando (es la Fase 1 de la hoja de ruta)
1. **Hecho el 30/09**: `datos/materiales.json` es **la única** biblioteca. Las tablas que
   estaban adentro de `00` se mudaron ahí con una `clave` por material, y `00` las lee.
   Falta elegir **un valor** donde hay dos (está anotado en la sección `_revisar`).
2. `datos/cargas.json` (nuevo) guarda **qué compone cada cosa**: paño de losa, cubierta,
   ancho tributario `b`, uso y si lleva viento. Eso hoy está adentro de `00` como
   `FORJADOS` / `CUBIERTAS`.
3. `calc/cargas.py` (nuevo) hace **la única cuenta**: superficie → carga lineal,
   combinaciones CIRSOC y reparto a cada pórtico.
4. **El pórtico y la losa leen el mismo resultado.** La losa deja de rearmar las
   combinaciones y de tener el 1,881 fijo; el pórtico deja de depender de un `.txt`
   leído con expresiones regulares (pasa a leer `datos/cargas.json`).

---

## 5. Decisión de fondo: ¿2D pórtico a pórtico, o 3D?

> **DECIDIDO (02/10/2026):** se adopta **Pynite como motor único**. La comparación con
> anaStruct ya se hizo (`tools/comparar_motores.py`) y para pórticos planos dan lo mismo
> (Pórticos 1 y 3, <0,1 %). Pynite se usa **en modo pórtico plano**, pero queda el 3D,
> las placas y los resortes disponibles para cuando la estructura se complique.
> anaStruct queda solo como control cruzado hasta terminar la migración de `P01`.
>
> Nota: los pórticos que hoy están cargados en `datos/estructura.json` son **datos de
> prueba**, no obras concretas. Por eso al Pórtico 2 le pueden faltar campos (p. ej.
> `nivel`): cuando se carguen obras reales, los datos van a estar completos y consistentes.

**El pórtico se sigue pensando en 2D, pórtico a pórtico.** Es como se revisa a mano y
como se calcula hoy (`P01` usa anaStruct 2D, `SystemElements`). Eso no se cambia.

La propuesta, para dejarla escrita de una vez:

- **Un solo motor: Pynite**, pero armado como **pórtico plano** (se bloquean los
  movimientos fuera del plano, o sea "2D con el motor de 3D"). Así no quedan dos
  programas de cálculo vivos dando números distintos.
- **El 3D queda disponible, no obligatorio.** Se usa solo donde aporta algo que el 2D no
  puede dar: losas con placas, terreno con resortes, viento en dos direcciones, o un
  modelo de todo el edificio cuando haya que mirar torsión.
- **anaStruct se retira** cuando los Pórticos 1, 2 y 3 den igual en los dos motores
  (±1-2 %). Hasta entonces se mantiene como control cruzado para la memoria.

Así el ingeniero sigue trabajando pórtico por pórtico (que es como razona), pero el
programa tiene un solo motor y una sola forma de decir las cosas.

## 6. Hoja de ruta

Se marca a medida que avanza. Cada casilla es un trabajo de una sesión, más o menos.

### Fase 0 · Cimientos — HECHO
- [x] Rutas en un solo lugar (`calc/rutas.py`)
- [x] Las 10 etapas con semáforo (`calc/pipeline.py`) y `estado.bat`
- [x] Ventana propia (`app/`, PySide6) y `calculador.bat`
- [x] Red de seguridad de resultados (`tools/regresion.py`)
- [x] `requirements.txt` + `.gitignore` (para no volver a perder la configuración)
- [x] Guía de git en `PUESTA_EN_MARCHA.md` (sección 7)

### Fase 1 · Una sola fuente de datos — EN CURSO
- [x] Mover `GAMMA` / `SISTEMAS` / `SOBRECARGAS` / `VIENTO` de `00_Analisis_cargas.py` a `datos/materiales.json` — 30/09: con `clave` por material; el informe quedó con los mismos números
- [x] Crear `calc/materiales.py` (la biblioteca, para que la lea todo el programa) — 30/09
- [x] Pasar los elementos (`MUROS` / `FORJADOS` / `CUBIERTAS` / `ENCADENADOS`) y sus anchos tributarios a `datos/cargas.json` — 30/09
- [x] Crear `calc/cargas.py` (la única cuenta: superficie → lineal + combinaciones) y dejar `00_…` como envoltorio — 30/09: informe idéntico byte a byte
- [x] Que se pueda calcular **un elemento solo** (una losa, un muro) sin armar un pórtico — 30/09: `py -m calc.cargas "<elemento>"`
- [ ] Unificar los valores repetidos (ver `_revisar` en `datos/materiales.json`): cascotes 17,0 vs 16,0 · cubierta 0,15 vs 0,07 · teja 0,65 vs 0,90 · losa 1,81 vs 1,881
- [ ] Que la losa use la cuenta única (`calc/cargas.py`): sacar el `1,881` fijo y las combinaciones propias (1.2D+1.6L / 1.4D) y que vea el viento — **desde el 01/10 el lugar es `calc/losas.py`** (el `1,881` es `PESO_PROPIO_DEFECTO`)
- [ ] Que `P00` lea `datos/cargas.json` (adiós al `.txt` leído con expresiones regulares)

### Fase 2 · Cerrar la obra que está abierta (Pórtico 4) — PENDIENTE
- [ ] Pórtico 4: cálculo del pórtico (etapa 4)
- [ ] Pórtico 4: vigas y planilla de vigas (etapas 5 y 6)
- [ ] Pórtico 4: columnas (etapa 7)
- [ ] Pórtico 4: bases (etapa 8)
- [ ] Pórtico 4: plano lateral DXF (etapa 9) — ojo: `P06` tiene `PORTICO = "Portico 3"` fijo adentro
- [ ] Guardar en git lo que quedó suelto

### Fase 3 · Migrar los scripts a `calc/`, de a uno — PENDIENTE
Cada migración termina con el script viejo llamando a la función nueva, así se compara
en el momento. Se empieza por la que más molesta: losas.
- [x] `L00_Losas_alivianadas.py` → `calc/losas.py` — 01/10: `calcular(datos)` + `datos/losas.json` (clásico → tabla de viguetas); `L00` quedó de envoltorio y la memoria salió idéntica
- [ ] `P00_Ingresar_datos_estructura.py` → geometría cargada desde datos
- [ ] `P02_Viga_portico.py` → `calc/vigas.py` — conectado a la app por `calc/diseno_vigas.py`,
      pero la lógica de diseño sigue en el script legacy.
- [ ] Crear/editar geometría de pórticos desde la app; actualmente `datos/crear_estructura.py`
      sigue siendo por consola y el selector superior solo elige pórticos existentes.
- [ ] `P04_Columnas_portico.py` → `calc/columnas.py`
- [ ] `P05_Bases_portico.py` → `calc/bases.py`
- [ ] `P06_Portico_dxf.py` → `calc/planos.py` (y sacarle el pórtico fijo)
- [ ] Que la ventana pueda ejecutar esas etapas sin abrir la consola
- [ ] Al terminar cada migración, mover el script viejo a `legacy\` y actualizarlo en `calc/pipeline.py`

### Fase 4 · Motor de cálculo — EN CURSO
- [x] Comparar anaStruct vs Pynite con los Pórticos 1 y 3 (reacciones, momentos, cortantes, axiales) — 02/10: `tools/comparar_motores.py`; coinciden **<0,1 %** → **Pynite validado**
- [ ] (opcional) Comparar el Pórtico 2 cuando sus datos estén completos — hoy es **dato de prueba** y le falta `nivel` en `C1-a`/`C1-b`; se revisa al cargar una obra real
- [ ] Pasar `P01` a Pynite en modo **pórtico plano**
- [ ] Retirar anaStruct y dejar anotado en la memoria por qué se cambió

### Fase 5 · Datos que faltan — PENDIENTE
- [ ] `datos/terreno.json` (capas, nivel freático, `q_adm`, módulo de balasto) → etapa 2
- [ ] `datos/tipos_losa.json` (vigueta / maciza / casetonada) → etapa 10
- [ ] Diagramas de interacción como dato canónico (hoy son archivos sueltos)

### Fase 6 · Salidas y memoria de cálculo — PENDIENTE
- [ ] Memoria de cálculo del proyecto, armada sola con lo que ya está en `salidas/`
- [ ] Cómputo y listado de planos consolidado
- [ ] Empaquetar el `.exe` (se puede: las licencias son MIT / LGPL)

### Fase 7 · Ordenar la casa (detalle en la sección 10) — PENDIENTE
- [ ] Mover los 3 `.txt` explicativos a `docs\`
- [ ] Sacarle a cada script viejo su propio `__file__` (que use `calc.rutas`): recién ahí se pueden mover sin romperse
- [ ] `datos\` y `salidas\` por obra (`obras\Casa Mercedes\...`): cada obra con su archivo y sus resultados, sin pisarse
- [ ] `legacy\` con los scripts ya migrados, solo de referencia

---

## 7. Cómo seguimos cada sesión (para no perder el hilo)

1. Abrir este README y mirar la **hoja de ruta**: elegir **un** ítem de la fase en curso.
2. Doble clic en **`estado.bat`** para ver el semáforo antes de tocar nada.
3. Mirar **`ver_cambios.bat`** para saber qué quedó sin guardar de la vez anterior.
4. Terminado el ítem: probarlo, marcarlo acá con `[x]` y anotar la fecha si hace falta.
5. Si se tocó algo de cálculo: `py tools\regresion.py comparar`. Después, guardar en git
   con un mensaje corto (ej. `se unifica biblioteca de materiales`).

**Regla de oro para no volver a los "scripts improvisados":** ningún número nuevo se
escribe dentro de un script. Va a `datos/*.json` y el script lo lee.

---

## 8. Mapa del proyecto

| Ruta | Qué es |
|---|---|
| `calculador.bat` / `app/` | La ventana (PySide6): semáforo, tablas y abrir archivos |
| `estado.bat` / `calc/pipeline.py` | El semáforo: las 10 etapas y sus dependencias |
| `calc/rutas.py` | Todas las rutas y el guardado seguro de los JSON |
| `calc/materiales.py` / `calc/cargas.py` | La biblioteca de materiales y la única cuenta de cargas (funciona con todo el conjunto o con un elemento solo) |
| `calc/losas.py` | La única cuenta de losa alivianada (clásico → tabla de viguetas); calcula una losa sola desde `datos/losas.json` |
| `datos/*.json` | **Los datos**: materiales, geometría, cargas, losas, coeficientes, viguetas |
| `salidas/` | Todo lo calculado: vigas, columnas, bases, losas, cargas, **reacciones**, planos |
| `tools/` | Utilidades de trabajo (regresión, informe de cambios, comparación de motores) |
| `tools/regresion.py` · `tools/comparar_motores.py` | Red de seguridad de resultados y validación cruzada anaStruct ↔ Pynite |
| `tests/golden/` | Copias congeladas de resultados de referencia (no se suben a git) |
| `00_…`, `P00_…` a `P06_…`, `L00_…`, `C00_…`, `V0x_…` | Los scripts de cálculo de hoy; se migran de a uno a `calc/` |

## 9. Documentos hermanos

| Documento | Para qué |
|---|---|
| `README.md` | Este archivo: el tablero de trabajo (qué falta y en qué orden) |
| `ARQUITECTURA.md` | Cómo está pensado el programa: las 4 capas, los 7 pasos y dónde vive cada dato |
| `PUESTA_EN_MARCHA.md` | Decisiones tomadas, cómo reinstalar todo y guía de git explicada |
| `00_Readme_analisis_cargas_py.txt` | Cómo se cargan las cubiertas y el viento en el análisis |
| `C00_Readme_columnas_py.txt` | Notas del cálculo de columnas |
| `10_Metalicosreticulaejem.txt` | Ejemplo de reticulado metálico (para más adelante) |

---

## 10. Cómo queremos que quede la estructura

### Por qué los scripts viejos todavía están en la raíz

No es por dejados: **9 de ellos averiguan dónde están los datos a partir de su propia
ubicación** (`__file__`). Si se los mueve a `legacy\`, van a buscar
`legacy\datos\estructura.json` y se rompen. Verificado uno por uno:

| Se rompen si se mueven (usan `__file__`) | Se rompen si se corre desde otra carpeta (rutas relativas) |
|---|---|
| `P00`, `P01`, `P02`, `P03`, `P04`, `V01`, `V02`, `V03` | `P05`, `P06`, `L00`, `10_metalicos_correas.py` |

Por eso el orden correcto es: **primero migrar la lógica a `calc\` (Fase 3) y recién
después archivar el script viejo.** Al revés se rompe todo junto.

### El árbol al que apuntamos

```
Calculador\
├─ calculador.bat        <- doble clic: abre la ventana (queda siempre en la raíz)
├─ estado.bat            <- doble clic: el semáforo
├─ ver_cambios.bat       <- doble clic: qué cambió y qué falta guardar
├─ README.md  PUESTA_EN_MARCHA.md  requirements.txt  .gitignore
├─ app\                  <- la ventana (no calcula nada)
│   ├─ principal.py
│   └─ paginas\          (una pantalla por etapa, cuando haga falta)
├─ calc\                 <- el cálculo (sin input(), sin rutas relativas)
│   ├─ rutas.py  pipeline.py  materiales.py
│   └─ cargas.py  portico.py  vigas.py  columnas.py  bases.py  losas.py  planos.py
├─ datos\                <- los datos: la única fuente de verdad
│   ├─ estructura.json  cargas.json  materiales.json
│   ├─ terreno.json  tipos_losa.json  viguetas.json
│   ├─ coeficientes_kd.json  perfiles_metalicos.json  moments_input.json
│   └─ diagramas_interaccion\
├─ salidas\              <- resultados (no se editan a mano)
│   ├─ analisis_cargas\  vigas\  columnas\  bases\  losas\  reacciones\  dxf\
├─ docs\                 <- los .txt explicativos y la memoria de cálculo
├─ legacy\               <- los scripts viejos ya migrados, solo de referencia
├─ tools\  tests\  imagenes\
```

### La regla de cada carpeta

| Carpeta | Qué va | Qué NO va |
|---|---|---|
| Raíz | los `.bat` y los documentos | cálculo, datos ni resultados |
| `app\` | pantalla, colores, tablas | ninguna fórmula |
| `calc\` | las fórmulas: una función `calcular(datos)` por etapa | `input()`, rutas relativas, reportes |
| `datos\` | todo número que se elige (materiales, geometría, coeficientes) | resultados |
| `salidas\` | todo lo que produce el cálculo | datos de entrada |
| `docs\` | explicaciones para humanos | código |
| `legacy\` | los scripts viejos ya migrados | **nada de `calc\` ni de `app\` los importa** |

### Lo que se puede ordenar ya, sin riesgo

- [ ] Mover los 3 `.txt` explicativos a `docs\` (ningún programa los lee)
- [ ] Que cada script viejo use `calc.rutas` en vez de su propio `__file__` (así se puede
      mover sin romperse; es el mismo trabajo que migrarlo)
- [ ] Separar `salidas\` por obra → resuelve el `[OJO!]` cruzado del semáforo
- [ ] Datos por obra (`obras\Casa Mercedes\datos\...`): cada obra con su `estructura.json`,
      así un cambio en un pórtico no marca a los otros como "a recalcular"

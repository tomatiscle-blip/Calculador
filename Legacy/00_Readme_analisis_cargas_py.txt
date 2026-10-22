========================================
README – analisis_cargas.py
========================================

DESCRIPCIÓN:
-------------
Este script permite generar un análisis de cargas para distintos elementos estructurales (cubiertas, losas, muros, encadenados, pórticos) y calcular combinaciones de diseño según CIRSOC.  
Permite incluir cargas permanentes (D), sobrecargas (L) y viento (W), tanto global como específico de cubiertas. 

----------------------------------------
1. NOMBRE DEL ANÁLISIS
----------------------------------------
Definir en la sección de activación de elementos:

NOMBRE_ANALISIS = "Análisis completo – pórtico con cargas generales"

Este nombre aparecerá en la cabecera del resumen y en el archivo de salida.

----------------------------------------
2. CUBIERTAS
----------------------------------------
Formato de cada cubierta:

CUBIERTAS = {
    1: {
        "nombre": "Cubierta liviana de chapa",
        "activo": 1,                # 1=incluida, 0=omitida
        "b": 2.025,                 # ancho tributario (m)
        "componentes": [             # lista de componentes permanentes
            (SISTEMAS["Cubiertas"]["chapa_ondulada"]["nombre"], SISTEMAS["Cubiertas"]["chapa_ondulada"]["q"]),
            (SISTEMAS["Cielorrasos"]["yeso_suspendido"]["nombre"], SISTEMAS["Cielorrasos"]["Cielorraso"]["q"])
        ],
        "sobrecarga": SOBRECARGAS["cubierta_acceso_poco_frecuente"],
        "viento_activo": 1,         # 1=se aplica succión de cubierta, 0=no
        "pendiente": 6,             # grados de inclinación
        "altura_vertical": 1.20     # altura vertical de la cubierta (m)
    }
}

- Si `viento_activo=0`, no se suma la succión de esta cubierta.
- La altura tributaria para viento se calcula automáticamente: `altura_vertical * cos(pendiente)`.

----------------------------------------
3. FORJADOS
----------------------------------------
Formato similar a cubiertas: incluir componentes permanentes y sobrecarga.

----------------------------------------
4. MUROS
----------------------------------------
Definidos con función `muro(tipo, e, h)`:

MUROS = {
    "muro_comun_18": {
        "activo": 1,
        "tipo": "ladrillo_comun",
        "e": 0.18,
        "h": 1.33
    }
}

- `activo` permite incluir o no el muro en el análisis.

----------------------------------------
5. ENCADENADOS
----------------------------------------
Definidos con función `encadenado(b, h)`:

ENCADENADOS = {
    "enc_20x20": {
        "activo": 1,
        "b": 0.20,
        "h": 0.20,
        "cantidad": 1
    }
}

- Se calcula carga lineal según volumen y peso específico.

----------------------------------------
6. VIENTO
----------------------------------------
- **Viento global**: se aplica a pórticos, columnas o muros.  

```python
viento_global_activo = True   # True=se calcula, False=no
altura_global = 3.0           # altura de referencia para columnas/muros/pórticos
Cd_global = VIENTO["Cd"]      # coeficiente (Santa Fe = 1.3)

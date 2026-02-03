from datetime import datetime
from pathlib import Path
"""
ANALISIS DE CARGAS – Edificación
Autor: Claudio
Unidades:
- Peso específico: kN/m3
- Carga superficial: kN/m2
- Carga lineal: kN/m
"""

# ======================================================
# PESOS ESPECIFICOS (kN/m3)
# ======================================================
GAMMA = {
    "Hormigon": {
        "armado": {"nombre": "Hormigón armado", "gamma": 25.0},
        "sin_armar": {"nombre": "Hormigón sin armar", "gamma": 23.5},
        "alivianado_eps_alta": {"nombre": "Hormigón alivianado (EPS alta densidad)", "gamma": 10.0},
        "alivianado_eps_250": {"nombre": "Hormigón alivianado (EPS 250 kg/m³)", "gamma": 2.45},
        "cascotes_cal": {"nombre": "Contrapiso de cascotes y cal", "gamma": 17.0}
    },
    "Madera": {
        "pino": {"nombre": "Pino", "gamma": 6.0},
        "eucalipto": {"nombre": "Eucalipto", "gamma": 8.0},
        "cedro": {"nombre": "Cedro", "gamma": 5.5}
    },
    "Mamposteria": {
        "ladrillo_hueco_portante": {"nombre": "Ladrillo hueco portante", "gamma": 12.0},
        "ladrillo_hueco": {"nombre": "Ladrillo hueco", "gamma": 10.5},
        "ladrillo_comun": {"nombre": "Ladrillo cerámico común", "gamma": 17.0}
    },
    "Morteros": {
        "cemento_arena": {"nombre": "Mortero cemento–arena", "gamma": 21.0},
        "cal_arena": {"nombre": "Mortero cal–arena", "gamma": 17.0}
    }
}



# ======================================================
# SISTEMAS CON CARGA DIRECTA (kN/m2)
# ======================================================

SISTEMAS = {
    "Pisos": {
        "porcelanato": {
            "nombre": "Porcelanato 10mm + pegamento",
            "q": 0.20
        },
        "mosaico_granitico": {
            "nombre": "Mosaico granítico",
            "q": 0.60
        },
        "ceramica": {
            "nombre": "Cerámica 12 mm + pegamento",
            "q": 0.28
        }
    },

    "Cubiertas": {
        "chapa_ondulada": {
            "nombre": "Cubierta liviana (chapa ondulada + fijaciones + correas)",
            "q": 0.15
        },
        "teja": {
            "nombre": "Teja cerámica con entablonado",
            "q": 0.65
        }
    },
        "EstructuraCubierta": {
        "correas_metalicas": {
            "nombre": "Correas metálicas livianas",
            "q": 0.05
        },
        "correas_madera": {
            "nombre": "Correas de madera liviana",
            "q": 0.08
        }
    },

    "Cielorrasos": {
        "yeso_suspendido": {
            "nombre": "Cielorraso de yeso suspendido",
            "q": 0.33
        },
        "yeso_adherido": {
            "nombre": "Cielorraso de yeso adherido",
            "q": 0.18
        }
    },
    "Forjados": {
        "Losa_alivianada": {
            "nombre": "Losa alivianada EPS P.P estimado",
            "q": 1.81
        }
    }
}


# ======================================================
# SOBRECARGAS (kN/m2)
# ======================================================

SOBRECARGAS = {
    "vivienda": 2.0,
    "oficina": 3.0,
    "pasillo": 4.0,
    "garaje": 5.0,
    "terraza": 5.0,
    "cubierta_acceso_poco_frecuente": 1.0
}

# ======================================================
# VIENTO – CIRSOC 102 (Santa Fe Capital)
# ======================================================

VIENTO = {
    "ciudad": "Santa Fe",
    "V": 51.0,          # m/s
    "rho": 1.25,        # kg/m3
    "Cd": 1.3,          # coeficiente global típico
}


# ======================================================
# FUNCIONES DE CONVERSION
# ======================================================

def carga_superficial(gamma, espesor):
    """kN/m3 * m = kN/m2"""
    return gamma * espesor

def carga_lineal_muro(gamma, espesor, altura):
    """kN/m3 * m * m = kN/m"""
    return gamma * espesor * altura

def carga_lineal_desde_superficie(q_superficial, ancho_tributario):
    """kN/m2 * m = kN/m"""
    return q_superficial * ancho_tributario

def carpeta_nivelacion(espesor_m):
    gamma = GAMMA["Morteros"]["cemento_arena"]["gamma"]  # 21.0 kN/m3
    nombre = f"Carpeta nivelación {espesor_m*100:.1f} cm"
    valor = carga_superficial(gamma, espesor_m)
    return (nombre, valor)

def contrapiso_cascotes(espesor_m):
    gamma = GAMMA["Hormigon"]["cascotes_cal"]["gamma"]  # 16.0 kN/m3
    nombre = f"Contrapiso cascotes/cal {espesor_m*100:.0f} cm"
    valor = carga_superficial(gamma, espesor_m)
    return (nombre, valor)

def presion_viento(V, rho=1.25):
    return 0.5 * rho * (V**2) / 1000  # kN/m2



# ======================================================
# CLASE ANALISIS DE CARGAS
# ======================================================

class AnalisisCargas:
    def __init__(self, nombre):
        self.nombre = nombre
        self.items = []

    def agregar(self, descripcion, tipo, valor, unidad, composicion=None):
        self.items.append({
            "descripcion": descripcion,
            "tipo": tipo,          # "D" o "L"
            "valor": valor,
            "unidad": unidad,
            "composicion": composicion
        })

    def carga_lineal_viento(V, ancho_tributario, Cd=1.3, rho=1.25):
        """Carga lineal de viento (kN/m)"""
        q = presion_viento(V, rho)
        return q * Cd * ancho_tributario


    def resumen_texto(self):
        lineas = []
        D = {}
        L = {}
        W = {}  # <-- nuevo

        lineas.append(f"ANALISIS DE CARGAS: {self.nombre}")
        lineas.append("-" * 60)

        for i in self.items:
            linea = f"{i['descripcion']:<35} {i['valor']:>7.2f} {i['unidad']} ({i['tipo']})"
            lineas.append(linea)

            if i["composicion"]:
                if isinstance(i["composicion"], list):
                    for linea_comp in i["composicion"]:
                        lineas.append(f"  {linea_comp}")
                else:
                    lineas.append(f"  {i['composicion']}")
            if i["tipo"] == "D":
                D[i["unidad"]] = D.get(i["unidad"], 0) + i["valor"]
            elif i["tipo"] == "L":
                L[i["unidad"]] = L.get(i["unidad"], 0) + i["valor"]
            elif i["tipo"] == "W":   # <-- nuevo
                W[i["unidad"]] = W.get(i["unidad"], 0) + i["valor"]

            lineas.append("")

        lineas.append("-" * 60)

        for u, v in D.items():
            lineas.append(f"D total = {v:.2f} {u}")
        for u, v in L.items():
            lineas.append(f"L total = {v:.2f} {u}")
        for u, v in W.items():      # <-- nuevo
            lineas.append(f"W total = {v:.2f} {u}")

        return "\n".join(lineas)

       # -------------------------------
 
    # combinaciones CIRSOC
    def combinaciones_CIRSOC(self):
        D_total = sum(i["valor"] for i in self.items if i["tipo"] == "D")
        L_total = sum(i["valor"] for i in self.items if i["tipo"] == "L")
        W_total = sum(i["valor"] for i in self.items if i["tipo"] == "W")

        combos = {}
        combos["1.4D"] = 1.4 * D_total
        combos["1.2D+1.6L"] = 1.2 * D_total + 1.6 * L_total
        combos["1.2D+0.5L+1.6W"] = 1.2 * D_total + 0.5 * L_total + 1.6 * W_total
        combos["0.9D+1.6W"] = 0.9 * D_total + 1.6 * W_total

        return combos

        #U=1.4 (D+F)
        #U=1.2 (D+F+T) + 1.6 (L+H) + (f1 Lr o 0.5S o 0.5R)
        #U=1.2 D + 1.6 (Lr o S o R) + (f1 Lr o 0.5S o 0.5R)
        #U=1.2 D + 1.6 W + f1 L + (f1 Lr o 0.5S o 0.5R)
        #U=1.2 D + 1.0 E + f1 (L + Lr) + f2 S
        #U= 0.9 D + 1.6 W + 1.6 H
        #U= 0.9 D + 1.0 E + 1.6 H   
    
# ===============================
# ACTIVACIÓN DE ELEMENTOS
# ===============================
# Nombre general del análisis
NOMBRE_ANALISIS = "Cargas Portico 2"
# -------------------------------
# Cubiertas completas con componentes, sobrecarga y viento opcional
# -------------------------------
# NOTA: Si se analiza una carga global para un pórtico, columnas o muros,
# desactivar 'viento_activo' aquí para que no se sume al viento general.

CUBIERTAS = {
    1: {
        "nombre": "Cubierta liviana de chapa c/estructura y cielorraso suspendido",
        "activo": 1,
        "b": 2.10, # ancho tributario para calcular Portico1.
        "componentes": [
            (SISTEMAS["Cubiertas"]["chapa_ondulada"]["nombre"], SISTEMAS["Cubiertas"]["chapa_ondulada"]["q"]),
            (SISTEMAS["EstructuraCubierta"]["correas_metalicas"]["nombre"], SISTEMAS["EstructuraCubierta"]["correas_metalicas"]["q"]),
            (SISTEMAS["Cielorrasos"]["yeso_suspendido"]["nombre"], SISTEMAS["Cielorrasos"]["yeso_suspendido"]["q"])
        ],
        "sobrecarga": SOBRECARGAS["cubierta_acceso_poco_frecuente"],
        "viento_activo": 0,     # activa succión
        "pendiente": 6,         # grados de la cubierta
        "altura_vertical": None  # se completa después
    },
    2: {
        "nombre": "Cubierta de teja cerámica con cielorraso incorporado",
        "activo": 0,
        "b": 3.0,
        "componentes": [
            ("Teja cerámica tipo Marsella con entablonado", SISTEMAS["Cubiertas"]["teja"]["q"]),
            ("Cabios de madera (espesor equivalente 0.02 m)", carga_superficial(GAMMA["Madera"]["pino"]["gamma"], 0.02))
        ],
        "sobrecarga": SOBRECARGAS["cubierta_acceso_poco_frecuente"],
        "viento_activo": 0,
        "pendiente": 6,
        "altura_vertical": None  # se completa después
    }
}

# Forjados completos con componentes y sobrecarga
FORJADOS = {
    1: {
        "nombre": "Losa Maciza completa",
        "activo": 0,# 1si, 0 no
        "b": 2.90, # ancho tributario
        "componentes": [
            (SISTEMAS["Pisos"]["ceramica"]["nombre"], SISTEMAS["Pisos"]["ceramica"]["q"]),
            carpeta_nivelacion(0.025),  # 2.5 cm
            #contrapiso_cascotes(0.05),  # 5 cm, elegir contrapiso de cascotes o alivianado
            (GAMMA["Hormigon"]["alivianado_eps_250"]["nombre"] + " 0.08 m", carga_superficial(GAMMA["Hormigon"]["alivianado_eps_250"]["gamma"], 0.08)),
            (GAMMA["Hormigon"]["armado"]["nombre"] + " 0.12 m", carga_superficial(GAMMA["Hormigon"]["armado"]["gamma"], 0.12)),
            (SISTEMAS["Cielorrasos"]["yeso_suspendido"]["nombre"], SISTEMAS["Cielorrasos"]["yeso_suspendido"]["q"])
        ],
        "sobrecarga": SOBRECARGAS["vivienda"]
    },
    2: {
        "nombre": "Losa Alivianada L0-1",
        "activo": 1, #1si, 0no
        "b": 2.99,  # ancho tributario
        "componentes": [
            (SISTEMAS["Pisos"]["porcelanato"]["nombre"], SISTEMAS["Pisos"]["porcelanato"]["q"]),
            #(SISTEMAS["Pisos"]["ceramica"]["nombre"], SISTEMAS["Pisos"]["ceramica"]["q"]),
            carpeta_nivelacion(0.03),  # 3 cm
            contrapiso_cascotes(0.08),  # 8cm
            #GAMMA["Hormigon"]["alivianado_eps_250"]["nombre"] + " 0.08 m", carga_superficial(GAMMA["Hormigon"]["alivianado_eps_250"]["gamma"], 0.08)),
            (SISTEMAS["Forjados"]["Losa_alivianada"]["nombre"], SISTEMAS["Forjados"]["Losa_alivianada"]["q"]),
            (SISTEMAS["Cielorrasos"]["yeso_suspendido"]["nombre"], SISTEMAS["Cielorrasos"]["yeso_suspendido"]["q"])
        ],
        "sobrecarga": SOBRECARGAS["vivienda"]
    },
    3: {
        "nombre": "Losa Alivianada completa",
        "activo": 0, #1si, 0no
        "b": 1.10,  # ancho tributario
        "componentes": [
            (SISTEMAS["Pisos"]["porcelanato"]["nombre"], SISTEMAS["Pisos"]["porcelanato"]["q"]),
            #(SISTEMAS["Pisos"]["ceramica"]["nombre"], SISTEMAS["Pisos"]["ceramica"]["q"]),
            carpeta_nivelacion(0.03),  # 3 cm
            contrapiso_cascotes(0.08),  # 8cm
            #GAMMA["Hormigon"]["alivianado_eps_250"]["nombre"] + " 0.08 m", carga_superficial(GAMMA["Hormigon"]["alivianado_eps_250"]["gamma"], 0.08)),
            (SISTEMAS["Forjados"]["Losa_alivianada"]["nombre"], SISTEMAS["Forjados"]["Losa_alivianada"]["q"]),
            (SISTEMAS["Cielorrasos"]["yeso_suspendido"]["nombre"], SISTEMAS["Cielorrasos"]["yeso_suspendido"]["q"])
        ],
        "sobrecarga": SOBRECARGAS["vivienda"]
    }
}

# Muros reutilizables (plantillas)
def muro(tipo, e, h):
    gamma = GAMMA["Mamposteria"][tipo]["gamma"]
    nombre_base = GAMMA["Mamposteria"][tipo]["nombre"]

    return {
        "tipo": tipo,
        "nombre": f"Muro de {nombre_base} (e={int(e*100)} cm · h={h:.2f} m)",
        "valor": carga_lineal_muro(gamma, e, h),
        "unidad": "kN/m",
        "detalle": f"γ={gamma} kN/m3 · e={e} m · h={h} m"
    }
# Lista de muros del proyecto
MUROS = {
    "muro_comun_24": {
        "activo": 1,
        "tipo": "ladrillo_comun",
        "e": 0.24,
        "h": 1.33
    },

    "muro_comun_18": {
        "activo": 0,
        "tipo": "ladrillo_comun",
        "e": 0.18,
        "h": 1.33
    },
    # Muro comun 12 cm pared doble exterior
    "muro_comun_12": {
        "activo": 1,
        "tipo": "ladrillo_comun",
        "e": 0.12,
        "h": 5.00 #2.65+2.35
    },

    "muro_hueco_port_20": {
        "activo": 0,
        "tipo": "ladrillo_hueco_portante",
        "e": 0.20,
        "h": 2.60
    },
    # Muro hueco 8 cm pared doble
    "muro_hueco_np_8": {
        "activo": 1,
        "tipo": "ladrillo_hueco",
        "e": 0.08,
        "h": 5.00 #2.65+2.35
    },
    "muro_hueco_np_18": {
        "activo": 0,
        "tipo": "ladrillo_hueco",
        "e": 0.18,
        "h": 2.80
    },
}
# Encadenados (carga lineal)
def encadenado(b, h):
    gamma = GAMMA["Hormigon"]["armado"]["gamma"]
    nombre = GAMMA["Hormigon"]["armado"]["nombre"]

    return {
        "nombre": f"Encadenado de {nombre} ({int(b*100)}x{int(h*100)} cm)",
        "valor": gamma * b * h,
        "unidad": "kN/m",
        "detalle": f"γ={gamma} kN/m3 · b={b} m · h={h} m"
    }
ENCADENADOS = {
    "enc_20x30": {
        "activo": 0,
        "b": 0.20,
        "h": 0.30,
        "cantidad": 1
    },

    "enc_20x20": {
        "activo": 1,
        "b": 0.20,
        "h": 0.20,
        "cantidad": 1
    }
}

# -------------------------------
# VIENTO GLOBAL – dato general para otros análisis
# -------------------------------
# Este viento se aplica solo en análisis generales (pórticos, muros, encadenados),
# y no se suma al viento de cubiertas locales si 'viento_activo' de la cubierta está activo.

viento_global_activo = 1   # 1 = activo, 0 = desactivado # True = se calcula; False = se desactiva para no sumar al análisis
ancho_tributario_global = 2.85   # Lado transversal del pórtico expuesto al viento (NO el largo de la viga)
Cd_global = VIENTO["Cd"]        # Coeficiente del viento
# Inicializamos el valor, se calculará solo si viento_global_activo=True
w_global = None # carga lineal de viento (kN/m) = presión * Cd * ancho_tributario




# ===============================
# MAIN
# ===============================
if __name__ == "__main__":

    analisis = AnalisisCargas(NOMBRE_ANALISIS) 
    # -------------------------------
    # Cubiertas
    # -------------------------------
    for idx, c in CUBIERTAS.items():
        if c["activo"]:
            # Calcula q total
            q_d_total = sum(q for _, q in c["componentes"])
            analisis.agregar(
                f"{c['nombre']} – cargas permanentes",
                "D",
                carga_lineal_desde_superficie(q_d_total, c["b"]),
                "kN/m",
                [f"{name:<35} {q:.2f} kN/m2" for name, q in c["componentes"]] +
                [f"q total = {q_d_total:.2f} kN/m2 · b = {c['b']:.2f} m"]
            )
            # Sobrecarga
            q_l = c["sobrecarga"]
            analisis.agregar(
                f"{c['nombre']} – sobrecarga",
                "L",
                carga_lineal_desde_superficie(q_l, c["b"]),
                "kN/m",
                [f"q = {q_l:.2f} kN/m2 · b = {c['b']:.2f} m"]
            )
            # -------------------------------
            # Viento / succión (opcional)
            # -------------------------------
            if c.get("viento_activo", 0):
                pendiente = c.get("pendiente", 0)  # grados, opcional
                b = c["b"]
                # Si no está definida la altura_vertical, calcularla automáticamente
                if c.get("altura_vertical") is None:
                        from math import tan, radians
                        altura_vertical = b * tan(radians(pendiente))
                else:
                    altura_vertical = c["altura_vertical"]

                from math import cos, radians
                # Altura tributaria ajustada por pendiente
                altura_tributaria = altura_vertical * cos(radians(pendiente))

                Cd_succion = 0.9  # coeficiente de succión
                w_succion = AnalisisCargas.carga_lineal_viento(
                    VIENTO["V"],
                    altura_tributaria,
                    Cd_succion
             )

                analisis.agregar(
                    f"Viento – succión cubierta {c['nombre']}",
                    "W",
                    w_succion,
                    "kN/m",
                    f"h = {altura_tributaria:.2f} m · Cd_succion = {Cd_succion}"
                )
    # -------------------------------
    # Forjados
    # -------------------------------
    for idx, f in FORJADOS.items():
        if f["activo"]:
            # Cargas permanentes
            q_d_total = sum(q for _, q in f["componentes"])
            analisis.agregar(
                f"{f['nombre']} – cargas permanentes",
                "D",
                carga_lineal_desde_superficie(q_d_total, f["b"]),
                "kN/m",
                [f"{name:<35} {q:.2f} kN/m2" for name, q in f["componentes"]] +
                [f"q total = {q_d_total:.2f} kN/m2 · b = {f['b']:.2f} m"]
            )

            # Sobrecarga
            q_l = f["sobrecarga"]
            analisis.agregar(
                f"{f['nombre']} – sobrecarga",
                "L",
                carga_lineal_desde_superficie(q_l, f["b"]),
                "kN/m",
                [f"q = {q_l:.2f} kN/m2 · b = {f['b']:.2f} m"]
            )


    # -------------------------------
    # Muros
    # -------------------------------
    for key, m in MUROS.items():
        if not m["activo"]:
            continue
        muro_calc = muro(m["tipo"], m["e"], m["h"])
        analisis.agregar(
            muro_calc["nombre"],
            "D",
            muro_calc["valor"],
            muro_calc["unidad"],
            muro_calc["detalle"]
        )
    # -------------------------------
    # Encadenados
    # -------------------------------
    for key, e in ENCADENADOS.items():
        if not e["activo"]:
            continue

        enc = encadenado(e["b"], e["h"])
        valor_total = enc["valor"] * e.get("cantidad", 1)
        analisis.agregar(
            f'{enc["nombre"]} × {e.get("cantidad",1)}',
            "D",
            valor_total,
            enc["unidad"],
            f'{enc["detalle"]} · cantidad={e.get("cantidad",1)}'
        )
    # -------------------------------
    # Viento general para el proyecto
    # -------------------------------
    if viento_global_activo:
        w_global = AnalisisCargas.carga_lineal_viento(
            VIENTO["V"],
            ancho_tributario_global,   # lado transversal del pórtico expuesto al viento
        Cd_global
        )

        analisis.agregar(
            f"Viento diseño Santa Fe (V=51 m/s) – general",
            "W",
            w_global,
            "kN/m",
            f"b = {ancho_tributario_global:.2f} m · Cd = {Cd_global}"
        )


    # -------------------------------
    # Resumen en pantalla
    # -------------------------------
    resumen = analisis.resumen_texto()
    
    # Mostrar combinaciones CIRSOC
    combos = analisis.combinaciones_CIRSOC()    
    resumen_completo = resumen + "\n\nCOMBINACIONES CIRSOC:\n"
    for nombre, valor in combos.items():
        resumen_completo += f"{nombre}: {valor:.2f} kN/m\n"

    print(resumen_completo)
    # -------------------------------
    # Guardar en archivo txt
    # -------------------------------
    SALIDAS_DIR = Path("salidas")/ "analisis_cargas"
    SALIDAS_DIR.mkdir(parents=True, exist_ok=True)

    fecha_hora = datetime.now().strftime("%d-%m-%Y_%H%M")
    archivo_salida = SALIDAS_DIR / f"{analisis.nombre.replace(' ', '_')}_{fecha_hora}.txt"

    with open(archivo_salida, "w", encoding="utf-8") as f:
        f.write(resumen_completo)

    print(f"\nResumen guardado en: {archivo_salida}")



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
        "cascotes_cal": {"nombre": "Contrapiso de cascotes y cal", "gamma": 16.0}
    },
    "Madera": {
        "pino": {"nombre": "Pino", "gamma": 6.0},
        "eucalipto": {"nombre": "Eucalipto", "gamma": 8.0},
        "cedro": {"nombre": "Cedro", "gamma": 5.5}
    },
    "Mamposteria": {
        "ladrillo_hueco_portante": {"nombre": "Ladrillo hueco portante", "gamma": 10.0},
        "ladrillo_comun": {"nombre": "Ladrillo común", "gamma": 16.0}
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
            "nombre": "Chapa ondulada + fijaciones",
            "q": 0.10
        },
        "teja": {
            "nombre": "Teja cerámica con entablonado",
            "q": 1.00
        }
    },

    "Cielorrasos": {
        "yeso_suspendido": {
            "nombre": "Cielorraso suspendido de yeso",
            "q": 0.35
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

    def resumen_texto(self):
        lineas = []
        D = {}
        L = {}

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
            else:
                L[i["unidad"]] = L.get(i["unidad"], 0) + i["valor"]

            lineas.append("")

        lineas.append("-" * 60)

        for u, v in D.items():
            lineas.append(f"D total = {v:.2f} {u}")
        for u, v in L.items():
            lineas.append(f"L total = {v:.2f} {u}")

        return "\n".join(lineas)
    
    # -------------------------------
    # NUEVO MÉTODO: combinaciones CIRSOC
    # -------------------------------
    def combinaciones_CIRSOC(self):
        """Calcula combinaciones típicas 1.4D y 1.2D+1.6L"""
        D_total = sum(i["valor"] for i in self.items if i["tipo"] == "D")
        L_total = sum(i["valor"] for i in self.items if i["tipo"] == "L")
        
        combos = {}
        combos["1.4D"] = 1.4 * D_total
        combos["1.2D+1.6L"] = 1.2 * D_total + 1.6 * L_total
        
        return combos


# ===============================
# ACTIVACIÓN DE ELEMENTOS
# ===============================

# Cubiertas completas con componentes y sobrecarga
CUBIERTAS = {
    1: {
        "nombre": "Cubierta liviana",
        "activo": 1,
        "b": 3.0,
        "componentes": [
            (SISTEMAS["Cubiertas"]["chapa_ondulada"]["nombre"], SISTEMAS["Cubiertas"]["chapa_ondulada"]["q"]),
            (SISTEMAS["Cielorrasos"]["yeso_suspendido"]["nombre"], SISTEMAS["Cielorrasos"]["yeso_suspendido"]["q"]),
            (GAMMA["Madera"]["pino"]["nombre"] + " 0.04 m", carga_superficial(GAMMA["Madera"]["pino"]["gamma"], 0.04))
        ],
        "sobrecarga": SOBRECARGAS["cubierta_acceso_poco_frecuente"]
    },
    2: {
        "nombre": "Cubierta pesada",
        "activo": 0,
        "b": 3.0,
        "componentes": [
            ("Teja cerámica con entablonado", SISTEMAS["Cubiertas"]["teja"]["q"]),
            ("Estructura de madera 0.05 m", carga_superficial(GAMMA["Madera"]["pino"]["gamma"], 0.05))
        ],
        "sobrecarga": SOBRECARGAS["cubierta_acceso_poco_frecuente"]
    },
    3: {
        "nombre": "Cubierta terraza",
        "activo": 0,
        "b": 3.0,
        "componentes": [
            ("Loseta cerámica 0.02 m", carga_superficial(GAMMA["Hormigon"]["armado"]["gamma"], 0.02)),
            ("Contrapiso 0.05 m", carga_superficial(GAMMA["Hormigon"]["armado"]["gamma"], 0.05))
        ],
        "sobrecarga": SOBRECARGAS["terraza"]
    }
}

# Forjados completos con componentes y sobrecarga
FORJADOS = {
    1: {
        "nombre": "Forjado completo",
        "activo": 1,
        "b": 3.0,
        "componentes": [
            (SISTEMAS["Pisos"]["ceramica"]["nombre"], SISTEMAS["Pisos"]["ceramica"]["q"]),
            carpeta_nivelacion(0.025),  # 2.5 cm
            contrapiso_cascotes(0.05),  # 5 cm
            (GAMMA["Hormigon"]["armado"]["nombre"] + " 0.08 m",
             carga_superficial(GAMMA["Hormigon"]["armado"]["gamma"], 0.08))
        ],
        "sobrecarga": SOBRECARGAS["vivienda"]
    },
    2: {
        "nombre": "Losa simple",
        "activo": 0, #1si, 0no
        "b": 4.0,
        "componentes": [
            ("Losa hormigón armado 0.10 m", carga_superficial(GAMMA["Hormigon"]["armado"]["gamma"], 0.10))
        ],
        "sobrecarga": SOBRECARGAS["vivienda"]
    }
}


# Muros reutilizables (plantillas)
def M1(e, h):
    gamma = GAMMA["Mamposteria"]["ladrillo_comun"]["gamma"]
    return {
        "nombre": "Muro de ladrillo común",
        "valor": carga_lineal_muro(gamma, e, h),
        "unidad": "kN/m",
        "detalle": f"γ={gamma} kN/m3 · e={e} m · h={h} m"
    }

def M2(e, h):
    gamma = GAMMA["Mamposteria"]["ladrillo_hueco_portante"]["gamma"]
    return {
        "nombre": "Muro de ladrillo hueco",
        "valor": carga_lineal_muro(gamma, e, h),
        "unidad": "kN/m",
        "detalle": f"γ={gamma} kN/m3 · e={e} m · h={h} m"
    }

# Lista de muros del proyecto
MUROS = {
    "M1": {"func": M1, "activo": 1, "e": 0.24, "h": 2.60},
    "M2": {"func": M2, "activo": 1, "e": 0.20, "h": 2.60},
}


# ===============================
# MAIN
# ===============================
if __name__ == "__main__":

    analisis = AnalisisCargas("Cargas Viga 01") #cambiar nombre proyecto
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
        if m["activo"]:
            muro = m["func"](m["e"], m["h"])
            analisis.agregar(muro["nombre"], "D", muro["valor"], muro["unidad"], muro["detalle"])

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
    SALIDAS_DIR = Path("salidas")
    SALIDAS_DIR.mkdir(exist_ok=True)
    fecha_hora = datetime.now().strftime("%d-%m-%Y_%H%M")
    archivo_salida = SALIDAS_DIR / f"{analisis.nombre.replace(' ', '_')}_{fecha_hora}.txt"

    with open(archivo_salida, "w", encoding="utf-8") as f:
        f.write(resumen_completo)

    print(f"\nResumen guardado en: {archivo_salida}")



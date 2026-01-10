import json
import math
from pathlib import Path


class DisenadorViga:
    def __init__(self, nombre, b_cm, h_cm, rec_cm=None, fc_MPa=25, fy_MPa=420):
        self.nombre = nombre
        self.b = b_cm / 100.0
        self.h = h_cm / 100.0
        # recubrimiento fijo de 40 mm si no se pasa otro
        self.rec = (rec_cm / 100.0) if rec_cm is not None else 0.04
        self.fc = fc_MPa
        self.fy = fy_MPa
        self.d = self.h - self.rec - 0.01

        # φ
        self.phi = {
            "flexion_traccion": 0.90,
            "compresion_espiral": 0.70,
            "compresion_otro": 0.65,
            "corte": 0.75,
            "torsion": 0.75,
            "aplastamiento": 0.65,
            "anclaje_postesado": 0.85,
            "bielas_puntales": 0.75,
            "pretesado_transferencia": 0.75,
            "pretesado_anclaje": (0.75, 0.90)
        }

        self.barras_comerciales = {
            6: 0.222, 8: 0.395, 10: 0.617, 12: 0.888,
            16: 1.578, 20: 2.466, 25: 3.853
        }

    def _fc_key(self):
        fc = self.fc
        if fc < 25:
            return "fc_20"
        elif fc < 30:
            return "fc_25"
        else:
            return "fc_30"

    def leer_tramos(viga_json):
        tramos = []
        for tramo in viga_json["tramos"]:
            if tramo.get("es_voladizo", False):
                pi = None
            else:
                pi = tramo.get("puntos_inflexion_m", [])
            tramos.append({
                "id": tramo["id"],
                "longitud": tramo["longitud_m"],
                "es_voladizo": tramo.get("es_voladizo", False),
                "puntos_inflexion_m": pi
            })
        return tramos

    def _kd_por_tu_tabla(self, Mn_kNm):
        d_m = self.d
        bw_m = self.b
        Mn_MNm = abs(Mn_kNm) / 1000.0
        base = Mn_MNm / bw_m
        if base <= 0:
            return float('inf')
        return d_m / math.sqrt(base)

    def calcular_as_necesaria(self, Mu_kNm, coef_kd):
        phi = self.phi["flexion_traccion"]
        Mn_kNm = abs(Mu_kNm) / phi

        fc_key = self._fc_key()
        tabla = coef_kd["coeficientes_flexion"][fc_key]
        kd_in = self._kd_por_tu_tabla(Mn_kNm)
        claves = sorted([float(k) for k in tabla.keys()])
        kd_bajo = max([k for k in claves if k <= kd_in], default=claves[0])
        kd_alto = min([k for k in claves if k >= kd_in], default=claves[-1])

        if kd_bajo == kd_alto:
            fila = tabla[str(kd_bajo)]
            kd_sel = str(kd_bajo)
        else:
            datos_bajo = tabla[str(kd_bajo)]
            datos_alto = tabla[str(kd_alto)]
            t = (kd_in - kd_bajo) / (kd_alto - kd_bajo)
            def interp(campo):
                return datos_bajo[campo] + (datos_alto[campo] - datos_bajo[campo]) * t
            fila = {"Ke": interp("Ke"), "Kc": interp("Kc"), "Kz": interp("Kz")}
            kd_sel = f"{kd_bajo:.3f}-{kd_alto:.3f}"

        Mn_MNm = abs(Mn_kNm) / 1000.0
        Ke = fila["Ke"]
        As_req_cm2 = Ke * (Mn_MNm / self.d)

        As_min_1 = (math.sqrt(self.fc) / (4.0 * self.fy)) * self.b * self.d * 10000.0
        As_min_2 = (1.4 / self.fy) * self.b * self.d * 10000.0
        As_min_cm2 = max(As_min_1, As_min_2)

        beta1 = 0.85 if self.fc <= 28 else 0.80
        eps_cu = 0.003
        Es = 200000.0
        eps_y = self.fy / Es
        rho_bal = 0.85 * beta1 * (self.fc / self.fy) * (eps_cu / (eps_cu + eps_y))
        As_bal_cm2 = rho_bal * self.b * self.d * 10000.0
        As_max_cm2 = 0.75 * As_bal_cm2
        As_final_cm2 = min(max(As_req_cm2, As_min_cm2), As_max_cm2)

        if As_final_cm2 < As_min_cm2:
            estado = "❌ Cuantía insuficiente (menor que mínima)"
        elif As_final_cm2 > As_max_cm2:
            estado = "❌ Cuantía excesiva (mayor que máxima)"
        else:
            estado = "✅ Cuantía dentro del rango permitido"
        return {
            "Mu_kNm": Mu_kNm,
            "phi_flexion": phi,
            "Mn_kNm": Mn_kNm,
            "As_req_cm2": As_req_cm2,
            "As_min_cm2": As_min_cm2,
            "As_bal_cm2": As_bal_cm2,
            "As_max_cm2": As_max_cm2,
            "As_final_cm2": As_final_cm2,
            "estado": estado,
            "kd_in": kd_in,
            "kd_sel": kd_sel,
            "fc_key": fc_key,
            "fila": fila
        }

    def seleccionar_armadura(self, as_cm2):
        opciones = [6, 8, 10, 12, 16, 20, 25]
        mejor_opcion = None
        area_min_exceso = float('inf')

        sep_min = 2.0
        b_util_cm = self.b * 100 - 2 * self.rec
        max_barras_simple = 8

        for diam in opciones:
            area_barra = math.pi * (diam / 10)**2 / 4
            n_barras = max(2, math.ceil(as_cm2 / area_barra))
            area_total = n_barras * area_barra
            exceso = area_total - as_cm2

            ancho_requerido = n_barras * (diam / 10) + (n_barras - 1) * sep_min
            entra_simple = (ancho_requerido <= b_util_cm) and (n_barras <= max_barras_simple)
            capa = "simple" if entra_simple else "doble"

            penalizacion = 0.0 if entra_simple else 1e6
            exceso_eval = exceso + penalizacion
            rec_geom_mm = self.rec * 1000                # recubrimiento geométrico en mm
            rec_ef_mm = rec_geom_mm + diam / 2.0         # recubrimiento efectivo
            Ld_mm = max(12 * diam, self.d * 1000)        # longitud de desarrollo

            if (exceso_eval < area_min_exceso) or (
                abs(exceso_eval - area_min_exceso) < 1e-6 and 
                (mejor_opcion is None or diam > mejor_opcion["diam"])
            ):
                area_min_exceso = exceso_eval
                mejor_opcion = {
                    "n": n_barras,
                    "diam": diam,
                    "area_total_cm2": area_total,
                    "exceso_cm2": exceso,
                    "capa": capa,
                    "ancho_requerido_cm": ancho_requerido,
                    "b_util_cm": b_util_cm,
                    "recubrimiento_efectivo_mm": rec_ef_mm,
                    "Ld_mm": Ld_mm      
                }
        return mejor_opcion
    
    def mapear_armaduras(self, viga_json, coef_kd):
        resultados = []
        for tramo in viga_json["tramos"]:
            if tramo.get("es_voladizo", False):
                Mu = tramo.get("Mu_kNm", 0.0)
                datos = self.calcular_as_necesaria(Mu, coef_kd)
                arm_sup = self.seleccionar_armadura(datos["As_final_cm2"])

                rec_geom_mm = self.rec * 1000
                rec_ef_mm = rec_geom_mm + arm_sup["diam"] / 2.0
                Ld_mm = max(12 * arm_sup["diam"], self.d * 1000)

                resultados.append({
                    "tramo": tramo["id"],
                    "tipo": "voladizo",
                    "armadura_superior": arm_sup,
                    "armadura_inferior": "mínima (2Ø10)",
                    "recubrimiento_efectivo_mm": rec_ef_mm,
                    "Ld_mm": Ld_mm,
                    "nota": "Armadura superior prolongada hasta extremo libre"
                })
            else:
                Mu_campo = tramo.get("Mu_kNm_campo", 0.0)
                Mu_apoyo_izq = tramo.get("Mu_kNm_apoyo_izq", 0.0)
                Mu_apoyo_der = tramo.get("Mu_kNm_apoyo_der", 0.0)

                datos_inf = self.calcular_as_necesaria(Mu_campo, coef_kd)
                arm_inf = self.seleccionar_armadura(datos_inf["As_final_cm2"])

                Mu_sup = max(abs(Mu_apoyo_izq), abs(Mu_apoyo_der))
                datos_sup = self.calcular_as_necesaria(Mu_sup, coef_kd)
                arm_sup = self.seleccionar_armadura(datos_sup["As_final_cm2"])

                rec_geom_mm = self.rec * 1000
                rec_ef_mm = rec_geom_mm + arm_sup["diam"] / 2.0
                Ld_mm = max(12 * arm_sup["diam"], self.d * 1000)

                resultados.append({
                    "tramo": tramo["id"],
                    "tipo": "continuo",
                    "armadura_inferior": arm_inf,
                    "armadura_superior": arm_sup,
                    "puntos_inflexion_m": tramo.get("puntos_inflexion_m", []),
                    "recubrimiento_efectivo_mm": rec_ef_mm,
                    "Ld_mm": Ld_mm,
                    "nota": "Inferior prolongada más allá del PI; superior hasta el PI"
                })
        return resultados

    def calcular_estribos(self, Vu_kN, diam_mm=6, ramas=2):
    # Conversión
        Vu_N = Vu_kN * 1000
        b_mm = self.b * 1000
        d_mm = self.d * 1000
        phi_corte = self.phi["corte"]

        # Hormigón
        sqrt_fc = min(math.sqrt(self.fc), 8.3)
        Vc_N = 0.17 * sqrt_fc * b_mm * d_mm

        # Mínimo normativo si el corte es bajo
        if Vu_N < (phi_corte * Vc_N / 2):
            return {
                "diam_mm": 6,
                "ramas": 2,
                "s_cm": 20,
                "modo": "mínimo (d/2 y 20 cm)",
                "Vc_kN": Vc_N/1000,
                "Vs_req_kN": 0.0,
                "φVn_kN": phi_corte * Vc_N / 1000
            }

        # Vs requerido
        Vs_N = (Vu_N / phi_corte) - Vc_N
        Vs_N = max(Vs_N, 0.0)

        # Área de ramas
        area_barra_mm2 = math.pi * (diam_mm**2) / 4
        Asv_mm2 = ramas * area_barra_mm2

        # Separación teórica
        if Vs_N > 0:
            s_mm = (Asv_mm2 * self.fy * d_mm) / Vs_N
        else:
            s_mm = d_mm / 2  # máximo permitido

        # Límites normativos
        s_mm = min(s_mm, d_mm / 2, 300)  # d/2 y 30 cm
        s_cm = s_mm / 10
        s_cm = max(6, min(round(s_cm / 2) * 2, 30))  # redondeo a múltiplos de 2 cm

        # Chequeo de cuantía mínima
        area_barra_cm2 = math.pi * (6 / 10)**2 / 4
        Asv_min_cm2 = 2 * area_barra_cm2
        av_s_min = Asv_min_cm2 / 20.0

        Asv_adopt_cm2 = ramas * (math.pi * (diam_mm / 10)**2 / 4)
        av_s_adopt = Asv_adopt_cm2 / s_cm

        if av_s_adopt < av_s_min:
            diam_mm = 6
            ramas = 2
            s_cm = 20
            modo = "forzado a mínimo Ø6 c/20 (2 ramas)"
        else:
            modo = "resistencia requerida"

        return {
            "diam_mm": diam_mm,
            "ramas": ramas,
            "s_cm": int(s_cm),
            "modo": modo,
            "Vc_kN": Vc_N/1000,
            "Vs_req_kN": Vs_N/1000,
            "φVn_kN": phi_corte * (Vc_N + Vs_N)/1000
        }


    def _Vc_N(self):
        b_mm = self.b * 1000
        d_mm = self.d * 1000
        sqrt_fc = min(math.sqrt(self.fc), 8.3)
        return 0.17 * sqrt_fc * b_mm * d_mm

    def _xcrit_continuo(self, Vu_izq_kN, Vu_der_kN, L_m):
        Vc_N = self._Vc_N()
        phiVc_kN = self.phi["corte"] * Vc_N / 1000.0

        denom = (Vu_der_kN - Vu_izq_kN)
        if abs(denom) < 1e-9:
            # Corte casi constante: si supera umbral, toda la viga es densa; si no, mínima
            if Vu_izq_kN <= phiVc_kN:
                xcrit_izq = 0.0
            else:
                xcrit_izq = L_m
        else:
            xcrit_izq = L_m * (phiVc_kN - Vu_izq_kN) / denom

        xcrit_izq = max(0.0, min(L_m, xcrit_izq))
        xcrit_der = L_m - xcrit_izq
        return xcrit_izq, xcrit_der, phiVc_kN


    def _xcrit_voladizo(self, Vu_emp_kN, L_m):
        Vc_N = self._Vc_N()
        phiVc_kN = self.phi["corte"] * Vc_N / 1000.0
        if Vu_emp_kN <= phiVc_kN:
            return 0.0, phiVc_kN
        xcrit = L_m * (1.0 - phiVc_kN / Vu_emp_kN)
        xcrit = max(0.0, min(L_m, xcrit))
        return xcrit, phiVc_kN


    def zonificar_estribos_continuo(self, L_m, Vu_izq_kN, Vu_der_kN, diam_mm=6, ramas=2):
        # Distancias críticas
        xcrit_izq, xcrit_der, phiVc_kN = self._xcrit_continuo(Vu_izq_kN, Vu_der_kN, L_m)

        # Longitud mínima densa
        lmin = min(self.d, 0.5)  # m

        # Zonas densas en apoyos
        zA_ini, zA_fin = 0.0, max(lmin, xcrit_izq)
        zC_ini, zC_fin = L_m - max(lmin, xcrit_der), L_m

        # Corrección si se solapan (apoyos se comen la central)
        if zA_fin > zC_ini:
            # Reducir densos a lmin y recalcular central
            zA_fin = lmin
            zC_ini = L_m - lmin
        # Zona central
        zB_ini, zB_fin = zA_fin, zC_ini

        zonas = []

        def Vu_x(x):
            return Vu_izq_kN + (Vu_der_kN - Vu_izq_kN) * (x / L_m)

        # Zona A (apoyo izq)
        xA = (zA_ini + zA_fin) / 2.0
        zonas.append({
            "zona": "apoyo_izq",
            "desde_m": zA_ini, "hasta_m": zA_fin,
            "Vu_kN_ref": Vu_x(xA),
            "estribos": self.calcular_estribos(Vu_x(xA), diam_mm, ramas)
        })
        # Zona B (central) si existe
        if zB_fin > zB_ini + 1e-6:
            xB = (zB_ini + zB_fin) / 2.0
            zonas.append({
                "zona": "central",
                "desde_m": zB_ini, "hasta_m": zB_fin,
                "Vu_kN_ref": Vu_x(xB),
                "estribos": self.calcular_estribos(Vu_x(xB), diam_mm, ramas)
            })
        # Zona C (apoyo der)
        xC = (zC_ini + zC_fin) / 2.0
        zonas.append({
            "zona": "apoyo_der",
            "desde_m": zC_ini, "hasta_m": zC_fin,
            "Vu_kN_ref": Vu_x(xC),
            "estribos": self.calcular_estribos(Vu_x(xC), diam_mm, ramas)
        })
        return {
            "phiVc_kN": phiVc_kN,
            "zonas": zonas
        }



    def zonificar_estribos_voladizo(self, L_m, Vu_emp_kN, diam_mm=6, ramas=2):
        xcrit, phiVc_kN = self._xcrit_voladizo(Vu_emp_kN, L_m)
        lmin = min(self.d, 0.5)

        zA_ini, zA_fin = 0.0, max(lmin, xcrit)   # denso desde empotramiento
        zB_ini, zB_fin = zA_fin, L_m            # resto del voladizo

        zonas = []

        def Vu_x(x):
            return Vu_emp_kN * (1.0 - x / L_m)

        # Zona A (empotramiento)
        xA = (zA_ini + zA_fin) / 2.0
        zonas.append({
            "zona": "empotramiento",
            "desde_m": zA_ini, "hasta_m": zA_fin,
            "Vu_kN_ref": Vu_x(xA),
            "estribos": self.calcular_estribos(Vu_x(xA), diam_mm, ramas)
        })
        # Zona B (resto)
        if zB_fin > zB_ini + 1e-6:
            xB = (zB_ini + zB_fin) / 2.0
            zonas.append({
                "zona": "extremo_libre",
                "desde_m": zB_ini, "hasta_m": zB_fin,
                "Vu_kN_ref": Vu_x(xB),
                "estribos": self.calcular_estribos(Vu_x(xB), diam_mm, ramas)
            })
        return {
            "phiVc_kN": phiVc_kN,
            "zonas": zonas
        }

# =========================================================
# FUNCIÓN PLANILLA (VOLADIZO + INTERIOR)
# =========================================================
def generar_planilla(
    viga_data,
    nombre_viga,
    id_tramo,
    coef_kd,
    L_viga,
    b_cm,
    h_cm,
    rec_cm,
    fc_MPa,
    fy_MPa,
    es_voladizo=False,
    puntos_inflexion_m=None
):

    lineas = []
    zonas_corte = {"zonas": []}

    v = DisenadorViga(nombre_viga, b_cm, h_cm, rec_cm, fc_MPa, fy_MPa)

    # ---- Cuantías ----
    as_izq = v.calcular_as_necesaria(viga_data["m_izq"], coef_kd)
    as_tra = v.calcular_as_necesaria(viga_data["m_tra"], coef_kd)
    as_der = v.calcular_as_necesaria(viga_data["m_der"], coef_kd)

    arm_sup_izq = v.seleccionar_armadura(as_izq["As_final_cm2"])
    arm_inf_tra = v.seleccionar_armadura(as_tra["As_final_cm2"])
    arm_sup_der = v.seleccionar_armadura(as_der["As_final_cm2"])

    def lon_gancho(d): return 15 * d / 1000
    def lon_anclaje(d): return 40 * d / 1000

    # ---- Encabezado ----
    lineas.append("=" * 80)
    tipo = "VOLADIZO" if es_voladizo else "TRAMO INTERIOR"
    lineas.append(f"PLANILLA DE DOBLADO – {nombre_viga} – {tipo} – TRAMO {id_tramo}")
    lineas.append("=" * 80)
    lineas.append(f"{'Pos':<4} | {'Cant':<5} | {'Ø':<4} | {'Detalle':<22} | {'L (m)':<7} | {'Peso (kg)':<10}")
    lineas.append("-" * 80)

    pos = 1
    total_peso = 0.0

    # =====================================================
    # VOLADIZO
    # =====================================================
    if es_voladizo:

        # Superior principal
        n, d = arm_sup_izq["n"], arm_sup_izq["diam"]
        L = L_viga + lon_anclaje(d)
        P = L * n * v.barras_comerciales[d]
        lineas.append(f"{pos:<4} | {n:<5} | {d:<4} | Sup. voladizo       | {L:<7.2f} | {P:<10.2f}")
        total_peso += P; pos += 1

        # Inferior mínima
        n, d = 2, 10
        L = L_viga
        P = L * n * v.barras_comerciales[d]
        lineas.append(f"{pos:<4} | {n:<5} | {d:<4} | Inf. mínima         | {L:<7.2f} | {P:<10.2f}")
        total_peso += P; pos += 1

        # Estribos
        zonas = v.zonificar_estribos_voladizo(L_viga, viga_data["v_izq"])
        for z in zonas["zonas"]:
            e = z["estribos"]
            n_est = math.ceil((z["hasta_m"] - z["desde_m"]) / (e["s_cm"] / 100))
            L = 2*(v.b-0.06) + 2*(v.h-0.06) + 0.10
            P = n_est * L * v.barras_comerciales[e["diam_mm"]]
            lineas.append(
                f"{pos:<4} | {n_est:<5} | {e['diam_mm']:<4} | Estribo {z['zona']:<12} | {L:<7.2f} | {P:<10.2f}"
            )
            total_peso += P; pos += 1

        # ---- Resumen de zonificación de voladizo ----
        z_emp = next((z for z in zonas["zonas"] if "empot" in z["zona"]), None)
        z_lib = next((z for z in zonas["zonas"] if "libre" in z["zona"]), None)

        if z_emp:
            if z_lib:
                lineas.append(
                    f"Nota: Zonificación de estribos (voladizo) → "
                    f"Empotramiento={z_emp['desde_m']:.2f}-{z_emp['hasta_m']:.2f} m | "
                    f"Extremo libre={z_lib['desde_m']:.2f}-{z_lib['hasta_m']:.2f} m"
                )
            else:
                lineas.append(
                    f"Nota: Zonificación de estribos (voladizo) → "
                    f"Empotramiento={z_emp['desde_m']:.2f}-{z_emp['hasta_m']:.2f} m"
                )


    # =====================================================
    # TRAMO INTERIOR
    # =====================================================
    else:

        if puntos_inflexion_m and len(puntos_inflexion_m) == 2:
            PI_izq, PI_der = puntos_inflexion_m
        else:
            PI_izq, PI_der = 0.15 * L_viga, 0.85 * L_viga

        # Superior apoyo izq
        n, d = arm_sup_izq["n"], arm_sup_izq["diam"]
        L = PI_izq + lon_anclaje(d)
        P = L * n * v.barras_comerciales[d]
        lineas.append(f"{pos:<4} | {n:<5} | {d:<4} | Sup. apoyo izq      | {L:<7.2f} | {P:<10.2f}")
        total_peso += P; pos += 1

        # Inferior tramo
        n, d = arm_inf_tra["n"], arm_inf_tra["diam"]
        L = (PI_der - PI_izq) + 2 * (12 * d / 1000)
        P = L * n * v.barras_comerciales[d]
        lineas.append(f"{pos:<4} | {n:<5} | {d:<4} | Inf. tramo          | {L:<7.2f} | {P:<10.2f}")
        total_peso += P; pos += 1

        # Superior apoyo der
        n, d = arm_sup_der["n"], arm_sup_der["diam"]
        L = (L_viga - PI_der) + lon_anclaje(d)
        P = L * n * v.barras_comerciales[d]
        lineas.append(f"{pos:<4} | {n:<5} | {d:<4} | Sup. apoyo der      | {L:<7.2f} | {P:<10.2f}")
        total_peso += P; pos += 1

        # Auxiliar superior
        n, d = 2, 8
        L = L_viga + 2 * lon_gancho(d)
        P = L * n * v.barras_comerciales[d]
        lineas.append(f"{pos:<4} | {n:<5} | {d:<4} | Auxiliar superior   | {L:<7.2f} | {P:<10.2f}")
        total_peso += P; pos += 1

        # Estribos
        zonas = v.zonificar_estribos_continuo(L_viga, viga_data["v_izq"], viga_data["v_der"])
        for z in zonas["zonas"]:
            e = z["estribos"]
            n_est = math.ceil((z["hasta_m"] - z["desde_m"]) / (e["s_cm"] / 100))
            L = 2*(v.b-0.06) + 2*(v.h-0.06) + 0.10
            P = n_est * L * v.barras_comerciales[e["diam_mm"]]
            lineas.append(
                f"{pos:<4} | {n_est:<5} | {e['diam_mm']:<4} | Estribo {z['zona']:<12} | {L:<7.2f} | {P:<10.2f}"
            )
            total_peso += P; pos += 1
        # ---- RESUMEN DE ZONIFICACIÓN (ACÁ, NO EN OTRO LADO) ----
        z_izq = next((z for z in zonas["zonas"] if "izq" in z["zona"]), None)
        z_cen = next((z for z in zonas["zonas"] if z["zona"] == "central"), None)
        z_der = next((z for z in zonas["zonas"] if "der" in z["zona"]), None)

        if z_izq and z_der:
            if z_cen:
                lineas.append(
                    f"Nota: Zonificación de estribos → "
                    f"A={z_izq['desde_m']:.2f}-{z_izq['hasta_m']:.2f} m | "
                    f"B={z_cen['desde_m']:.2f}-{z_cen['hasta_m']:.2f} m | "
                    f"C={z_der['desde_m']:.2f}-{z_der['hasta_m']:.2f} m"
                )
            else:
                lineas.append(
                    f"Nota: Zonificación de estribos → "
                    f"A={z_izq['desde_m']:.2f}-{z_izq['hasta_m']:.2f} m | "
                    f"C={z_der['desde_m']:.2f}-{z_der['hasta_m']:.2f} m (sin zona central)"
        )
    # ---- Cierre ----
    lineas.append("-" * 80)
    lineas.append(f"{'PESO TOTAL ACERO (kg):':<55} {total_peso:.2f}")
    lineas.append("=" * 80)
    lineas.append("")
    lineas.append("NOTAS TÉCNICAS:")
    lineas.append(f"- Cuantía requerida (As req): {as_tra['As_req_cm2']:.2f} cm²")
    lineas.append(f"- Cuantía mínima (As min): {as_tra['As_min_cm2']:.2f} cm²")
    lineas.append(f"- Cuantía máxima (As max): {as_tra['As_max_cm2']:.2f} cm²")
    lineas.append(f"- Cuantía adoptada: {arm_inf_tra['area_total_cm2']:.2f} cm²")
    lineas.append(f"- Estado de diseño: {as_tra['estado']}")

    if 'kd_sel' in as_tra:
        lineas.append(f"- kd utilizado: {as_tra['kd_sel']} (fila {as_tra.get('fila','-')})")

    return "\n".join(lineas)


# =========================================================
# PROGRAMA PRINCIPAL
# =========================================================
BASE = Path(__file__).parent

with open(BASE / "datos" / "moments_input.json", encoding="utf-8") as f:
    datos_vigas = json.load(f)

with open(BASE / "datos" / "coeficientes_kd.json", encoding="utf-8") as f:
    coef_kd = json.load(f)

salidas = BASE / "salidas"
salidas.mkdir(exist_ok=True)

for viga_id, viga in datos_vigas.items():

    b = viga["b_cm"]
    h = viga["h_cm"]
    rec = viga["recubrimiento_cm"]
    fc = viga.get("fc_MPa", 21.0)
    fy = viga.get("fy_MPa", 420.0)

    for tramo in viga["tramos"]:

        if tramo["es_voladizo"]:
            datos_tramo = {
                "m_izq": tramo["Mu_kNm"],
                "m_tra": tramo["Mu_kNm"],
                "m_der": 0.0,
                "v_izq": tramo["Vu_kN_emp"],
                "v_der": 0.0
            }
        else:
            datos_tramo = {
                "m_izq": tramo["Mu_kNm_apoyo_izq"],
                "m_tra": tramo["Mu_kNm_campo"],
                "m_der": tramo["Mu_kNm_apoyo_der"],
                "v_izq": tramo["Vu_kN_izq"],
                "v_der": tramo["Vu_kN_der"]
            }

        texto = generar_planilla(
            viga_data=datos_tramo,
            nombre_viga=viga_id,
            id_tramo=tramo["id"],
            coef_kd=coef_kd,
            L_viga=tramo["longitud_m"],
            b_cm=b,
            h_cm=h,
            rec_cm=rec,
            fc_MPa=fc,
            fy_MPa=fy,
            es_voladizo=tramo["es_voladizo"],
            puntos_inflexion_m=tramo.get("puntos_inflexion_m")
        )

        with open(salidas / f"planilla_{viga_id}_{tramo['id']}.txt", "w", encoding="utf-8") as f:
            f.write(texto)


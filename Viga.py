import json
import math
from pathlib import Path


class DisenadorViga:
    def __init__(self, nombre, b_cm, h_cm, rec_cm=None, fc_MPa=25, fy_MPa=420,rec_sup_cm=None):
        self.nombre = nombre
        self.b = b_cm / 100.0
        self.h = h_cm / 100.0
        # recubrimiento fijo de 40 mm si no se pasa otro
        self.rec = (rec_cm / 100.0) if rec_cm is not None else 0.04
        self.recubrimiento_sup_cm = rec_sup_cm if rec_sup_cm is not None else rec_cm
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

    def _kd_fisico(self, Mn_kNm):
        # Convertir a MN·m
        Mn_MNm = Mn_kNm / 1000.0
        kd = self.d / math.sqrt(Mn_MNm / self.b)
        return kd

    def _kd_lim_por_fc(self, fc):
        # kd límite según resistencia del hormigón
        # ejemplo: tabla simplificada
        if fc <= 25: return 0.45
        elif fc <= 30: return 0.50
        else: return 0.55
   
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

    def _fila_kd_interpolada(self, kd_in, coef_kd, dato_externo=None):
        fc_key = self._fc_key()
        tabla = coef_kd["coeficientes_flexion"][fc_key]
        claves = sorted([float(k) for k in tabla.keys()])
        kd_min, kd_max = claves[0], claves[-1]
        # Saturación en bordes
        if kd_in <= kd_min:
            fila = tabla[f"{kd_min:.3f}"]
            modo = "saturado_inferior"
            info = {}
            # Límite inferior: aplicar kd × dato si corresponde
            if dato_externo is not None:
                info["kd_por_dato"] = kd_in * dato_externo
            return fila, kd_min, modo, info

        if kd_in >= kd_max:
            fila = tabla[f"{kd_max:.3f}"]
            return fila, kd_max, "saturado_superior", {}

        # Interpolación lineal
        kd_bajo = max([k for k in claves if k <= kd_in], default=claves[0])
        kd_alto = min([k for k in claves if k >= kd_in], default=claves[-1])

        if abs(kd_bajo - kd_alto) < 1e-12:
            fila = tabla[str(kd_bajo)]
            return fila, kd_bajo, "exacto", {}

        datos_bajo = tabla[str(kd_bajo)]
        datos_alto = tabla[str(kd_alto)]
        t = (kd_in - kd_bajo) / (kd_alto - kd_bajo)

        def interp(campo):
            return datos_bajo[campo] + (datos_alto[campo] - datos_bajo[campo]) * t

        fila = {"Ke": interp("Ke"), "Kc": interp("Kc"), "Kz": interp("Kz")}
        return fila, kd_in, "interpolado", {"intervalo": [kd_bajo, kd_alto], "t": t}
    
    def _fs_prima_cirsoc(self, k, k_star, x, eps_c=0.003, eps_c_star=0.003, Es_MPa=200000.0):
        # f_s' = Es * [ (1 - k* x) * eps_c* - (1 - k x) * eps_c ]
        eps_sp = (1 - k_star * x) * eps_c_star - (1 - k * x) * eps_c
        fs_p = Es_MPa * eps_sp
        fs_p = min(fs_p, self.fy)  # tope por fy
        return {"eps_s_prima": eps_sp, "fs_prima_MPa": fs_p}
    
    def _lookup_ke_kep(self, tabla_comp, k_ratio, x_rel, fs_target):
        rows = tabla_comp["rows"]

        # Buscar la fila más cercana en kd_ratio y x_rel
        row = min(
            rows,
            key=lambda r: abs(r["kd_kd"] - k_ratio) + abs(r["x"] - x_rel)
        )

        # Localizar columna más cercana en fs_cols
        fs_cols = tabla_comp["fs_cols"]
        col_idx = min(range(len(fs_cols)), key=lambda i: abs(fs_cols[i] - fs_target))

        # Proteger índice por si la fila tiene menos columnas
        if col_idx >= len(row["ke"]):
            col_idx = len(row["ke"]) - 1

        ke_val = row["ke"][col_idx]      # coeficiente Ke
        kep_val = row["ke_p"][col_idx]   # coeficiente Ke'

        return ke_val, kep_val, {"kd_ratio": row["kd_kd"], "x": row["x"], "col_idx": col_idx}
    
    def _as_por_mn(self, Mu_kNm, phi, ke_cm2_por_MN_por_m, kep_cm2_por_MN_por_m):
        Mn_MNm = abs(Mu_kNm) / phi / 1000.0  # kNm -> MN·m
        As_cm2  = ke_cm2_por_MN_por_m  * (Mn_MNm / self.d)
        Asp_cm2 = kep_cm2_por_MN_por_m * (Mn_MNm / self.d)
        return round(As_cm2, 2), round(Asp_cm2, 2)
    
    def calcular_as_necesaria(
        self, Mu_kNm, coef_kd, tabla_comp=None,
        dato_externo=None, x_rel=None, k_star=None, fs_politica=None, debug_compresion=False
    ):
        import math

        # 0) Parámetros base
        phi = self.phi["flexion_traccion"]
        Mn_kNm = abs(Mu_kNm) / phi
        Mn_MNm = Mn_kNm / 1000.0
        fs_obj = fs_politica or 360  # política/nearest para columna fs

        # 1) kd físico y kd límite (kd*)
        kd_fisico = self._kd_fisico(Mn_kNm)
        kd_lim    = self._kd_lim_por_fc(self.fc)

        # 2) Decidir lookup: simple (kd_físico) vs doble (kd_ratio)
        if kd_fisico >= kd_lim:
            kd_lookup = kd_fisico
            doble_armadura = False
        else:
            kd_lookup = kd_fisico / kd_lim
            doble_armadura = True

        # 3) x_rel para el lookup (geométrico simple: recubrimiento superior / d)
        if x_rel is None:
            try:
                x_rel = (self.recubrimiento_sup_cm / 100.0) / self.d
            except Exception:
                x_rel = None  # si no hay datos, se deja None
        # Inicialización de coeficientes
        Ke = None
        Kep = None
        Kc = None
        Kz = None
        k_used = None
        modo_k = None
        info_k = {}

        # 4) Obtener coeficientes desde tabla de flexión (fc_XX): SOLO Ke, Kc, Kz
        # Si hay tabla_comp y x_rel, usamos lookup bilineal para Ke/Kep (pero Kc/Kz siguen viniendo de flexión)
        if tabla_comp is not None and x_rel is not None:
            # Lookup bilineal (kd_lookup + x_rel) para Ke/Kep
            Ke, Kep_tab, info_k = self._lookup_ke_kep(
                k_ratio=kd_lookup,
                x_rel=x_rel,
                tabla_comp=tabla_comp,
                fs_target=fs_obj
            )
            Kep = Kep_tab
            k_used = kd_lookup
            modo_k = info_k.get("modo", "bilineal")

            # Kc/Kz deben venir de la tabla de flexión (fc_key)
            try:
                fc_key = self._fc_key()  # p.ej. "fc_30"
                fila_flex = self._fila_kd_interpolada(kd_lookup, coef_kd, dato_externo=dato_externo)[0]
                Kc = fila_flex.get("Kc", None)
                Kz = fila_flex.get("Kz", None)
                # Si no existen, dejar None (no 0.0) y documentar
                if Kc is None or Kz is None:
                    info_k["nota"] = "Kc/Kz no presentes en flexión para este kd; se dejan como None."
            except Exception as e:
                info_k["error_flexion"] = f"No se pudo leer Kc/Kz de flexión: {e}"
                Kc, Kz = None, None

        else:
            # Fallback al interpolador anterior (solo kd_lookup) desde flexión
            fila, k_used, modo_k, info_k = self._fila_kd_interpolada(
                kd_lookup, coef_kd, dato_externo=dato_externo
            )
            Ke = fila.get("Ke")
            # Kc/Kz SOLO desde flexión; si faltan, dejar None (no 0.4/1.0 por defecto)
            Kc = fila.get("Kc", None)
            Kz = fila.get("Kz", None)
            if Kc is None or Kz is None:
                info_k["nota"] = "Kc/Kz no presentes en flexión para este kd; se dejan como None."

        # 5) Cálculo de As por tracción y compresión (primera pasada)
        As_req_cm2  = Ke * (Mn_MNm / self.d) if Ke is not None else None
        Asp_req_cm2 = None
        if Kep is not None and doble_armadura:
            Asp_req_cm2 = Kep * (Mn_MNm / self.d)

        # 6) Mínima / balanceada / máxima
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

        # Si As_req_cm2 es None, no podemos evaluar rango; mantenemos None y estado informativo
        # Si As_req_cm2 es None, no podemos evaluar rango; mantenemos None y estado informativo
        if As_req_cm2 is None:
            As_final_cm2 = None
            estado = "⚠️ Ke no disponible; no se puede evaluar cuantía."
        else:
            if doble_armadura:
                # En doble armadura, la adoptada coincide con la requerida
                As_final_cm2 = As_req_cm2
                estado = "ℹ️ Doble armadura activa; menor ductilidad (Kc límite aplicado)"
            else:
                # En flexión simple, sí se controla mínima/máxima
                As_final_cm2 = min(max(As_req_cm2, As_min_cm2), As_max_cm2)
                if As_final_cm2 < As_min_cm2:
                    estado = "❌ Cuantía insuficiente (menor que mínima)"
                elif As_final_cm2 > As_max_cm2:
                    estado = "❌ Cuantía excesiva (mayor que máxima)"
                else:
                    estado = "✅ Cuantía dentro del rango permitido"

        # 8) Salida clara y trazable (Kc/Kz solo de flexión; nunca cero por defecto)
        salida = {
            "Mu_kNm": Mu_kNm,
            "phi_flexion": phi,
            "Mn_kNm": Mn_kNm,
            "As_req_cm2": round(As_req_cm2, 2) if As_req_cm2 is not None else None,
            "Asp_req_cm2": round(Asp_req_cm2, 2) if Asp_req_cm2 is not None else 0.0,
            "As_min_cm2": round(As_min_cm2, 2),
            "As_bal_cm2": round(As_bal_cm2, 2),
            "As_max_cm2": round(As_max_cm2, 2),
            "As_final_cm2": round(As_final_cm2, 2) if As_final_cm2 is not None else None,
            "estado": estado,
            "kd_fisico": round(kd_fisico, 3),
            "kd_lim": round(kd_lim, 3),
            "kd_lookup": round(kd_lookup, 3),
            "doble_armadura": doble_armadura,
            "x_rel": round(x_rel, 4) if x_rel is not None else None,
            "k_usado": k_used,
            "modo_k": modo_k,
            "info_k": info_k,
            "fc_key": self._fc_key(),
            # Kc/Kz vienen de flexión; si faltan, se dejan como None (no 0.0)
            "fila": {"Ke": Ke, "Kc": Kc, "Kz": Kz}
        }
        
        # 9) Activación de compresión (si hay tabla y x_rel) — NO pisa Kc/Kz
        if tabla_comp is not None and x_rel is not None:    
            # Compatibilidad para f_s' (si no hay política fija)
            if fs_politica is None:
                k_star = k_star if k_star is not None else k_used
                comp = self._fs_prima_cirsoc(k_used, k_star, x_rel, eps_c=eps_cu, eps_c_star=eps_cu, Es_MPa=Es)
                fs_target = comp["fs_prima_MPa"]
            else:
                fs_target = fs_politica  # p. ej., 360 MPa por ductilidad

            # Selección de fila más cercana en tabla_compresion
            try:
                fila_match = min(
                    tabla_comp["rows"],
                    key=lambda r: abs(r["kd_kd"] - k_used) + abs(r["x"] - x_rel)
                )
                ke_tab = float(fila_match["ke"][0])
                kep_tab = float(fila_match["ke_p"][0])
            except Exception as e:
                if debug_compresion:
                    print(f"[Compresión] Error en lookup de tabla: {e}")
                ke_tab, kep_tab, fila_match = None, None, None

            # Calcular As desde compresión si hay datos
            As_tab_cm2 = ke_tab * (Mn_MNm / self.d) if ke_tab is not None else None
            Asp_tab_cm2 = kep_tab * (Mn_MNm / self.d) if (kep_tab is not None and doble_armadura) else 0.0

            if debug_compresion:
                print("\n=== DEBUG TRAMO ===")
                print(f"Mu_kNm        : {Mu_kNm:.2f}")
                print(f"phi           : {phi:.3f}")
                print(f"Mn_kNm        : {Mn_kNm:.2f}")
                print(f"Mn_MNm        : {Mn_MNm:.3f}")
                print(f"fc            : {self.fc} MPa")
                print(f"kd_fisico     : {kd_fisico:.4f}")
                print(f"kd_lim        : {kd_lim:.4f}")
                print(f"kd_lookup     : {kd_lookup:.4f}")
                print(f"x_rel         : {x_rel}")
                print(f"Ke            : {Ke}")
                print(f"Kep           : {Kep}")
                print(f"As_req_cm2    : {As_req_cm2 if As_req_cm2 is not None else 'None'}")
                print(f"Asp_req_cm2   : {Asp_req_cm2 if Asp_req_cm2 is not None else 'None'}")
                print(f"As_min_cm2    : {As_min_cm2:.2f}")
                print(f"As_bal_cm2    : {As_bal_cm2:.2f}")
                print(f"As_max_cm2    : {As_max_cm2:.2f}")
                print(f"As_final_cm2  : {As_final_cm2 if As_final_cm2 is not None else 'None'}")
                print(f"Estado        : {estado}")

            # Actualizar salida principal con Asp_tab_cm2 si corresponde
            if Asp_tab_cm2 is not None:
                salida["Asp_req_cm2"] = round(Asp_tab_cm2, 2)

            salida.update({
                "compresion": {
                    "fs_target_MPa": round(fs_target, 1),
                    "ke_cm2_por_MN": round(ke_tab, 3) if ke_tab is not None else None,
                    "ke_p_cm2_por_MN": round(kep_tab, 3) if kep_tab is not None else None,
                    "As_tab_cm2": round(As_tab_cm2, 2) if As_tab_cm2 is not None else None,
                    "As_prima_tab_cm2": round(Asp_tab_cm2, 2) if Asp_tab_cm2 is not None else None,
                    "fila_match": fila_match
                }
            })

        return salida
    
    def Ec_MPa(self):
        """
        Módulo de elasticidad del hormigón según CIRSOC / ACI
        """
        return 4700.0 * math.sqrt(self.fc)
    
    def propiedades_seccion(self):
        """
        Propiedades geométricas brutas de la sección
        """
        b = self.b
        h = self.h
        rec_inf = self.rec              # recubrimiento inferior en metros
        rec_sup = getattr(self, "rec_sup", self.rec)  # recubrimiento superior (si no está definido, usa el mismo)

        # Distancias útiles
        d_trac = h - rec_inf
        d_comp = rec_sup  # si querés más exacto: rec_sup + Ø/2

        return {
            "b_m": b,
            "h_m": h,
            "rec_inf_m": rec_inf,
            "rec_sup_m": rec_sup,
            "d_m": d_trac,       # útil a tracción
            "d_prime_m": d_comp  # útil a compresión
        }
    
    def momento_fisuracion(self):
        """
        Momento de fisuración del hormigón
        """
        fr_MPa = 0.625 * math.sqrt(self.fc)
        Ig = self.b * self.h**3 / 12.0
        yt = self.h / 2.0

        Mcr_Nm = fr_MPa * 1e6 * Ig / yt
        return Mcr_Nm / 1000.0  # kNm
     
    def inercias_seccion(self, As_cm2, Kc=0.40, M_kNm=None, factor_servicio=1.0, Es_MPa=210000):
        """
        Devuelve:
        Ig  : inercia bruta
        Icr : inercia fisurada (Steiner)
        Mcr : momento de fisuración
        Ie  : inercia efectiva (Branson) si se pasa M
        """
        # --- Inercia bruta ---
        Ig = self.b * self.h**3 / 12.0

        # --- Momento de fisuración ---
        fr_MPa = 0.625 * math.sqrt(self.fc)
        Mcr = fr_MPa * 1e6 * Ig / (self.h / 2.0) / 1000.0  # kNm

        # --- Icr estilo Steiner ---
        As = As_cm2 / 10000.0  # cm² → m²
        Ec = 4700.0 * math.sqrt(self.fc)  # MPa
        n = Es_MPa / Ec  # módulo relativo

        # c (eje neutro) ajustado a servicio
        c = Kc * self.d * factor_servicio
        # Distancias centroides
        y_hormigon = c / 2.0           # centroide bloque hormigón
        y_acero = self.d - c            # centroide acero (aprox desde el eje neutro)

        # Inercia del bloque de hormigón
        I_h = (self.b * c**3) / 12.0 + (self.b * c) * y_hormigon**2

        # Inercia del acero transformado a hormigón
        I_s = n * As * y_acero**2
        Icr = I_h + I_s
        # --- Inercia efectiva ---
        Ie = None
        if M_kNm is not None and M_kNm > Mcr:
            ratio = (Mcr / M_kNm) ** 3
            Ie = ratio * Ig + (1 - ratio) * Icr
        else:
            Ie = Ig
        return Ig, Icr, Mcr, Ie

    def hormigon_tramo(self, L_m, gamma_kg_m3=2500):
        """
        Volumen y peso propio del hormigón
        """
        V_m3 = self.b * self.h * L_m
        peso_kg = V_m3 * gamma_kg_m3

        return {
            "V_m3": V_m3,
            "peso_kg": peso_kg
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
    as_izq = v.calcular_as_necesaria(
        viga_data["m_izq"], coef_kd,
        tabla_comp=tabla_compresion,
        x_rel=(v.recubrimiento_sup_cm / 100.0) / v.d,
        debug_compresion=False
    )
    as_tra = v.calcular_as_necesaria(
        viga_data["m_tra"], coef_kd,
        tabla_comp=tabla_compresion,
        x_rel=(v.recubrimiento_sup_cm / 100.0) / v.d,
        debug_compresion=False
    )

    as_der = v.calcular_as_necesaria(
        viga_data["m_der"], coef_kd,
        tabla_comp=tabla_compresion,
        x_rel=(v.recubrimiento_sup_cm / 100.0) / v.d,
        debug_compresion=False
    )

    # ---- Armaduras adoptadas ----
    arm_sup_izq = v.seleccionar_armadura(as_izq["As_final_cm2"])
    arm_inf_tra = v.seleccionar_armadura(as_tra["As_final_cm2"])
    arm_sup_der = v.seleccionar_armadura(as_der["As_final_cm2"])
    
    # Armadura superior (compresión) -> sólo si doble armadura
    arm_sup_comp = {}
    if as_tra.get("doble_armadura", False):
        arm_sup_comp = v.seleccionar_armadura(as_tra.get("Asp_req_cm2", 0.0))

    
    # ---- Funciones auxiliares ----
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
        # --- Inferior tramo ---
        n, d = arm_inf_tra["n"], arm_inf_tra["diam"]
        L = (PI_der - PI_izq) + 2 * (12 * d / 1000)
        P = L * n * v.barras_comerciales[d]
        lineas.append(f"{pos:<4} | {n:<5} | {d:<4} | Inf. tramo          | {L:<7.2f} | {P:<10.2f}")
        total_peso += P; pos += 1
        # --- Armadura comprimida según tabla (si corresponde) ---
        if as_tra["Asp_req_cm2"] > 0:
            arm_comp = v.seleccionar_armadura(as_tra["Asp_req_cm2"])
            n, d = arm_comp["n"], arm_comp["diam"]
            L = (PI_der - PI_izq)  # misma longitud de tramo
            P = L * n * v.barras_comerciales[d]
            lineas.append(f"{pos:<4} | {n:<5} | {d:<4} | Arm. comprimida    | {L:<7.2f} | {P:<10.2f}")
            total_peso += P
            pos += 1

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
    # =====================================================
    # DATOS GLOBALES PARA NOTAS TÉCNICAS
    # =====================================================

    # =======================================
    # Cálculo de inercias usando Kc y kd
    # =======================================
    

    # 1) Decidir Kc_usado
    if as_tra["fila"].get("Kc") is not None and as_tra["fila"].get("Kc") != 0.0:
        # Caso flexión: usar el Kc de la fila
        Kc_usado = as_tra["fila"]["Kc"]
    else:
        # Caso compresión: no hay Kc en la fila
        try:
            fila_flex = as_tra.get("info_k", {}).get("fila_flexion")
            Kc_usado = fila_flex.get("Kc") if fila_flex else None
        except Exception:
            Kc_usado = None

        # Si sigue siendo None o cero, calcular dinámicamente con triángulos semejantes
        if Kc_usado is None or Kc_usado == 0.0:
            eps_c = 0.003   # 3 por mil
            eps_s = 0.005   # 5 por mil
            Kc_usado = eps_c / (eps_c + eps_s)  # da 0.375

    # 2) Calcular inercias con el Kc elegido
    Ig, Icr, Mcr, Ie = v.inercias_seccion(
        As_cm2 = arm_inf_tra["area_total_cm2"],
        Kc = Kc_usado if Kc_usado is not None else 0.40,  # fallback si querés un default
        M_kNm = as_tra.get("Mu_kNm")
    )

    # 3) Guardar salida
    inercias = {
        "Ig_m4": Ig,
        "Icr_m4": Icr,
        "Ie_m4": Ie,
        "Kc_usado": Kc_usado,
        "kd_lookup": as_tra["kd_lookup"]
    }

    # Hormigón
    V_horm = v.b * v.h * L_viga
    peso_horm = V_horm * 2500  # kg/m³

    horm = {
        "V_m3": V_horm,
        "peso_kg": peso_horm
    }

    # Notas técnicas (UNIFICADO)
    notas_tecnicas = {
        "flexion": {
            "As_req_cm2": as_tra["As_req_cm2"],
            "Asp_req_cm2": as_tra.get("Asp_req_cm2", 0.0),
            "As_min_cm2": as_tra["As_min_cm2"],
            "As_max_cm2": as_tra["As_max_cm2"],
            "As_adopt_cm2": arm_inf_tra["area_total_cm2"],
            "Asp_adopt_cm2": arm_sup_comp.get("area_total_cm2", 0.0),
            "estado": as_tra["estado"],
            "kd": as_tra.get("kd_sel"),
            "fila": as_tra.get("fila"),
            "doble_armadura": as_tra.get("doble_armadura", False),# <- bandera clara
            "Ke": as_tra["fila"].get("Ke"),
            "Kep": as_tra.get("compresion", {}).get("ke_p_cm2_por_MN")
        },
        "materiales": {
            "fc_MPa": v.fc,
            "fy_MPa": v.fy,
            "Ec_MPa": v.Ec_MPa()
        },
        "seccion": {
            "b_cm": v.b * 100,
            "h_cm": v.h * 100,
            "d_cm": v.d * 100,             # útil a tracción
            "rec_inf_cm": v.rec * 100,     # recubrimiento inferior
            "rec_sup_cm": getattr(v, "rec_sup", v.rec) * 100,  # recubrimiento superior si lo definiste
            "d_prime_cm": getattr(v, "d_prime", None)          # útil a compresión si lo calculás
        },
        "inercias": inercias,
        "hormigon": horm
    }  

    # ---- Cierre ----
    lineas.append("")
    lineas.append("NOTAS TÉCNICAS")
    lineas.append("-" * 80)
    lineas.append("1 FLEXIÓN")

    if notas_tecnicas['flexion'].get('doble_armadura', False):
        # Mostrar coeficientes de compresión
        Ke_val  = notas_tecnicas['flexion'].get('Ke')
        Kep_val = notas_tecnicas['flexion'].get('Kep')

        lineas.append(f"   Ke (tracción)      = {Ke_val:.3f} cm²/MN" if Ke_val else "   Ke (tracción)      = —")
        lineas.append(f"   K’e (compresión)   = {Kep_val:.3f} cm²/MN" if Kep_val else "   K’e (compresión)   = —")

        # Mostrar As y As’ requeridas/adoptadas (reales)
        As_req   = notas_tecnicas['flexion']['As_req_cm2']
        Asp_req  = notas_tecnicas['flexion']['Asp_req_cm2']
        As_adopt = notas_tecnicas['flexion'].get('As_adopt_cm2', 0.0)   # área real inferior
        Asp_adopt= notas_tecnicas['flexion'].get('Asp_adopt_cm2', 0.0)  # área real superior

        lineas.append(f"   As requerida       = {As_req:.2f} cm²")
        lineas.append(f"   As adoptada        = {As_adopt:.2f} cm²")
        lineas.append(f"   ΔAs (inf)          = {As_adopt - As_req:+.2f} cm²")

        lineas.append(f"   As’ requerida      = {Asp_req:.2f} cm²")
        lineas.append(f"   As’ adoptada       = {Asp_adopt:.2f} cm²")
        lineas.append(f"   ΔAs’ (sup)         = {Asp_adopt - Asp_req:+.2f} cm²")

        # Estado específico
        lineas.append("   Estado             = ℹ️ Doble armadura activa; menor ductilidad (Kc límite aplicado)")

    else:
        # Caso flexión simple: mantener lógica tradicional
        lineas.append(f"   As requerida      = {notas_tecnicas['flexion']['As_req_cm2']:.2f} cm²")
        lineas.append(f"   As mínima         = {notas_tecnicas['flexion']['As_min_cm2']:.2f} cm²")
        lineas.append(f"   As máxima         = {notas_tecnicas['flexion']['As_max_cm2']:.2f} cm²")
        lineas.append(f"   As adoptada       = {notas_tecnicas['flexion']['As_adopt_cm2']:.2f} cm²")

        if notas_tecnicas['flexion'].get('Asp_req_cm2', 0.0) > 0:
            lineas.append(f"   As’ (tabla K’e)   = {notas_tecnicas['flexion']['Asp_req_cm2']:.2f} cm²")

        exceso = (notas_tecnicas['flexion']['As_adopt_cm2'] +
                notas_tecnicas['flexion'].get('Asp_req_cm2', 0.0) -
                notas_tecnicas['flexion']['As_max_cm2'])

        if exceso > 0:
            lineas.append(f"   ⚠️ Exceso de As      = {exceso:.2f} cm²")
            lineas.append("   Nota: Se pasó de la cuantía máxima; considerar ajuste de sección o peralte")

        lineas.append(f"   Estado            = {notas_tecnicas['flexion']['estado']}")

    # kd y fila usados
    if notas_tecnicas["flexion"]["kd"] is not None:
        lineas.append(
            f"   kd utilizado      = {notas_tecnicas['flexion']['kd']} "
            f"(fila {notas_tecnicas['flexion']['fila']})"
        )
    lineas.append("")

    lineas.append("2 SECCIÓN")
    lineas.append(f"   b = {notas_tecnicas['seccion']['b_cm']:.1f} cm")
    lineas.append(f"   h = {notas_tecnicas['seccion']['h_cm']:.1f} cm")
    lineas.append(f"   d = {notas_tecnicas['seccion']['d_cm']:.1f} cm")
    # Caso flexión simple: solo recubrimiento inferior
    if not notas_tecnicas["flexion"].get("doble_armadura", False):
        lineas.append(f"   recubrimiento = {notas_tecnicas['seccion']['rec_inf_cm']:.1f} cm")
    else:
        # Caso doble armadura: mostrar ambos
        lineas.append(f"   recubrimiento inf = {notas_tecnicas['seccion']['rec_inf_cm']:.1f} cm")
        lineas.append(f"   recubrimiento sup = {notas_tecnicas['seccion']['rec_sup_cm']:.1f} cm")
        if notas_tecnicas["seccion"].get("d_prime_cm") is not None:
            lineas.append(f"   d’ = {notas_tecnicas['seccion']['d_prime_cm']:.1f} cm (compresión)")

    lineas.append("")
    lineas.append("3 MATERIALES")
    fc = notas_tecnicas["materiales"]["fc_MPa"]
    fy = notas_tecnicas["materiales"]["fy_MPa"]
    Ec = notas_tecnicas["materiales"]["Ec_MPa"]
    lineas.append(f"   f'c = {fc} MPa")
    lineas.append(f"   fy  = {fy} MPa")
    lineas.append(f"   Ec  = {Ec:.0f} MPa")

    lineas.append("")
    lineas.append("4 INERCIAS")
    # Conversión a cm⁴
    Ig_cm4  = Ig * 1e8
    Icr_cm4 = Icr * 1e8
    Ie_cm4  = Ie * 1e8

    lineas.append(f"- Inercia bruta Ig: {Ig:.6f} m⁴ | {Ig_cm4:,.0f} cm⁴")
    lineas.append(f"- Inercia fisurada Icr: {Icr:.6f} m⁴ | {Icr_cm4:,.0f} cm⁴")
    lineas.append(f"- Inercia efectiva Ie (Branson): {Ie:.6f} m⁴ | {Ie_cm4:,.0f} cm⁴")
    lineas.append(f"- Momento de fisuración Mcr: {Mcr:.2f} kNm")

    # Añadir trazabilidad de Kc y kd_lookup
    lineas.append(f"- Kc usado para sección: {inercias.get('Kc_usado', 0.0):.3f}")
    lineas.append(f"- kd_lookup: {inercias.get('kd_lookup', 0.0):.3f}")


    lineas.append("")
    lineas.append("5 HORMIGÓN")
    lineas.append(f"   Volumen = {notas_tecnicas['hormigon']['V_m3']:.3f} m³")
    lineas.append(f"   Peso ≈ {notas_tecnicas['hormigon']['peso_kg']:.0f} kg")
    # Puntos de inflexión
    if puntos_inflexion_m and len(puntos_inflexion_m) == 2:
        PI_izq, PI_der = puntos_inflexion_m
    else:
        PI_izq, PI_der = 0.15 * L_viga, 0.85 * L_viga

    lineas.append("")
    lineas.append("6 PUNTOS DE INFLEXIÓN")
    lineas.append(f"   PI izquierdo = {PI_izq:.2f} m")
    lineas.append(f"   PI derecho   = {PI_der:.2f} m")
    # 🔎 Debug en planilla
    
    
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
tabla_flexion = coef_kd["coeficientes_flexion"]
tabla_compresion = coef_kd["tabla_compresion"]["H20_H25_H30_fy420"]


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


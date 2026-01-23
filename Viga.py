import json
import math
from pathlib import Path

RESULTADOS = {
    "tramos": [],
    "flexion": [],
    "estado": [],
    "materiales": [],
    "inercias": [],
    "armaduras": [],
    "estribos": [],
    "flecha": [],
    "corte": []
}

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
        Mn_MNm = Mn_kNm / 1000.0
        if self.b <= 0:
            raise ValueError("El ancho b de la viga no puede ser cero.")
        if Mn_MNm <= 0:
            # Caso voladizo o momento nulo → kd físico no aplica
            return 0.0
        return self.d / math.sqrt(Mn_MNm / self.b)

    def _kd_lim_por_fc(self, fc: float) -> float:
        if fc <= 20:
            return 0.469
        elif fc <= 25:
            return 0.419
        elif fc <= 30:
            return 0.383
        else:
            # fallback: si fc > 30, usar el límite de H30
         return 0.383
   
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
            fila = tabla[f"{kd_bajo:.3f}"]
            return fila, kd_bajo, "exacto", {}

        # 🔑 Aquí el cambio: usar f"{:.3f}" en vez de str()
        datos_bajo = tabla[f"{kd_bajo:.3f}"]
        datos_alto = tabla[f"{kd_alto:.3f}"]

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
        kd_vals = [r["kd_kd"] for r in rows]
        x_vals = [r["x"] for r in rows]

        notas = []

        # Ajuste kd fuera de dominio
        if k_ratio < min(kd_vals):
            k_ratio = min(kd_vals)
            notas.append(f"fuera de dominio kd; se usa límite inferior ({min(kd_vals)})")
        elif k_ratio > max(kd_vals):
            k_ratio = max(kd_vals)
            notas.append(f"fuera de dominio kd; se usa límite superior ({max(kd_vals)})")

        # Ajuste x_rel fuera de dominio
        if x_rel < min(x_vals):
            x_rel = min(x_vals)
            notas.append(f"fuera de dominio x_rel; se usa límite inferior ({min(x_vals)})")
        elif x_rel > max(x_vals):
            x_rel = max(x_vals)
            notas.append(f"fuera de dominio x_rel; se usa límite superior ({max(x_vals)})")

        # Buscar la fila más cercana en kd_ratio y x_rel
        row = min(
            rows,
            key=lambda r: abs(r["kd_kd"] - k_ratio) + abs(r["x"] - x_rel)
        )

        # Localizar columna más cercana en fs_cols
        fs_cols = tabla_comp["fs_cols"]
        col_idx = min(range(len(fs_cols)), key=lambda i: abs(fs_cols[i] - fs_target))

        if col_idx >= len(row["ke"]):
            col_idx = len(row["ke"]) - 1

        ke_val = row["ke"][col_idx]
        kep_val = row["ke_p"][col_idx]

        return ke_val, kep_val, {
            "kd_ratio": row["kd_kd"],
            "x": row["x"],
            "col_idx": col_idx,
            "nota": "; ".join(notas) if notas else "dentro de dominio"
        }
    
    def _as_por_mn(self, Mu_kNm, phi, ke_cm2_por_MN_por_m, kep_cm2_por_MN_por_m):
        Mn_MNm = abs(Mu_kNm) / phi / 1000.0  # kNm -> MN·m
        As_cm2  = ke_cm2_por_MN_por_m  * (Mn_MNm / self.d)
        Asp_cm2 = kep_cm2_por_MN_por_m * (Mn_MNm / self.d)
        return round(As_cm2, 2), round(Asp_cm2, 2)
    
    def calcular_as_necesaria(self, Mu_kNm, coef_kd, tabla_comp=None,
                          dato_externo=None, x_rel=None, k_star=None,
                          fs_politica=None, debug_compresion=False,
                          es_voladizo=False):

        # 0) Parámetros base
        phi = self.phi["flexion_traccion"]
        Mn_kNm = abs(Mu_kNm) / phi
        Mn_MNm = Mn_kNm / 1000.0
        fs_obj = fs_politica or 360  # política/nearest para columna fs

        # 1) kd físico y kd límite (kd*)
        kd_fisico = self._kd_fisico(Mn_kNm)
        kd_lim    = self._kd_lim_por_fc(self.fc)

        # 2) Decidir lookup: simple (kd_físico) vs doble (kd_ratio)
        if kd_fisico >= kd_lim or es_voladizo:
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

        # 4) Obtener coeficientes según criterio kd_lim
        if doble_armadura and tabla_comp is not None and x_rel is not None:
            # Caso compresión activa → SOLO tabla de compresión
            Ke, Kep_tab, info_k = self._lookup_ke_kep(
                k_ratio=kd_lookup,
                x_rel=x_rel,
                tabla_comp=tabla_comp,
                fs_target=fs_obj
            )
            Kep = Kep_tab
            Kc = 0.375   # fijo por deformaciones
            Kz = None
            k_used = kd_lookup
            modo_k = info_k.get("modo", "bilineal")

        else:
            # Caso flexión simple → SOLO tabla de flexión
            fila, k_used, modo_k, info_k = self._fila_kd_interpolada(
                kd_lookup, coef_kd, dato_externo=dato_externo
            )
            Ke = fila.get("Ke")
            Kc = fila.get("Kc", None)
            Kz = fila.get("Kz", None)
            Kep = None   # no aplica compresión
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
            "fila": {"Ke": Ke, "Kc": Kc, "Kz": Kz},
            # 🔎 Nota de auditoría del lookup (kd/x_rel clamp)
            "nota": info_k.get("nota", "")        
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
     
    def inercias_seccion(self, As_cm2, c=None, M_kNm=None, factor_servicio=0.73, Es_MPa=210000):
        """
        Devuelve:
        Ig   : inercia bruta (sección completa)
        Jh   : inercia transformada completa (Steiner con Yg)
        Jhf  : inercia fisurada (Steiner con bloque comprimido c)
        Mcr  : momento de fisuración (usando Jh)
        Ie   : inercia efectiva (Branson, usando Jh si M < Mcr)
        """

        # --- Inercia bruta ---
        Ig = self.b * self.h**3 / 12.0

        # --- Datos de acero y módulo relativo ---
        As = As_cm2 / 10000.0  # cm² → m²
        Ec = 4700.0 * math.sqrt(self.fc)  # MPa
        n = Es_MPa / Ec  # módulo relativo

        # --- Jh: inercia transformada completa (Steiner con Yg) ---
        Ah = self.b * self.h + (n - 1) * As
        y_horm = self.h / 2.0
        y_s = self.d
        Yg = (self.b * self.h * y_horm + (n - 1) * As * y_s) / Ah

        Jh = (self.b * self.h**3) / 12.0 \
            + (self.b * self.h) * (y_horm - Yg)**2 \
            + n * As * (y_s - Yg)**2

        # --- Momento de fisuración (según libro) ---
        fr_MPa = 0.625 * math.sqrt(self.fc)
        Mcr = fr_MPa * 1e6 * Jh / (self.h - Yg) / 1000.0  # kNm

        # --- Jhf: inercia fisurada (bloque comprimido c) ---
        if c is None:
            c = 0.40 * self.d * factor_servicio
        else:
            c = c * factor_servicio

        y_hormigon = c / 2.0
        y_acero = self.d - c

        I_h = (self.b * c**3) / 12.0 + (self.b * c) * y_hormigon**2
        I_s = n * As * y_acero**2
        Jhf = I_h + I_s

        # --- Ie: inercia efectiva (Branson corregida) ---
        if M_kNm is not None:
            M_servicio = M_kNm * factor_servicio
            if M_servicio > Mcr:
                ratio = (Mcr / M_servicio) ** 3
                Ie = ratio * Jh + (1 - ratio) * Jhf
            else:
                Ie = Jh
        else:
            Ie = Jh
        
        return Ig, Jh, Jhf, Mcr, Ie

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

        sep_min = 2.5  # separación mínima normativa en cm (≈25 mm)
        b_util_cm = self.b * 100 - 2 * self.rec
        max_barras_simple = 8

        for diam in opciones:
            area_barra = math.pi * (diam / 10)**2 / 4
            n_barras = max(2, math.ceil(as_cm2 / area_barra))  # mínimo 2 barras

            area_total = n_barras * area_barra
            exceso = area_total - as_cm2

            # chequeo de ancho requerido
            ancho_requerido = n_barras * (diam / 10) + (n_barras - 1) * sep_min
            entra_simple = (ancho_requerido <= b_util_cm) and (n_barras <= max_barras_simple)
            capa = "simple" if entra_simple else "doble"

            # penalización: doble capa o demasiadas barras finas
            penalizacion = 0.0
            if not entra_simple:
                penalizacion += 1e6
            if n_barras > 6 and diam < 12:
                penalizacion += 5000  # penaliza soluciones tipo "7 Ø8"

            exceso_eval = exceso + penalizacion

            # recubrimientos
            rec_geom_mm = self.rec * 1000
            rec_ef_mm = rec_geom_mm + diam / 2.0
            rec_sup_mm = (self.recubrimiento_sup_cm or self.rec * 100) * 10
            rec_sup_ef_mm = rec_sup_mm + diam / 2.0

            # longitud de desarrollo
            Ld_mm = max(12 * diam, self.d * 1000)

            # selección de mejor opción
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
                    "recubrimiento_inf_ef_mm": rec_ef_mm,
                    "recubrimiento_sup_ef_mm": rec_sup_ef_mm,
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

    def calcular_estribos(self, Vu_kN, diam_mm=6, ramas=2, zona="central"):
        # ---------- 1. Conversión ----------
        Vu_N = Vu_kN * 1000               # kN -> N
        b_m = self.b                       # ancho en m
        d_m = self.d                       # altura efectiva en m
        b_mm = b_m * 1000
        d_mm = d_m * 1000
        phi_corte = self.phi["corte"]      # factor phi
        # ---------- 2. Cortante nominal ----------
        Vn_N = Vu_N / phi_corte            # Vn = Vu / φ

        # ---------- 3. Tau sección ----------
        tau_n = Vn_N / (b_m * d_m * 1e6)  # N / m² -> Pa -> MPa
        tau_limite = (5/6) * math.sqrt(self.fc)  # MPa

        if tau_n > tau_limite:
            print(f"[{zona}] ¡Alerta! Tau sección = {tau_n:.2f} MPa > Tau límite = {tau_limite:.2f} MPa")

        # ---------- 4. Cortante que resiste el concreto ----------
        sqrt_fc = min(math.sqrt(self.fc), 8.3)
        Vc_N = 0.1666 * sqrt_fc * b_mm * d_mm  # N

        # ---------- 5. Cortante que deben resistir los estribos ----------
        Vs_N = max(Vn_N - Vc_N, 0.0)  # N

        # ---------- 6. Área de ramas ----------
        area_barra_mm2 = math.pi * (diam_mm**2) / 4
        Asv_mm2 = ramas * area_barra_mm2

        # ---------- 7. Separación teórica ----------
        if Vs_N > 0:
            s_mm_teor = (Asv_mm2 * self.fy * d_mm) / Vs_N
        else:
            s_mm_teor = d_mm / 2  # zona no crítica

        # ---------- 8. Límites normativos ----------
        if zona == "apoyo":
            limite = min(d_mm/2, 150)
        else:
            limite = min(d_mm/2, 200)

        # ---------- 9. Separación final adoptada ----------
        s_mm = min(s_mm_teor, limite)
        s_cm = max(6, min(round((s_mm/10)/2)*2, 30))

        print(
            f"[{zona}] Vu={Vu_kN:.2f} kN → "
            f"Vn={Vn_N/1000:.2f} kN → "
            f"τ_sección={tau_n:.2f} MPa → "
            f"Vc={Vc_N/1000:.2f} kN → "
            f"Vs={Vs_N/1000:.2f} kN → "
            f"s_teor={s_mm_teor:.1f} mm → "
            f"s_final={s_cm} cm"
        )

        return {
            "diam_mm": diam_mm,
            "ramas": ramas,
            "s_cm": int(s_cm),
            "s_mm_teor": s_mm_teor,
            "s_mm_final": s_mm,
            "modo": "resistencia requerida",
            "Vc_kN": Vc_N/1000,
            "Vs_req_kN": Vs_N/1000,
            "Vn_kN": Vn_N/1000,
            "tau_seccion_MPa": tau_n,
            "tau_limite_MPa": tau_limite
    }

    def calcular_estribos_alt(self, Vu_kN, diam_mm=8, ramas=2, zona="central"):
        datos = self.calcular_estribos(Vu_kN, diam_mm=diam_mm, ramas=ramas, zona=zona)
        print(f"[{zona}] opción Ø{diam_mm} → sep={datos['s_cm']}cm")
        return datos

    def _Vc_N(self):
        b_mm = self.b * 1000
        d_mm = self.d * 1000
        sqrt_fc = min(math.sqrt(self.fc), 8.3)
        return 0.17 * sqrt_fc * b_mm * d_mm

    def _xcrit_continuo(self, Vu_izq_kN, Vu_der_kN, L_m, q_kNm):
        Vc_N = self._Vc_N()
        phiVc_kN = self.phi["corte"] * Vc_N / 1000.0

        # Distancia crítica desde cada apoyo
        xcrit_izq = max(0.0, min(L_m, (Vu_izq_kN - phiVc_kN) / q_kNm))
        xcrit_der = max(0.0, min(L_m, (Vu_der_kN - phiVc_kN) / q_kNm))

        # Debug prints
        print("\n--- DEBUG _xcrit_continuo ---")
        print(f"L = {L_m:.2f} m")
        print(f"Vu_izq = {Vu_izq_kN:.2f} kN, Vu_der = {Vu_der_kN:.2f} kN")
        print(f"phiVc = {phiVc_kN:.2f} kN")
        print(f"xcrit_izq = {xcrit_izq:.2f} m desde apoyo izq")
        print(f"xcrit_der = {xcrit_der:.2f} m desde apoyo der")
        print(f"q = {q_kNm:.2f} kN/m")
        print("-----------------------------\n")

        return xcrit_izq, xcrit_der, phiVc_kN

    def _xcrit_voladizo(self, Vu_emp_kN, L_m):
        Vc_N = self._Vc_N()
        phiVc_kN = self.phi["corte"] * Vc_N / 1000.0
        if Vu_emp_kN <= phiVc_kN:
            return 0.0, phiVc_kN
        xcrit = L_m * (1.0 - phiVc_kN / Vu_emp_kN)
        xcrit = max(0.0, min(L_m, xcrit))
        return xcrit, phiVc_kN

    def zonificar_estribos_continuo(self, L_m, Vu_izq_kN, Vu_der_kN, q_kN_m, diam_mm=6, ramas=2):
        xcrit_izq, xcrit_der, phiVc_kN = self._xcrit_continuo(Vu_izq_kN, Vu_der_kN, L_m, q_kN_m)

        # Zonas basadas en los puntos críticos
        zA_ini, zA_fin = 0.0, xcrit_izq
        zC_ini, zC_fin = L_m - xcrit_der, L_m
        zB_ini, zB_fin = zA_fin, zC_ini

        def Vu_x(x):
            return Vu_izq_kN + (Vu_der_kN - Vu_izq_kN) * (x / L_m)

        zonas = []
        if zA_fin > zA_ini:
            zonas.append({
                "zona": "apoyo_izq",
                "desde_m": zA_ini,
                "hasta_m": zA_fin,
                "Vu_kN_ref": Vu_x((zA_ini+zA_fin)/2),
                "estribos": {
                    "Ø6": self.calcular_estribos(
                        Vu_x((zA_ini + zA_fin) / 2),
                        diam_mm=6,
                        ramas=ramas,
                        zona="apoyo"
                    ),
                    "Ø8": self.calcular_estribos(
                    Vu_x((zA_ini + zA_fin) / 2),
                     diam_mm=8,
                    ramas=ramas,
                     zona="apoyo"
                    )
                }
            })
        if zB_fin > zB_ini:
            zonas.append({
                "zona": "central",
                "desde_m": zB_ini,
                "hasta_m": zB_fin,
                "Vu_kN_ref": Vu_x((zB_ini+zB_fin)/2),
                "estribos": {
                    "Ø6": self.calcular_estribos(
                        phiVc_kN,
                        diam_mm=6,
                        ramas=ramas,
                        zona="minimo"
                    ),
                    "Ø8": self.calcular_estribos(
                        phiVc_kN,
                        diam_mm=8,
                         ramas=ramas,
                         zona="minimo"
                    )
                }

            })

        if zC_fin > zC_ini:
            zonas.append({
                "zona": "apoyo_der",
                "desde_m": zC_ini,
                "hasta_m": zC_fin,
                "Vu_kN_ref": Vu_x((zC_ini+zC_fin)/2),
                "estribos": {
                    "Ø6": self.calcular_estribos(
                        Vu_x((zC_ini + zC_fin) / 2),
                        diam_mm=6,
                        ramas=ramas,
                        zona="apoyo"
                    ),
                    "Ø8": self.calcular_estribos(
                        Vu_x((zC_ini + zC_fin) / 2),
                        diam_mm=8,
                        ramas=ramas,
                        zona="apoyo"
                    )
                }

            })

        return {"phiVc_kN": phiVc_kN, "zonas": zonas}

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
            "estribos": self.calcular_estribos(Vu_x(xA), diam_mm, ramas, zona="apoyo")
        })
        # Zona B (resto)
        if zB_fin > zB_ini + 1e-6:
            xB = (zB_ini + zB_fin) / 2.0
            zonas.append({
                "zona": "extremo_libre",
                "desde_m": zB_ini, "hasta_m": zB_fin,
                "Vu_kN_ref": Vu_x(xB),
                "estribos": self.calcular_estribos(Vu_x(xB), diam_mm, ramas, zona="central")
            })
        return {
            "phiVc_kN": phiVc_kN,
            "zonas": zonas
        }

# Función para calcular k y c
# =======================================
def calcular_k_c(b_cm, d_cm, As_cm2, fc_MPa, Es_MPa=210000):
    Ec = 4700.0 * math.sqrt(fc_MPa)
    n = Es_MPa / Ec

    A = (b_cm * d_cm**2) / 2.0
    B = n * As_cm2 * d_cm       # positivo
    C = - n * As_cm2 * d_cm     # negativo

    disc = B**2 - 4*A*C
    print(f"[DEBUG] A={A:.6f}, B={B:.6f}, C={C:.6f}, disc={disc:.6f}")
    if disc < 0:
        print(f"[DEBUG] Discriminante negativo: {disc:.6f}")
        # fallback: usar c aproximado del libro
        k = 0.40
        c = 0.40 * d_cm
        return k, c

    k1 = (-B + math.sqrt(disc)) / (2*A)
    k2 = (-B - math.sqrt(disc)) / (2*A)
    k = k1 if 0 < k1 < 1 else k2
    c = k * d_cm
    print(f"[DEBUG] k={k:.6f}, c={c:.6f}")
    return k, c

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
    puntos_inflexion_m=None,
    cargas=None 
):

    lineas = []
    v = DisenadorViga(nombre_viga, b_cm, h_cm, rec_cm, fc_MPa, fy_MPa)
    tramo_id = f"{nombre_viga}.{id_tramo}"
    # 🔹 Definir las tablas acá, usando coef_kd
    tabla_compresion = coef_kd["tabla_compresion"]["H20_H25_H30_fy420"]
    tabla_flexion = coef_kd["coeficientes_flexion"]


    # -----------------------------------------------------
    # Funciones auxiliares
    # -----------------------------------------------------
    def lon_gancho(d_mm):
        return 15 * d_mm / 1000

    def lon_anclaje(d_mm):
        return 40 * d_mm / 1000

    # -----------------------------------------------------
    # Encabezado
    # -----------------------------------------------------
    tipo = "VOLADIZO" if es_voladizo else "TRAMO INTERIOR"
    lineas.extend([
        "=" * 80,
        f"PLANILLA DE DOBLADO – {nombre_viga} – {tipo} – TRAMO {id_tramo}",
        "=" * 80,
        f"{'Pos':<4} | {'Cant':<5} | {'Ø':<4} | {'Detalle':<21} | {'L (m)':<7} | {'Peso (kg)':<10}",
        "-" * 80
    ])

    pos = 1
    total_peso = 0.0

    # =====================================================
    # VOLADIZO
    # =====================================================
    if es_voladizo:

        Mu = viga_data.get("m_tra", 0.0)

        as_emp = v.calcular_as_necesaria(
            Mu, coef_kd,
            tabla_comp=tabla_compresion,
            x_rel=(v.recubrimiento_sup_cm / 100) / v.d,
            debug_compresion=False,
            es_voladizo=True
        )

        arm_sup = v.seleccionar_armadura(as_emp["As_final_cm2"])

        # Superior principal
        n, d = arm_sup["n"], arm_sup["diam"]
        L = L_viga + lon_anclaje(d)
        P = L * n * v.barras_comerciales[d]
        lineas.append(f"{pos:<4} | {n:<5} | {d:<4} | Sup. voladizo         | {L:<7.2f} | {P:<10.2f}")
        RESULTADOS["armaduras"].append({
            "tramo": tramo_id,
            "pos": pos,
            "cantidad": n,
            "diam_mm": d,
            "detalle": "Sup. voladizo",
            "longitud_m": L,
            "peso_kg": P
        })

        total_peso += P
        pos += 1

        # Inferior mínima
        n, d = 2, 10
        L = L_viga
        P = L * n * v.barras_comerciales[d]
        lineas.append(f"{pos:<4} | {n:<5} | {d:<4} | Inf. mínima           | {L:<7.2f} | {P:<10.2f}")
        RESULTADOS["armaduras"].append({
            "tramo": tramo_id,
            "pos": pos,
            "cantidad": n,
            "diam_mm": d,
            "detalle": "Inf. mínima",
            "longitud_m": L,
            "peso_kg": P
        })
        total_peso += P
        pos += 1

        # Estribos
        zonas = v.zonificar_estribos_voladizo(L_viga, viga_data["Vu_emp"])

        for z in zonas["zonas"]:
            e = z["estribos"]
            n_est = math.ceil((z["hasta_m"] - z["desde_m"]) / (e["s_cm"] / 100))
            L_est = 2 * (v.b - 0.06) + 2 * (v.h - 0.06) + 0.10
            P = n_est * L_est * v.barras_comerciales[e["diam_mm"]]

            lineas.append(
                f"{pos:<4} | {n_est:<5} | {e['diam_mm']:<4} | Estribo {z['zona']:<12} "
                f"| {L_est:<7.2f} | {P:<10.2f} | c/{e['s_cm']}cm"
            )
            RESULTADOS["estribos"].append({
                "tramo": tramo_id,
                "zona": z["zona"],
                "diam_mm": e["diam_mm"],
                "separacion_cm": e["s_cm"],
                "cantidad": n_est,
                "longitud_m": L_est,
                "peso_kg": P
            })
            total_peso += P
            pos += 1

        lineas.append("-" * 80)
        lineas.append(f"{'TOTAL PESO':<43} | {total_peso:<10.2f} kg")

        # -----------------------------------------------------
        # Inercias para voladizo (usando misma lógica que tramo)
        # -----------------------------------------------------
        if es_voladizo:
            # --- Calcular k y c con la función cuadrática ---
            k_emp, c_emp = calcular_k_c(
                b_cm=b_cm,
                d_cm=v.d*100,                   # altura útil en cm
                As_cm2=arm_sup["area_total_cm2"],
                fc_MPa=fc_MPa
            )

            # --- Calcular inercias solo del empotramiento ---
            Ig_emp, Jh_emp, Jhf_emp, Mcr_emp, Ie_emp = v.inercias_seccion(
                As_cm2=arm_sup["area_total_cm2"],
                c=c_emp/ 100.0,  # ✅ convertir cm → m
                M_kNm=Mu         # momento máximo en el empotramiento
            )

            # --- Guardar resultados en diccionario ---
            viga_data["inercias"] = {
                "empotramiento": {
                    "Ig": Ig_emp,
                    "Jh": Jh_emp,
                    "Jhf": Jhf_emp,
                    "Icr": Jhf_emp,
                    "Ie": Ie_emp,
                    "Mcr": Mcr_emp,
                    "k": k_emp,
                    "c": c_emp
                },
                "Ief": Ie_emp  # para la flecha
            }

            # --- Flecha de servicio ---
            q_serv = cargas["servicio"]             # kN/m
            q_serv = q_serv * 1e3 / 1000           # convertir a N/mm si lo tenés así
            L = L_viga * 1000                       # mm
            Ec = 4700 * math.sqrt(fc_MPa)           # N/mm²

            delta = q_serv * L**4 / (8 * Ec * (Ie_emp * 1e12))  # mm

        # Volumen y peso de hormigón (válido para cualquier tramo)
        A_cm2 = b_cm * h_cm
        V_cm3 = A_cm2 * (L_viga * 100)   # L en cm
        V_m3 = V_cm3 / 1e6
        peso_hormigon = V_m3 * 2400  # kg

        lineas.append("")
        lineas.append("SECCIÓN Y VOLUMEN (Voladizo)")
        lineas.append(f"   b = {b_cm} cm | h = {h_cm} cm | Recubrimiento = {rec_cm} cm")
        lineas.append(f"   Área sección = {A_cm2:.0f} cm²")
        lineas.append(f"   Longitud L = {L_viga:.2f} m")
        lineas.append(f"   Volumen hormigón = {V_m3:.3f} m³")
        lineas.append(f"   Peso aprox. = {peso_hormigon:.0f} kg")

        lineas.append("")
        lineas.append("CUANTÍA DE ARMADURA (Voladizo)")
        lineas.append(f"   As requerido = {as_emp['As_req_cm2']:.2f} cm²")
        lineas.append(f"   As mínimo    = {as_emp['As_min_cm2']:.2f} cm²")
        lineas.append(f"   As balance   = {as_emp['As_bal_cm2']:.2f} cm²")
        lineas.append(f"   As máximo    = {as_emp['As_max_cm2']:.2f} cm²")
        lineas.append(f"   As adoptado  = {as_emp['As_final_cm2']:.2f} cm²")

        lineas.append("")
        lineas.append("FLECHA DE SERVICIO (Voladizo)")
        lineas.append(f"   δmax = {delta:.2f} mm")
        lineas.append(f"   Flecha adm L/180 (Voladizos) = {L/180:.2f} mm | q_serv = {cargas['servicio']:.2f} kN/m → {'✅ Cumple' if delta <= L/180 else '❌ No cumple'}")

        lineas.append("")
        lineas.append("INERCIAS Y MOMENTO DE FISURACIÓN")
        lineas.append(f"   Ig  (bruta)       = {Ig_emp:.3e} m⁴")
        lineas.append(f"   Icr (fisurada)    = {Jhf_emp:.3e} m⁴")
        lineas.append(f"   Ie  (efectiva)    = {Ie_emp:.3e} m⁴")
        lineas.append(f"   Mcr (fisuración)  = {Mcr_emp:.2f} kNm")
        lineas.append(f"   Relación Ie/Ig    = {Ie_emp/Ig_emp:.3f}")
        lineas.append(f"   k usado           = {k_emp:.3f} | c usado = {c_emp:.2f} cm")
        lineas.append(f"   Materiales: fc = {fc_MPa} MPa | fy = {fy_MPa} MPa | Ec = {Ec:.0f} MPa")  

        # Guardado homogéneo en RESULTADOS
        tramo_id = f"{nombre_viga}.{id_tramo}"

        RESULTADOS["tramos"].append({
            "id": tramo_id,
            "viga": nombre_viga,
            "tipo": "voladizo",
            "L_m": L_viga
        })
        RESULTADOS["flexion"].append({
            "tramo": tramo_id,
            "As_req_cm2": as_emp["As_req_cm2"],
            "As_adop_cm2": as_emp["As_final_cm2"],
            "cumple": as_emp["As_final_cm2"] >= as_emp["As_req_cm2"]
        })

        RESULTADOS["flecha"].append({
            "tramo": tramo_id,
            "tipo_tramo": "voladizo",
            "L_m": L_viga,
            "q_serv_kN_m": cargas["servicio"],
            "delta_mm": delta,
            "lim_L180_mm": L/180,
            "cumple_L180": delta <= L/180,
            "Ie_m4": Ie_emp
        })

        RESULTADOS["materiales"].append({
            "tramo": tramo_id,
            "fc_MPa": fc_MPa,
            "fy_MPa": fy_MPa,
            "Ec_MPa": Ec,
            "V_m3": V_m3,
            "peso_kg": peso_hormigon
        })
        
        # ✅ Guardar corte aquí solo una vez
        RESULTADOS["corte"].append({
            "tramo": tramo_id,
            "V_kN": viga_data["Vu_emp"],
            "cumple": True
        })

        return "\n".join(lineas)

    # =====================================================
    # TRAMO INTERIOR
    # =====================================================

    # Puntos de inflexión
    if puntos_inflexion_m and len(puntos_inflexion_m) == 2:
        PI_izq, PI_der = puntos_inflexion_m
    else:
        PI_izq, PI_der = 0.15 * L_viga, 0.85 * L_viga
    RESULTADOS["tramos"].append({
        "id": tramo_id,
        "viga": nombre_viga,
        "tipo": "voladizo" if es_voladizo else "interior",
        "L_m": L_viga
    })
    RESULTADOS["tramos"][-1].update({
        "PI_izq_m": PI_izq,
        "PI_der_m": PI_der
    })

    # Materiales y hormigón
    A_cm2 = b_cm * h_cm
    V_cm3 = A_cm2 * (L_viga * 100)   # L en cm
    V_m3 = V_cm3 / 1e6
    peso_hormigon = V_m3 * 2400

    RESULTADOS["materiales"].append({
        "tramo": tramo_id,
        "fc_MPa": v.fc,
        "fy_MPa": v.fy,
        "Ec_MPa": v.Ec_MPa(),
        "V_m3": V_m3,
        "peso_kg": peso_hormigon
    })

    # -----------------------------------------------------
    # CÁLCULO REAL DE ARMADURAS
    # -----------------------------------------------------

    ## 1) Campo (armadura inferior)
    as_tra = v.calcular_as_necesaria(
        Mu_kNm=viga_data["m_tra"],
        coef_kd=coef_kd,
        tabla_comp=tabla_compresion
    )
    arm_inf_tra = v.seleccionar_armadura(as_tra["As_final_cm2"])

    # 2) Apoyo izquierdo (armadura superior izq)
    as_izq = v.calcular_as_necesaria(
        Mu_kNm=viga_data["m_izq"],
        coef_kd=coef_kd,
        tabla_comp=tabla_compresion
    )
    arm_sup_izq = v.seleccionar_armadura(as_izq["As_final_cm2"])

    # 3) Apoyo derecho (armadura superior der)
    as_der = v.calcular_as_necesaria(
        Mu_kNm=viga_data["m_der"],
        coef_kd=coef_kd,
        tabla_comp=tabla_compresion
    )
    arm_sup_der = v.seleccionar_armadura(as_der["As_final_cm2"])

    # 4) Si el tramo pide doble armadura, también compresión
    arm_sup_comp = {}
    if as_tra["doble_armadura"] and as_tra["Asp_req_cm2"] > 0:
        arm_sup_comp = v.seleccionar_armadura(as_tra["Asp_req_cm2"])

    nota_tra = "Verificar armado según normativa"

    # -----------------------------------------------------
    # Superiores e inferiores
    # -----------------------------------------------------
    # Superior apoyo izquierdo
    n, d = arm_sup_izq["n"], arm_sup_izq["diam"]
    L = PI_izq + 40 * d / 1000  # lon_anclaje
    P = L * n * v.barras_comerciales[d]
    lineas.append(f"{pos:<4} | {n:<5} | {d:<4} | Sup. apoyo izq        | {L:<7.2f} | {P:<10.2f}")
    RESULTADOS["armaduras"].append({
        "tramo": tramo_id,
        "pos": pos,
        "cantidad": n,
        "diam_mm": d,
        "detalle": "Sup. apoyo izq",
        "longitud_m": L,
        "peso_kg": P
    })
    total_peso += P
    pos += 1

    # Armadura inferior tramo (solo vano)
    n_total, d = arm_inf_tra["n"], arm_inf_tra["diam"]
    n_prol = max(2, n_total//2)   # mínimo 2 prolongadas
    n_tramo = n_total - n_prol    # las que quedan solo en vano

    if n_tramo > 0:
        L_tramo = (PI_der - PI_izq) + 2 * (12 * d / 1000)
        P_tramo = L_tramo * n_tramo * v.barras_comerciales[d]
        lineas.append(f"{pos:<4} | {n_tramo:<5} | {d:<4} | Inf. tramo            | {L_tramo:<7.2f} | {P_tramo:<10.2f}")
        RESULTADOS["armaduras"].append({
            "tramo": tramo_id,
            "pos": pos,
            "cantidad": n_tramo,
            "diam_mm": d,
            "detalle": "Inf. tramo",
            "longitud_m": L_tramo,
            "peso_kg": P_tramo
        })
        total_peso += P_tramo
        pos += 1

    # Armadura inferior tramo completo (prolongada a apoyos)
    L_completo = L_viga + 2 * (12 * d / 1000)  # recorre todo el vano + ganchos
    P_completo = L_completo * n_prol * v.barras_comerciales[d]
    lineas.append(f"{pos:<4} | {n_prol:<5} | {d:<4} | Inf. tramo completo   | {L_completo:<7.2f} | {P_completo:<10.2f}")
    RESULTADOS["armaduras"].append({
        "tramo": tramo_id,
        "pos": pos,
        "cantidad": n_prol,
        "diam_mm": d,
        "detalle": "Inf. tramo completo",
        "longitud_m": L_completo,
        "peso_kg": P_completo
    })
    total_peso += P_completo
    pos += 1



    # Armadura comprimida según tabla (si corresponde)
    if as_tra["Asp_req_cm2"] > 0:
        arm_comp = v.seleccionar_armadura(as_tra["Asp_req_cm2"])
        n, d = arm_comp["n"], arm_comp["diam"]
        L = (PI_der - PI_izq)
        P = L * n * v.barras_comerciales[d]
        lineas.append(f"{pos:<4} | {n:<5} | {d:<4} | Arm. comprimida   | {L:<7.2f} | {P:<10.2f}")
        RESULTADOS["armaduras"].append({
            "tramo": tramo_id,
            "pos": pos,
            "cantidad": n,
            "diam_mm": d,
            "detalle": "Arm. comprimida",
            "longitud_m": L,
            "peso_kg": P
        })
        total_peso += P
        pos += 1

    # Superior apoyo derecho
    n, d = arm_sup_der["n"], arm_sup_der["diam"]
    L = (L_viga - PI_der) + 40 * d / 1000
    P = L * n * v.barras_comerciales[d]
    lineas.append(f"{pos:<4} | {n:<5} | {d:<4} | Sup. apoyo der        | {L:<7.2f} | {P:<10.2f}")
    RESULTADOS["armaduras"].append({
        "tramo": tramo_id,
        "pos": pos,
        "cantidad": n,
        "diam_mm": d,
        "detalle": "Sup. apoyo der",
        "longitud_m": L,
        "peso_kg": P
    })
    total_peso += P
    pos += 1

    # Auxiliar superior
    n, d = 2, 8
    L = L_viga + 2 * 15 * d / 1000  # lon_gancho
    P = L * n * v.barras_comerciales[d]
    lineas.append(f"{pos:<4} | {n:<5} | {d:<4} | Auxiliar superior     | {L:<7.2f} | {P:<10.2f}")
    RESULTADOS["armaduras"].append({
        "tramo": tramo_id,
        "pos": pos,
        "cantidad": n,
        "diam_mm": d,
        "detalle": "Auxiliar superior",
        "longitud_m": L,
        "peso_kg": P
    })
    total_peso += P
    pos += 1
    # Estribos
    zonas = v.zonificar_estribos_continuo(
        L_viga,
        viga_data["Vu_izq"],
        viga_data["Vu_der"],
        cargas["1.2D+1.6L"]
    )

    total_peso_6 = 0
    total_peso_8 = 0

    for z in zonas["zonas"]:
        e = z["estribos"]["Ø6"]   # solo la opción normativa Ø6
        n_est = math.ceil((z["hasta_m"] - z["desde_m"]) / (e["s_cm"] / 100))
        L_est = 2 * (v.b - 0.06) + 2 * (v.h - 0.06) + 0.10
        P = n_est * L_est * v.barras_comerciales[e["diam_mm"]]

        lineas.append(
            f"{pos:<4} | {n_est:<5} | {e['diam_mm']:<4} | Estribo {z['zona']:<13} "
            f"| {L_est:<7.2f} | {P:<10.2f} | c/{e['s_cm']}cm (Ø6)"
        )
        RESULTADOS["estribos"].append({
            "tramo": tramo_id,
            "zona": z["zona"],
            "diam_mm": e["diam_mm"],
            "separacion_cm": e["s_cm"],
            "cantidad": n_est,
            "longitud_m": L_est,
            "peso_kg": P
        })
        
        total_peso += P
        pos += 1
    
    # después de terminar de agregar todas las barras
    lineas.append("--------------------------------------------------------------------------------")
    lineas.append(f"TOTAL PESO                                            | {total_peso:.2f} kg")

    #=============================
    # Recopilación de As para inercias

    # Tramo central
    # Armadura inferior total (ya la tenés en arm_inf_tra)
    As_total_inf = arm_inf_tra["area_total_cm2"]

    # Armadura inferior tramo completo (prolongadas)
    arm_inf_completo = {
        "n": n_prol,
        "diam": d,
        "area_total_cm2": (n_prol * (math.pi * d**2 / 4)) / 100
    }
    # Armadura inferior tramo (solo vano)
    As_tramo = As_total_inf - arm_inf_completo["area_total_cm2"]

    # Ahora sí podés usar:
    As_tramo_completo = arm_inf_completo["area_total_cm2"]

    As_aux_comp = arm_sup_comp.get("area_total_cm2", 0.0)
    arm_comp = {}
    if as_tra["Asp_req_cm2"] > 0:
        arm_comp = v.seleccionar_armadura(as_tra["Asp_req_cm2"])

    As_comp = arm_comp.get("area_total_cm2", 0.0)

    As_tramo_total = As_tramo + As_tramo_completo + As_aux_comp + As_comp

    # Apoyo izquierdo
    As_apoyo_izq = arm_sup_izq["area_total_cm2"] + As_tramo_completo

    # Apoyo derecho
    As_apoyo_der = arm_sup_der["area_total_cm2"] + As_tramo_completo
       
    # --- Calculo de inercias ---
    # Tramo central
    k_tramo, c_tramo = calcular_k_c(b_cm, v.d*100, As_tramo_total, fc_MPa)
    Ig_t, Jh_t, Jhf_t, Mcr_t, Ie_t = v.inercias_seccion(
     As_cm2=As_tramo_total,
        c=c_tramo / 100.0,   # ✅ cm → m
        M_kNm=viga_data["m_tra"]
    )

    # Apoyo izquierdo
    k_izq, c_izq = calcular_k_c(b_cm, v.d*100, As_apoyo_izq, fc_MPa)
    Ig_i, Jh_i, Jhf_i, Mcr_i, Ie_i = v.inercias_seccion(
        As_cm2=As_apoyo_izq,
        c=c_tramo / 100.0,   # ✅ cm → m
        M_kNm=viga_data["m_izq"]
    )

    # Apoyo derecho
    k_der, c_der = calcular_k_c(b_cm, v.d*100, As_apoyo_der, fc_MPa)
    Ig_d, Jh_d, Jhf_d, Mcr_d, Ie_d = v.inercias_seccion(
        As_cm2=As_apoyo_der,
        c=c_tramo / 100.0,   # ✅ cm → m
        M_kNm=viga_data["m_der"]
    )   

    # Promedio de apoyos y efectiva
    Ie_apoyo = (Ie_i + Ie_d)
    Ief = 0.50 * Ie_t + 0.25 * Ie_apoyo

    # --- Guardado ---
    inercias = {
        "tramo": {
            "Ig": Ig_t, "Jh": Jh_t, "Jhf": Jhf_t,
            "Icr": Jhf_t, "Ie": Ie_t, "Mcr": Mcr_t,
            "k": k_tramo, "c": c_tramo
        },
        "apoyo_izq": {
            "Ig": Ig_i, "Jh": Jh_i, "Jhf": Jhf_i,
            "Icr": Jhf_i, "Ie": Ie_i, "Mcr": Mcr_i,
            "k": k_izq, "c": c_izq
        },
        "apoyo_der": {
            "Ig": Ig_d, "Jh": Jh_d, "Jhf": Jhf_d,
            "Icr": Jhf_d, "Ie": Ie_d, "Mcr": Mcr_d,
            "k": k_der, "c": c_der
        },
        "Ie_apoyo_prom": Ie_apoyo,
        "Ief": Ief
    }
    viga_data["inercias"] = inercias


    # Hormigón
    V_horm = v.b * v.h * L_viga
    peso_horm = V_horm * 2500  # kg/m³
    horm = {"V_m3": V_horm, "peso_kg": peso_horm}

    # Notas técnicas
    notas_tecnicas = {
        "flexion": {
            "As_req_cm2": as_tra["As_req_cm2"],
            "Asp_req_cm2": as_tra.get("Asp_req_cm2", 0.0),
            "As_min_cm2": as_tra["As_min_cm2"],
            "As_max_cm2": as_tra["As_max_cm2"],
            "As_adopt_cm2": arm_inf_tra["area_total_cm2"],
            "Asp_adopt_cm2": arm_sup_comp["area_total_cm2"] if arm_sup_comp else 0.0,
            "estado": as_tra["estado"],
            "nota": nota_tra,
            "kd": as_tra.get("kd_sel"),
            "fila": as_tra.get("fila"),
            "doble_armadura": as_tra.get("doble_armadura", False),
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
            "d_cm": v.d * 100,
            "rec_inf_cm": v.rec * 100,
            "rec_sup_cm": getattr(v, "rec_sup", v.rec) * 100,
            "d_prime_cm": getattr(v, "d_prime", None)
        },
        "inercias": inercias,
        "hormigon": horm
    }

    lineas.append("")
    lineas.append("NOTAS TÉCNICAS")
    lineas.append("-" * 80)
    lineas.append("1 FLEXIÓN")

    if notas_tecnicas['flexion'].get('doble_armadura', False):
        Ke_val  = notas_tecnicas['flexion'].get('Ke')
        Kep_val = notas_tecnicas['flexion'].get('Kep')

        lineas.append(f"   Ke (tracción)      = {Ke_val:.3f} cm²/MN" if Ke_val else "   Ke (tracción)      = —")
        lineas.append(f"   K’e (compresión)   = {Kep_val:.3f} cm²/MN" if Kep_val else "   K’e (compresión)   = —")

        As_req    = notas_tecnicas['flexion']['As_req_cm2']
        Asp_req   = notas_tecnicas['flexion']['Asp_req_cm2']
        As_adopt  = notas_tecnicas['flexion'].get('As_adopt_cm2', 0.0)
        Asp_adopt = notas_tecnicas['flexion'].get('Asp_adopt_cm2', 0.0)

        lineas.append(f"   As requerida       = {As_req:.2f} cm²")
        lineas.append(f"   As adoptada        = {As_adopt:.2f} cm²")
        lineas.append(f"   ΔAs (inf)          = {As_adopt - As_req:+.2f} cm²")

        lineas.append(f"   As’ requerida      = {Asp_req:.2f} cm²")
        lineas.append(f"   As’ adoptada       = {Asp_adopt:.2f} cm²")
        lineas.append(f"   ΔAs’ (sup)         = {Asp_adopt - Asp_req:+.2f} cm²")

        estado_flex = notas_tecnicas['flexion'].get('estado', "")
        lineas.append(f"   Estado             = {estado_flex}")

        nota_flex = notas_tecnicas['flexion'].get('nota')
        if nota_flex:
            lineas.append(f"   Nota               = {nota_flex}")

        lineas.append(f"   Luz libre = {L_viga:.2f} m")

    else:
        # Caso flexión simple: mantener lógica tradicional
        lineas.append(f"   As requerida       = {notas_tecnicas['flexion']['As_req_cm2']:.2f} cm²")
        lineas.append(f"   As mínima          = {notas_tecnicas['flexion']['As_min_cm2']:.2f} cm²")
        lineas.append(f"   As máxima          = {notas_tecnicas['flexion']['As_max_cm2']:.2f} cm²")
        lineas.append(f"   As adoptada        = {notas_tecnicas['flexion']['As_adopt_cm2']:.2f} cm²")

        if notas_tecnicas['flexion'].get('Asp_req_cm2', 0.0) > 0:
            lineas.append(f"   As’ (tabla K’e)    = {notas_tecnicas['flexion']['Asp_req_cm2']:.2f} cm²")

        exceso = (notas_tecnicas['flexion']['As_adopt_cm2'] +
                notas_tecnicas['flexion'].get('Asp_req_cm2', 0.0) -
                notas_tecnicas['flexion']['As_max_cm2'])

        if exceso > 0:
            lineas.append(f"   ⚠️ Exceso de As      = {exceso:.2f} cm²")
            lineas.append("   Nota: Se pasó de la cuantía máxima; considerar ajuste de sección o peralte")

        lineas.append(f"   Estado             = {notas_tecnicas['flexion']['estado']}")
        lineas.append(f"   Luz libre = {L_viga:.2f} m")

    # kd y fila usados
    if notas_tecnicas["flexion"]["kd"] is not None:
        lineas.append(
            f"   kd utilizado      = {notas_tecnicas['flexion']['kd']} "
            f"(fila {notas_tecnicas['flexion']['fila']})"
        )

    RESULTADOS["flexion"].append({
        "viga": nombre_viga,
        "tramo": id_tramo,
        "As_req_cm2": notas_tecnicas['flexion']['As_req_cm2'],
        "As_adop_cm2": notas_tecnicas['flexion'].get('As_adopt_cm2', 0.0),
        "cumple": notas_tecnicas['flexion'].get('As_adopt_cm2', 0.0) >= notas_tecnicas['flexion']['As_req_cm2']
    })
    RESULTADOS["estado"].append({
        "tramo": tramo_id,
        "flexion": notas_tecnicas["flexion"]["estado"],
        "nota": notas_tecnicas["flexion"]["nota"]
    })


    lineas.append("")
    lineas.append("2 SECCIÓN")
    lineas.append(f"   b = {notas_tecnicas['seccion']['b_cm']:.1f} cm")
    lineas.append(f"   h = {notas_tecnicas['seccion']['h_cm']:.1f} cm")
    lineas.append(f"   d = {notas_tecnicas['seccion']['d_cm']:.1f} cm")

    if not notas_tecnicas["flexion"].get("doble_armadura", False):
        lineas.append(f"   recubrimiento = {notas_tecnicas['seccion']['rec_inf_cm']:.1f} cm")
    else:
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
    lineas.append("4 INERCIAS (Rigidez global de la viga)")

    lineas.append(f"- Inercia bruta Ig              = {Ig_t:.6f} m⁴")
    lineas.append(f"- Inercia efectiva global Ief   = {Ief:.6f} m⁴")
    ratio = Ief / Ig_t
    if ratio > 0.95:
        nota = "sección poco fisurada"
    elif ratio > 0.85:
        nota = "buen comportamiento en servicio"
    elif ratio > 0.70:
        nota = "fisuración normal (esperable)"
    else:
        nota = "fisuración elevada – revisar flechas"

    lineas.append(f"- Relación Ief / Ig             = {ratio:.3f}  ({nota})")
    lineas.append(f"- Momento de fisuración Mcr     = {Mcr_t:.2f} kNm (tramo)")
    lineas.append(f"- c usado                       = {c_tramo:.2f} cm")  # Mostrar c usado


    RESULTADOS["inercias"].append({
        "tramo": tramo_id,
        "Ig_m4": Ig_t,
        "Ief_m4": Ief,
        "ratio_Ief_Ig": Ief / Ig_t,
        "Mcr_tramo_kNm": Mcr_t,
        "c_usado_cm": c_tramo,
        "criterio_Ief": "0.50·tramo + 0.25·apoyo_izq + 0.25·apoyo_der"
    })



    lineas.append("")
    lineas.append("5 HORMIGÓN")
    lineas.append(f"   Volumen = {notas_tecnicas['hormigon']['V_m3']:.3f} m³")
    lineas.append(f"   Peso ≈ {notas_tecnicas['hormigon']['peso_kg']:.0f} kg")

    if puntos_inflexion_m and len(puntos_inflexion_m) == 2:
        PI_izq, PI_der = puntos_inflexion_m
    else:
        PI_izq, PI_der = 0.15 * L_viga, 0.85 * L_viga

    lineas.append("")
    lineas.append("6 PUNTOS DE INFLEXIÓN")
    lineas.append(f"   PI izquierdo = {PI_izq:.2f} m")
    lineas.append(f"   PI derecho   = {PI_der:.2f} m")

    lineas.append("")
    lineas.append("7 ARMADURAS (As calculadas)")
    lineas.append(f"   As apoyo izq = {arm_sup_izq['area_total_cm2']:.2f} cm²")
    lineas.append(f"   As tramo     = {arm_inf_tra['area_total_cm2']:.2f} cm²")
    lineas.append(f"   As apoyo der = {arm_sup_der['area_total_cm2']:.2f} cm²")

    if as_tra.get("Asp_req_cm2", 0.0) > 0:
        lineas.append(f"   As’ (compresión) = {as_tra['Asp_req_cm2']:.2f} cm²")

    # Estribos: totales y zonificación

    total_estribos_6 = sum(
        math.ceil((z["hasta_m"] - z["desde_m"]) / (z["estribos"]["Ø6"]["s_cm"] / 100))
        for z in zonas["zonas"]
    )

    total_estribos_8 = sum(
        math.ceil((z["hasta_m"] - z["desde_m"]) / (z["estribos"]["Ø8"]["s_cm"] / 100))
        for z in zonas["zonas"]
    )

    # tomás los valores de la primera zona como referencia
    e6 = zonas["zonas"][0]["estribos"]["Ø6"] if zonas["zonas"] else {"s_mm_teor":0,"s_cm":0}
    e8 = zonas["zonas"][0]["estribos"]["Ø8"] if zonas["zonas"] else {"s_mm_teor":0,"s_cm":0}

    lineas.append("")
    lineas.append("8 ESTRIBOS")
    lineas.append(
        f"   Cantidad total Ø6 = {total_estribos_6} unidades "
        f"(s_teor≈{e6['s_mm_teor']:.1f} mm → s_final={e6['s_cm']} cm)"
    )
    lineas.append(
        f"   Cantidad total Ø8 = {total_estribos_8} unidades "
        f"(s_teor≈{e8['s_mm_teor']:.1f} mm → s_final={e8['s_cm']} cm)"
    )

    zonas_txt = " | ".join(
        f"{z['zona']}={z['desde_m']:.2f}-{z['hasta_m']:.2f} m"
        for z in zonas["zonas"]
    )
    lineas.append(f"   Zonificación   = {zonas_txt}")
    # Guardar estribos en RESULTADOS
    for z in zonas["zonas"]:
        for diam in ["Ø6", "Ø8"]:
            e = z["estribos"].get(diam)
            if e:
                n_est = math.ceil((z["hasta_m"] - z["desde_m"]) / (e["s_cm"] / 100))

    # =============================================================================
    # 9 FLECHA DE SERVICIO – MÉTODO DE CARGA VIRTUAL
    # =============================================================================

    # Inercias ya calculadas para ESTA viga
    inercias = viga_data.get("inercias", None)
    if inercias is None:
        raise ValueError("Inercias no calculadas para este tramo")
    
    # Inercias ya calculadas para ESTA viga
    Ec = 4700 * math.sqrt(fc_MPa)  # N/mm²
    L = L_viga * 1000              # mm
    q = cargas["servicio"]         # N/mm

    Ie_t = inercias["tramo"]["Ie"] * 1e12
    Ie_i = inercias["apoyo_izq"]["Ie"] * 1e12
    Ie_d = inercias["apoyo_der"]["Ie"] * 1e12

    # Coeficientes del método de carga virtual
    alpha_tramo = 5 / 12
    alpha_apoyo = 1 / 4

    inv_Ie_eq = (
        alpha_tramo / Ie_t +
        alpha_apoyo / Ie_i +
        alpha_apoyo / Ie_d
    )
    # --- Flecha máxima de servicio ---
    coef_global = 5 / 384

    delta = coef_global * q * L**4 / Ec * inv_Ie_eq   # mm


    lineas.append("")
    lineas.append("9 FLECHA DE SERVICIO")
    lineas.append("   Método: Carga virtual con rigidez variable (tramo + apoyos)")
    lineas.append(f"   δmax = {delta:.2f} mm")

    lineas.append(
        f"   Flecha adm L/250 = {L/250:.2f} mm | q_serv = {cargas['servicio']:.2f} kN/m → "
        f"{'✅ Cumple' if delta <= L/250 else '❌ No cumple'}"
    )
    lineas.append(
        f"   Flecha adm L/360 = {L/360:.2f} mm | q_serv = {cargas['servicio']:.2f} kN/m → "
        f"{'✅ Cumple' if delta <= L/360 else '❌ No cumple'}"
    )
    lineas.append(
        f"   Flecha adm L/480 = {L/480:.2f} mm | q_serv = {cargas['servicio']:.2f} kN/m → "
        f"{'✅ Cumple' if delta <= L/480 else '❌ No cumple'}"
    )
    lineas.append("")
    lineas.append("   Inercias usadas (tramo + apoyos)")
    lineas.append(f"   Ie tramo       = {inercias['tramo']['Ie']:.3e} m⁴")
    lineas.append(f"   Ie apoyo izq   = {inercias['apoyo_izq']['Ie']:.3e} m⁴")
    lineas.append(f"   Ie apoyo der   = {inercias['apoyo_der']['Ie']:.3e} m⁴")
    lineas.append(f"   c usado        = {c_tramo:.2f} cm")
    RESULTADOS["flecha"].append({
        "viga": nombre_viga,
        "tramo": id_tramo,
        "tipo_tramo": "voladizo" if es_voladizo else "interior",
        "L_m": L_viga,
        "q_serv_kN_m": cargas["servicio"],
        "delta_mm": delta,
        "lim_L250_mm": L / 250,
        "lim_L360_mm": L / 360,
        "lim_L480_mm": L / 480,
        "cumple_L250": delta <= L / 250,
        "cumple_L360": delta <= L / 360,
        "cumple_L480": delta <= L / 480,
        "Ie_tramo_m4": inercias["tramo"]["Ie"],
        "Ie_apoyo_izq_m4": inercias["apoyo_izq"]["Ie"],
        "Ie_apoyo_der_m4": inercias["apoyo_der"]["Ie"],
        "criterio": "Carga virtual – rigidez variable"
    })

    # Corte (interior)
    RESULTADOS["corte"].append({
        "tramo": tramo_id,
        "V_izq_kN": viga_data.get("Vu_izq", 0.0),
        "V_der_kN": viga_data.get("Vu_der", 0.0),
        "cumple": True
    })
    # -----------------------------
    # 10 ABERTURA DE FISURAS
    # -----------------------------
    beta = 1.25  # viga, 1.35 losa
    fy = notas_tecnicas["materiales"]["fy_MPa"]
    fs = 2/3 * fy  # MPa
    # recubrimiento hasta centroide de la barra más cercana [mm]
    diametro_barra_mm = arm_inf_tra["diam"]  # tomamos la barra inferior del tramo
    dc = notas_tecnicas["seccion"]["rec_inf_cm"]*10 + (diametro_barra_mm/2)

    # número de barras y área efectiva
    n_barras = arm_inf_tra["n"]  # total de barras inferiores
    A_barra_cm2 = arm_inf_tra["area_total_cm2"] / n_barras
    A_barra_mm2 = A_barra_cm2 * 100  # convertir a mm²

    # cálculo aproximado de apertura de fisuras
    w_um = 0.011 * beta * fs * (dc * A_barra_mm2) ** (1/3)  # µm
    w_mm = w_um / 1000  # convertir a mm

    lineas.append("")
    lineas.append("10 ABERTURA DE FISURAS")
    lineas.append(f"   dc = {dc:.1f} mm | fs = {fs:.1f} MPa | A_barra = {A_barra_mm2:.1f} mm²")
    lineas.append(f"   w ≈ {w_um:.1f} μm / {w_mm:.3f} mm")
    # límites según CIRSOC 201
    limites_w_mm = {
        "aire_seco": 0.41,
        "aire_humedo": 0.30,
        "contencion_agua": 0.10
    }

    # elegimos el límite más restrictivo para comparar
    w_lim_mm = min(limites_w_mm.values())  # 0.10 mm
    cumple = w_mm <= w_lim_mm

    lineas.append(f"   Limite CIRSOC = {w_lim_mm:.2f} mm → {'✅ Cumple' if cumple else '❌ No cumple'}")

    return "\n".join(lineas)

def procesar_vigas(datos_vigas, coef_kd, carpeta_salida):
    resultados = {}  # diccionario para acumular resultados

    for viga_id, viga in datos_vigas.items():
        b = viga["b_cm"]
        h = viga["h_cm"]
        rec = viga["recubrimiento_cm"]
        fc = viga.get("fc_MPa", 21.0)
        fy = viga.get("fy_MPa", 420.0)

        resultados[viga_id] = []  # cada viga tendrá una lista de tramos

        for tramo in viga["tramos"]:
            if tramo["es_voladizo"]:
                datos_tramo = {
                    "m_izq": tramo["Mu_kNm"],
                    "m_tra": tramo["Mu_kNm"],
                    "m_der": 0.0,
                    "Vu_emp": tramo["Vu_kN_emp"],
                    "Vu_izq": tramo["Vu_kN_emp"],  # opcional
                    "Vu_der": 0.0
                }
            else:
                datos_tramo = {
                    "m_izq": tramo["Mu_kNm_apoyo_izq"],
                    "m_tra": tramo["Mu_kNm_campo"],
                    "m_der": tramo["Mu_kNm_apoyo_der"],
                    "Vu_izq": tramo["Vu_kN_izq"],
                    "Vu_der": tramo["Vu_kN_der"]
                }

            # Generar texto de planilla
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
                puntos_inflexion_m=tramo.get("puntos_inflexion_m"),
                cargas=viga["cargas"]
            )

            # Guardar TXT en carpeta de salida
            with open(carpeta_salida / f"planilla_{viga_id}_{tramo['id']}.txt", "w", encoding="utf-8") as f:
                f.write(texto)

            print(f"✅ Guardado: {carpeta_salida / f'planilla_{viga_id}_{tramo['id']}.txt'}")

            # Acumular resultados en el diccionario
            resultados[viga_id].append({
                "tramo_id": tramo["id"],
                "longitud_m": tramo["longitud_m"],
                "es_voladizo": tramo["es_voladizo"],
                "datos_tramo": datos_tramo
            })

    return resultados

def guardar_resultados_json(ruta, resultados):
    with open(ruta, "w", encoding="utf-8") as f:
        json.dump(resultados, f, indent=2, ensure_ascii=False)


# =========================================================
# PROGRAMA PRINCIPAL
# =========================================================
BASE = Path(__file__).parent

with open(BASE / "datos" / "moments_input.json", encoding="utf-8") as f:
    datos_vigas = json.load(f)

with open(BASE / "datos" / "coeficientes_kd.json", encoding="utf-8") as f:
    coef_kd = json.load(f)

salidas = BASE / "salidas"
vigas_dir = salidas / "vigas"
vigas_dir.mkdir(parents=True, exist_ok=True)

# 🔹 solo llamás a procesar_vigas para que llene RESULTADOS global
procesar_vigas(datos_vigas, coef_kd, vigas_dir)

# 🔹 guardás el RESULTADOS global que ya se fue llenando
guardar_resultados_json(vigas_dir / "resultados_vigas.json", RESULTADOS)

print(f"✅ Resultados guardados en: {vigas_dir / 'resultados_vigas.json'}")



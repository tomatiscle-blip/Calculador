import math

class DisenadorViga:

    # =====================================================
    # MATERIALES
    # =====================================================
    def Ec_MPa(self):
        return 4700.0 * math.sqrt(self.fc)

    # =====================================================
    # PROPIEDADES DE SECCIÓN
    # =====================================================
    def propiedades_seccion(self):
        Ig = self.b * self.h**3 / 12.0
        return {
            "b_m": self.b,
            "h_m": self.h,
            "d_m": self.d,
            "rec_m": self.rec,
            "Ig_m4": Ig
        }

    # =====================================================
    # INERCIAS
    # =====================================================
    def inercias_seccion(self, As_cm2, Kc=0.40, M_kNm=None):
        """
        Ig  : inercia bruta
        Icr : inercia fisurada
        Mcr : momento de fisuración
        Ie  : inercia efectiva (Branson)
        """
        # --- Inercia bruta ---
        Ig = self.b * self.h**3 / 12.0

        # --- Momento de fisuración ---
        fr_MPa = 0.62 * math.sqrt(self.fc)
        yt = self.h / 2.0
        Mcr = fr_MPa * 1e6 * Ig / yt / 1000.0  # kNm

        # --- Inercia fisurada ---
        As = As_cm2 / 10000.0
        x = Kc * self.d
        Icr = (self.b * x**3) / 3.0 + As * (self.d - x)**2

        # --- Inercia efectiva (Branson) ---
        if M_kNm is not None and M_kNm > Mcr:
            ratio = (Mcr / M_kNm) ** 3
            Ie = ratio * Ig + (1.0 - ratio) * Icr
        else:
            Ie = Ig

        return {
            "Ig_m4": Ig,
            "Icr_m4": Icr,
            "Ie_m4": Ie,
            "Mcr_kNm": Mcr
        }

    # =====================================================
    # FLECHAS
    # =====================================================
    def flecha_instantanea(self, L_m, q_kN_m, Ie_m4):
        """
        Viga biapoyada – carga distribuida
        """
        q = q_kN_m * 1000      # N/m
        E = self.Ec_MPa() * 1e6
        L = L_m

        delta = (5 * q * L**4) / (384 * E * Ie_m4)
        return delta  # m

    def flecha_final(self, delta_inst_m, lambda_dif=2.0):
        return delta_inst_m * (1.0 + lambda_dif)

    def chequeo_flecha(self, delta_m, L_m, limite=250):
        adm = L_m / limite
        return {
            "delta_mm": delta_m * 1000,
            "adm_mm": adm * 1000,
            "cumple": delta_m <= adm
        }

    # =====================================================
    # FISURACIÓN
    # =====================================================
    def ancho_fisura(self, As_cm2, diam_mm, Ma_kNm):
        """
        Modelo simplificado ACI / CIRSOC
        """
        Es = 200000.0  # MPa
        As = As_cm2 / 10000.0
        z = 0.9 * self.d

        sigma_s = (Ma_kNm * 1e6) / (As * z) / 1e6  # MPa
        sr_mm = 150 + 0.25 * diam_mm
        w_mm = (sigma_s / Es) * sr_mm

        return {
            "sigma_s_MPa": sigma_s,
            "w_mm": w_mm
        }

    # =====================================================
    # HORMIGÓN
    # =====================================================
    def hormigon_tramo(self, L_m, gamma_kg_m3=2500):
        V_m3 = self.b * self.h * L_m
        peso_kg = V_m3 * gamma_kg_m3
        return {
            "V_m3": V_m3,
            "peso_kg": peso_kg
        }

    # =====================================================
    # CÓMPUTO DE ACERO
    # =====================================================
    def computo_acero(self, barras):
        """
        barras = [
            { "diam": mm, "L_m": m, "n": cant }
        ]
        """
        resumen = {}

        for b in barras:
            d = b["diam"]
            peso_ml = self.barras_comerciales[d]
            kg = b["L_m"] * b["n"] * peso_ml

            if d not in resumen:
                resumen[d] = 0.0
            resumen[d] += kg

        total = sum(resumen.values())

        return {
            "por_diametro": resumen,
            "total_kg": total
        }

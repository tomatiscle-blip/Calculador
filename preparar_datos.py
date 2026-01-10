import json
import os
import math

def d_efectivo(h_cm, rec_cm, diam_mm):
    h = h_cm / 100.0
    rec = rec_cm / 100.0
    phi = diam_mm / 1000.0
    return h - rec - phi/2.0

def as_min(bw_cm, d_m, fc_MPa, fy_MPa):
    bw = bw_cm / 100.0
    area_bd = bw * d_m
    criterio_1 = (math.sqrt(fc_MPa) / (4.0 * fy_MPa)) * area_bd
    criterio_2 = (1.4 * area_bd) / fy_MPa
    return {
        "criterio_1_m2": criterio_1,
        "criterio_2_m2": criterio_2,
        "as_min_m2": max(criterio_1, criterio_2)
    }

def as_max(bw_cm, d_m, fc_MPa, fy_MPa, beta1=0.85):
    bw = bw_cm / 100.0
    # cuantía balanceada según tu reglamento
    rho_bal = 0.85 * beta1 * (fc_MPa / fy_MPa) * (600.0 / (600.0 + fy_MPa))
    # límite máximo: 0.75 * rho_bal
    rho_max = 0.75 * rho_bal
    as_max_val = rho_max * bw * d_m
    return {
        "rho_bal": rho_bal,
        "rho_max": rho_max,
        "as_max_m2": as_max_val
    }

def preparar_datos_y_guardar():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    carpeta = os.path.join(base_dir, "datos")
    os.makedirs(carpeta, exist_ok=True)

    archivo_in = os.path.join(carpeta, "moments_input.json")
    archivo_out = os.path.join(carpeta, "moments_preparado.json")

    if not os.path.exists(archivo_in):
        print("No se encontró moments_input.json")
        return

    with open(archivo_in, "r") as f:
        datos = json.load(f)

    diam_mm = 16
    phi_m = diam_mm / 1000.0
    Ld_m = 50.0 * phi_m

    preparado = {}

    for viga_id, viga in datos.items():
        b = viga["b_cm"] / 100.0
        h = viga["h_cm"] / 100.0
        rec = viga["recubrimiento_cm"] / 100.0
        fc = viga["fc_MPa"]
        fy = viga["fy_MPa"]

        Ec = 4700 * math.sqrt(fc)
        I = b * (h**3) / 12
        d = d_efectivo(viga["h_cm"], viga["recubrimiento_cm"], diam_mm)

        as_min_vals = as_min(viga["b_cm"], d, fc, fy)
        as_max_vals = as_max(viga["b_cm"], d, fc, fy)

        viga_out = {
            "seccion": {
                "b_cm": viga["b_cm"],
                "h_cm": viga["h_cm"],
                "recubrimiento_cm": viga["recubrimiento_cm"],
                "Ec_MPa": round(Ec),
                "I_m4": float(f"{I:.6f}"),
                "d_m": round(d, 3),
                "As_min": as_min_vals,
                "As_max": as_max_vals
            },
            "tramos": []
        }

        for tramo in viga["tramos"]:
            L = tramo["longitud_m"]
            es_vol = tramo["es_voladizo"]
            tramo_out = {
                "id": tramo["id"],
                "longitud_m": L,
                "es_voladizo": es_vol
            }

            if es_vol:
                tramo_out.update({
                    "momentos": {"voladizo_kNm": tramo["Mu_kNm"]},
                    "cortantes": {"emp_kN": tramo["Vu_kN_emp"]},
                    "validacion": {
                        "voladizo_sin_pi_ok": ("puntos_inflexion_m" not in tramo or not tramo.get("puntos_inflexion_m")),
                        "observaciones": []
                    },
                    "prolongaciones": {
                        "diam_mm": diam_mm,
                        "Ld_m": Ld_m,
                        "hasta_extremo_m": L,
                        "anclaje_apoyo_m": Ld_m
                    }
                })
            else:
                puntos = tramo.get("puntos_inflexion_m", [])
                pi_rango_ok = bool(puntos) and all(0 < p < L for p in puntos)
                pi_orden_ok = bool(puntos) and all(puntos[i] < puntos[i+1] for i in range(len(puntos)-1))
                obs = []
                if not puntos:
                    obs.append("Tramo no voladizo sin puntos de inflexión.")
                if puntos and not pi_rango_ok:
                    obs.append("Hay puntos de inflexión fuera de rango.")
                if puntos and not pi_orden_ok:
                    obs.append("Los puntos de inflexión no están en orden estricto.")

                tramo_out.update({
                    "momentos": {
                        "campo_kNm": tramo["Mu_kNm_campo"],
                        "apoyo_izq_kNm": tramo["Mu_kNm_apoyo_izq"],
                        "apoyo_der_kNm": tramo["Mu_kNm_apoyo_der"]
                    },
                    "cortantes": {
                        "izq_kN": tramo["Vu_kN_izq"],
                        "der_kN": tramo["Vu_kN_der"]
                    },
                    "puntos_inflexion_m": puntos,
                    "validacion": {
                        "pi_rango_ok": bool(puntos) and pi_rango_ok,
                        "pi_orden_ok": bool(puntos) and pi_orden_ok,
                        "voladizo_sin_pi_ok": True,
                        "observaciones": obs
                    }
                })

                # Prolongaciones en función de Pi y Ld
                if puntos and pi_rango_ok and pi_orden_ok:
                    prolong_desde_pi = [round(p + Ld_m, 3) for p in puntos]
                    pi_izq = puntos[0]
                    pi_der = puntos[-1]
                    tramo_out["prolongaciones"] = {
                        "diam_mm": diam_mm,
                        "Ld_m": Ld_m,
                        "desde_pi_m": prolong_desde_pi,
                        "apoyo_izq_m": round(pi_izq + Ld_m, 3),
                        "apoyo_der_m": round(pi_der + Ld_m, 3)
                    }
                else:
                    ext = max(0.3 * L, Ld_m)
                    tramo_out["prolongaciones"] = {
                        "diam_mm": diam_mm,
                        "Ld_m": Ld_m,
                        "criterio_conservador_m": round(ext, 3)
                    }

            viga_out["tramos"].append(tramo_out)
        preparado[viga_id] = viga_out

    with open(archivo_out, "w") as f:
        json.dump(preparado, f, indent=2)
    print(f"\nDatos preparados guardados en {archivo_out}")

if __name__ == "__main__":
    preparar_datos_y_guardar()
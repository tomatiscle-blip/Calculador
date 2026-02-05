import json
import os
import math

def ingresar_datos():
    viga_id = input("ID de la viga: ")

    b_cm = float(input("Ancho b [cm]: "))
    h_cm = float(input("Altura h [cm]: "))
    rec_cm = float(input("Recubrimiento [cm]: "))
    rec_sup_cm = float(input("Recubrimiento superior [cm]: "))
    fc = float(input("f'c [MPa]: "))
    fy = float(input("fy [MPa]: "))
    # 🔹 Cargas distribuidas (nuevo bloque)
    D = float(input("Carga muerta D [kN/m]: "))
    L = float(input("Carga viva L [kN/m]: "))

    q_14D = 1.4 * D
    q_12D16L = 1.2 * D + 1.6 * L
    q_servicio = D + L

    cargas = {
        "D": D,
        "L": L,
        "1.4D": q_14D,
        "1.2D+1.6L": q_12D16L,
        "servicio": q_servicio
    }


    tramos = []
    n_tramos = int(input("Cantidad de tramos: "))
    for i in range(n_tramos):
        print(f"\n--- Tramo {i+1} ---")
        tramo_id = input("ID tramo (ej. V01-T1, V01-T2_v): ")
        longitud_m = float(input("Longitud del tramo [m]: "))
        es_voladizo = input("¿Es voladizo? (s/n): ").lower() == "s"

        if es_voladizo:
            Mu = float(input("Momento Mu [kN·m] (con signo): "))
            Vu = float(input("Cortante Vu en empotramiento [kN]: "))
            tramos.append({
                "id": tramo_id,
                "longitud_m": longitud_m,
                "es_voladizo": True,
                "Mu_kNm": Mu,
                "Vu_kN_emp": Vu
            })
        else:
            Mu_campo = float(input("Momento en campo Mu [kN·m]: "))
            Mu_ap_izq = float(input("Momento en apoyo izq Mu [kN·m]: "))
            Mu_ap_der = float(input("Momento en apoyo der Mu [kN·m]: "))
            Vu_izq = float(input("Cortante en apoyo izq Vu [kN]: "))
            Vu_der = float(input("Cortante en apoyo der Vu [kN]: "))
            
            # Puntos de inflexión
            texto = input("Puntos de inflexión [m], separados por coma: ").strip()
            puntos = [float(x.strip()) for x in texto.split(",")] if texto else []

            tramos.append({
                "id": tramo_id,
                "longitud_m": longitud_m,
                "es_voladizo": False,
                "Mu_kNm_campo": Mu_campo,
                "Mu_kNm_apoyo_izq": Mu_ap_izq,
                "Mu_kNm_apoyo_der": Mu_ap_der,
                "Vu_kN_izq": Vu_izq,
                "Vu_kN_der": Vu_der,
                "puntos_inflexion_m": puntos
            })

    nueva_viga = {
        viga_id: {
            "b_cm": b_cm,
            "h_cm": h_cm,
            "recubrimiento_cm": rec_cm,
            "recubrimiento_sup_cm": rec_sup_cm,
            "fc_MPa": fc,
            "fy_MPa": fy,
            "cargas": cargas,
            "tramos": tramos
        }
    }

    base_dir = os.path.dirname(os.path.abspath(__file__))
    carpeta = os.path.join(base_dir, "datos")
    os.makedirs(carpeta, exist_ok=True)
    archivo = os.path.join(carpeta, "moments_input.json")

    if os.path.exists(archivo):
        with open(archivo, "r") as f:
            datos = json.load(f)
    else:
        datos = {}

    datos.update(nueva_viga)

    with open(archivo, "w") as f:
        json.dump(datos, f, indent=2)

    print(f"\nDatos guardados en {archivo}")

    # 🔹 Preparación de datos: resumen técnico
    preparar_datos(datos)

def preparar_datos(datos):
    for viga_id, viga in datos.items():
        print(f"\n=== Viga {viga_id} ===")
        b = viga["b_cm"] / 100.0
        h = viga["h_cm"] / 100.0
        rec = viga["recubrimiento_cm"] / 100.0
        fc = viga["fc_MPa"]
        fy = viga["fy_MPa"]

        Ec = 4700 * math.sqrt(fc)  # MPa
        I = b * (h**3) / 12        # m^4
        rho_min = 0.7 * math.sqrt(fc) / fy

        print(f"Sección: {b*100:.1f} x {h*100:.1f} cm, recubrimiento {rec*100:.1f} cm")
        print(f"Ec = {Ec:.0f} MPa, I = {I:.6f} m^4, ρmin = {rho_min:.4f}")

        for tramo in viga["tramos"]:
            print(f"\n--- Tramo {tramo['id']} ---")
            L = tramo["longitud_m"]

            if tramo["es_voladizo"]:
                print(f"Voladizo de {L:.2f} m")
                print(f"Mu = {tramo['Mu_kNm']} kN·m, Vu_emp = {tramo['Vu_kN_emp']} kN")
            else:
                print(f"Tramo normal de {L:.2f} m")
                print(f"Mu campo = {tramo['Mu_kNm_campo']} kN·m")
                print(f"Mu apoyo izq = {tramo['Mu_kNm_apoyo_izq']} kN·m")
                print(f"Mu apoyo der = {tramo['Mu_kNm_apoyo_der']} kN·m")
                print(f"Vu izq = {tramo['Vu_kN_izq']} kN, Vu der = {tramo['Vu_kN_der']} kN")

                puntos = tramo.get("puntos_inflexion_m", [])
                if not puntos:
                    print("⚠ Error: tramo no voladizo sin puntos de inflexión")
                else:
                    if any(p <= 0 or p >= L for p in puntos):
                        print("⚠ Error: Pi fuera de rango")
                    elif any(puntos[i] >= puntos[i+1] for i in range(len(puntos)-1)):
                        print("⚠ Error: Pi no ordenados")
                    else:
                        print(f"Puntos de inflexión: {puntos}")

if __name__ == "__main__":
    ingresar_datos()
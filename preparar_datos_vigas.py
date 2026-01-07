import json
import os

def ingresar_datos():
    viga_id = input("ID de la viga: ")

    b_cm = float(input("Ancho b [cm]: "))
    h_cm = float(input("Altura h [cm]: "))
    rec_cm = float(input("Recubrimiento [cm]: "))
    fc = float(input("f'c [MPa]: "))
    fy = float(input("fy [MPa]: "))

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
            "fc_MPa": fc,
            "fy_MPa": fy,
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

if __name__ == "__main__":
    ingresar_datos()
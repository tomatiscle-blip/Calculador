import pandas as pd
import json
import os
import re

def exportar_vigas_excel(json_dir, excel_path):
    acumulado = {}

    for archivo in os.listdir(json_dir):
        if archivo.startswith("resultados_Portico") and archivo.endswith("_vigas.json"):
            json_path = os.path.join(json_dir, archivo)
            with open(json_path, "r", encoding="utf-8") as f:
                resultados = json.load(f)

            # Identificar pórtico
            match = re.search(r"Portico_(\d+)", archivo)
            if match:
                portico_num = match.group(1)
                clave_portico = f"Portico_{portico_num}"
            else:
                clave_portico = os.path.splitext(archivo)[0]

            filas = []
            for tramo in resultados.get("tramos", []):
                tramo_id = tramo["id"]
                tramo_id_simple = tramo_id.split(".")[-1]

                fila = {
                    "Piso": tramo["viga"].split("-")[0],
                    "Tramo_ID": tramo_id,
                    "Tipo": tramo.get("tipo", ""),
                    "Longitud_m": tramo.get("L_m", 0)
                }

                # Materiales
                mat = next((m for m in resultados.get("materiales", []) if m.get("tramo") in [tramo_id, tramo_id_simple]), {})
                fila.update({
                    "fc_MPa": mat.get("fc_MPa", 0),
                    "fy_MPa": mat.get("fy_MPa", 0),
                    "Ec_MPa": mat.get("Ec_MPa", 0),
                    "b_cm": mat.get("b_cm", 0),
                    "h_cm": mat.get("h_cm", 0),
                    "recubrimiento_cm": mat.get("recubrimiento_cm", 0),
                    "V_m3": mat.get("V_m3", 0),
                    "peso_kg": mat.get("peso_kg", 0)
                })

                # Flexión
                flex = next((f for f in resultados.get("flexion", []) if f.get("tramo") in [tramo_id, tramo_id_simple]), {})
                fila.update({
                    "As_req_cm2": flex.get("As_req_cm2", 0),
                    "As_adop_cm2": flex.get("As_adop_cm2", 0),
                    "Cumple_flexion": flex.get("cumple", "")
                })

                # Inercias
                ine = next((i for i in resultados.get("inercias", []) if i.get("tramo") in [tramo_id, tramo_id_simple]), {})
                fila.update({
                    "Ig_cm4": ine.get("Ig_cm4", 0),
                    "Ief_cm4": ine.get("Ief_cm4", 0),
                    "ratio_Ief_Ig": ine.get("ratio_Ief_Ig", 0),
                    "Mcr_kNm": ine.get("Mcr_tramo_kNm", 0),
                    "c_usado_cm": ine.get("c_usado_cm", 0),
                    "criterio_Ief": ine.get("criterio_Ief", "")
                })

                # Flecha
                fle = next((f for f in resultados.get("flecha", []) if f.get("tramo") in [tramo_id, tramo_id_simple]), {})
                fila.update({
                    "delta_mm": fle.get("delta_mm", 0),
                    "limite_mm": fle.get("lim_L360_mm", fle.get("lim_L180_mm", 0)),
                    "cumple_flecha": fle.get("cumple_L360", fle.get("cumple_L180", "")),
                    "Ie_tramo_cm4": fle.get("Ie_tramo_cm4", fle.get("Ie_cm4", 0)),
                    "Ie_apoyo_izq_cm4": fle.get("Ie_apoyo_izq_cm4", 0),
                    "Ie_apoyo_der_cm4": fle.get("Ie_apoyo_der_cm4", 0),
                    "criterio_flecha": fle.get("criterio", "")
                })

                # Corte
                cor = next((c for c in resultados.get("corte", []) if c.get("tramo") in [tramo_id, tramo_id_simple]), {})
                fila.update({
                    "V_izq_kN": cor.get("V_izq_kN", 0),
                    "V_der_kN": cor.get("V_der_kN", 0),
                    "V_kN": cor.get("V_kN", 0),
                    "cumple_corte": cor.get("cumple", "")
                })

                # Fisuración
                fis = next((fi for fi in resultados.get("fisuracion", []) if fi.get("tramo") in [tramo_id, tramo_id_simple]), {})
                fila.update({
                    "dc_mm": fis.get("dc_mm", 0),
                    "fs_MPa": fis.get("fs_MPa", 0),
                    "A_barra_mm2": fis.get("A_barra_mm2", 0),
                    "w_mm": fis.get("w_mm", 0),
                    "w_lim_mm": fis.get("w_lim_mm", 0),
                    "cumple_fisuracion": fis.get("cumple", "")
                })

                # Estado
                est = next((e for e in resultados.get("estado", []) if e.get("tramo") in [tramo_id, tramo_id_simple]), {})
                fila.update({
                    "Estado": est.get("flexion", ""),
                    "Nota": est.get("nota", "")
                })

                # Armaduras → columnas separadas
                armaduras = [a for a in resultados.get("armaduras", []) if a.get("tramo") in [tramo_id, tramo_id_simple]]
                for i, a in enumerate(armaduras, start=1):
                    fila[f"Armadura_{i}_Pos"] = a.get("pos", 0)
                    fila[f"Armadura_{i}_Detalle"] = a.get("detalle", "")
                    fila[f"Armadura_{i}_Diam_mm"] = a.get("diam_mm", 0)
                    fila[f"Armadura_{i}_Cantidad"] = a.get("cantidad", 0)
                    fila[f"Armadura_{i}_Longitud_m"] = a.get("longitud_m", 0)
                    fila[f"Armadura_{i}_Peso_kg"] = a.get("peso_kg", 0)

                # Estribos → columnas separadas
                estribos = [e for e in resultados.get("estribos", []) if e.get("tramo") in [tramo_id, tramo_id_simple]]
                for i, e in enumerate(estribos, start=1):
                    fila[f"Estribo_{i}_Zona"] = e.get("zona", "")
                    fila[f"Estribo_{i}_Diam_mm"] = e.get("diam_mm", 0)
                    fila[f"Estribo_{i}_Cantidad"] = e.get("cantidad", 0)
                    fila[f"Estribo_{i}_Separacion_cm"] = e.get("separacion_cm", 0)
                    fila[f"Estribo_{i}_Longitud_m"] = e.get("longitud_m", 0)
                    fila[f"Estribo_{i}_Peso_kg"] = e.get("peso_kg", 0)

                filas.append(fila)

            acumulado[clave_portico] = filas

    # Exportar a Excel, cada pórtico en hoja propia
    with pd.ExcelWriter(excel_path, engine="xlsxwriter") as writer:
        for portico, filas in acumulado.items():
            df = pd.DataFrame(filas)
            df.to_excel(writer, sheet_name=portico[:31], index=False)

if __name__ == "__main__":
    base_dir = os.path.dirname(__file__)
    vigas_dir = os.path.join(base_dir, "salidas", "vigas")
    os.makedirs(vigas_dir, exist_ok=True)

    excel_file = os.path.join(vigas_dir, "Planilla_Vigas_Portico.xlsx")

    exportar_vigas_excel(vigas_dir, excel_file)
    print(f"✅ Excel generado en: {excel_file}")



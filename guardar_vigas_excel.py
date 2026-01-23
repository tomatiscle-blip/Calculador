import pandas as pd
import json
import os

def exportar_vigas_excel(json_path, excel_path):
    # Leer el JSON
    with open(json_path, "r", encoding="utf-8") as f:
        resultados = json.load(f)

    # Crear Excel con cada bloque en una hoja
    with pd.ExcelWriter(excel_path, engine="xlsxwriter") as writer:
        for clave, lista in resultados.items():
            if lista:  # solo si hay datos
                df = pd.DataFrame(lista)
                df.to_excel(writer, sheet_name=clave, index=False)

if __name__ == "__main__":
    base_dir = os.path.dirname(__file__)
    vigas_dir = os.path.join(base_dir, "salidas", "vigas")
    os.makedirs(vigas_dir, exist_ok=True)

    json_file = os.path.join(vigas_dir, "resultados_vigas.json")
    excel_file = os.path.join(vigas_dir, "planilla_vigas.xlsx")

    exportar_vigas_excel(json_file, excel_file)
    print(f"✅ Excel generado en: {excel_file}")

# tests/probar_varios_muros.py
"""
Script de pruebas y análisis rápido para diferentes tipos de muros de mampostería
según el módulo calc.mamposteria (CIRSOC 501).
"""

import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent
if str(RAIZ) not in sys.path:
    sys.path.insert(0, str(RAIZ))

from calc import mamposteria as mamp  # noqa: E402
import pprint  # noqa: E402


def ejecutar_pruebas_muros():
    print("=" * 60)
    print(" ANÁLISIS Y VERIFICACIÓN DE DISTINTOS TIPOS DE MUROS")
    print("=" * 60)

    # Definimos casos de estudio con diferentes espesores y componentes
    casos = [
        {
            "nombre": "Muro Portante de 18 cm con Mortero Normal (Carga moderada)",
            "capas": [
                {"clave": "cemento_arena", "grupo": "Morteros", "espesor_m": 0.015},
                {"clave": "ladrillo_hueco_portante", "grupo": "Mamposteria", "espesor_m": 0.18},
                {"clave": "cal_arena", "grupo": "Morteros", "espesor_m": 0.015},
            ],
            "altura_m": 2.60,
            "mortero": "mortero_normal_n",
            "armada": False,
            "carga_axial_kN_m": 35.0
        },
        {
            "nombre": "Muro Portante de 12 cm con Mortero Normal (Más esbelto)",
            "capas": [
                {"clave": "ladrillo_hueco_portante", "grupo": "Mamposteria", "espesor_m": 0.12},
                {"clave": "cal_arena", "grupo": "Morteros", "espesor_m": 0.015},
            ],
            "altura_m": 2.80,
            "mortero": "mortero_normal_n",
            "armada": False,
            "carga_axial_kN_m": 20.0
        },
        {
            "nombre": "Muro Portante de 20 cm con Mortero Elevado y Carga Alta",
            "capas": [
                {"clave": "ladrillo_hueco_portante", "grupo": "Mamposteria", "espesor_m": 0.20},
            ],
            "altura_m": 3.20,
            "mortero": "mortero_elevado_e",
            "armada": True,
            "carga_axial_kN_m": 75.0
        }
    ]

    for i, caso in enumerate(casos, 1):
        print(f"\n--- CASO {i}: {caso['nombre']} ---")
        
        # Ejecutamos la verificación usando tu módulo
        res = mamp.verificar_muro(
            capas=caso["capas"],
            altura_m=caso["altura_m"],
            mortero=caso["mortero"],
            armada=caso["armada"],
            carga_axial_kN_m=caso["carga_axial_kN_m"]
        )

        if not res["ok"]:
            print(f"  ❌ No se pudo verificar: {res['pendiente']}")
            continue

        # Mostramos los resultados clave de forma limpia
        portante = res["hoja_portante"]
        comp = res["compresion"]
        
        print(f"  • Material hoja portante : {portante['clave']}")
        print(f"  • Resistencia fm         : {portante['fm_mpa']} MPa")
        print(f"  • Espesor                : {portante['espesor_m'] * 100:.1f} cm")
        print(f"  • Esbeltez (h/r)         : {res['esbeltez_hr']:.2f} (¿Esbelta?: {'Sí ⚠' > '' if res['es_esbelta'] else 'No'})")
        print(f"  • Capacidad axial φ·Pn   : {comp['phi_Pn_kN']:.2f} kN/m")
        print(f"  • Carga axial aplicada   : {comp['carga_axial_kN_m']} kN/m")
        
        if comp["pasa"] is not None:
            estado = "✅ ¡PASA!" if comp["pasa"] else "❌ NO PASA (Sobrecargado)"
            print(f"  • Estado compresión      : {estado}")
            
        print("  • Recomendaciones técnicas:")
        for rec in res["recomendaciones"]:
            print(f"    - {rec}")

    print("\n" + "=" * 60)


if __name__ == "__main__":
    ejecutar_pruebas_muros()
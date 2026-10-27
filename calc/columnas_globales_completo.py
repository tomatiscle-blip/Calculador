from __future__ import annotations
import sys
from pathlib import Path
_RAIZ = Path(__file__).resolve().parent.parent
if str(_RAIZ) not in sys.path:
    sys.path.insert(0, str(_RAIZ))

def main():
    from calc import columnas_persistencia, columnas_unificacion, columnas_biaxial_global, columnas_planilla_final
    print("Fase 1: persistencia...")
    r1 = columnas_persistencia.generar_dimensionados_json()
    print(f"  {len(r1)} JSON generados")
    print("Fase 2: unificacion...")
    r2 = columnas_unificacion.unificar()
    print(f"  {len(r2)} columnas unificadas")
    print("Fase 4: biaxial...")
    r4 = columnas_biaxial_global.ejecutar()
    print(f"  {len(r4)} resultados finales")
    print("Planilla final...")
    columnas_planilla_final.generar()
    print("Listo.")

if __name__ == "__main__":
    main()

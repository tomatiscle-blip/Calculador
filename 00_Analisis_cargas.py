"""
00_Analisis_cargas.py — Análisis de cargas (envoltorio de la etapa 1).

La cuenta ya NO vive en este archivo: vive en `calc/cargas.py`, y los datos en
`datos/cargas.json` (los elementos que reciben carga) y `datos/materiales.json`
(la biblioteca de materiales). Acá solo se llama, se imprime y se guarda el
informe, así la etapa sigue funcionando igual que siempre (el .txt que lee P00
se guarda con el mismo nombre y el mismo texto).

    Para ver las cargas sin guardar nada:  py -m calc.cargas
    Para ver UN elemento solo:             py -m calc.cargas "Losa Alivianada L0-1"
    Para ver qué elementos hay:            py -m calc.cargas --lista
"""

import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parent
if str(RAIZ) not in sys.path:
    sys.path.insert(0, str(RAIZ))

from calc import cargas  # noqa: E402


def main() -> int:
    analisis = cargas.cargas_del_conjunto()
    texto = cargas.informe_completo(analisis)
    print(texto)
    archivo = cargas.guardar_informe(analisis, texto)
    print(f"\nResumen guardado en: {archivo}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

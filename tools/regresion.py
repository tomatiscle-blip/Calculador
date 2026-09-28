"""
tools.regresion — Red de seguridad: congelar resultados y detectar cambios.

Para qué sirve
--------------
Vamos a tocar el motor de cálculo (anaStruct 2D -> Pynite 3D) y a reordenar
los scripts. Si no hay una foto de "cómo daban los números antes", cualquier
cambio es a ciegas. Esto saca esa foto y después compara.

Uso (desde la carpeta del proyecto o desde cualquier lado):

    py tools\\regresion.py congelar          -> guarda la referencia actual
    py tools\\regresion.py comparar          -> avisa qué cambió respecto de la referencia
    py tools\\regresion.py comparar --diff   -> además muestra las diferencias
    py tools\\regresion.py ver <archivo>      -> muestra la referencia guardada

La referencia vive en tests/golden/ (no se versiona en git).
"""

from __future__ import annotations

import difflib
import json
import shutil
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent
if str(RAIZ) not in sys.path:
    sys.path.insert(0, str(RAIZ))

from calc import rutas  # noqa: E402  (después de arreglar sys.path)

GOLDEN = RAIZ / "tests" / "golden"
MANIFIESTO = GOLDEN / "manifiesto.json"

# Qué se congela: los resultados que hoy damos por buenos.
OBJETIVOS = (
    "datos/estructura.json",
    "datos/moments_input.json",
    "salidas/analisis_cargas",
    "salidas/vigas",
    "salidas/columnas",
    "salidas/bases",
    "salidas/losas",
    "salidas/dxf",
)

EXTENSIONES_BINARIAS = {".xlsx", ".xls", ".png", ".jpg", ".pdf"}


def _archivos() -> list[Path]:
    """Lista de archivos a congelar/comparar, en orden."""
    encontrados: list[Path] = []
    for relativo in OBJETIVOS:
        camino = RAIZ / relativo
        if camino.is_file():
            encontrados.append(camino)
        elif camino.is_dir():
            encontrados.extend(sorted(f for f in camino.rglob("*") if f.is_file()))
    return encontrados


def congelar() -> int:
    """Copia los resultados actuales como referencia y guarda el manifiesto."""
    archivos = _archivos()
    if not archivos:
        print("No se encontró ningún archivo de resultados para congelar.")
        return 1

    if GOLDEN.exists():
        shutil.rmtree(GOLDEN)
    GOLDEN.mkdir(parents=True, exist_ok=True)

    manifiesto = {}
    for archivo in archivos:
        relativo = archivo.relative_to(RAIZ).as_posix()
        destino = GOLDEN / relativo
        destino.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(archivo, destino)
        manifiesto[relativo] = {
            "huella": rutas.huella(archivo),
            "bytes": archivo.stat().st_size,
        }

    rutas.guardar_json(MANIFIESTO, manifiesto, indent=2)
    print(f"Referencia congelada: {len(archivos)} archivos en tests/golden/")
    print("A partir de ahora, 'py tools\\regresion.py comparar' avisa si cambian.")
    return 0


def _diff_de(actual: Path, referencia: Path, max_lineas: int = 40) -> str:
    if actual.suffix.lower() in EXTENSIONES_BINARIAS:
        return "  (archivo binario: cambió el contenido)"
    try:
        a = actual.read_text(encoding="utf-8", errors="replace").splitlines()
        b = referencia.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError as exc:
        return f"  (no se pudo leer: {exc})"
    lineas = list(
        difflib.unified_diff(b, a, fromfile="antes", tofile="ahora", lineterm="", n=1)
    )
    if len(lineas) > max_lineas:
        lineas = lineas[:max_lineas] + [f"... ({len(lineas) - max_lineas} líneas más)"]
    return "\n".join("  " + l for l in lineas)


def comparar(mostrar_diff: bool = False) -> int:
    """Compara los resultados actuales contra la referencia congelada."""
    manifiesto = rutas.leer_json(MANIFIESTO)
    if not manifiesto:
        print("Todavía no hay referencia. Corré primero:  py tools\\regresion.py congelar")
        return 1

    iguales, cambiados, faltantes = [], [], []
    for relativo in sorted(manifiesto):
        referencia = GOLDEN / relativo
        actual = RAIZ / relativo
        if not actual.exists():
            faltantes.append(relativo)
        elif rutas.huella(actual) == manifiesto[relativo]["huella"]:
            iguales.append(relativo)
        else:
            cambiados.append(relativo)

    actuales = {a.relative_to(RAIZ).as_posix() for a in _archivos()}
    nuevos = sorted(actuales - set(manifiesto))

    print("COMPARACIÓN CONTRA LA REFERENCIA")
    print("=" * 78)
    print(f"Sin cambios : {len(iguales)}")
    print(f"Cambiados   : {len(cambiados)}")
    print(f"Faltantes   : {len(faltantes)}")
    print(f"Nuevos      : {len(nuevos)}")

    for titulo, lista in (("CAMBIARON", cambiados), ("FALTAN", faltantes), ("NUEVOS", nuevos)):
        if not lista:
            continue
        print(f"\n{titulo}:")
        for relativo in lista:
            print(f"  - {relativo}")
            if mostrar_diff and titulo == "CAMBIARON":
                print(_diff_de(RAIZ / relativo, GOLDEN / relativo))

    return 0 if not (cambiados or faltantes) else 2


def ver(relativo: str) -> int:
    """Muestra la referencia guardada de un archivo."""
    archivo = GOLDEN / relativo.replace("\\", "/")
    if not archivo.exists():
        print(f"No hay referencia para {relativo}")
        return 1
    print(archivo.read_text(encoding="utf-8", errors="replace"))
    return 0


def _main(argv: list[str]) -> int:
    accion = argv[0] if argv else "comparar"
    if accion == "congelar":
        return congelar()
    if accion == "comparar":
        return comparar(mostrar_diff="--diff" in argv)
    if accion == "ver" and len(argv) > 1:
        return ver(argv[1])
    print(__doc__)
    return 1


if __name__ == "__main__":
    raise SystemExit(_main(sys.argv[1:]))

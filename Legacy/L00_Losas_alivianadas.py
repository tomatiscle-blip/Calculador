"""
L00_Losas_alivianadas.py — Cálculo de losa alivianada (envoltorio de consola).

La cuenta ya NO vive en este archivo: vive en `calc/losas.py` (la única cuenta
de losa). Acá solo se piden los datos por teclado, como siempre, y se llama a
`calc.losas.calcular`. Así el script sigue funcionando igual que antes, pero el
cálculo se puede usar también sin este archivo:

    py -m calc.losas --lista        qué losas hay (datos/losas.json)
    py -m calc.losas L00            calcula la losa L00 sin preguntar nada
    py -m calc.losas L00 --guardar  y la guarda en salidas/losas
"""

import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parent
if str(RAIZ) not in sys.path:
    sys.path.insert(0, str(RAIZ))

from calc import losas, materiales, rutas  # noqa: E402


def _elegir_materiales(biblioteca: dict, categoria: str) -> list[dict]:
    """Pide (por número) los materiales de una categoría. Enter = ninguno."""
    opciones = biblioteca.get(categoria, [])
    print(f"\n{categoria}:")
    for i, material in enumerate(opciones, start=1):
        print(f"{i}. {material['nombre']}")
    elegidos = input(
        f"Seleccione los números del {categoria} separados por coma "
        f"(o Enter para ninguno): "
    )
    seleccion = []
    if elegidos.strip():
        for x in elegidos.split(","):
            indice = int(x.strip()) - 1
            seleccion.append({"categoria": categoria, "nombre": opciones[indice]["nombre"]})
    return seleccion


def main() -> int:
    biblioteca = materiales.cargar()

    nombre = input("Nombre de la losa: ")

    seleccion: list[dict] = []
    for categoria in ("Pisos", "Contrapiso", "Cielorrasos"):
        seleccion.extend(_elegir_materiales(biblioteca, categoria))

    print("\nTipos de sobrecarga disponibles:")
    for i, sobrecarga in enumerate(biblioteca["Sobrecargas"], start=1):
        print(f"{i}. {sobrecarga['uso']} → {sobrecarga['valor_kNm2']} kN/m²")
    indice = int(input("Seleccione el número de sobrecarga: ")) - 1
    sobrecarga = biblioteca["Sobrecargas"][indice]
    print(f"Sobrecarga seleccionada: {sobrecarga['uso']}")

    luz_libre = float(input("Ingrese luz libre de la losa (m): "))
    ancho_losa = float(input("Ingrese ancho de la losa (m): "))

    datos = {
        "nombre": nombre,
        "luz_libre_m": luz_libre,
        "ancho_losa_m": ancho_losa,
        "materiales": seleccion,
        "sobrecarga": sobrecarga.get("clave") or sobrecarga["uso"],
    }

    resultado = losas.calcular(datos)
    print()
    print(resultado.memoria)

    # Mostrar la imagen de la vigueta elegida, como antes (si está disponible)
    try:
        from PIL import Image

        imagen = rutas.RAIZ / resultado.imagen
        if imagen.exists():
            Image.open(imagen).show()
            print(f"Imagen seleccionada: {resultado.imagen}")
    except Exception as exc:  # la imagen es un extra: nunca debe frenar el cálculo
        print(f"(no se pudo mostrar la imagen: {exc})")

    memoria, _computo = losas.guardar_resultado(resultado)
    print(f"\n📄 Memoria técnica generada en '{memoria}'")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


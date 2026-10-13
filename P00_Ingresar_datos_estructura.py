"""Lanzador compatible para cargar interactivamente un pórtico o un muro."""

from datos.crear_estructura import ingresar_datos_estructura, ingresar_muros

if __name__ == "__main__":
    # Por compatibilidad, podés llamar a pórticos directamente
    ingresar_datos_estructura()

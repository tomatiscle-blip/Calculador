"""
calc.rutas — Única fuente de verdad para las rutas del proyecto.

¿Por qué existe este archivo?
-----------------------------
Varios scripts viejos abren archivos con rutas "relativas" (por ejemplo
`open("datos/materiales.json")`). Eso funciona SOLO si el programa se lanza
parado en la carpeta del proyecto. Si mañana lo lanza la app, un acceso
directo del escritorio o un .exe, la ruta se rompe y el cálculo falla.

Acá todas las rutas se arman a partir de la ubicación REAL de este archivo,
así que funcionan desde cualquier lado.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Raíz del proyecto: ...\Calculador
# (rutas.py está en ...\Calculador\calc\rutas.py  ->  subir dos niveles)
# ---------------------------------------------------------------------------
RAIZ: Path = Path(__file__).resolve().parent.parent

DATOS: Path = RAIZ / "datos"
SALIDAS: Path = RAIZ / "salidas"
DIAGRAMAS_INTERACCION: Path = DATOS / "diagramas_interaccion"

# Carpetas de salida por etapa (mismos nombres que ya usa el proyecto)
SAL_ANALISIS_CARGAS: Path = SALIDAS / "analisis_cargas"
SAL_VIGAS: Path = SALIDAS / "vigas"
SAL_COLUMNAS: Path = SALIDAS / "columnas"
SAL_BASES: Path = SALIDAS / "bases"
SAL_DXF: Path = SALIDAS / "dxf"
SAL_LOSAS: Path = SALIDAS / "losas"
SAL_REACCIONES: Path = SALIDAS / "reacciones"
SAL_SOLICITACIONES: Path = SALIDAS / "solicitaciones"

CARPETAS_SALIDA = (
    SALIDAS,
    SAL_ANALISIS_CARGAS,
    SAL_VIGAS,
    SAL_COLUMNAS,
    SAL_BASES,
    SAL_DXF,
    SAL_LOSAS,
    SAL_REACCIONES,
    SAL_SOLICITACIONES,
)

# ---------------------------------------------------------------------------
# Archivos de datos (entradas del cálculo)
# ---------------------------------------------------------------------------
ESTRUCTURA: Path = DATOS / "estructura.json"
MATERIALES: Path = DATOS / "materiales.json"
VIGUETAS: Path = DATOS / "viguetas.json"
PERFILES_METALICOS: Path = DATOS / "perfiles_metalicos.json"
COEFICIENTES_KD: Path = DATOS / "coeficientes_kd.json"
MOMENTS_INPUT: Path = DATOS / "moments_input.json"

# Datos que la app todavía no tiene (se crean en la Fase D del plan)
TERRENO: Path = DATOS / "terreno.json"
CARGAS: Path = DATOS / "cargas.json"
TIPOS_LOSA: Path = DATOS / "tipos_losa.json"
LOSAS: Path = DATOS / "losas.json"

# ---------------------------------------------------------------------------
# Archivos de salida con nombre fijo
# ---------------------------------------------------------------------------
PLANILLA_COLUMNAS: Path = SAL_COLUMNAS / "planilla_columnas.csv"
COMPUTO_LOSAS: Path = SAL_LOSAS / "computo_losas.csv"


def ruta_salida(*partes: str) -> Path:
    """Devuelve una ruta dentro de `salidas/` (sin crearla)."""
    return SALIDAS.joinpath(*partes)


def asegurar_directorios() -> list[Path]:
    """Crea las carpetas de salida si no existen. Devuelve las carpetas usadas."""
    for carpeta in CARPETAS_SALIDA:
        carpeta.mkdir(parents=True, exist_ok=True)
    return list(CARPETAS_SALIDA)


# ---------------------------------------------------------------------------
# Nombres de archivo seguros
# ---------------------------------------------------------------------------
_CARACTERES_INVALIDOS = re.compile(r'[<>:"/\\|?*\x00-\x1f]')


def nombre_seguro(nombre: str) -> str:
    """
    Convierte un nombre (ej. "Portico 3(mercedes) ") en algo válido como
    nombre de archivo en Windows: sin caracteres prohibidos ni espacios o
    puntos al final (Windows los recorta y después el archivo "desaparece").
    """
    limpio = _CARACTERES_INVALIDOS.sub("_", str(nombre)).strip()
    return limpio.rstrip(". ") or "sin_nombre"


# ---------------------------------------------------------------------------
# Lectura y escritura de archivos
# ---------------------------------------------------------------------------
def existe(ruta: os.PathLike | str) -> bool:
    return Path(ruta).exists()


def leer_json(ruta: os.PathLike | str, por_defecto: Any = None) -> Any:
    """
    Lee un JSON. Si el archivo no existe o está roto devuelve `por_defecto`,
    así la app nunca se cae por un archivo faltante o mal guardado.
    """
    camino = Path(ruta)
    if not camino.exists():
        return por_defecto
    try:
        with open(camino, "r", encoding="utf-8") as f:
            return json.load(f)
    except (json.JSONDecodeError, OSError, UnicodeDecodeError):
        return por_defecto


def guardar_json(ruta: os.PathLike | str, datos: Any, indent: int = 2) -> Path:
    """
    Guarda un JSON de forma SEGURA: primero escribe un archivo temporal y
    recién al final reemplaza el original. Si se corta la luz a mitad del
    guardado, el archivo bueno anterior queda intacto.
    """
    camino = Path(ruta)
    camino.parent.mkdir(parents=True, exist_ok=True)
    temporal = camino.with_name(camino.name + ".tmp")
    with open(temporal, "w", encoding="utf-8") as f:
        json.dump(datos, f, indent=indent, ensure_ascii=False)
    os.replace(temporal, camino)
    return camino


def leer_texto(ruta: os.PathLike | str, por_defecto: str = "") -> str:
    camino = Path(ruta)
    if not camino.exists():
        return por_defecto
    try:
        return camino.read_text(encoding="utf-8")
    except OSError:
        return por_defecto


def guardar_texto(ruta: os.PathLike | str, texto: str) -> Path:
    camino = Path(ruta)
    camino.parent.mkdir(parents=True, exist_ok=True)
    temporal = camino.with_name(camino.name + ".tmp")
    temporal.write_text(texto, encoding="utf-8")
    os.replace(temporal, camino)
    return camino


def listar(ruta: os.PathLike | str, patron: str = "*") -> list[Path]:
    """Lista archivos dentro de una carpeta (lista vacía si no existe)."""
    camino = Path(ruta)
    if not camino.exists():
        return []
    return sorted(p for p in camino.glob(patron) if p.is_file())


# ---------------------------------------------------------------------------
# Fechas y huellas (para saber si un resultado quedó viejo)
# ---------------------------------------------------------------------------
def fecha_modificacion(ruta: os.PathLike | str) -> float | None:
    """Momento de la última modificación (segundos). None si no existe."""
    camino = Path(ruta)
    try:
        return camino.stat().st_mtime
    except OSError:
        return None


def fecha_mas_reciente(rutas) -> float | None:
    """La fecha más nueva entre varios archivos o carpetas."""
    fechas: list[float] = []
    for ruta in rutas:
        camino = Path(ruta)
        if camino.is_dir():
            fechas.extend(f.stat().st_mtime for f in camino.rglob("*") if f.is_file())
        elif camino.exists():
            fechas.append(camino.stat().st_mtime)
    return max(fechas) if fechas else None


def huella(ruta: os.PathLike | str) -> str:
    """
    Huella (hash) del contenido de un archivo o carpeta, en 12 caracteres.
    Sirve para detectar si un resultado cambió contra la referencia guardada.
    """
    camino = Path(ruta)
    h = hashlib.sha256()
    if camino.is_file():
        h.update(camino.read_bytes())
    elif camino.is_dir():
        for f in sorted(p for p in camino.rglob("*") if p.is_file()):
            h.update(str(f.relative_to(camino)).encode("utf-8"))
            h.update(f.read_bytes())
    return h.hexdigest()[:12]


# ---------------------------------------------------------------------------
# Estructura (datos/estructura.json): el archivo central del proyecto
# ---------------------------------------------------------------------------
def cargar_estructura() -> dict:
    """Lee datos/estructura.json (diccionario vacío si no existe)."""
    return leer_json(ESTRUCTURA, por_defecto={}) or {}


def guardar_estructura(estructura: dict) -> Path:
    return guardar_json(ESTRUCTURA, estructura)


def listar_porticos() -> list[str]:
    """Nombres de los pórticos cargados ('Portico 1', 'Portico 2', ...)."""
    return list(cargar_estructura().keys())


def proximo_nombre_portico() -> str:
    """
    Próximo nombre libre de pórtico. Respeta el criterio viejo de P00
    (cantidad de pórticos + 1) sin pisar nombres ya usados.
    """
    existentes = set(listar_porticos())
    numero = len(existentes) + 1
    while f"Portico {numero}" in existentes:
        numero += 1
    return f"Portico {numero}"

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
import shutil
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Raíz del proyecto: ...\Calculador
# (rutas.py está en ...\Calculador\calc\rutas.py  ->  subir dos niveles)
# ---------------------------------------------------------------------------
RAIZ: Path = Path(__file__).resolve().parent.parent
DATOS_GLOBAL: Path = RAIZ / "datos"
OBRAS: Path = RAIZ / "Obras"
CONFIG_OBRA_ACTIVA: Path = DATOS_GLOBAL / "obra_activa.json"
CARPETA_OBRA: Path | None = None
DATOS: Path = DATOS_GLOBAL
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

_NOMBRES_CARPETAS_SALIDA = (
    "analisis_cargas",
    "vigas",
    "columnas",
    "bases",
    "dxf",
    "losas",
    "reacciones",
    "solicitaciones",
)


def _carpetas_salida(raiz: Path) -> tuple[Path, ...]:
    return (raiz, *(raiz / nombre for nombre in _NOMBRES_CARPETAS_SALIDA))


CARPETAS_SALIDA = _carpetas_salida(SALIDAS)

# ---------------------------------------------------------------------------
# Archivos de datos (entradas del cálculo)
# ---------------------------------------------------------------------------
ESTRUCTURA: Path = DATOS / "estructura.json"
MATERIALES: Path = DATOS_GLOBAL / "materiales.json"
VIGUETAS: Path = DATOS_GLOBAL / "viguetas.json"
PERFILES_METALICOS: Path = DATOS_GLOBAL / "perfiles_metalicos.json"
COEFICIENTES_KD: Path = DATOS_GLOBAL / "coeficientes_kd.json"
MOMENTS_INPUT: Path = DATOS / "moments_input.json"

# Datos de entrada propios de cada obra, gestionados desde la app.
TERRENO: Path = DATOS / "terreno.json"
CARGAS: Path = DATOS / "cargas.json"
TIPOS_LOSA: Path = DATOS / "tipos_losa.json"
LOSAS: Path = DATOS / "losas.json"

# ---------------------------------------------------------------------------
# Archivos de salida con nombre fijo
# ---------------------------------------------------------------------------
PLANILLA_COLUMNAS: Path = SAL_COLUMNAS / "planilla_columnas.csv"
COMPUTO_LOSAS: Path = SAL_LOSAS / "computo_losas.csv"


def _leer_obra_activa() -> Path | None:
    """Lee la selección persistida sin depender de helpers definidos más abajo."""
    try:
        with CONFIG_OBRA_ACTIVA.open("r", encoding="utf-8") as archivo:
            nombre = str(json.load(archivo).get("obra", "")).strip()
    except (OSError, json.JSONDecodeError, AttributeError):
        return None
    if not nombre:
        return None
    carpeta = (OBRAS / nombre).resolve()
    try:
        carpeta.relative_to(OBRAS.resolve())
    except ValueError:
        return None
    return carpeta if carpeta.is_dir() else None


def _actualizar_rutas_obra(carpeta: Path | None) -> None:
    """Actualiza las rutas de los datos y resultados propios de la obra activa."""
    global CARPETA_OBRA, DATOS, SALIDAS, DIAGRAMAS_INTERACCION
    global SAL_ANALISIS_CARGAS, SAL_VIGAS, SAL_COLUMNAS, SAL_BASES, SAL_DXF
    global SAL_LOSAS, SAL_REACCIONES, SAL_SOLICITACIONES, CARPETAS_SALIDA
    global ESTRUCTURA, CARGAS, TERRENO, TIPOS_LOSA, LOSAS, MOMENTS_INPUT
    global PLANILLA_COLUMNAS, COMPUTO_LOSAS
    CARPETA_OBRA = carpeta
    DATOS = carpeta / "datos" if carpeta else DATOS_GLOBAL
    SALIDAS = carpeta / "salidas" if carpeta else RAIZ / "salidas"
    DIAGRAMAS_INTERACCION = DATOS_GLOBAL / "diagramas_interaccion"
    SAL_ANALISIS_CARGAS = SALIDAS / "analisis_cargas"
    SAL_VIGAS = SALIDAS / "vigas"
    SAL_COLUMNAS = SALIDAS / "columnas"
    SAL_BASES = SALIDAS / "bases"
    SAL_DXF = SALIDAS / "dxf"
    SAL_LOSAS = SALIDAS / "losas"
    SAL_REACCIONES = SALIDAS / "reacciones"
    SAL_SOLICITACIONES = SALIDAS / "solicitaciones"
    CARPETAS_SALIDA = _carpetas_salida(SALIDAS)
    ESTRUCTURA = DATOS / "estructura.json"
    CARGAS = DATOS / "cargas.json"
    TERRENO = DATOS / "terreno.json"
    TIPOS_LOSA = DATOS / "tipos_losa.json"
    LOSAS = DATOS / "losas.json"
    MOMENTS_INPUT = DATOS / "moments_input.json"
    PLANILLA_COLUMNAS = SAL_COLUMNAS / "planilla_columnas.csv"
    COMPUTO_LOSAS = SAL_LOSAS / "computo_losas.csv"


def listar_obras() -> list[Path]:
    OBRAS.mkdir(parents=True, exist_ok=True)
    return sorted((p for p in OBRAS.iterdir() if p.is_dir() and not p.name.startswith("_")),
                  key=lambda p: p.name.casefold())


def crear_obra(nombre: str, copiar_actual: bool = False) -> Path:
    """Crea una obra; la migración conserva los resultados antiguos en su archivo."""
    limpio = nombre_seguro(nombre).strip()
    if not limpio or limpio in (".", ".."):
        raise ValueError("Ingresá un nombre válido para la obra.")
    carpeta = OBRAS / limpio
    if carpeta.exists():
        raise FileExistsError(f"Ya existe la obra {limpio}.")
    (carpeta / "datos").mkdir(parents=True)
    for directorio in _carpetas_salida(carpeta / "salidas"):
        directorio.mkdir(parents=True, exist_ok=True)
    if copiar_actual:
        origen_datos = DATOS
        for archivo in ("cargas.json", "estructura.json", "losas.json", "terreno.json",
                        "tipos_losa.json", "moments_input.json"):
            origen = origen_datos / archivo
            if origen.is_file():
                shutil.copy2(origen, carpeta / "datos" / archivo)
        if SALIDAS.exists():
            shutil.copytree(SALIDAS, carpeta / "archivo_migracion" / "salidas", dirs_exist_ok=True)
    else:
        guardar_json(carpeta / "datos" / "cargas.json", {
            "obra": limpio, "id_proyecto": limpio.lower().replace(" ", "_"),
            "notas": "", "viento": {"activo": 0, "ancho_tributario_m": 0.0},
            "elementos": {}, "aplicaciones": [],
        })
        guardar_json(carpeta / "datos" / "estructura.json", {})
    cargas = leer_json(carpeta / "datos" / "cargas.json", {}) or {}
    ficha = {
        "nombre": str(cargas.get("obra") or limpio),
        "id": str(cargas.get("id_proyecto") or limpio.lower().replace(" ", "_")),
        "ubicacion": cargas.get("viento", {}).get("referencia_cirsoc_102_25", {}).get("ubicacion", {}),
    }
    guardar_json(carpeta / "obra.json", ficha)
    return carpeta


def seleccionar_obra(carpeta: Path | str) -> Path:
    destino = Path(carpeta).resolve()
    try:
        destino.relative_to(OBRAS.resolve())
    except ValueError as exc:
        raise ValueError("La obra debe estar dentro de la carpeta Obras.") from exc
    if not destino.is_dir():
        raise FileNotFoundError(f"No existe la carpeta de obra: {destino}")
    _actualizar_rutas_obra(destino)
    asegurar_directorios()
    guardar_json(CONFIG_OBRA_ACTIVA, {"obra": destino.name})
    return destino


def inicializar_obras() -> Path:
    """Prepara Obras; copia datos fuente y archiva salidas previas la primera vez."""
    OBRAS.mkdir(parents=True, exist_ok=True)
    activa = leer_json(CONFIG_OBRA_ACTIVA, {}) or {}
    seleccion = OBRAS / str(activa.get("obra", "")) if activa.get("obra") else None
    if seleccion and seleccion.is_dir():
        return seleccionar_obra(seleccion)
    existentes = listar_obras()
    if not existentes:
        legado = leer_json(DATOS_GLOBAL / "cargas.json", {}) or {}
        nombre = str(legado.get("obra") or "Obra nueva").replace("_", " ").strip()
        existente = OBRAS / nombre_seguro(nombre)
        if existente.exists():
            existente = OBRAS / (nombre_seguro(nombre) + " (2)")
        carpeta = crear_obra(existente.name, copiar_actual=bool(legado or (DATOS_GLOBAL / "estructura.json").exists()))
        existentes = [carpeta]
    return seleccionar_obra(existentes[0])


def numero_portico(nombre: str, estructura: dict | None = None) -> int:
    """Obtiene el número del nombre o propone el siguiente número acumulativo."""
    usados = set()
    for portico_id, datos in (estructura or {}).items():
        m = re.search(r"(\d+)\s*$", str(portico_id))
        if m:
            usados.add(int(m.group(1)))
        for viga_id in datos.get("vigas", {}):
            m = re.match(r"V\d+-(\d+)$", str(viga_id))
            if m:
                usados.add(int(m.group(1)))
    coincidencia = re.search(r"(\d+)\s*$", str(nombre))
    if coincidencia:
        numero = int(coincidencia.group(1))
        if numero in usados:
            otro = str(numero)
            for portico_id, datos in (estructura or {}).items():
                m_nombre = re.search(r"(\d+)\s*$", str(portico_id))
                coincide_nombre = bool(m_nombre and int(m_nombre.group(1)) == numero)
                coincide_viga = any(
                    re.match(rf"V\d+-{numero}$", str(viga_id))
                    for viga_id in datos.get("vigas", {})
                )
                if portico_id != nombre and (coincide_nombre or coincide_viga):
                    otro = str(portico_id)
                    break
            raise ValueError(f"El número de pórtico {numero} ya está usado por {otro}.")
        return numero
    return max(usados, default=0) + 1


def siguiente_numero_columna(estructura: dict | None = None) -> int:
    """Devuelve el próximo número de columna libre en toda la obra."""
    mayor = 0
    cantidad_columnas_existentes = 0
    for datos in (estructura or {}).values():
        por_nivel: dict[str, set[str]] = {}
        for columna_id in (datos.get("columnas", {}) or {}):
            identificador = str(columna_id)
            coincidencia = re.fullmatch(
                r"C(\d+)-(.+)", identificador, re.IGNORECASE
            )
            if not coincidencia:
                continue
            nivel, sufijo = coincidencia.groups()
            por_nivel.setdefault(nivel, set()).add(sufijo)
            numero = re.fullmatch(r"(?:\d+-)?(\d+)", sufijo)
            if numero:
                mayor = max(mayor, int(numero.group(1)))
        cantidad_columnas_existentes += max(
            (len(columnas) for columnas in por_nivel.values()),
            default=0,
        )
    return max(mayor, cantidad_columnas_existentes) + 1


def resolver_portico(nombre: str, estructura: dict | None = None) -> str:
    """Acepta el nombre completo del pórtico o su número final."""
    estructura = estructura if estructura is not None else cargar_estructura()
    texto = str(nombre).strip()
    if texto in estructura:
        return texto
    m = re.search(r"(\d+)\s*$", texto)
    if not m:
        raise KeyError(f"No existe el pórtico '{texto}'.")
    numero = int(m.group(1))
    candidatos = []
    for portico_id in estructura:
        n = re.search(r"(\d+)\s*$", str(portico_id))
        if n and int(n.group(1)) == numero:
            candidatos.append(str(portico_id))
    if texto.isdigit() and texto in candidatos:
        return texto
    if len(candidatos) == 1:
        return candidatos[0]
    if candidatos:
        raise KeyError(f"El número {numero} coincide con varios pórticos: {', '.join(candidatos)}.")
    raise KeyError(f"No existe el pórtico número {numero}.")


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


_actualizar_rutas_obra(_leer_obra_activa())


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

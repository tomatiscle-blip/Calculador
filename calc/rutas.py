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

import copy
import hashlib
import json
import os
import re
import shutil
import stat
import time
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
SOLICITACIONES_COLUMNAS: Path = SALIDAS / "solicitaciones_columnas.json"
SOLICITACIONES_BASES: Path = SALIDAS / "solicitaciones_bases.json"
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
MAMPOSTERIA: Path = DATOS_GLOBAL / "mamposteria.json"
VIGUETAS: Path = DATOS_GLOBAL / "viguetas.json"
PERFILES_METALICOS: Path = DATOS_GLOBAL / "perfiles_metalicos.json"
COEFICIENTES_KD: Path = DATOS_GLOBAL / "coeficientes_kd.json"
MOMENTS_INPUT: Path = DATOS / "moments_input.json"

# Datos de entrada propios de cada obra, gestionados desde la app.
TERRENO: Path = DATOS / "terreno.json"
CARGAS: Path = DATOS / "cargas.json"
TIPOS_LOSA: Path = DATOS / "tipos_losa.json"
LOSAS: Path = DATOS / "losas.json"
PARAMETROS_MOTOR: Path = DATOS / "parametros_motor.json"

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
    global SAL_LOSAS, SAL_REACCIONES, SAL_SOLICITACIONES, SOLICITACIONES_COLUMNAS, SOLICITACIONES_BASES, CARPETAS_SALIDA
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
    SOLICITACIONES_COLUMNAS = SALIDAS / "solicitaciones_columnas.json"
    SOLICITACIONES_BASES = SALIDAS / "solicitaciones_bases.json"
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
                        "tipos_losa.json", "moments_input.json", "ejes.json", "niveles.json"):
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
        # Ejes y niveles empiezan vacíos para la nueva obra
        guardar_json(carpeta / "datos" / "ejes.json", {"x": [], "y": []})
        guardar_json(carpeta / "datos" / "niveles.json", {"niveles": []})
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
    anterior = CARPETA_OBRA
    _actualizar_rutas_obra(destino)
    if anterior != destino:
        # Cambió la obra: las cachés de estructura, listados y fechas ya no sirven.
        _invalidar_cache_estructura(ESTRUCTURA)
        _cache_listados.clear()
        invalidar_cache_firmas()
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
    for portico_id, datos in (estructura or {}).get("pórticos", {}).items():
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
            for portico_id, datos in (estructura or {}).get("pórticos", {}).items():
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
    global _generacion_datos
    camino = Path(ruta)
    camino.parent.mkdir(parents=True, exist_ok=True)
    temporal = camino.with_name(camino.name + ".tmp")
    with open(temporal, "w", encoding="utf-8") as f:
        json.dump(datos, f, indent=indent, ensure_ascii=False)
    os.replace(temporal, camino)
    _generacion_datos += 1
    _invalidar_cache_estructura(camino)
    invalidar_cache_firmas()
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
    invalidar_cache_firmas()
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
# Caché de los stat() con vida corta: durante un refresco la app consulta las
# mismas rutas muchas veces. TTL chico para no perder cambios externos.
_TTL_FIRMAS = 0.5  # segundos
_cache_firmas: dict[str, tuple[float, int | None, float | None, bool]] = {}


def _stat_cache(ruta: os.PathLike | str) -> tuple[int | None, float | None, bool]:
    """(mtime_ns, mtime, es_carpeta) de una ruta, con caché de un instante.

    El semáforo consulta el estado de las mismas decenas de rutas varias
    veces por refresco. Consultar el disco en cada consulta es lo más caro
    de la app, así que el resultado se guarda un instante y se invalida
    explícitamente cuando el programa escribe algo.
    """
    camino = os.fspath(ruta)
    ahora = time.monotonic()
    guardada = _cache_firmas.get(camino)
    if guardada is not None and ahora - guardada[0] <= _TTL_FIRMAS:
        return guardada[1], guardada[2], guardada[3]
    try:
        info = os.stat(camino)
    except OSError:
        resultado: tuple[int | None, float | None, bool] = (None, None, False)
    else:
        resultado = (info.st_mtime_ns, info.st_mtime, stat.S_ISDIR(info.st_mode))
    _cache_firmas[camino] = (ahora, *resultado)
    return resultado


def invalidar_cache_firmas() -> None:
    """Olvida las fechas cacheadas (tras escribir archivos o correr un script)."""
    _cache_firmas.clear()


def fecha_modificacion(ruta: os.PathLike | str) -> float | None:
    """Momento de la última modificación (segundos). None si no existe."""
    return _stat_cache(ruta)[1]


def fecha_mas_reciente(rutas) -> float | None:
    """La fecha más nueva entre varios archivos o carpetas."""
    fechas: list[float] = []
    for ruta in rutas:
        mtime_ns, mtime, es_carpeta = _stat_cache(ruta)
        if es_carpeta:
            for archivo in _listar_contenido(os.fspath(ruta)):
                marca = fecha_modificacion(archivo)
                if marca is not None:
                    fechas.append(marca)
        elif mtime is not None:
            fechas.append(mtime)
    return max(fechas) if fechas else None


def firma_modificacion(ruta: os.PathLike | str) -> int | None:
    """mtime en nanosegundos (más preciso que los segundos). None si no existe."""
    return _stat_cache(ruta)[0]


# Listado de carpetas cacheado: se vuelve a listar solo si la carpeta cambió
# (si se agregó o borró un archivo, cambia su mtime).
_cache_listados: dict[str, tuple[int | None, tuple[str, ...]]] = {}


def _listar_contenido(carpeta: str) -> tuple[str, ...]:
    """Archivos (recursivos) de una carpeta, cacheados por su mtime."""
    firma = firma_modificacion(carpeta)
    guardada = _cache_listados.get(carpeta)
    if guardada is not None and guardada[0] == firma:
        return guardada[1]
    archivos: list[str] = []
    for raiz, _dirs, nombres in os.walk(carpeta):
        for nombre in nombres:
            archivos.append(os.path.join(raiz, nombre))
    archivos.sort()
    resultado = tuple(archivos)
    _cache_listados[carpeta] = (firma, resultado)
    return resultado


def firma_ruta(ruta: os.PathLike | str) -> int | None:
    """mtime_ns de un archivo; de una carpeta, el más nuevo de su contenido."""
    camino = os.fspath(ruta)
    mtime_ns, _mtime, es_carpeta = _stat_cache(camino)
    if mtime_ns is None:
        return None
    if not es_carpeta:
        return mtime_ns
    firmas = [
        f for f in (firma_modificacion(a) for a in _listar_contenido(camino)) if f is not None
    ]
    return max(firmas) if firmas else mtime_ns


def firma_mas_reciente(rutas) -> int | None:
    """La firma (mtime_ns) más nueva entre varios archivos o carpetas."""
    firmas = [f for f in (firma_ruta(r) for r in rutas) if f is not None]
    return max(firmas) if firmas else None


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
# Caché de la estructura: se relee solo si el archivo cambió (mtime/tamaño).
# Antes se leía 16 veces por refresco de la app.
_cache_estructura: tuple[Path, tuple[int, int] | None, dict] | None = None

# Contador que cambia cada vez que se guarda un JSON del proyecto.
# Lo usan las cachés de otros módulos (el semáforo) para invalidarse.
_generacion_datos: int = 0


def _firma_archivo(camino: Path) -> tuple[int, int] | None:
    """(mtime_ns, tamaño) del archivo; None si no existe."""
    try:
        info = camino.stat()
    except OSError:
        return None
    return (info.st_mtime_ns, info.st_size)


def generacion_datos() -> int:
    """Número de guardados JSON realizados desde que arrancó el programa.

    Sirve para invalidar cachés de otros módulos (por ejemplo, el semáforo
    de etapas) sin tener que adivinar qué archivos cambiaron.
    """
    return _generacion_datos


def _invalidar_cache_estructura(camino: Path) -> None:
    """Si se acaba de guardar la estructura, la caché vieja ya no sirve."""
    global _cache_estructura
    if _cache_estructura is not None and _cache_estructura[0] == Path(camino):
        _cache_estructura = None


def cargar_estructura() -> dict:
    """Lee datos/estructura.json (diccionario vacío si no existe).

    El contenido se cachea y se vuelve a leer solo cuando el archivo cambia
    (o cuando se cambia de obra). Se devuelve una copia para que quien la
    modifique no altere la caché por accidente.
    """
    global _cache_estructura
    camino = Path(ESTRUCTURA)
    firma = _firma_archivo(camino)
    if (
        _cache_estructura is not None
        and _cache_estructura[0] == camino
        and _cache_estructura[1] == firma
    ):
        return copy.deepcopy(_cache_estructura[2])
    datos = leer_json(camino, por_defecto={}) or {}
    _cache_estructura = (camino, firma, datos)
    return copy.deepcopy(datos)


def guardar_estructura(estructura: dict) -> Path:
    global _cache_estructura
    camino = guardar_json(ESTRUCTURA, estructura)
    _cache_estructura = (camino, _firma_archivo(camino), copy.deepcopy(estructura))
    return camino


def listar_porticos() -> list[str]:
    """Nombres de los pórticos cargados ('Pórtico 1', 'Pórtico 2', ...).

    Solo devuelve elementos que son realmente pórticos (tienen vigas y/o
    columnas). Los muros estructurales (tipo == "muro") se excluyen porque
    no se resuelven con el motor Pynite (pórtico plano): se verifican con
    CIRSOC 501 (calc.mamposteria).
    """
    estructura = cargar_estructura()
    return [
        nombre
        for nombre, datos in estructura.items()
        if isinstance(datos, dict)
        and str(datos.get("tipo", "")).lower() != "muro"
        and (datos.get("vigas") or datos.get("columnas"))
    ]


def listar_muros() -> list[str]:
    """Nombres de los muros estructurales cargados (tipo == "muro")."""
    estructura = cargar_estructura()
    return [
        nombre
        for nombre, datos in estructura.items()
        if isinstance(datos, dict)
        and str(datos.get("tipo", "")).lower() == "muro"
        and nombre
        not in (
            "columnas",
            "nudos",
            "bases",
            "pórticos",
            "nodos",
        )
    ]


def proximo_nombre_portico() -> str:
    """
    Próximo nombre libre de pórtico. Respeta el criterio viejo de P00
    (cantidad de pórticos + 1) sin pisar nombres ya usados.

    Usa "Pórtico N" con tilde: es la convención del resto del proyecto
    (la app crea "Pórtico N" y los datos existentes usan esa forma).
    """
    existentes = set(listar_porticos())
    numero = len(existentes) + 1
    while f"Pórtico {numero}" in existentes or f"Portico {numero}" in existentes:
        numero += 1
    return f"Pórtico {numero}"

def partes(estructura: dict) -> tuple[dict, dict, dict, dict]:
    """Separa la estructura en sus cuatro partes:
    (columnas globales, nudos globales, bases globales, pórticos).

    En este proyecto los pórticos se guardan como claves de nivel raíz de
    estructura.json (no bajo una clave "pórticos"), así que la cuarta parte
    son las claves de nivel raíz que parecen pórticos.
    """
    return (
        estructura.get("columnas", {}) or {},
        estructura.get("nudos", {}) or {},
        estructura.get("bases", {}) or {},
        {
            nombre: datos
            for nombre, datos in estructura.items()
            if isinstance(datos, dict)
            and str(datos.get("tipo", "")).lower() != "muro"
            and (datos.get("vigas") or datos.get("columnas"))
        },
    )


def columnas_globales(estructura: dict) -> dict:
    """Catálogo global de columnas definido a nivel de obra."""
    return partes(estructura)[0]


def nudos_globales(estructura: dict) -> dict:
    return partes(estructura)[1]


def bases_globales(estructura: dict) -> dict:
    return partes(estructura)[2]


def porticos(estructura: dict) -> dict:
    return partes(estructura)[3]


def pórtico(estructura: dict, nombre: str):
    """Devuelve el diccionario del pórtico o None si no existe."""
    return porticos(estructura).get(nombre)


def vigas_por_portico(nombre: str, estructura: dict | None = None) -> dict:
    """Vigas de un pórtico, sin importar cómo esté guardado el resto."""
    p = pórtico(estructura, nombre) if estructura is not None else {}
    return p.get("vigas", {}) if isinstance(p, dict) else {}


def columnas_referencia(nombre: str, estructura: dict | None = None) -> list[str]:
    """IDs globales de las columnas que utiliza un pórtico."""
    p = pórtico(estructura, nombre) if estructura is not None else {}
    return p.get("columnas", []) if isinstance(p, dict) else []


def bases_referencia(nombre: str, estructura: dict | None = None) -> list[str]:
    """IDs globales de las columnas que tienen base en un pórtico."""
    p = pórtico(estructura, nombre) if estructura is not None else {}
    return p.get("bases", []) if isinstance(p, dict) else []


def columna_global(id_columna: str, estructura: dict | None = None) -> dict | None:
    """Encuentra una columna por su ID global."""
    return columnas_globales(estructura).get(id_columna)


def base_global(id_columna: str, estructura: dict | None = None) -> dict | None:
    """Encuentra la base de una columna por su ID global."""
    return bases_globales(estructura).get(id_columna)


def eliminar_portico(nombre: str, estructura: dict | None = None) -> bool:
    """Elimina un pórtico del proyecto, su traza analítica y dimensionados,
    pero conserva las columnas/bases globales que todavía usan otros
    elementos de la obra.

    Devuelve True si se eliminó, False si el pórtico no existía.
    """
    if estructura is None:
        estructura = cargar_estructura()
    if not isinstance(estructura, dict):
        return False

    datos = estructura.get(nombre)
    if not isinstance(datos, dict):
        return False
    if str(datos.get("tipo", "")).lower() == "muro":
        return False
    if not (datos.get("vigas") or datos.get("columnas")):
        return False

    del estructura[nombre]

    # Columnas y bases que siguen siendo útiles para el resto de la obra.
    siguen_usadas = set()
    for nombre_p, datos in estructura.items():
        if not isinstance(datos, dict):
            continue
        if str(datos.get("tipo", "")).lower() == "muro":
            continue
        if not (datos.get("vigas") or datos.get("columnas")):
            continue
        siguen_usadas.update(datos.get("columnas") or [])
        siguen_usadas.update(datos.get("bases") or [])

    columnas = estructura.get("columnas", {}) or {}
    bases = estructura.get("bases", {}) or {}
    for gid in list(columnas.keys()):
        if gid not in siguen_usadas:
            del columnas[gid]
    for gid in list(bases.keys()):
        if gid not in siguen_usadas:
            del bases[gid]

    estructura["columnas"] = columnas
    estructura["nudos"] = estructura.get("nudos", {}) or {}
    estructura["bases"] = bases

    guardar_estructura(estructura)
    _borrar_salidas_portico(nombre)
    return True


def _borrar_salidas_portico(nombre: str) -> None:
    """Borra los archivos de cálculo de un pórtico."""
    clave = nombre_seguro(nombre)
    for carpeta in (SAL_SOLICITACIONES, SAL_VIGAS, SAL_COLUMNAS, SAL_BASES):
        patron = f"{clave}.json"
        for archivo in listar(carpeta, patron):
            try:
                archivo.unlink()
            except OSError:
                pass


# ---------------------------------------------------------------------------
# Principio 4 — Acumulación de solicitudes biaxial por columna global
# ---------------------------------------------------------------------------

def _clave_columna_global(datos_col: dict) -> tuple | None:
    """Clave única (x, nivel) para identificar una columna física."""
    try:
        x = round(float(datos_col["x"]), 4)
        nivel = round(float(datos_col.get("nivel", datos_col.get("nivel_inicio", 0.0))), 4)
        return (x, nivel)
    except (KeyError, TypeError, ValueError):
        return None


def _catalogo_columnas_globales(estructura: dict) -> dict[str, dict]:
    """Construye un catálogo de columnas globales a partir de la estructura."""
    columnas_raiz = estructura.get("columnas", {}) or {}
    if isinstance(columnas_raiz, dict) and columnas_raiz:
        return {gid: dict(d) for gid, d in columnas_raiz.items()}
    catalogo: dict[str, dict] = {}
    contador = 1
    for nombre, datos in estructura.items():
        if not isinstance(datos, dict):
            continue
        if str(datos.get("tipo", "")).lower() == "muro":
            continue
        if not (datos.get("vigas") or datos.get("columnas")):
            continue
        columnas_p = datos.get("columnas", {}) or {}
        if isinstance(columnas_p, list):
            continue
        for cid, cd in columnas_p.items():
            if not isinstance(cd, dict):
                continue
            clave = _clave_columna_global(cd)
            if clave is None:
                continue
            existente = next(
                (gid for gid, d in catalogo.items()
                 if _clave_columna_global(d) == clave),
                None,
            )
            if existente is None:
                gid = f"COL-{contador:03d}"
                contador += 1
                catalogo[gid] = {
                    "x": cd.get("x"),
                    "nivel": cd.get("nivel", cd.get("nivel_inicio", 0.0)),
                    "altura_m": cd.get("altura_m"),
                    "eje_id": cd.get("eje_id"),
                    "eje_nombre": cd.get("eje_nombre"),
                }
    return catalogo


def siguiente_id_columna_global(estructura: dict | None = None) -> str:
    """Devuelve el próximo ID de columna global libre (COL-001, COL-002, ...)."""
    estructura = estructura if estructura is not None else cargar_estructura()
    columnas = estructura.get("columnas", {}) or {}
    mayor = 0
    for gid in columnas:
        m = re.fullmatch(r"COL-(\d+)", str(gid))
        if m:
            mayor = max(mayor, int(m.group(1)))
    return f"COL-{mayor + 1:03d}"


def obtener_o_crear_columna_global(
    x: float,
    nivel: float,
    altura_m: float,
    eje_id: str,
    eje_nombre: str,
    nivel_inicio: str,
    estructura: dict | None = None,
) -> str:
    """Devuelve el ID global de una columna por coordenadas (x, nivel).

    Si ya existe una columna global con esas coordenadas en
    estructura["columnas"], devuelve su ID. Si no existe, crea una nueva
    con el próximo ID COL-XXX disponible y la registra en el catálogo
    global de la obra.

    Este es el paso clave del principio 4: las columnas son entidades
    globales de la obra, compartidas por todos los pórticos que las usan.
    """
    if estructura is None:
        estructura = cargar_estructura()

    columnas = estructura.get("columnas")
    if not isinstance(columnas, dict):
        columnas = {}
        estructura["columnas"] = columnas

    clave_buscar = (round(float(x), 4), round(float(nivel), 4))

    # 1. Buscar si ya existe una columna global con esas coordenadas
    for gid, datos in columnas.items():
        if not isinstance(datos, dict):
            continue
        clave_existente = _clave_columna_global(datos)
        if clave_existente == clave_buscar:
            return gid

    # 2. Crear una nueva columna global
    gid = siguiente_id_columna_global(estructura)
    columnas[gid] = {
        "x": float(x),
        "nivel": float(nivel),
        "altura_m": float(altura_m),
        "eje_id": eje_id,
        "eje_nombre": eje_nombre,
        "nivel_inicio": nivel_inicio,
    }
    return gid


def siguiente_id_base_global(estructura: dict | None = None) -> str:
    """Devuelve el próximo ID de base global libre (BASE-001, BASE-002, ...)."""
    estructura = estructura if estructura is not None else cargar_estructura()
    bases = estructura.get("bases", {}) or {}
    mayor = 0
    for gid in bases:
        m = re.fullmatch(r"BASE-(\d+)", str(gid))
        if m:
            mayor = max(mayor, int(m.group(1)))
    return f"BASE-{mayor + 1:03d}"


def obtener_o_crear_base_global(
    x: float,
    nivel: float,
    eje_id: str,
    eje_nombre: str,
    estructura: dict | None = None,
) -> str:
    """Devuelve el ID global de una base por coordenadas (x, nivel).

    Si ya existe una base global con esas coordenadas en
    estructura["bases"], devuelve su ID. Si no existe, crea una nueva
    con el próximo ID BASE-XXX disponible y la registra en el catálogo
    global de la obra.
    """
    if estructura is None:
        estructura = cargar_estructura()

    bases = estructura.get("bases")
    if not isinstance(bases, dict):
        bases = {}
        estructura["bases"] = bases

    clave_buscar = (round(float(x), 4), round(float(nivel), 4))

    for gid, datos in bases.items():
        if not isinstance(datos, dict):
            continue
        clave_existente = (round(float(datos.get("x", 0.0)), 4),
                           round(float(datos.get("nivel", 0.0)), 4))
        if clave_existente == clave_buscar:
            return gid

    gid = siguiente_id_base_global(estructura)
    bases[gid] = {
        "x": float(x),
        "nivel": float(nivel),
        "tipo": "empotramiento",
        "eje_id": eje_id,
        "eje_nombre": eje_nombre,
    }
    return gid


def _mapear_cids_a_globales(nombre_portico: str, estructura: dict,
                            catalogo: dict[str, dict]) -> dict[str, str]:
    """Devuelve {cid_local: id_global} para un pórtico."""
    p = pórtico(estructura, nombre_portico) or {}
    columnas_p = p.get("columnas", {}) or {}
    if isinstance(columnas_p, list):
        return {gid: gid for gid in columnas_p}
    mapeo: dict[str, str] = {}
    for cid, cd in columnas_p.items():
        if not isinstance(cd, dict):
            continue
        clave = _clave_columna_global(cd)
        if clave is None:
            continue
        for gid, d in catalogo.items():
            if _clave_columna_global(d) == clave:
                mapeo[cid] = gid
                break
    return mapeo


def acumular_solicitaciones_columnas(estructura: dict | None = None) -> dict[str, str]:
    """Lee todos los JSON de salidas/solicitaciones y acumula por columna global.

    Post-procesamiento del principio 4 de ARQUITECTURA.md: vincula las
    solicitudes de cada pórtico al ID global de la columna para que los
    dimensionadores evalúen el efecto bidireccional (flexión esviada).

    Devuelve {id_global: ruta_archivo_guardado}.
    """
    if estructura is None:
        estructura = cargar_estructura()

    catalogo = _catalogo_columnas_globales(estructura)
    if not catalogo:
        return {}

    SAL_COLUMNAS.mkdir(parents=True, exist_ok=True)

    resultados: dict[str, str] = {}
    archivos = listar(SAL_SOLICITACIONES, "*.json")
    for archivo in archivos:
        datos = leer_json(archivo, {}) or {}
        nombre_portico = datos.get("portico") or archivo.stem
        mapeo = _mapear_cids_a_globales(nombre_portico, estructura, catalogo)
        if not mapeo:
            continue

        envolvente = datos.get("envolvente", {}).get("columnas", {})
        solicitudes = datos.get("solicitaciones", {})

        for cid, valores in envolvente.items():
            if not valores:
                continue
            gid = mapeo.get(cid)
            if gid is None:
                continue

            entrada_por_combo: dict[str, dict] = {}
            for combo, bloque in solicitudes.items():
                col_bloque = bloque.get("columnas", {}).get(cid, {})
                if col_bloque:
                    entrada_por_combo[combo] = {
                        "M_inf": col_bloque.get("M_inf", 0.0),
                        "M_sup": col_bloque.get("M_sup", 0.0),
                        "N": col_bloque.get("N", 0.0),
                        "deriva": col_bloque.get("deriva", 0.0),
                    }

            m_inf = valores.get("M_inf", 0.0)
            m_sup = valores.get("M_sup", 0.0)
            envol_entry = {
                "M_inf": m_inf,
                "M_sup": m_sup,
                "N": valores.get("N", 0.0),
                "deriva": valores.get("deriva", 0.0),
                "|M|_máx": max(abs(m_inf), abs(m_sup)),
                "|N|_máx": abs(valores.get("N", 0.0)),
            }

            salida = {
                "columna_global": gid,
                "pórtico": nombre_portico,
                "combinaciones": entrada_por_combo,
                "envolvente_pórtico": envol_entry,
            }

            archivo_salida = SAL_COLUMNAS / f"{gid}.json"
            previo = leer_json(archivo_salida, {}) or {}
            porticos = previo.get("pórticos", {})
            porticos[nombre_portico] = salida

            env_global: dict[str, float] = {}
            for p_datos in porticos.values():
                env_p = p_datos.get("envolvente_pórtico", {})
                for clave2 in ("M_inf", "M_sup", "N", "deriva", "|M|_máx", "|N|_máx"):
                    val = env_p.get(clave2, 0.0)
                    if clave2 in ("M_inf", "M_sup"):
                        if clave2 not in env_global or abs(val) > abs(env_global[clave2]):
                            env_global[clave2] = val
                    else:
                        env_global[clave2] = max(env_global.get(clave2, 0.0), abs(val))

            salida_global = {
                "columna_global": gid,
                "catalogo": catalogo.get(gid, {}),
                "pórticos": porticos,
                "envolvente_global": env_global,
                "n_pórticos": len(porticos),
            }
            guardar_json(archivo_salida, salida_global)
            resultados[gid] = str(archivo_salida)

    return resultados



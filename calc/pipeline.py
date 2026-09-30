"""
calc.pipeline — Orden de las etapas, dependencias y ESTADO (semáforo).

Este módulo es el "plano de obra" del cálculo: dice qué etapa va antes de cuál
y con qué archivos. Con eso la app puede mostrar un semáforo por etapa:

    OK              -> la salida existe y es más nueva que sus entradas
    DESACTUALIZADA  -> cambió una entrada después de calcularse la salida
    PENDIENTE       -> nunca se calculó (falta la salida)
    SIN DATOS       -> faltan las entradas (hay que cargar datos primero)
    A DESARROLLAR   -> etapa planificada, todavía no programada

Nada de esto calcula: solo mira archivos y contenido. Es lo que permite que la
app avise "ojo, cambiaste la geometría: hay que recalcular vigas y columnas".
"""

from __future__ import annotations

import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from . import rutas

# ---------------------------------------------------------------------------
# Definición de una etapa
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Etapa:
    clave: str                      # nombre corto, sin espacios (ej. "vigas")
    nombre: str                     # nombre para mostrar en pantalla
    descripcion: str                # qué hace, en castellano
    script: str | None              # script de consola que hoy la ejecuta
    entradas: tuple = ()            # de qué archivos depende
    salidas: tuple = ()             # qué archivos produce ("{portico}" se reemplaza)
    depende_de: tuple = ()          # claves de otras etapas
    interactiva: bool = False       # todavía pide datos por teclado (input)
    implementada: bool = True       # False = planificada, sin programar
    marca: Callable[..., tuple[bool, str]] | None = None  # chequeo por contenido
    nota: str = ""


# ---------------------------------------------------------------------------
# Marcas: chequeos por CONTENIDO (más confiables que la fecha del archivo)
# ---------------------------------------------------------------------------
def _marca_cargas(_portico: str = "") -> tuple[bool, str]:
    txts = rutas.listar(rutas.SAL_ANALISIS_CARGAS, "*.txt")
    if txts:
        return True, f"{len(txts)} análisis en salidas/analisis_cargas"
    return False, "sin análisis de cargas guardados"


def _marca_geometria(portico: str = "") -> tuple[bool, str]:
    est = rutas.cargar_estructura()
    if not est:
        return False, "datos/estructura.json no existe o está vacío"
    nombre = portico or sorted(est)[0]
    p = est.get(nombre)
    if not p:
        return False, f"el pórtico '{nombre}' no está en estructura.json"
    n_col = len(p.get("columnas", {}))
    n_tra = sum(len(v.get("tramos", [])) for v in p.get("vigas", {}).values())
    return True, f"{len(est)} pórtico(s); {nombre}: {n_col} columnas, {n_tra} tramos"


def _marca_portico(portico: str = "") -> tuple[bool, str]:
    est = rutas.cargar_estructura()
    if not est or not portico or portico not in est:
        return False, "sin geometría para analizar"
    cols = est[portico].get("columnas", {})
    con_momentos = [c for c in cols.values() if "Mu_kNm_inf" in c]
    if cols and len(con_momentos) == len(cols):
        return True, f"{len(cols)} columnas con momentos calculados"
    return False, f"{len(con_momentos)}/{len(cols)} columnas calculadas"


def _marca_vigas(portico: str = "") -> tuple[bool, str]:
    if not portico:
        return False, "sin pórtico seleccionado"
    archivo = rutas.SAL_VIGAS / f"resultados_{portico}_vigas.json"
    datos = rutas.leer_json(archivo)
    if not datos:
        return False, f"falta {archivo.name}"
    n = len(datos.get("tramos", []))
    return True, f"{n} tramo{'s' if n != 1 else ''} en {archivo.name}"


def _marca_columnas(portico: str = "") -> tuple[bool, str]:
    """
    Lee salidas/columnas/planilla_columnas.csv (separado por ';').
    Cuenta las columnas del pórtico elegido y avisa si alguna quedó con el
    estribo fuera de norma (columna 'cumple_estribo' en False).
    """
    texto = rutas.leer_texto(rutas.PLANILLA_COLUMNAS)
    if not texto:
        return False, "falta salidas/columnas/planilla_columnas.csv"
    lineas = [l for l in texto.splitlines()[1:] if l.strip()]
    propias = [l for l in lineas if l.split(";")[0].strip() == portico]
    if not propias:
        return False, f"el pórtico {portico} todavía no tiene columnas calculadas"
    revisar = [l for l in propias if len(l.split(";")) > 17 and l.split(";")[17].strip() == "False"]
    detalle = f"{len(propias)} columnas de {portico}"
    if revisar:
        detalle += f" · OJO: {len(revisar)} con estribo fuera de norma"
    memoria = rutas.SAL_COLUMNAS / f"memoria_{rutas.nombre_seguro(portico)}.txt"
    if memoria.exists():
        detalle += " + memoria de cálculo"
    return True, detalle


def _marca_bases(portico: str = "") -> tuple[bool, str]:
    archivo = rutas.SAL_BASES / f"bases_{rutas.nombre_seguro(portico)}.json"
    datos = rutas.leer_json(archivo)
    if not datos:
        return False, f"falta salidas/bases/bases_{portico}.json"
    bases = [k for k in datos if k != "vigas_fundacion"]
    return True, f"{len(bases)} bases calculadas"


def _marca_planos(portico: str = "") -> tuple[bool, str]:
    archivo = rutas.SAL_DXF / f"{rutas.nombre_seguro(portico)}_lateral.dxf"
    if not archivo.exists():
        return False, f"falta salidas/dxf/{portico}_lateral.dxf"
    return True, f"{archivo.name} ({archivo.stat().st_size // 1024} kB)"


def _marca_losas(_portico: str = "") -> tuple[bool, str]:
    memorias = rutas.listar(rutas.SAL_LOSAS, "memoria_losa_*.txt")
    if memorias:
        return True, f"{len(memorias)} memorias de losa"
    return False, "sin memorias de losa"


def _marca_terreno(_portico: str = "") -> tuple[bool, str]:
    if rutas.TERRENO.exists():
        return True, "datos/terreno.json cargado"
    return False, "falta datos/terreno.json (a desarrollar)"


# ---------------------------------------------------------------------------
# LAS ETAPAS, EN ORDEN DE EJECUCIÓN
# ---------------------------------------------------------------------------

ETAPAS: tuple[Etapa, ...] = (
    Etapa(
        clave="cargas",
        nombre="1. Análisis de cargas",
        descripcion="Peso propio, sobrecargas, viento y combinaciones CIRSOC.",
        script="00_Analisis_cargas.py",
        entradas=(),
        salidas=(rutas.SAL_ANALISIS_CARGAS,),
        marca=_marca_cargas,
        nota="Hoy la configuración está dentro del script; pasará a datos/cargas.json.",
    ),
    Etapa(
        clave="terreno",
        nombre="2. Datos del terreno",
        descripcion="Capas, nivel freático, q_adm y módulo de balasto.",
        script=None,
        entradas=(),
        salidas=(rutas.TERRENO,),
        implementada=False,
        marca=_marca_terreno,
        nota="A desarrollar: hoy q_adm y profundidad se piden a mano en P05.",
    ),
    Etapa(
        clave="geometria",
        nombre="3. Geometría del pórtico",
        descripcion="Columnas, tramos, voladizos y cargas puntuales.",
        script="P00_Ingresar_datos_estructura.py",
        entradas=(rutas.SAL_ANALISIS_CARGAS,),
        salidas=(rutas.ESTRUCTURA,),
        depende_de=("cargas",),
        interactiva=True,
        marca=_marca_geometria,
        nota="Todavía se carga por teclado; después pasa a tabla en pantalla.",
    ),
    Etapa(
        clave="portico",
        nombre="4. Cálculo del pórtico",
        descripcion="Solicitaciones en vigas, columnas y bases (motor Pynite).",
        script="P01_dimensionado_portico_hormigon.py",
        entradas=(rutas.ESTRUCTURA,),
        salidas=(rutas.ESTRUCTURA,),
        depende_de=("geometria",),
        interactiva=True,
        marca=_marca_portico,
        nota="A migrar de anaStruct 2D a Pynite 3D, con comparación previa.",
    ),
    Etapa(
        clave="vigas",
        nombre="5. Vigas de hormigón",
        descripcion="Flexión, corte, flecha, fisuración y planilla de armado.",
        script="P02_Viga_portico.py",
        entradas=(rutas.ESTRUCTURA, rutas.COEFICIENTES_KD),
        salidas=(rutas.SAL_VIGAS / "resultados_{portico}_vigas.json",),
        depende_de=("portico",),
        interactiva=True,
        marca=_marca_vigas,
    ),
    Etapa(
        clave="vigas_excel",
        nombre="6. Planilla de vigas (Excel)",
        descripcion="Junta los resultados de vigas en un Excel por pórtico.",
        script="P03_Viga_portico_guardar_excel.py",
        entradas=(rutas.ruta_salida("vigas", "resultados_*_vigas.json"),),
        salidas=(rutas.SAL_VIGAS / "Planilla_Vigas_Portico.xlsx",),
        depende_de=("vigas",),
    ),
    Etapa(
        clave="columnas",
        nombre="7. Columnas",
        descripcion="Esbeltez, cuantías y diagramas de interacción.",
        script="P04_Columnas_portico.py",
        entradas=(rutas.ESTRUCTURA, rutas.DIAGRAMAS_INTERACCION),
        salidas=(
            rutas.PLANILLA_COLUMNAS,
            rutas.SAL_COLUMNAS / "memoria_{portico}.txt",
        ),
        depende_de=("portico",),
        interactiva=True,
        marca=_marca_columnas,
    ),
    Etapa(
        clave="bases",
        nombre="8. Bases (zapatas)",
        descripcion="Dimensionado, tensiones del suelo, punzonado y armadura.",
        script="P05_Bases_portico.py",
        entradas=(rutas.ESTRUCTURA, rutas.PLANILLA_COLUMNAS),
        salidas=(
            rutas.SAL_BASES / "bases_{portico}.json",
            rutas.SAL_BASES / "bases_completas_{portico}.txt",
        ),
        depende_de=("columnas",),
        interactiva=True,
        marca=_marca_bases,
        nota="A desarrollar: tomar q_adm y profundidad desde datos/terreno.json.",
    ),
    Etapa(
        clave="planos",
        nombre="9. Plano lateral (DXF)",
        descripcion="Plano de armaduras de vigas, columnas y bases.",
        script="P06_Portico_dxf.py",
        entradas=(
            rutas.ESTRUCTURA,
            rutas.PLANILLA_COLUMNAS,
            rutas.ruta_salida("vigas", "resultados_{portico}_vigas.json"),
            rutas.ruta_salida("vigas", "planilla_{portico}_*"),
            rutas.ruta_salida("bases", "bases_{portico}.json"),
        ),
        salidas=(rutas.SAL_DXF / "{portico}_lateral.dxf",),
        depende_de=("vigas", "columnas", "bases"),
        marca=_marca_planos,
        nota="Hoy el pórtico está fijo dentro del script (PORTICO = 'Portico 3').",
    ),
    Etapa(
        clave="losas",
        nombre="10. Losas",
        descripcion="Losas alivianadas, macizas y casetonadas (por tipología).",
        script="L00_Losas_alivianadas.py",
        entradas=(rutas.MATERIALES, rutas.VIGUETAS),
        salidas=(rutas.COMPUTO_LOSAS, rutas.SAL_LOSAS),
        interactiva=True,
        marca=_marca_losas,
        nota="A desarrollar: separar por tipología (vigueta / maciza / casetonada).",
    ),
)

ETAPAS_POR_CLAVE: dict[str, Etapa] = {e.clave: e for e in ETAPAS}

# ---------------------------------------------------------------------------
# Estado de cada etapa (el semáforo)
# ---------------------------------------------------------------------------

ICONOS = {
    "ok": "[ OK ]",
    "desactualizada": "[OJO!]",
    "pendiente": "[  -  ]",
    "sin_datos": "[FALTA]",
    "a_desarrollar": "[TODO ]",
}


def _marca_excel(_portico: str = "") -> tuple[bool, str]:
    """El Excel de vigas está al día si es más nuevo que los JSON que lo alimentan."""
    excel = rutas.SAL_VIGAS / "Planilla_Vigas_Portico.xlsx"
    if not excel.exists():
        return False, "falta salidas/vigas/Planilla_Vigas_Portico.xlsx"
    fuentes = rutas.listar(rutas.SAL_VIGAS, "resultados_*_vigas.json")
    if not fuentes:
        return False, "no hay resultados de vigas todavía"
    f_excel = rutas.fecha_modificacion(excel) or 0.0
    f_fuente = rutas.fecha_mas_reciente(fuentes) or 0.0
    if f_fuente > f_excel + 1:
        return False, "hay resultados de vigas más nuevos que la planilla"
    return True, f"{excel.name} con {len(fuentes)} pórtico(s)"


# Marcas que no se pueden declarar dentro de ETAPAS (se definen después)
MARCAS_ESPECIALES: dict[str, Callable[..., tuple[bool, str]]] = {
    "vigas_excel": _marca_excel,
}


def _resolver(coleccion, portico: str) -> tuple[Path, ...]:
    """
    Reemplaza {portico} en las rutas y expande los comodines (*).
    Si falta el pórtico, descarta las rutas que lo necesitan.
    """
    resueltas: list[Path] = []
    for elemento in coleccion:
        texto = str(elemento)
        if "{portico}" in texto:
            if not portico:
                continue
            texto = texto.replace("{portico}", rutas.nombre_seguro(portico))
        ruta = Path(texto)
        if "*" in ruta.name or "?" in ruta.name:
            if ruta.parent.exists():
                coincidencias = sorted(ruta.parent.glob(ruta.name))
                if coincidencias:
                    resueltas.extend(coincidencias)
                    continue
            resueltas.append(ruta)  # comodín sin resultados -> cuenta como faltante
        else:
            resueltas.append(ruta)
    return tuple(resueltas)


def portico_por_defecto(portico: str = "") -> str:
    """Si no se indica pórtico, toma el último cargado (el más nuevo)."""
    if portico:
        return portico
    disponibles = rutas.listar_porticos()
    return disponibles[-1] if disponibles else ""


def estado_etapa(clave: str, portico: str = "") -> dict:
    """Estado de una etapa: {'estado': 'ok'|'pendiente'|..., 'detalle': ...}."""
    etapa = ETAPAS_POR_CLAVE[clave]
    portico = portico_por_defecto(portico)
    entradas = _resolver(etapa.entradas, portico)
    salidas = _resolver(etapa.salidas, portico)

    dato = {
        "clave": etapa.clave,
        "nombre": etapa.nombre,
        "descripcion": etapa.descripcion,
        "script": etapa.script,
        "depende_de": list(etapa.depende_de),
        "interactiva": etapa.interactiva,
        "nota": etapa.nota,
        "portico": portico,
        "entradas": [str(p) for p in entradas],
        "salidas": [str(p) for p in salidas],
    }

    if not etapa.implementada:
        return {**dato, "estado": "a_desarrollar", "detalle": etapa.nota or "planificada"}

    faltan = [Path(p).name for p in entradas if not Path(p).exists()]
    if faltan:
        return {**dato, "estado": "sin_datos", "detalle": "falta: " + ", ".join(faltan)}

    detalle = ""
    marca = etapa.marca if etapa.marca is not None else MARCAS_ESPECIALES.get(clave)
    if marca is not None:
        try:
            listo, detalle = marca(portico)
        except Exception as exc:  # una marca NUNCA debe tumbar la app
            listo, detalle = False, f"no se pudo verificar ({exc})"
        if not listo:
            return {**dato, "estado": "pendiente", "detalle": detalle}

    f_entrada = rutas.fecha_mas_reciente(entradas)
    f_salida = rutas.fecha_mas_reciente(salidas)
    if f_salida is None:
        return {**dato, "estado": "pendiente", "detalle": detalle or "todavía no se calculó"}
    if f_entrada and f_entrada > f_salida + 1:
        aviso = "hay datos más nuevos que este resultado: conviene recalcular"
        return {**dato, "estado": "desactualizada", "detalle": f"{detalle} · {aviso}".strip(" ·")}
    return {**dato, "estado": "ok", "detalle": detalle}


def semaforo(portico: str = "") -> list[dict]:
    """Estado de todas las etapas, en orden de ejecución."""
    return [estado_etapa(etapa.clave, portico) for etapa in ETAPAS]


def resumen_texto(portico: str = "") -> str:
    """Tabla de texto del semáforo (para consola o para copiar a la memoria)."""
    portico = portico_por_defecto(portico)
    lineas = [
        f"ESTADO DEL PROYECTO — Pórtico: {portico or '(ninguno cargado)'}",
        "=" * 78,
    ]
    for dato in semaforo(portico):
        icono = ICONOS.get(dato["estado"], "[ ??? ]")
        lineas.append(f"{icono} {dato['nombre']:<34} {dato['detalle']}")
    return "\n".join(lineas)


def pendientes(portico: str = "") -> list[str]:
    """Claves de las etapas que hay que resolver (sin datos, viejas o sin calcular)."""
    a_hacer = []
    for dato in semaforo(portico):
        if dato["estado"] in ("pendiente", "desactualizada", "sin_datos"):
            a_hacer.append(dato["clave"])
    return a_hacer


# ---------------------------------------------------------------------------
# Ejecutar una etapa (puente hacia los scripts actuales)
# ---------------------------------------------------------------------------
def ejecutar(
    clave: str,
    portico: str = "",
    respuestas: str | None = None,
    timeout: int = 1800,
) -> dict:
    """
    Corre el programa de la etapa y devuelve {'ok', 'codigo', 'salida', 'error'}.

    `respuestas` es el texto que se le entrega por teclado (una respuesta por
    línea), para las etapas que todavía piden datos con input(). Ejemplo:
        respuestas="0\\n1\\n4.85\\n3.0\\n"

    Se ejecuta con la carpeta del proyecto como directorio de trabajo, por eso
    las rutas relativas de los scripts viejos siguen funcionando.
    """
    etapa = ETAPAS_POR_CLAVE.get(clave)
    if etapa is None:
        return {"ok": False, "codigo": None, "salida": "", "error": f"etapa desconocida: {clave}"}
    if not etapa.script:
        return {
            "ok": False,
            "codigo": None,
            "salida": "",
            "error": f"{etapa.nombre}: todavía no está programada",
        }

    script = rutas.RAIZ / etapa.script
    if not script.exists():
        return {"ok": False, "codigo": None, "salida": "", "error": f"no se encontró {script.name}"}

    comando = [sys.executable, str(script)]
    if portico:
        comando.append(portico)  # los scripts migrados podrán tomarlo como argumento

    # Los scripts viejos imprimen símbolos (γ, ·, →). Cuando su salida se captura,
    # Windows usa cp1252 y el programa se cae con UnicodeEncodeError; por eso se les
    # fuerza UTF-8. Es lo que permite ejecutar una etapa desde la ventana.
    entorno = {**os.environ, "PYTHONIOENCODING": "utf-8"}

    try:
        proceso = subprocess.run(
            comando,
            cwd=str(rutas.RAIZ),
            env=entorno,
            input=respuestas,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return {
            "ok": False,
            "codigo": None,
            "salida": "",
            "error": f"{etapa.nombre}: tardó más de {timeout} segundos",
        }
    except OSError as exc:
        return {
            "ok": False,
            "codigo": None,
            "salida": "",
            "error": f"no se pudo lanzar {script.name}: {exc}",
        }

    return {
        "ok": proceso.returncode == 0,
        "codigo": proceso.returncode,
        "salida": proceso.stdout,
        "error": proceso.stderr,
    }


def _main(argv: list[str]) -> int:
    portico = argv[0] if argv else ""
    print(resumen_texto(portico))
    a_hacer = pendientes(portico)
    if a_hacer:
        print("\nPara resolver: " + ", ".join(a_hacer))
    return 0


if __name__ == "__main__":
    raise SystemExit(_main(sys.argv[1:]))


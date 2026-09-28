"""
tools.ver_cambios — Muestra, en castellano, QUÉ archivos se agregaron o cambiaron.

Para qué sirve
--------------
No hace falta saber git. Este programa le pregunta a git qué cambió respecto
del último guardado y te lo explica en una lista, con una descripción de para
qué sirve cada archivo nuevo.

Uso (desde la carpeta del proyecto):

    py tools\\ver_cambios.py                  -> informe en pantalla
    py tools\\ver_cambios.py --abrir          -> además abre los archivos nuevos en VS Code
    py tools\\ver_cambios.py --diff P06_Portico_dxf.py
                                             -> abre VS Code comparando antes/después

También hay un `ver_cambios.bat` para hacer doble clic.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent

# ---------------------------------------------------------------------------
# Para qué sirve cada archivo nuevo (en castellano, sin vueltas)
# ---------------------------------------------------------------------------
DESCRIPCIONES = {
    ".gitignore": "Le dice a git qué carpetas NO guardar (.venv, __pycache__, .vs, copias de prueba).",
    "requirements.txt": "Lista de librerías del proyecto con su versión. Permite reinstalar todo después de actualizar Windows o Python.",
    "PUESTA_EN_MARCHA.md": "Documento de puesta en marcha: decisiones (Pynite / PySide6), cómo reinstalar, cómo leer el semáforo y qué sigue.",
    "estado.bat": "Doble clic: muestra el semáforo del proyecto (qué está calculado, qué falta, qué quedó viejo).",
    "ver_cambios.bat": "Doble clic: muestra este mismo informe.",
    "calc/": "Carpeta nueva: el núcleo de cálculo (por ahora son los cimientos, sin fórmulas todavía).",
    "calc/__init__.py": "Presentación del paquete calc: deja escrita la regla de oro (el cálculo no conoce la pantalla).",
    "calc/rutas.py": "TODAS las rutas del proyecto en un solo lugar. Arregla las rutas relativas que se rompían al abrir el programa desde otra carpeta. Guarda JSON en forma segura (si se corta la luz no se corrompe estructura.json).",
    "calc/pipeline.py": "Las 10 etapas del cálculo en orden, con sus entradas y salidas, el semáforo de estado y la función que ejecuta una etapa.",
    "tools/": "Carpeta nueva: utilidades de trabajo (no son cálculo).",
    "tools/regresion.py": "Red de seguridad: congela los resultados actuales y avisa si cambian. Se usa para validar el cambio de motor a Pynite.",
    "tools/ver_cambios.py": "Este informe.",
    "tests/": "Carpeta nueva: pruebas y copias de referencia.",
    "tests/golden/": "Copias congeladas de los 45 resultados de referencia (no se suben a git).",
}

# Los que conviene leer primero: explican todo lo demás
DESTACADOS = (
    "PUESTA_EN_MARCHA.md",
    "calc/rutas.py",
    "calc/pipeline.py",
    "tools/regresion.py",
)


def _git(argumentos: list[str]) -> str | None:
    """Ejecuta un comando de git. Devuelve None si git no está disponible."""
    try:
        proceso = subprocess.run(
            ["git", "-c", "core.quotepath=false", *argumentos],
            cwd=str(RAIZ),
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if proceso.returncode != 0:
        return None
    return proceso.stdout


def _estado() -> tuple[list[str], list[str]] | None:
    """Devuelve (nuevos, modificados) según git. None si no hay git."""
    salida = _git(["status", "--porcelain", "-uall"])
    if salida is None:
        return None
    nuevos: list[str] = []
    modificados: list[str] = []
    for linea in salida.splitlines():
        if len(linea) < 4:
            continue
        codigo = linea[:2]
        ruta = linea[3:].strip().strip('"')
        if codigo.strip() == "??":
            nuevos.append(ruta)
        elif "D" in codigo:
            modificados.append(f"{ruta}   (BORRADO)")
        else:
            modificados.append(ruta)
    return sorted(nuevos), sorted(modificados)


def _descripcion(ruta: str) -> str:
    if ruta in DESCRIPCIONES:
        return DESCRIPCIONES[ruta]
    partes = ruta.split("/")
    for i in range(len(partes) - 1, 0, -1):
        carpeta = "/".join(partes[:i]) + "/"
        if carpeta in DESCRIPCIONES:
            return DESCRIPCIONES[carpeta]
    return "(sin descripción anotada)"


def _lineas(ruta: str) -> str:
    camino = RAIZ / ruta
    extensiones = {".py", ".md", ".txt", ".bat", ".csv", ".json"}
    if camino.is_file() and camino.suffix.lower() in extensiones:
        try:
            n = len(camino.read_text(encoding="utf-8", errors="replace").splitlines())
            return f"   [{n} líneas, {max(1, camino.stat().st_size // 1024)} kB]"
        except OSError:
            return ""
    return ""


def informe() -> str:
    """Arma el texto del informe."""
    estado = _estado()
    lineas = ["QUÉ CAMBIÓ EN LA CARPETA DEL PROYECTO", "=" * 78]

    if estado is None:
        nuevos = [r for r in DESCRIPCIONES if "/" not in r or r.endswith("/")]
        modificados: list[str] = []
        lineas.append("No se pudo consultar git (¿está instalado?).")
        lineas.append("Igual, todo lo que agregué son archivos NUEVOS:")
    else:
        nuevos, modificados = estado

    lineas.append("")
    de_trabajo = [r for r in nuevos if not r.startswith("salidas/")]
    de_resultados = [r for r in nuevos if r.startswith("salidas/")]

    lineas.append(f"ARCHIVOS QUE AGREGUÉ ({len(de_trabajo)}):")
    for ruta in de_trabajo:
        lineas.append(f"  + {ruta}{_lineas(ruta)}")
        lineas.append(f"      {_descripcion(ruta)}")

    lineas.append("")
    lineas.append(f"RESULTADOS NUEVOS QUE TODAVÍA NO ESTÁN GUARDADOS EN GIT ({len(de_resultados)}):")
    lineas.append("  (los generó el propio cálculo cuando corriste P02 / P04 / P05 / P06 / L00)")
    for ruta in de_resultados:
        lineas.append(f"    - {ruta}")

    lineas.append("")
    lineas.append(f"ARCHIVOS YA EXISTENTES QUE FIGURAN MODIFICADOS ({len(modificados)}):")
    if modificados:
        for ruta in modificados:
            lineas.append(f"  * {ruta}")
        lineas.append("")
        lineas.append("  Son cambios de tus corridas anteriores. En esta sesión no se modificó")
        lineas.append("  ningún script de cálculo. Para ver antes/después de uno de ellos:")
        lineas.append("      py tools\\ver_cambios.py --diff P06_Portico_dxf.py")
    else:
        lineas.append("  (ninguno)")

    lineas.append("")
    lineas.append("PARA VER EL CONTENIDO:")
    lineas.append("  py tools\\ver_cambios.py --abrir        abre los archivos nuevos en VS Code")
    lineas.append("  py tools\\ver_cambios.py --diff <arch>  compara antes/después en VS Code")
    lineas.append("")
    lineas.append("Los scripts viejos siguen funcionando igual. Todo lo nuevo se puede")
    lineas.append("borrar sin consecuencias si no te convence.")
    return "\n".join(lineas)


def _code(argumentos: list[str]) -> bool:
    """Abre VS Code. Devuelve False si no se pudo."""
    try:
        subprocess.run(["code", *argumentos], cwd=str(RAIZ), shell=True)
        return True
    except (OSError, subprocess.SubprocessError):
        return False


def abrir_destacados() -> int:
    """Abre en VS Code los archivos que explican el cambio."""
    existentes = [d for d in DESTACADOS if (RAIZ / d).exists()]
    if not existentes:
        print("No se encontraron los archivos para abrir.")
        return 1
    if not _code(existentes):
        print("No se pudo abrir VS Code con el comando 'code'. Abrilos a mano:")
        for ruta in existentes:
            print("   " + ruta)
        return 1
    print("Abriendo en VS Code:")
    for ruta in existentes:
        print("   " + ruta)
    return 0


def abrir_diff(relativo: str) -> int:
    """Compara la versión guardada en git contra la actual, en VS Code."""
    relativo = relativo.replace("\\", "/")
    actual = RAIZ / relativo
    if not actual.exists():
        print(f"No existe el archivo: {relativo}")
        return 1

    viejo_texto = _git(["show", f"HEAD:{relativo}"])
    if viejo_texto is None:
        print(f"'{relativo}' es un archivo NUEVO: no hay versión anterior para comparar.")
        print("Se abre tal cual: todo su contenido es lo que se agregó.")
        _code([relativo])
        return 0

    temporal = RAIZ / "tests" / "_tmp"
    temporal.mkdir(parents=True, exist_ok=True)
    viejo = temporal / ("antes__" + Path(relativo).name)
    viejo.write_text(viejo_texto, encoding="utf-8")

    print(f"Comparando antes/después de {relativo} en VS Code...")
    print("   izquierda = lo último guardado en git  |  derecha = lo que hay ahora")
    if not _code(["--diff", str(viejo), str(actual)]):
        print("No se pudo abrir VS Code. La versión anterior quedó guardada en:")
        print("   " + str(viejo))
        return 1
    return 0


def _main(argv: list[str]) -> int:
    if "--diff" in argv:
        i = argv.index("--diff")
        if i + 1 >= len(argv):
            print("Falta el nombre del archivo. Ejemplo:")
            print("   py tools\\ver_cambios.py --diff P06_Portico_dxf.py")
            return 1
        return abrir_diff(argv[i + 1])
    if "--abrir" in argv:
        print(informe())
        print()
        return abrir_destacados()
    print(informe())
    return 0


if __name__ == "__main__":
    raise SystemExit(_main(sys.argv[1:]))

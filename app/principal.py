"""
app.principal — Ventana principal del Calculador (PySide6).

Esta primera versión hace:
  * Semáforo de las 10 etapas (verde / ojo / falta / a desarrollar).
  * Ver los resultados que ya existen: vigas, columnas, bases y losas.
  * Abrir cualquier archivo o carpeta de salida con doble clic.
  * Ejecutar las etapas que todavía no piden datos por teclado.

No calcula nada por su cuenta: el cálculo sigue estando en los scripts y en
`calc/`. Si acá ves un número, es porque está escrito en un archivo de
`salidas/` (y se muestra de dónde salió).
"""

from __future__ import annotations

import csv
import importlib
import os
import re
import subprocess
import sys
from pathlib import Path

# Permite abrirla tanto con `py -m app` como con `py app\principal.py`
RAIZ = Path(__file__).resolve().parent.parent
if str(RAIZ) not in sys.path:
    sys.path.insert(0, str(RAIZ))

from PySide6.QtCore import Qt  # noqa: E402
from PySide6.QtGui import QColor  # noqa: E402
from PySide6.QtWidgets import (  # noqa: E402
    QAbstractItemView,
    QApplication,
    QComboBox,
    QDoubleSpinBox,
    QDialog,
    QFormLayout,
    QGridLayout,
    QDialogButtonBox,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QPlainTextEdit,
    QScrollArea,
    QSpinBox,
    QPushButton,
    QInputDialog,
    QSplitter,
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from calc import cargas, columnas, bases, pipeline, rutas, vigas_globales  # noqa: E402
from app.cargas_proyecto import PaginaCargas  # noqa: E402
from app.muros_proyecto import PaginaMuros  # noqa: E402
from app.diseno_columnas_bases import datos_bases as pedir_datos_bases  # noqa: E402
from app.diseno_columnas_bases import respuestas_columnas  # noqa: E402
from app.diseno_vigas import dimensionar_desde_app  # noqa: E402
from app.ejes import (  # noqa: E402
    EditorAsociacionPorticos,
    EditorEjes,
    EditorPorticoPorNiveles,
    EditorMuroPorNiveles,
    EditorNiveles,
    cargar_ejes,
    cargar_niveles,
    resumen_ejes,
    resumen_niveles,
)
from app.visualizacion import VistaDiagramaMomentos, VistaPortico, VistaPlantaPorticos  # noqa: E402

# ---------------------------------------------------------------------------
# Colores del semáforo: (fondo, texto)
# ---------------------------------------------------------------------------
COLORES_ESTADO = {
    "ok": (QColor("#d6f2d6"), QColor("#14532d")),
    "desactualizada": (QColor("#fff0c2"), QColor("#7a4d00")),
    "pendiente": (QColor("#e6e6e6"), QColor("#3a3a3a")),
    "sin_datos": (QColor("#ffd9d9"), QColor("#7f1d1d")),
    "a_desarrollar": (QColor("#dbe8ff"), QColor("#1e3a8a")),
}
COLOR_MAL = QColor("#ffd9d9")     # no cumple / requiere atención
COLOR_BIEN = QColor("#eaf7ea")    # cumple

# Advertencias antes de ejecutar una etapa desde la ventana
ADVERTENCIAS = {
    "planos": (
        "P06_Portico_dxf.py todavía tiene el pórtico fijo adentro del programa\n"
        "(PORTICO = \"Portico 3\").\n\n"
        "Se va a generar el plano de 'Portico 3', sin importar qué pórtico\n"
        "tengas elegido en la ventana."
    ),
    "vigas_excel": (
        "El Excel junta los resultados de TODOS los pórticos que encuentre\n"
        "en salidas/vigas, no solo el elegido."
    ),
}


# ---------------------------------------------------------------------------
# Ayudas de formato
# ---------------------------------------------------------------------------
def numero(valor, decimales: int = 2) -> str:
    """Formatea 1234.5 -> '1.234,50' (formato argentino). Vacío si no es número."""
    try:
        texto = f"{float(valor):,.{decimales}f}"
    except (TypeError, ValueError):
        return "" if valor is None else str(valor)
    return texto.replace(",", "@").replace(".", ",").replace("@", ".")


def si_no(valor) -> str:
    if valor is None or valor == "":
        return ""
    return "Sí" if valor else "No"


def a_float(valor):
    """Convierte '297.94' o '297,94' en número. None si no se puede."""
    if isinstance(valor, (int, float)):
        return float(valor)
    try:
        return float(str(valor).replace(",", ".").strip())
    except (TypeError, ValueError):
        return None


def abrir_con_windows(ruta: Path) -> str:
    """Abre un archivo o carpeta con el programa que Windows tenga asociado."""
    try:
        os.startfile(str(ruta))  # type: ignore[attr-defined]
        return ""
    except OSError as exc:
        return f"No se pudo abrir {ruta.name}: {exc}"


def buscar_archivo(carpeta: Path, plantilla: str, portico: str) -> Path | None:
    """
    Busca un archivo de salida probando el nombre exacto (como lo escriben los
    scripts) y también con el nombre saneado.
    """
    nombres = [plantilla.format(portico=portico)]
    saneado = rutas.nombre_seguro(portico)
    if saneado != portico:
        nombres.append(plantilla.format(portico=saneado))
    for nombre in nombres:
        camino = carpeta / nombre
        if camino.exists():
            return camino
    return None


def _buscar(lista, tramo_id: str, clave: str = "tramo") -> dict:
    """Busca un tramo en una de las listas del JSON de vigas."""
    simple = tramo_id.split(".")[-1]
    for elemento in lista or []:
        if elemento.get(clave) in (tramo_id, simple):
            return elemento
    return {}


# ---------------------------------------------------------------------------
# Lectura de resultados: devuelven (archivo, filas) listos para las tablas.
# Cada fila es (celdas, esta_mal): si esta_mal es True, se pinta en rojo.
# ---------------------------------------------------------------------------
ENC_VIGAS = (
    "Tramo", "Viga", "L (m)", "b×h (cm)", "As req (cm²)", "As adop (cm²)",
    "Flexión", "δ (mm)", "δ lím (mm)", "Flecha", "Corte", "Fisur.", "Nota",
)


def datos_vigas(portico: str):
    archivo = buscar_archivo(rutas.SAL_VIGAS, "resultados_{portico}_vigas.json", portico)
    if archivo is None:
        return None, []
    datos = rutas.leer_json(archivo) or {}
    filas = []
    for tramo in datos.get("tramos", []):
        tid = tramo.get("id", "")
        mat = _buscar(datos.get("materiales", []), tid)
        flex = _buscar(datos.get("flexion", []), tid)
        fle = _buscar(datos.get("flecha", []), tid)
        corte = _buscar(datos.get("corte", []), tid)
        fis = _buscar(datos.get("fisuracion", []), tid)
        estado = _buscar(datos.get("estado", []), tid)

        cumple_flex = flex.get("cumple")
        cumple_flecha = fle.get("cumple", fle.get("cumple_L360", fle.get("cumple_L180")))
        cumple_corte = corte.get("cumple")
        cumple_fis = fis.get("cumple")
        comprobaciones = [c for c in (cumple_flex, cumple_flecha, cumple_corte, cumple_fis) if c is not None]
        malo = any(c is False for c in comprobaciones)

        b = a_float(mat.get("b_cm"))
        h = a_float(mat.get("h_cm"))
        seccion = f"{numero(b, 0)}×{numero(h, 0)}" if b and h else ""

        filas.append((
            [
                tid, tramo.get("viga", ""), numero(tramo.get("L_m")), seccion,
                numero(flex.get("As_req_cm2")), numero(flex.get("As_adop_cm2")), si_no(cumple_flex),
                numero(fle.get("delta_mm"), 1),
                numero(fle.get("limite_mm", fle.get("lim_L360_mm", fle.get("lim_L180_mm"))), 1),
                si_no(cumple_flecha), si_no(cumple_corte), si_no(cumple_fis),
                estado.get("nota", "") or estado.get("flexion", ""),
            ],
            malo,
        ))
    return archivo, filas


ENC_COLUMNAS = (
    "Columna", "Tipo", "Sección", "Dimensiones", "Altura (m)", "Pu (kN)", "Mu (kNm)",
    "f'c (MPa)", "λ", "Clasificación", "Ast (cm²)", "Armadura", "Estribo Ø",
    "Paso (cm)", "Cumple", "Nota del diagrama",
)


def datos_columnas(portico: str):
    archivo = rutas.PLANILLA_COLUMNAS
    if not archivo.exists():
        return None, []
    filas = []
    with open(archivo, newline="", encoding="utf-8") as f:
        for fila in csv.DictReader(f, delimiter=";"):
            if (fila.get("Pórtico") or "").strip() != portico:
                continue
            cumple = (fila.get("cumple_estribo") or "").strip()
            malo = cumple == "False"
            n_barras = fila.get("n_barras")
            diam = fila.get("diam_long(cm)")
            armadura = f"{n_barras} Ø{diam}" if n_barras and diam else ""
            diam_estr = (fila.get("diam_estribo(cm)") or "").strip()
            filas.append((
                [
                    fila.get("Columna", ""), fila.get("Tipo", ""), fila.get("Tipo_seccion", ""),
                    fila.get("Dimensiones", ""), numero(fila.get("Altura_libre(m)")),
                    numero(fila.get("Pu(kN)"), 1), numero(fila.get("Mu(kNm)"), 1),
                    fila.get("f'c(MPa)", ""), numero(fila.get("λ")), fila.get("Clasificación", ""),
                    numero(fila.get("Ast(cm²)")), armadura,
                    f"Ø{diam_estr}" if diam_estr else "",
                    numero(fila.get("paso_estribo(cm)"), 1),
                    ("NO" if malo else "Sí") if cumple else "",
                    (fila.get("Nota_diagrama") or "").strip(),
                ],
                malo,
            ))
    return archivo, filas


ENC_BASES = (
    "Base", "Lado (m)", "Ãrea (m²)", "Excentricidad (m)", "Espesor (m)",
    "q máx (kPa)", "q mín (kPa)", "Viga de fundación", "Armadura", "Peso (kg)",
)


def datos_bases(portico: str):
    """Devuelve (archivo, filas, nota de vigas de fundación)."""
    archivo = buscar_archivo(rutas.SAL_BASES, "bases_{portico}.json", portico)
    if archivo is None:
        return None, [], ""
    datos = rutas.leer_json(archivo) or {}
    filas = []
    for nombre, base in datos.items():
        if nombre == "vigas_fundacion" or not isinstance(base, dict):
            continue
        geotecnia = base.get("geotecnia", {}) or {}
        armadura = base.get("armadura_zapata", {}) or {}
        q_min = a_float(geotecnia.get("q_min_kPa")) or 0.0
        requiere = bool(geotecnia.get("requiere_viga_fundacion"))
        malo = requiere or q_min < 0
        pasa = armadura.get("barras_por_sentido")
        filas.append((
            [
                nombre,
                numero(geotecnia.get("lado_m")),
                numero(geotecnia.get("area_m2"), 3),
                numero(geotecnia.get("excentricidad_m")),
                numero(base.get("espesor_m")),
                numero(geotecnia.get("q_max_kPa"), 1),
                numero(geotecnia.get("q_min_kPa"), 1),
                si_no(requiere),
                f"{armadura.get('diametro', '')} c/{numero(armadura.get('paso_cm'), 0)}"
                + (f" ({pasa} por sentido)" if pasa else ""),
                numero(armadura.get("peso_total_kg"), 1),
            ],
            malo,
        ))

    nota = ""
    vigas = datos.get("vigas_fundacion")
    if isinstance(vigas, dict) and vigas:
        nota = "Vigas de fundación: " + "  |  ".join(
            f"{clave}: {v.get('seccion_cm', '')}, {v.get('armadura_longitudinal', '')}, "
            f"{numero(v.get('longitud_m'))} m"
            for clave, v in vigas.items()
        )
    return archivo, filas, nota


# ---------------------------------------------------------------------------
# Solicitaciones del motor (salidas/solicitaciones/, lo que produce calc.portico)
# ---------------------------------------------------------------------------
ENC_ENVOLVENTE = ("Barra", "Tipo", "|M| máx (kN·m)", "|V| máx (kN)", "|N| máx (kN)")


def datos_solicitaciones(portico: str):
    """
    Lee salidas/solicitaciones/<portico>.json (la salida del MOTOR) y devuelve
    (archivo, datos, filas) con la envolvente lista para una tabla.
    """
    archivo = buscar_archivo(rutas.SAL_SOLICITACIONES, "{portico}.json", portico)
    if archivo is None:
        return None, None, []
    datos = rutas.leer_json(archivo) or {}
    env = datos.get("envolvente", {})
    filas = []
    for cid, d in env.get("columnas", {}).items():
        if not d:
            continue
        m = max(abs(d.get("M_inf", 0.0)), abs(d.get("M_sup", 0.0)))
        filas.append(([cid, "columna", numero(m), "", numero(d.get("N"))], False))
    for grupo, tipo in (("vigas", "viga"), ("voladizos", "voladizo")):
        for clave, d in env.get(grupo, {}).items():
            if not d:
                continue
            v = max(abs(d.get("V_izq", 0.0)), abs(d.get("V_der", 0.0)))
            filas.append(([clave, tipo, numero(d.get("M_campo")), numero(v), ""], False))
    return archivo, datos, filas


# ---------------------------------------------------------------------------
# Losas y resumen general
# ---------------------------------------------------------------------------
def memorias_losas() -> list[Path]:
    """Memorias de losa existentes, de la más nueva a la más vieja."""
    memorias = rutas.listar(rutas.SAL_LOSAS, "memoria_losa_*.txt")
    return sorted(memorias, key=lambda p: p.stat().st_mtime, reverse=True)


def resumen_etapas(portico: str) -> str:
    """Resumen ordenado del semaforo, priorizando lo que requiere accion."""
    conteo: dict[str, int] = {}
    for dato in pipeline.semaforo(portico):
        estado = dato["estado"]
        conteo[estado] = conteo.get(estado, 0) + 1
    nombres = {
        "sin_datos": "incompletas",
        "desactualizada": "por recalcular",
        "pendiente": "pendientes",
        "ok": "al dia",
        "a_desarrollar": "a desarrollar",
    }
    orden = ("sin_datos", "desactualizada", "pendiente", "ok", "a_desarrollar")
    partes = [f"{conteo[e]} {nombres[e]}" for e in orden if conteo.get(e)]
    return "Estado: " + (" \u00b7 ".join(partes) if partes else "sin etapas")


# ---------------------------------------------------------------------------
# Ventana principal
# ---------------------------------------------------------------------------
class VentanaPrincipal(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Calculador — Estructuras")
        self.resize(1340, 820)
        self.setMinimumSize(1050, 680)
        self.setStyleSheet("""
            QMainWindow, QWidget {
                background: #f3f6fb;
                color: #172033;
                font-family: "Segoe UI";
                font-size: 10pt;
            }
            QGroupBox {
                background: #ffffff;
                border: 1px solid #d9e1ec;
                border-radius: 9px;
                margin-top: 12px;
                padding: 12px 10px 10px;
                font-weight: 600;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 12px;
                padding: 0 5px;
                color: #334155;
            }
            QPushButton {
                background: #ffffff;
                color: #26354d;
                border: 1px solid #cbd5e1;
                border-radius: 7px;
                padding: 9px 16px;
                min-height: 32px;
            }
            QPushButton:hover {
                background: #eff6ff;
                border-color: #93b4e8;
            }
            QPushButton:pressed { background: #dbeafe; }
            QPushButton:disabled {
                color: #94a3b8;
                background: #f1f5f9;
                border-color: #e2e8f0;
            }
            QPushButton[role="primary"] {
                background: #1d4ed8;
                color: #ffffff;
                border: 1px solid #1d4ed8;
                font-weight: 600;
                min-height: 32px;
            }
            QPushButton[role="primary"]:hover {
                background: #1e40af;
                border-color: #1e40af;
            }
            QComboBox, QLineEdit, QPlainTextEdit, QTreeWidget, QListWidget,
            QTableWidget {
                background: #ffffff;
                color: #172033;
                border: 1px solid #d5deea;
                border-radius: 6px;
                selection-background-color: #dbeafe;
                selection-color: #172033;
            }
            QComboBox, QLineEdit { padding: 6px 9px; }
            QComboBox::drop-down { border: none; width: 20px; }
            QTabWidget::pane {
                background: #f8fafc;
                border: 1px solid #d9e1ec;
                border-radius: 7px;
                top: -1px;
            }
            QTabBar::tab {
                background: #e8edf5;
                color: #475569;
                border: 1px solid #d9e1ec;
                border-bottom: none;
                padding: 9px 14px;
                margin-right: 3px;
                border-top-left-radius: 6px;
                border-top-right-radius: 6px;
            }
            QTabBar::tab:selected {
                background: #ffffff;
                color: #1d4ed8;
                font-weight: 600;
            }
            QHeaderView::section {
                background: #edf2f8;
                color: #334155;
                border: none;
                border-bottom: 1px solid #cbd5e1;
                padding: 8px 6px;
                font-weight: 600;
            }
            QTableWidget { gridline-color: #e7edf5; }
            QStatusBar {
                background: #eaf0f8;
                color: #475569;
                border-top: 1px solid #d9e1ec;
            }

            /* Ventanas emergentes: diálogos y mensajes */
            QDialog {
                background: #f3f6fb;
                color: #172033;
                font-family: "Segoe UI";
                font-size: 10pt;
            }
            QMessageBox {
                background: #ffffff;
                color: #172033;
                min-width: 420px;
            }
            QMessageBox QLabel {
                color: #172033;
                font-size: 10pt;
                padding: 6px;
            }
            QMessageBox QPushButton {
                min-width: 90px;
                min-height: 30px;
                padding: 6px 14px;
                font-size: 10pt;
            }

            /* Barras de desplazamiento más discretas */
            QScrollBar:vertical {
                background: #f1f5f9;
                width: 12px;
                border-radius: 6px;
                margin: 2px;
            }
            QScrollBar::handle:vertical {
                background: #cbd5e1;
                border-radius: 6px;
                min-height: 40px;
            }
            QScrollBar::handle:vertical:hover { background: #94a3b8; }
            QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {
                background: none;
                border: none;
                height: 0;
            }
            QScrollBar:horizontal {
                background: #f1f5f9;
                height: 12px;
                border-radius: 6px;
                margin: 2px;
            }
            QScrollBar::handle:horizontal {
                background: #cbd5e1;
                border-radius: 6px;
                min-width: 40px;
            }
            QScrollBar::handle:horizontal:hover { background: #94a3b8; }
            QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {
                background: none;
                border: none;
                width: 0;
            }
        """)

        # --- barra de arriba ---
        self.titulo_aplicacion = QLabel("CALCULADOR")
        self.titulo_aplicacion.setObjectName("tituloAplicacion")
        self.combo_portico = QComboBox()
        self.combo_portico.setMinimumWidth(190)
        self.combo_obras = QComboBox()
        self.combo_obras.setMinimumWidth(190)
        self.combo_nivel_planta = QComboBox()
        self.boton_nueva_obra = QPushButton("Nueva obra")
        self.boton_eliminar_obra = QPushButton("Eliminar obra")
        self.boton_actualizar = QPushButton("Actualizar")
        self.etiqueta_resumen = QLabel()
        self.etiqueta_resumen.setWordWrap(True)

        # --- pestaña Etapas ---
        self.tabla_etapas = self._tabla(("Estado", "Etapa", "Detalle", "Programa"))
        self.detalle_etapa = QPlainTextEdit()
        self.detalle_etapa.setReadOnly(True)
        self.boton_ejecutar = QPushButton("Ejecutar esta etapa")
        self.boton_ejecutar.setEnabled(False)
        self.boton_abrir_salida = QPushButton("Abrir la salida de la etapa")
        self.boton_abrir_salida.setEnabled(False)

        # --- pestañas de resultados ---
        self.tabla_vigas = self._tabla(ENC_VIGAS)
        self.boton_dimensionar_vigas = QPushButton("Dimensionar vigas")
        self.advertencias_vigas: list[str] = []
        self.etiqueta_planillas_vigas = QLabel()
        self.lista_planillas_vigas = QListWidget()
        self.texto_planilla_viga = QPlainTextEdit()
        self.texto_planilla_viga.setReadOnly(True)
        self.texto_planilla_viga.setStyleSheet(
            "QPlainTextEdit { color: #111827; background: #ffffff; }"
        )
        self.tabla_columnas = self._tabla(ENC_COLUMNAS)
        self.boton_dimensionar_columnas = QPushButton("Dimensionar columnas")
        self.etiqueta_planillas_columnas = QLabel()
        self.lista_planillas_columnas = QListWidget()
        self.texto_planilla_columnas = QPlainTextEdit()
        self.texto_planilla_columnas.setReadOnly(True)
        self.texto_planilla_columnas.setStyleSheet(
            "QPlainTextEdit { color: #111827; background: #ffffff; }"
        )
        self.tabla_bases = self._tabla(ENC_BASES)
        self.boton_dimensionar_bases = QPushButton("Dimensionar bases")
        self.etiqueta_planillas_bases = QLabel()
        self.lista_planillas_bases = QListWidget()
        self.texto_planilla_bases = QPlainTextEdit()
        self.texto_planilla_bases.setReadOnly(True)
        self.texto_planilla_bases.setStyleSheet(
            "QPlainTextEdit { color: #111827; background: #ffffff; }"
        )
        self.lista_losas = QListWidget()
        self.texto_losa = QPlainTextEdit()
        self.texto_losa.setReadOnly(True)
        self.arbol = QTreeWidget()
        self.arbol.setHeaderLabels(("Archivo", "Tamaño"))

        # --- pestaña Inicio ---
        self.etiqueta_inicio = QLabel()
        # Inicializar widgets que se conectarán después
        self.boton_crear_muro = None
        self.etiqueta_inicio.setWordWrap(True)
        self.nombre_proyecto = QLineEdit()
        self.id_proyecto_inicio = QLabel()
        self.id_proyecto_inicio.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self.ubicacion_proyecto = QLineEdit()
        self.ubicacion_proyecto.setPlaceholderText("Localidad, provincia, país")
        self.comitente_proyecto = QLineEdit()
        self.comitente_proyecto.setPlaceholderText("Nombre del cliente o comitente")
        self.responsable_proyecto = QLineEdit()
        self.responsable_proyecto.setPlaceholderText("Profesional responsable")
        self.notas_proyecto = QPlainTextEdit()
        self.notas_proyecto.setPlaceholderText("Criterios, contexto y notas de la obra")
        self.notas_proyecto.setMaximumHeight(44)
        self.boton_guardar_proyecto = QPushButton("Guardar carátula")
        self.boton_crear_portico = QPushButton("Crear pórtico…")
        self.boton_crear_portico.setToolTip(
            "Crear la geometría de un pórtico y definir columnas y vigas por nivel."
        )
        self.boton_crear_muro = QPushButton("Crear muro…", self)
        self.boton_crear_muro.setToolTip(
            "Crear la geometría de un muro y definir sus parámetros."
        )
        self.boton_eliminar_portico = QPushButton("Eliminar pórtico / muro…")
        self.boton_eliminar_portico.setEnabled(False)
        self.boton_definir_ejes = QPushButton("Definir ejes X/Y…")
        self.boton_definir_niveles = QPushButton("Definir niveles Z…")
        self.boton_asociar_porticos = QPushButton("Asociar ejes…")
        self.boton_definir_ejes.setToolTip("Definir las coordenadas de los ejes X e Y.")
        self.boton_definir_niveles.setToolTip("Definir cotas de niveles y entrepisos.")
        self.boton_asociar_porticos.setToolTip(
            "Vincula cada pórtico existente a los ejes de planta (plano y primera columna). "
            "Sirve para decir sobre qué eje transversal y a qué distancia queda cada pórtico, "
            "sin desfases. Úsalo después de definir ejes y niveles; no es obligatorio si los "
            "pórticos ya tienen su referencia de planta cargada."
        )
        self.resumen_ejes_inicio = QLabel()
        self.resumen_ejes_inicio.setWordWrap(True)
        self.resumen_niveles_inicio = QLabel()
        self.resumen_niveles_inicio.setWordWrap(True)
        self.referencia_portico_inicio = QLabel()
        self.referencia_portico_inicio.setWordWrap(True)
        self.boton_inicio_cargas = QPushButton("Definir cargas del proyecto")
        self.boton_inicio_aplicaciones = QPushButton("Asignar cargas a barras")
        self.boton_inicio_losas = QPushButton("Nueva losa en ejes…")
        self.vista_portico = VistaPortico()
        self.vista_planta = VistaPlantaPorticos()
        self.scroll_vista_portico_p = QScrollArea()
        self.scroll_vista_portico_p.setWidgetResizable(True)
        self.scroll_vista_portico_p.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.scroll_vista_portico_p.setMinimumHeight(340)
        self.scroll_vista_portico_p.setWidget(self.vista_planta)
        self.boton_referencias_cargas = QPushButton("Referencias de cargas…")
        self.boton_referencias_cargas.setToolTip(
            "Abrir en otra ventana los nombres, tipos y categorías de las cargas representadas."
        )
        self.scroll_vista_portico = QScrollArea()
        self.scroll_vista_portico.setWidgetResizable(True)
        self.scroll_vista_portico.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.scroll_vista_portico.setMinimumHeight(340)
        self.scroll_vista_portico.setWidget(self.vista_portico)
        self.boton_referencias_cargas.setMaximumHeight(30)
        self.lbl_cargas = QLabel()
        self.lbl_cargas.setWordWrap(True)
        self.lbl_motor = QLabel()
        self.lbl_motor.setWordWrap(True)
        self.boton_ver_cargas = QPushButton("Abrir el último análisis")
        self.boton_resolver_motor = QPushButton("Resolver solicitaciones del pórtico")
        self.boton_ver_motor = QPushButton("Abrir el JSON del motor")
        self.boton_diagrama_momentos = QPushButton("Ver diagrama de momentos")
        self.boton_diagrama_momentos.setEnabled(False)
        self.tabla_envolvente = self._tabla(ENC_ENVOLVENTE)
        self.tabla_portico_resumen = self._tabla(("Elemento", "Cantidad", "Detalle"))
        self.ultimo_analisis: Path | None = None

        # --- inicializar estructura, ejes y niveles ---
        self.estructura = rutas.cargar_estructura()
        self.ejes = cargar_ejes()
        self.niveles = cargar_niveles()

        for boton in (
            self.boton_resolver_motor,
            self.boton_dimensionar_vigas,
            self.boton_dimensionar_columnas,
            self.boton_dimensionar_bases,
        ):
            boton.setProperty("role", "primary")

        self.pagina_cargas = PaginaCargas(self.refrescar, self._portico)
        self.pagina_muros = PaginaMuros(self)
        self._armar_interfaz()
        self._conectar()
        self.refrescar()

    # ------------------------------------------------------------------
    # Armado de la interfaz
    # ------------------------------------------------------------------
    def _tabla(self, encabezados) -> QTableWidget:
        tabla = QTableWidget(0, len(encabezados))
        tabla.setHorizontalHeaderLabels(list(encabezados))
        tabla.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        tabla.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        tabla.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        tabla.setAlternatingRowColors(True)
        tabla.verticalHeader().setVisible(False)
        tabla.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        return tabla

    def _llenar(self, tabla: QTableWidget, filas) -> None:
        """filas = [(celdas, esta_mal), ...]"""
        tabla.setRowCount(len(filas))
        for i, (celdas, malo) in enumerate(filas):
            for j, valor in enumerate(celdas):
                item = QTableWidgetItem(str(valor))
                item.setToolTip(str(valor))
                alineacion = Qt.AlignmentFlag.AlignLeft if j == 0 else Qt.AlignmentFlag.AlignCenter
                item.setTextAlignment(alineacion | Qt.AlignmentFlag.AlignVCenter)
                if malo:
                    item.setBackground(COLOR_MAL)
                tabla.setItem(i, j, item)

    def _pagina(self, texto_etiqueta: str, tabla: QTableWidget):
        pagina = QWidget()
        caja = QVBoxLayout(pagina)
        caja.setContentsMargins(14, 12, 14, 14)
        caja.setSpacing(10)
        etiqueta = QLabel(texto_etiqueta)
        etiqueta.setWordWrap(True)
        caja.addWidget(etiqueta)
        caja.addWidget(tabla)
        return pagina, etiqueta

    def _armar_interfaz(self) -> None:
        contenedor = QWidget()
        principal = QVBoxLayout(contenedor)
        principal.setContentsMargins(14, 12, 14, 10)
        principal.setSpacing(10)

        barra = QHBoxLayout()
        barra.setSpacing(9)
        barra.addWidget(self.titulo_aplicacion)
        barra.addWidget(QLabel("Obra:"))
        barra.addWidget(self.combo_obras)
        barra.addWidget(self.boton_nueva_obra)
        barra.addWidget(self.boton_eliminar_obra)
        barra.addWidget(self.boton_actualizar)
        barra.addWidget(self.etiqueta_resumen, 1)
        banda_superior = QWidget()
        banda_superior.setObjectName("bandaSuperior")
        banda_superior.setStyleSheet(
            "#bandaSuperior { background: #ffffff; border: 1px solid #d9e1ec; "
            "border-radius: 8px; }"
            "#tituloAplicacion { color: #17356f; font-size: 13pt; "
            "font-weight: 700; padding: 4px 8px; }"
        )
        banda_superior.setLayout(barra)
        principal.addWidget(banda_superior)

        pestanias = QTabWidget()
        pestanias.setDocumentMode(True)
        pestanias.addTab(self._pagina_inicio(), "Inicio")
        pestanias.addTab(self.pagina_cargas, "Cargas")
        pestanias.addTab(self._pagina_porticos(), "Pórticos")
        pestanias.addTab(self._pagina_estructura_global(), "Estructura global")
        pestanias.addTab(self._pagina_etapas(), "1 · Estado y etapas")


        pestanias.addTab(self._pagina_losas(), "5 · Losas")
        pestanias.addTab(self.pagina_muros, "6 · Muros")
        pestanias.addTab(self._pagina_archivos(), "7 · Archivos")

        principal.addWidget(pestanias, 1)
        self.setCentralWidget(contenedor)
        self.statusBar().showMessage(f"Carpeta del proyecto: {rutas.RAIZ}")

    def _pagina_vigas(self) -> QWidget:
        pagina = QWidget()
        caja_vigas = QVBoxLayout(pagina)
        self.etiqueta_vigas = QLabel()
        self.etiqueta_vigas.setWordWrap(True)
        caja_vigas.addWidget(self.etiqueta_vigas)
        fila_dimensionado = QHBoxLayout()
        fila_dimensionado.addWidget(self.boton_dimensionar_vigas)
        fila_dimensionado.addStretch(1)
        caja_vigas.addLayout(fila_dimensionado)
        contenido_vigas = QSplitter(Qt.Orientation.Vertical)
        contenido_vigas.addWidget(self.tabla_vigas)
        panel_planillas = QWidget()
        caja_planillas = QVBoxLayout(panel_planillas)
        caja_planillas.setContentsMargins(0, 0, 0, 0)
        caja_planillas.addWidget(self.etiqueta_planillas_vigas)
        division_planillas = QSplitter(Qt.Orientation.Horizontal)
        division_planillas.addWidget(self.lista_planillas_vigas)
        division_planillas.addWidget(self.texto_planilla_viga)
        division_planillas.setSizes([260, 700])
        caja_planillas.addWidget(division_planillas, 1)
        contenido_vigas.addWidget(panel_planillas)
        contenido_vigas.setSizes([300, 320])
        caja_vigas.addWidget(contenido_vigas, 1)
        self.etiqueta_vigas_fuente = self.etiqueta_vigas
        return pagina

    def _pagina_columnas(self) -> QWidget:
        pagina = QWidget()
        caja_columnas = QVBoxLayout(pagina)
        self.etiqueta_columnas = QLabel()
        self.etiqueta_columnas.setWordWrap(True)
        caja_columnas.addWidget(self.etiqueta_columnas)
        nota_columnas = QLabel(
            "Nota: esta pestaña calcula cada columna en forma plana según el pórtico seleccionado. "
            "La consolidación global (Mx/My por ID de columna compartida) se realiza en “Columnas globales”."
        )
        nota_columnas.setWordWrap(True)
        caja_columnas.addWidget(nota_columnas)
        caja_columnas.addWidget(self.boton_dimensionar_columnas, 0, Qt.AlignmentFlag.AlignLeft)
        contenido_columnas = QSplitter(Qt.Orientation.Vertical)
        contenido_columnas.addWidget(self.tabla_columnas)
        panel_planillas_columnas = QWidget()
        caja_planillas_columnas = QVBoxLayout(panel_planillas_columnas)
        caja_planillas_columnas.setContentsMargins(0, 0, 0, 0)
        caja_planillas_columnas.addWidget(self.etiqueta_planillas_columnas)
        division_planillas_columnas = QSplitter(Qt.Orientation.Horizontal)
        division_planillas_columnas.addWidget(self.lista_planillas_columnas)
        division_planillas_columnas.addWidget(self.texto_planilla_columnas)
        division_planillas_columnas.setSizes([260, 700])
        caja_planillas_columnas.addWidget(division_planillas_columnas, 1)
        contenido_columnas.addWidget(panel_planillas_columnas)
        contenido_columnas.setSizes([300, 320])
        caja_columnas.addWidget(contenido_columnas, 1)
        return pagina

    def _pagina_columnas_globales(self) -> QWidget:
        pagina = QWidget()
        caja = QVBoxLayout(pagina)
        etiqueta = QLabel("Columnas consolidadas por ID global (COL-XXX).")
        etiqueta.setWordWrap(True)
        caja.addWidget(etiqueta)
        boton = QPushButton("Calcular columnas globales")
        boton.clicked.connect(self._calcular_columnas_globales)
        caja.addWidget(boton, 0, Qt.AlignmentFlag.AlignLeft)
        self.tabla_columnas_globales = QTableWidget(0, 9)
        self.tabla_columnas_globales.setHorizontalHeaderLabels(
            ("ID", "Eje", "Nivel", "Mx", "My", "N", "Deriva", "Pórticos", "Estado")
        )
        caja.addWidget(self.tabla_columnas_globales, 1)
        return pagina

    def _calcular_columnas_globales(self) -> None:
        try:
            resultados = columnas.calcular_global()
            self.tabla_columnas_globales.setRowCount(0)
            for r in resultados:
                fila = self.tabla_columnas_globales.rowCount()
                self.tabla_columnas_globales.insertRow(fila)
                self.tabla_columnas_globales.setItem(fila, 0, QTableWidgetItem(str(r.get("id"))))
                self.tabla_columnas_globales.setItem(fila, 1, QTableWidgetItem(str(r.get("eje_nombre") or r.get("eje_id"))))
                self.tabla_columnas_globales.setItem(fila, 2, QTableWidgetItem(str(r.get("nivel"))))
                self.tabla_columnas_globales.setItem(fila, 3, QTableWidgetItem(str(r.get("Mx"))))
                self.tabla_columnas_globales.setItem(fila, 4, QTableWidgetItem(str(r.get("My"))))
                self.tabla_columnas_globales.setItem(fila, 5, QTableWidgetItem(str(r.get("N"))))
                self.tabla_columnas_globales.setItem(fila, 6, QTableWidgetItem(str(r.get("deriva"))))
                self.tabla_columnas_globales.setItem(fila, 7, QTableWidgetItem(", ".join(r.get("porticos", []))))
                self.tabla_columnas_globales.setItem(fila, 8, QTableWidgetItem(str(r.get("estado"))))
        except Exception as exc:
            QMessageBox.warning(self, "Columnas globales", str(exc))

    def _pagina_vigas_globales(self) -> QWidget:
        pagina = QWidget()
        caja = QVBoxLayout(pagina)
        etiqueta = QLabel("Vigas consolidadas por tramo de pórtico (VIG-P{pórtico}-{tramo}).")
        etiqueta.setWordWrap(True)
        caja.addWidget(etiqueta)
        boton = QPushButton("Calcular vigas globales")
        boton.clicked.connect(self._calcular_vigas_globales)
        caja.addWidget(boton, 0, Qt.AlignmentFlag.AlignLeft)
        self.tabla_vigas_globales = QTableWidget(0, 6)
        self.tabla_vigas_globales.setHorizontalHeaderLabels(
            ("ID", "Pórtico", "Nivel", "M", "V", "Estado")
        )
        caja.addWidget(self.tabla_vigas_globales, 1)
        return pagina

    def _calcular_vigas_globales(self) -> None:
        try:
            resultados = vigas_globales.calcular_global()
            self.tabla_vigas_globales.setRowCount(0)
            for r in resultados:
                fila = self.tabla_vigas_globales.rowCount()
                self.tabla_vigas_globales.insertRow(fila)
                self.tabla_vigas_globales.setItem(fila, 0, QTableWidgetItem(str(r.get("id"))))
                self.tabla_vigas_globales.setItem(fila, 1, QTableWidgetItem(str(r.get("portico"))))
                self.tabla_vigas_globales.setItem(fila, 2, QTableWidgetItem(str(r.get("nivel"))))
                self.tabla_vigas_globales.setItem(fila, 3, QTableWidgetItem(str(r.get("M"))))
                self.tabla_vigas_globales.setItem(fila, 4, QTableWidgetItem(str(r.get("V"))))
                self.tabla_vigas_globales.setItem(fila, 5, QTableWidgetItem(str(r.get("estado"))))
        except Exception as exc:
            QMessageBox.warning(self, "Vigas globales", str(exc))

    def _pagina_bases_globales(self) -> QWidget:
        pagina = QWidget()
        caja = QVBoxLayout(pagina)
        etiqueta = QLabel("Bases consolidadas por ID global (BASE-XXX), con Mx/My acumulados.")
        etiqueta.setWordWrap(True)
        caja.addWidget(etiqueta)
        boton = QPushButton("Calcular bases globales")
        boton.clicked.connect(self._calcular_bases_globales)
        caja.addWidget(boton, 0, Qt.AlignmentFlag.AlignLeft)
        self.tabla_bases_globales = QTableWidget(0, 8)
        self.tabla_bases_globales.setHorizontalHeaderLabels(
            ("ID", "Mx", "My", "M", "Fx", "Fy", "Pórticos", "Estado")
        )
        caja.addWidget(self.tabla_bases_globales, 1)
        return pagina

    def _calcular_bases_globales(self) -> None:
        try:
            resultados = bases.calcular_global()
            self.tabla_bases_globales.setRowCount(0)
            for r in resultados:
                fila = self.tabla_bases_globales.rowCount()
                self.tabla_bases_globales.insertRow(fila)
                self.tabla_bases_globales.setItem(fila, 0, QTableWidgetItem(str(r.get("id"))))
                self.tabla_bases_globales.setItem(fila, 1, QTableWidgetItem(str(r.get("Mx"))))
                self.tabla_bases_globales.setItem(fila, 2, QTableWidgetItem(str(r.get("My"))))
                self.tabla_bases_globales.setItem(fila, 3, QTableWidgetItem(str(r.get("M"))))
                self.tabla_bases_globales.setItem(fila, 4, QTableWidgetItem(str(r.get("Fx"))))
                self.tabla_bases_globales.setItem(fila, 5, QTableWidgetItem(str(r.get("Fy"))))
                self.tabla_bases_globales.setItem(fila, 6, QTableWidgetItem(", ".join(r.get("porticos", []))))
                self.tabla_bases_globales.setItem(fila, 7, QTableWidgetItem(str(r.get("estado"))))
        except Exception as exc:
            QMessageBox.warning(self, "Bases globales", str(exc))

    def _pagina_estructura_global(self) -> QWidget:
        contenedor = QWidget()
        caja = QVBoxLayout(contenedor)
        subpestanias = QTabWidget()
        subpestanias.setDocumentMode(True)
        subpestanias.addTab(self._pagina_columnas_globales(), "Columnas globales")
        subpestanias.addTab(self._pagina_vigas_globales(), "Vigas globales")
        subpestanias.addTab(self._pagina_bases_globales(), "Bases globales")
        caja.addWidget(subpestanias)
        return contenedor

    def _pagina_bases(self) -> QWidget:
        pagina = QWidget()
        caja_bases = QVBoxLayout(pagina)
        self.etiqueta_bases = QLabel()
        self.etiqueta_bases.setWordWrap(True)
        caja_bases.addWidget(self.etiqueta_bases)
        caja_bases.addWidget(self.boton_dimensionar_bases, 0, Qt.AlignmentFlag.AlignLeft)
        contenido_bases = QSplitter(Qt.Orientation.Vertical)
        contenido_bases.addWidget(self.tabla_bases)
        panel_planillas_bases = QWidget()
        caja_planillas_bases = QVBoxLayout(panel_planillas_bases)
        caja_planillas_bases.setContentsMargins(0, 0, 0, 0)
        caja_planillas_bases.addWidget(self.etiqueta_planillas_bases)
        division_planillas_bases = QSplitter(Qt.Orientation.Horizontal)
        division_planillas_bases.addWidget(self.lista_planillas_bases)
        division_planillas_bases.addWidget(self.texto_planilla_bases)
        division_planillas_bases.setSizes([260, 700])
        caja_planillas_bases.addWidget(division_planillas_bases, 1)
        contenido_bases.addWidget(panel_planillas_bases)
        contenido_bases.setSizes([300, 320])
        caja_bases.addWidget(contenido_bases, 1)
        return pagina

    def _pagina_etapas(self) -> QWidget:
        pagina = QWidget()
        caja = QVBoxLayout(pagina)
        caja.addWidget(self.tabla_etapas, 3)

        abajo = QHBoxLayout()
        abajo.addWidget(self.detalle_etapa, 3)
        botones = QVBoxLayout()
        botones.addWidget(self.boton_ejecutar)
        botones.addWidget(self.boton_abrir_salida)
        botones.addStretch(1)
        abajo.addLayout(botones, 1)
        caja.addLayout(abajo, 2)
        return pagina

    def _pagina_losas(self) -> QWidget:
        pagina = QWidget()
        caja = QVBoxLayout(pagina)
        self.etiqueta_losas = QLabel("Doble clic en una memoria para leerla acá al lado.")
        caja.addWidget(self.etiqueta_losas)
        division = QSplitter(Qt.Orientation.Horizontal)
        division.addWidget(self.lista_losas)
        division.addWidget(self.texto_losa)
        division.setSizes([280, 900])
        caja.addWidget(division, 1)
        return pagina

    def _pagina_terreno(self) -> QWidget:
        pagina = QWidget()
        caja = QVBoxLayout(pagina)
        grupo = QGroupBox("Datos de fundación / terreno (toda la obra)")
        form = QFormLayout(grupo)

        self.spin_profundidad_fundacion = QDoubleSpinBox()
        self.spin_profundidad_fundacion.setRange(0.1, 20.0)
        self.spin_profundidad_fundacion.setDecimals(2)
        self.spin_profundidad_fundacion.setSuffix(" m")

        self.spin_q_adm = QDoubleSpinBox()
        self.spin_q_adm.setRange(1.0, 10000.0)
        self.spin_q_adm.setDecimals(2)
        self.spin_q_adm.setSuffix(" kPa")

        form.addRow("Profundidad hasta el fondo de zapata:", self.spin_profundidad_fundacion)
        form.addRow("Presión admisible del suelo q_adm:", self.spin_q_adm)

        btn = QPushButton("Guardar datos de terreno")
        btn.clicked.connect(self._guardar_terreno)

        caja.addWidget(grupo)
        caja.addWidget(btn)
        caja.addStretch(1)
        self._cargar_terreno()
        return pagina

    def _cargar_terreno(self) -> None:
        datos = rutas.leer_json(rutas.TERRENO, {}) or {}
        self.spin_profundidad_fundacion.setValue(float(datos.get("profundidad_fundacion_m", 0.8)))
        self.spin_q_adm.setValue(float(datos.get("q_adm_kPa", 0.87 * 98.1)))

    def _guardar_terreno(self) -> None:
        datos = rutas.leer_json(rutas.TERRENO, {}) or {}
        datos.update({
            "profundidad_fundacion_m": float(self.spin_profundidad_fundacion.value()),
            "q_adm_kPa": float(self.spin_q_adm.value()),
        })
        rutas.guardar_json(rutas.TERRENO, datos)
        QMessageBox.information(self, "Terreno", "Datos de fundación guardados.")


    def _pagina_parametros_motor(self) -> QWidget:
        pagina = QWidget()
        caja = QVBoxLayout(pagina)
        grupo = QGroupBox("Parámetros del motor de pórticos")
        form = QFormLayout(grupo)

        def spin(minimo, maximo, decimales=6):
            c = QDoubleSpinBox()
            c.setRange(minimo, maximo)
            c.setDecimals(decimales)
            return c

        self.spin_e_mod = spin(0.001, 1e6, 6)
        self.spin_g_mod = spin(0.001, 1e6, 6)
        self.spin_ea = spin(0.0, 1e9, 3)
        self.spin_ei = spin(0.0, 1e9, 3)
        self.spin_poisson = spin(0.0, 0.5, 3)
        self.spin_peso_hormigon = spin(1.0, 100.0, 2)
        self.txt_reglamento = QLineEdit("CIRSOC 201/2024")

        form.addRow("Módulo de elasticidad E_MOD:", self.spin_e_mod)
        form.addRow("Módulo de corte G_MOD:", self.spin_g_mod)
        form.addRow("Rigidez axial EA:", self.spin_ea)
        form.addRow("Rigidez a flexión EI:", self.spin_ei)
        form.addRow("Coeficiente de Poisson ν:", self.spin_poisson)
        form.addRow("Peso específico hormigón (kN/m³):", self.spin_peso_hormigon)
        form.addRow("Reglamento aplicado:", self.txt_reglamento)

        btn = QPushButton("Guardar parámetros")
        btn.clicked.connect(self._guardar_parametros_motor)

        caja.addWidget(grupo)
        caja.addWidget(btn)
        caja.addStretch(1)
        self._cargar_parametros_motor()
        return pagina

    def _cargar_parametros_motor(self) -> None:
        datos = rutas.leer_json(rutas.PARAMETROS_MOTOR, {}) or {}
        self.spin_e_mod.setValue(float(datos.get("E_MOD", 1.0)))
        self.spin_g_mod.setValue(float(datos.get("G_MOD", 1.0)))
        self.spin_ea.setValue(float(datos.get("EA", 15000.0)))
        self.spin_ei.setValue(float(datos.get("EI", 5000.0)))
        self.spin_poisson.setValue(float(datos.get("POISSON", 0.3)))
        self.spin_peso_hormigon.setValue(float(datos.get("PESO_ESPECIFICO_HORMIGON_KN_M3", 25.0)))
        self.txt_reglamento.setText(str(datos.get("VERSION_MOTOR", "CIRSOC 201/2024")))

    def _guardar_parametros_motor(self) -> None:
        datos = {
            "E_MOD": float(self.spin_e_mod.value()),
            "G_MOD": float(self.spin_g_mod.value()),
            "EA": float(self.spin_ea.value()),
            "EI": float(self.spin_ei.value()),
            "POISSON": float(self.spin_poisson.value()),
            "PESO_ESPECIFICO_HORMIGON_KN_M3": float(self.spin_peso_hormigon.value()),
            "VERSION_MOTOR": self.txt_reglamento.text().strip(),
        }
        rutas.guardar_json(rutas.PARAMETROS_MOTOR, datos)
        import importlib
        import calc.portico as _cp
        importlib.reload(_cp)
        QMessageBox.information(self, "Parámetros", "Parámetros guardados.")


    def _pagina_porticos(self) -> QWidget:
        pagina = QWidget()
        caja = QVBoxLayout(pagina)

        grupo_portico = QGroupBox("Pórtico activo")
        form_portico = QFormLayout(grupo_portico)
        form_portico.addWidget(QLabel("Elegí el pórtico que mirás en este visor:"))
        form_portico.addWidget(self.combo_portico)

        panel_visor = QWidget()
        panel_visor.setMinimumSize(500, 400)
        caja_vista = QVBoxLayout(panel_visor)
        caja_vista.addWidget(self.scroll_vista_portico, 1)

        grupo_motor = QGroupBox("Solicitaciones del pórtico")
        caja_motor = QVBoxLayout(grupo_motor)
        caja_motor.addWidget(self.lbl_motor)
        fila_resolver = QHBoxLayout()
        fila_resolver.addWidget(self.boton_resolver_motor)
        fila_resolver.addStretch(1)
        caja_motor.addLayout(fila_resolver)

        fila_ver = QHBoxLayout()
        fila_ver.addWidget(self.boton_ver_motor)
        fila_ver.addStretch(1)
        caja_motor.addLayout(fila_ver)
        caja_motor.addWidget(self.boton_diagrama_momentos)
        caja_motor.addWidget(self.boton_eliminar_portico)

        # Panel lateral (a la derecha de la imagen): controles y botones
        panel_lateral_portico = QWidget()
        caja_lateral = QVBoxLayout(panel_lateral_portico)
        caja_lateral.setContentsMargins(0, 0, 0, 0)
        caja_lateral.addWidget(grupo_portico)
        caja_lateral.addWidget(grupo_motor)
        caja_lateral.addStretch(1)

        panel_lateral_portico.setMaximumWidth(340)
        # Parte superior: imagen a la izquierda (más espacio), controles a la derecha
        division_superior = QSplitter(Qt.Orientation.Horizontal)
        division_superior.addWidget(panel_visor)
        division_superior.addWidget(panel_lateral_portico)
        division_superior.setStretchFactor(0, 5)
        division_superior.setStretchFactor(1, 1)
        division_superior.setSizes([820, 340])
        division_superior.setChildrenCollapsible(False)

        # Parte inferior: datos del pórtico (resumen + envolvente) en pantalla dividida
        panel_datos = QWidget()
        caja_datos = QVBoxLayout(panel_datos)
        caja_datos.setContentsMargins(0, 0, 0, 0)
        fila_resumen = QHBoxLayout()
        fila_resumen.addWidget(self.tabla_portico_resumen)
        caja_datos.addLayout(fila_resumen)
        caja_datos.addWidget(QLabel("Envolvente del pórtico elegido (máximos en módulo, según el motor):"))
        caja_datos.addWidget(self.tabla_envolvente, 1)

        # La tabla de envolvente puede quedar contraída (solo encabezado) cuando el
        # motor aún no corrió; garantizamos una altura mínima y que la última
        # columna estire el ancho disponible para que la tabla no quede oculta.
        header = self.tabla_envolvente.horizontalHeader()
        header.setStretchLastSection(True)
        self.tabla_envolvente.setMinimumHeight(130)
        panel_datos.setMinimumHeight(160)

        division_vertical = QSplitter(Qt.Orientation.Vertical)
        division_vertical.addWidget(division_superior)
        division_vertical.addWidget(panel_datos)
        division_vertical.setStretchFactor(0, 3)
        division_vertical.setStretchFactor(1, 2)
        division_vertical.setSizes([700, 320])
        division_vertical.setChildrenCollapsible(False)
        caja.addWidget(division_vertical, 1)

        pagina_portico = pagina
        contenedor_porticos = QWidget()
        caja_porticos = QVBoxLayout(contenedor_porticos)
        subpestanias = QTabWidget()
        subpestanias.setDocumentMode(True)
        subpestanias.addTab(pagina_portico, "Pórtico")
        subpestanias.addTab(self._pagina_vigas(), "Vigas")
        subpestanias.addTab(self._pagina_columnas(), "Columnas")
        subpestanias.addTab(self._pagina_bases(), "Bases")
        subpestanias.addTab(self._pagina_terreno(), "Terreno")
        subpestanias.addTab(self._pagina_parametros_motor(), "Parámetros")
        caja_porticos.addWidget(subpestanias)
        return contenedor_porticos


    def _pagina_archivos(self) -> QWidget:
        pagina = QWidget()
        caja = QVBoxLayout(pagina)
        self.etiqueta_archivos = QLabel(
            "Doble clic en un archivo para abrirlo con el programa que Windows tenga asociado "
            "(Excel, LibreOffice, visor de DXF, bloc de notas…)."
        )
        self.etiqueta_archivos.setWordWrap(True)
        caja.addWidget(self.etiqueta_archivos)
        caja.addWidget(self.arbol, 1)
        return pagina

    def _conectar(self) -> None:
        self.boton_actualizar.clicked.connect(self.refrescar)
        self.combo_obras.currentIndexChanged.connect(self._seleccionar_obra_ui)
        self.boton_nueva_obra.clicked.connect(self._crear_obra_ui)
        self.boton_eliminar_obra.clicked.connect(self._eliminar_obra_ui)
        self.boton_guardar_proyecto.clicked.connect(self._guardar_ficha_proyecto)
        self.boton_crear_portico.clicked.connect(self._crear_portico)
        self.boton_eliminar_portico.clicked.connect(self._eliminar_portico)
        self.boton_crear_muro.clicked.connect(self._ingresar_muros)
        self.boton_definir_ejes.clicked.connect(self._editar_ejes)
        self.boton_definir_niveles.clicked.connect(self._editar_niveles)
        self.boton_asociar_porticos.clicked.connect(self._asociar_porticos)
        self.boton_inicio_cargas.clicked.connect(self._ir_a_definir_cargas)
        self.boton_inicio_aplicaciones.clicked.connect(self._ir_a_aplicaciones_cargas)
        self.boton_inicio_losas.clicked.connect(self._nueva_losa_desde_inicio)
        self.boton_referencias_cargas.clicked.connect(self._mostrar_referencias_cargas)
        self.boton_ver_cargas.clicked.connect(self._abrir_ultimo_analisis)
        self.boton_resolver_motor.clicked.connect(self._resolver_motor)
        self.boton_dimensionar_vigas.clicked.connect(self._dimensionar_vigas)
        self.boton_dimensionar_columnas.clicked.connect(self._dimensionar_columnas)
        self.boton_dimensionar_bases.clicked.connect(self._dimensionar_bases)
        self.boton_ver_motor.clicked.connect(self._abrir_json_motor)
        self.boton_diagrama_momentos.clicked.connect(self._mostrar_diagrama_momentos)
        self.combo_portico.currentTextChanged.connect(lambda _: self.refrescar())
        self.combo_portico.currentTextChanged.connect(self.pagina_cargas.actualizar_portico)
        self.combo_portico.currentTextChanged.connect(self.vista_planta.set_portico_activo)
        self.combo_nivel_planta.currentTextChanged.connect(self._cambiar_nivel_planta)
        self.tabla_etapas.itemSelectionChanged.connect(self._al_elegir_etapa)
        self.boton_ejecutar.clicked.connect(self._ejecutar_etapa)
        self.boton_abrir_salida.clicked.connect(self._abrir_salida_etapa)
        self.lista_losas.currentItemChanged.connect(self._al_elegir_losa)
        self.lista_planillas_vigas.currentItemChanged.connect(self._al_elegir_planilla_viga)
        self.lista_planillas_columnas.currentItemChanged.connect(
            self._al_elegir_planilla_columnas
        )
        self.lista_planillas_bases.currentItemChanged.connect(
            self._al_elegir_planilla_bases
        )
        self.arbol.itemDoubleClicked.connect(self._abrir_del_arbol)

    # ------------------------------------------------------------------
    # Refresco de la información
    # ------------------------------------------------------------------
    def _portico(self) -> str:
        return self.combo_portico.currentText().strip()


    def _ingresar_muros(self):
        dialogo = EditorMuroPorNiveles(self.ejes, self.niveles, self.estructura, self)
        if dialogo.exec() == QDialog.DialogCode.Accepted:
            datos = dialogo.datos_muro()
            if not datos:
                return

            # Generar nombre único contando solo muros existentes, para no
            # chocar con la numeración independiente de los pórticos.
            muros_existentes = {
                int(m.group(1))
                for nombre_e in self.estructura
                if (m := re.search(r"^Muro\s*(\d+)\s*$", str(nombre_e)))
            }
            siguiente = max(muros_existentes, default=0) + 1
            nombre = f"Muro {siguiente}"

            self.estructura[nombre] = datos
            rutas.guardar_estructura(self.estructura)
            self.refrescar()
            QMessageBox.information(self, "Muro guardado", f"Se guardó el muro {nombre}.")

    def _editar_ejes(self) -> None:
        dialogo = EditorEjes(cargar_ejes(), self)
        if dialogo.exec() == QDialog.DialogCode.Accepted:
            self.refrescar()

    def _editar_niveles(self) -> None:
        dialogo = EditorNiveles(cargar_niveles(), self)
        if dialogo.exec() == QDialog.DialogCode.Accepted:
            self.refrescar()

    def _asociar_porticos(self) -> None:
        estructura = rutas.cargar_estructura()
        if not estructura:
            self._aviso("Sin pórticos", "Creá al menos un pórtico antes de asociarlo a ejes.")
            return
        dialogo = EditorAsociacionPorticos(estructura, cargar_ejes(), self)
        if dialogo.exec() == QDialog.DialogCode.Accepted:
            self.refrescar()

    def _rel(self, ruta) -> str:
        """Muestra la ruta recortada respecto de la carpeta del proyecto."""
        try:
            return str(Path(ruta).relative_to(rutas.RAIZ))
        except ValueError:
            return str(ruta)

    def _cambiar_nivel_planta(self, nombre_nivel: str) -> None:
        """Filtra la vista de planta al nivel seleccionado."""
        if str(nombre_nivel).strip().casefold() == "todos los niveles":
            self.vista_planta.set_nivel_activo(None)
            return
        cota = None
        for n in self.niveles:
            if str(n.get("nombre", "")) == str(nombre_nivel):
                try:
                    cota = float(n.get("cota_m", 0.0))
                except (TypeError, ValueError):
                    cota = None
                break
        self.vista_planta.set_nivel_activo(cota)

    def _actualizar_tabla_portico_resumen(self) -> None:
        """Llena la tabla resumen con el número de elementos del pórtico activo."""
        portico = self._portico()
        if not portico:
            self.tabla_portico_resumen.setRowCount(0)
            return
        estructura = rutas.cargar_estructura() or {}
        datos = estructura.get(portico, {}) or {}
        columnas = datos.get("columnas", {}) or {}
        vigas = datos.get("vigas", {}) or {}
        bases = datos.get("bases", {}) or {}
        filas: list[tuple[str, str, str]] = [
            ("Columnas", str(len(columnas)), "N<nivel>-<columna>"),
            ("Vigas", str(len(vigas)), "V<nivel>-N T<tramo>"),
            ("Bases", str(len(bases)), "BASE-XXX (global)"),
        ]
        for nombre, cantidad, detalle in filas:
            fila = self.tabla_portico_resumen.rowCount()
            self.tabla_portico_resumen.insertRow(fila)
            self.tabla_portico_resumen.setItem(fila, 0, QTableWidgetItem(nombre))
            self.tabla_portico_resumen.setItem(fila, 1, QTableWidgetItem(cantidad))
            self.tabla_portico_resumen.setItem(fila, 2, QTableWidgetItem(detalle))

    def refrescar(self) -> None:
        obras = rutas.listar_obras()
        activa = rutas.CARPETA_OBRA
        # Recargar ejes y niveles de la obra activa (cambian al cambiar de obra)
        from app.ejes import cargar_ejes, cargar_niveles
        self.ejes = cargar_ejes()
        self.niveles = cargar_niveles()
        # Releer la estructura desde disco: si se creó, editó o eliminó un
        # pórtico/muro, la planta tiene que mostrar el estado real y no la
        # copia que quedó en memoria.
        self.estructura = rutas.cargar_estructura()
        self.vista_planta.set_datos(self.ejes, self.niveles, self.estructura)
        # Datos de viento para dibujar la flecha direccional en la planta
        try:
            proyecto = rutas.leer_json(rutas.CARGAS, {}) or {}
            configuracion_viento = proyecto.get("viento", {}) or {}
            direccion = str(configuracion_viento.get("direccion_viento", "x")).strip().lower()
            if direccion not in ("x", "y"):
                direccion = "x"
            activo = bool(configuracion_viento.get("activo"))
            lineal = 0.0
            ancho = 0.0
            nodos = []
            if activo:
                try:
                    from calc import materiales
                    items_viento = cargas.items_viento_general(proyecto, materiales.cargar())
                    lineal = sum(float(item["valor"]) for item in items_viento)
                    ancho = float(configuracion_viento.get("ancho_tributario_m", 0.0))
                    columnas_por_piso = {}
                    for _cid, _col in self.estructura.items():
                        if not isinstance(_col, dict):
                            continue
                        _clave = round(float(_col.get("nivel", 0.0)), 4)
                        columnas_por_piso.setdefault(_clave, []).append(_col)
                    for _columnas in columnas_por_piso.values():
                        if not _columnas:
                            continue
                        _altura = float(_columnas[0].get("altura_m", 0.0))
                        _fuerza = lineal * _altura / len(_columnas)
                        for _col in _columnas:
                            nodos.append({
                                "x": float(_col.get("x", 0.0)),
                                "y": float(_col.get("nivel", 0.0)) + float(_col.get("altura_m", 0.0)),
                                "valor_kN": _fuerza,
                            })
                except Exception:
                    nodos = []
            self.vista_planta.set_datos_viento(nodos, direccion, lineal, ancho)
        except Exception:
            pass
        # El selector de nivel se repuebla en _cargar_inicio (único lugar que
        # lo gobierna) para no duplicar la lista de niveles.
        self.combo_obras.blockSignals(True)
        self.combo_obras.clear()
        for carpeta in obras:
            self.combo_obras.addItem(carpeta.name, str(carpeta))
        if activa:
            self.combo_obras.setCurrentIndex(max(self.combo_obras.findData(str(activa)), 0))
        self.combo_obras.blockSignals(False)
        porticos = rutas.listar_porticos()
        elegido = self._portico()
        self.combo_portico.blockSignals(True)
        self.combo_portico.clear()
        self.combo_portico.addItems(porticos)
        if elegido in porticos:
            self.combo_portico.setCurrentText(elegido)
        self.combo_portico.blockSignals(False)
        self.boton_eliminar_portico.setEnabled(bool(porticos))

        portico = self._portico()

        self.vista_portico.establecer_niveles(self.niveles)
        self.vista_portico.actualizar(portico)
        self.pagina_cargas.actualizar_portico()
        if porticos:
            self.etiqueta_resumen.setText(resumen_etapas(portico))
        else:
            self.etiqueta_resumen.setText(
                "No hay pórticos cargados en datos/estructura.json (empezá por la etapa 3, geometría)."
            )

        self._cargar_inicio(portico)
        self._cargar_etapas(portico)
        self._cargar_vigas(portico)
        self._cargar_columnas(portico)
        self._cargar_bases(portico)
        self._cargar_losas()
        self.pagina_muros.refrescar()
        self._cargar_arbol()
        self.statusBar().showMessage(
            f"{len(porticos)} pórtico(s) en {rutas.CARPETA_OBRA or rutas.RAIZ}"
        )
        self._actualizar_tabla_portico_resumen()

    def _seleccionar_obra_ui(self, indice: int) -> None:
        ruta = self.combo_obras.itemData(indice)
        if not ruta or Path(ruta).resolve() == rutas.CARPETA_OBRA:
            return
        try:
            rutas.seleccionar_obra(ruta)
            importlib.reload(pipeline)
            self.pagina_cargas.recargar()
            self.refrescar()
        except (OSError, ValueError) as exc:
            QMessageBox.critical(self, "No se pudo abrir la obra", str(exc))

    def _crear_obra_ui(self) -> None:
        nombre, ok = QInputDialog.getText(self, "Nueva obra", "Nombre de la obra:")
        if not ok or not nombre.strip():
            return
        try:
            carpeta = rutas.crear_obra(nombre.strip())
            rutas.seleccionar_obra(carpeta)
            importlib.reload(pipeline)
            self.pagina_cargas.recargar()
            self.refrescar()
            self.statusBar().showMessage(f"Obra creada en {carpeta}", 7000)
        except (OSError, ValueError) as exc:
            QMessageBox.critical(self, "No se pudo crear la obra", str(exc))

    def _eliminar_obra_ui(self) -> None:
        carpeta = rutas.CARPETA_OBRA
        if carpeta is None:
            QMessageBox.information(self, "Eliminar obra", "No hay obra activa.")
            return
        try:
            carpeta = Path(carpeta).resolve()
            obras_raiz = rutas.OBRAS.resolve()
            if obras_raiz not in carpeta.parents:
                raise ValueError("La carpeta seleccionada no es una obra válida.")
        except Exception as exc:
            QMessageBox.critical(self, "Eliminar obra", str(exc))
            return
        respuesta = QMessageBox.question(
            self,
            "Eliminar obra",
            f"¿Eliminar definitivamente la obra '{carpeta.name}' y todos sus datos/salidas?\n\nEsta acción no se puede deshacer.",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        if respuesta != QMessageBox.StandardButton.Yes:
            return
        try:
            import shutil
            if carpeta.exists():
                shutil.rmtree(carpeta)
            # si era la activa, seleccionar otra o ninguna
            obras = rutas.listar_obras()
            if obras:
                rutas.seleccionar_obra(obras[0])
            else:
                rutas.seleccionar_obra(None)
            importlib.reload(pipeline)
            self.pagina_cargas.recargar()
            self.refrescar()
            self.statusBar().showMessage(f"Obra eliminada: {carpeta.name}", 7000)
        except Exception as exc:
            QMessageBox.critical(self, "No se pudo eliminar la obra", str(exc))

    # ------------------------------------------------------------------
    # Pestaña Inicio
    # ------------------------------------------------------------------
    def _pagina_inicio(self) -> QWidget:
        pagina = QWidget()
        caja = QVBoxLayout(pagina)

        grupo_proyecto = QGroupBox("Carátula de la obra")
        form_proyecto = QGridLayout(grupo_proyecto)
        form_proyecto.addWidget(QLabel("Nombre:"), 0, 0)
        form_proyecto.addWidget(self.nombre_proyecto, 0, 1, 1, 2)
        form_proyecto.addWidget(QLabel("ID:"), 0, 3)
        form_proyecto.addWidget(self.id_proyecto_inicio, 0, 4)
        form_proyecto.addWidget(self.boton_guardar_proyecto, 0, 5)
        form_proyecto.addWidget(QLabel("Lugar:"), 1, 0)
        form_proyecto.addWidget(self.ubicacion_proyecto, 1, 1, 1, 2)
        form_proyecto.addWidget(QLabel("Comitente:"), 1, 3)
        form_proyecto.addWidget(self.comitente_proyecto, 1, 4, 1, 2)
        form_proyecto.addWidget(QLabel("Responsable:"), 2, 0)
        form_proyecto.addWidget(self.responsable_proyecto, 2, 1, 1, 2)
        form_proyecto.addWidget(QLabel("Notas:"), 2, 3)
        form_proyecto.addWidget(self.notas_proyecto, 2, 4, 1, 2)
        caja.addWidget(grupo_proyecto)

        self.etiqueta_inicio.setStyleSheet("font-size: 11pt;")

        grupo_muros = QGroupBox("Muros de esta obra")
        caja_muros = QVBoxLayout(grupo_muros)
        self.combo_muro = QComboBox()
        self.combo_muro.setMinimumWidth(190)
        #caja_muros.addWidget(QLabel("Muro activo:"))
        #caja_muros.addWidget(self.combo_muro)
        caja_muros.addWidget(self.boton_crear_muro)
        # grupo_muros se reubica en el panel lateral, justo debajo de los ejes
        # caja_porticos.addStretch(1)         # selector de pórtico en la pestaña "Pórticos"
        # caja_porticos.addWidget(grupo_muros)  # selector de pórtico en la pestaña "Pórticos"
        # caja_muros.addWidget(self.boton_crear_muro)  # duplicado, eliminado
        # caja_porticos.addStretch(1)          # selector de pórtico en la pestaña "Pórticos"
        # caja_porticos.addWidget(grupo_muros)   # selector de pórtico en la pestaña "Pórticos"

        fila_accesos = QVBoxLayout()
        fila_accesos.setSpacing(4)
        fila_accesos.addWidget(self.boton_inicio_cargas)
        fila_accesos.addWidget(self.boton_inicio_aplicaciones)
        fila_accesos.addWidget(self.boton_inicio_losas)
        grupo_accesos = QGroupBox("Accesos rápidos")
        grupo_accesos.setLayout(fila_accesos)

        grupo_cargas = QGroupBox("Cargas del proyecto")
        caja_cargas = QVBoxLayout(grupo_cargas)
        caja_cargas.addWidget(self.lbl_cargas)
        fila_cargas = QHBoxLayout()
        fila_cargas.addWidget(self.boton_ver_cargas)
        fila_cargas.addStretch(1)
        caja_cargas.addLayout(fila_cargas)
        panel_lateral = QWidget()
        panel_lateral.setMaximumWidth(420)
        panel_lateral.setMinimumWidth(300)
        estados = QVBoxLayout(panel_lateral)
        estados.setContentsMargins(0, 0, 0, 0)
        estados.setSpacing(5)
        estados.addWidget(self.etiqueta_inicio)
        # estados.addWidget(grupo_porticos)  # Selector de pórtico en la pestaña "Pórticos"
        estados.addWidget(self.resumen_ejes_inicio)
        estados.addWidget(self.resumen_niveles_inicio)
        referencias = QVBoxLayout()
        referencias.setSpacing(4)
        for boton in (
            self.boton_definir_ejes,
            self.boton_definir_niveles,
            self.boton_asociar_porticos,
        ):
            boton.setMinimumHeight(32)
            boton.setStyleSheet("QPushButton { padding: 4px 8px; }")
            referencias.addWidget(boton)
        estados.addLayout(referencias)
        self.boton_crear_portico.setMinimumHeight(34)
        self.boton_crear_portico.setStyleSheet("QPushButton { padding: 5px 9px; font-weight: 600; }")
        estados.addWidget(self.boton_crear_portico)
        estados.addWidget(grupo_muros)   # "Crear muro", justo debajo de los ejes
        estados.addWidget(self.referencia_portico_inicio)
        estados.addWidget(grupo_accesos)
        estados.addWidget(grupo_cargas)
        estados.addStretch(1)
        panel_lateral_scroll = QScrollArea()
        panel_lateral_scroll.setWidgetResizable(True)
        panel_lateral_scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff
        )
        panel_lateral_scroll.setWidget(panel_lateral)
        panel_lateral_scroll.setMinimumWidth(320)
        panel_lateral_scroll.setMaximumWidth(440)

        self.boton_inicio_cargas.setText("Cargas")
        self.boton_inicio_cargas.setToolTip("Definir o editar las cargas del proyecto.")
        self.boton_inicio_aplicaciones.setText("Asignar cargas a barras")
        self.boton_inicio_aplicaciones.setToolTip(
            "Ir a las aplicaciones manuales de cargas dentro de Cargas."
        )
        self.boton_inicio_losas.setText("Losa desde ejes…")
        self.boton_inicio_losas.setToolTip(
            "Delimitar un paño por ejes, elegir las vigas de apoyo y definir sus cargas."
        )
        for boton in (
            self.boton_inicio_cargas,
            self.boton_inicio_aplicaciones,
            self.boton_inicio_losas,
        ):
            boton.setMinimumHeight(34)
            boton.setStyleSheet("QPushButton { padding: 5px 8px; font-weight: 600; }")

        self.grupo_planta_inicio = QGroupBox("Planta de la obra")
        caja_planta = QVBoxLayout(self.grupo_planta_inicio)

        # Barra de herramientas superior: solo el selector de nivel activo.
        fila_nivel = QHBoxLayout()
        fila_nivel.addWidget(QLabel("Nivel visible:"))
        self.combo_nivel_planta.addItems(["Todos los niveles"])
        fila_nivel.addWidget(self.combo_nivel_planta, 1)
        fila_nivel.addStretch(1)
        caja_planta.addLayout(fila_nivel)

        caja_planta.addWidget(self.scroll_vista_portico_p, 1)

        self.division_inicio = QSplitter(Qt.Orientation.Horizontal)
        self.division_inicio.addWidget(self.grupo_planta_inicio)
        self.division_inicio.addWidget(panel_lateral_scroll)
        self.division_inicio.setStretchFactor(0, 3)
        self.division_inicio.setStretchFactor(1, 1)
        self.division_inicio.setSizes([900, 360])
        self.division_inicio.setChildrenCollapsible(False)
        self.division_inicio.setMinimumHeight(340)
        caja.addWidget(self.division_inicio, 5)
        return pagina

    def _cargar_inicio(self, portico: str) -> None:
        import html

        self.grupo_planta_inicio.setTitle("Planta de la obra")
        self.ejes = cargar_ejes()
        self.niveles = cargar_niveles()
        self.resumen_ejes_inicio.setText(resumen_ejes(self.ejes))
        self.resumen_niveles_inicio.setText(resumen_niveles(self.niveles))
        proyecto = rutas.leer_json(rutas.CARGAS, {}) or {}
        ficha = (
            rutas.leer_json(rutas.CARPETA_OBRA / "obra.json", {}) or {}
            if rutas.CARPETA_OBRA else {}
        )
        self.nombre_proyecto.setText(str(proyecto.get("obra", "Obra")))
        self.id_proyecto_inicio.setText(PaginaCargas._id_proyecto(proyecto))
        ubicacion_viento = proyecto.get("viento", {}).get("referencia_cirsoc_102_25", {}).get("ubicacion", {})
        ubicacion_heredada = ", ".join(
            str(ubicacion_viento.get(k, "")).strip() for k in ("ciudad", "provincia", "pais")
            if str(ubicacion_viento.get(k, "")).strip()
        )
        self.ubicacion_proyecto.setText(str(
            ficha.get("ubicacion_obra", "")
            or ficha.get("lugar", "")
            or ubicacion_heredada
        ))
        self.comitente_proyecto.setText(str(ficha.get("comitente", "")))
        self.responsable_proyecto.setText(str(ficha.get("responsable", "")))
        self.notas_proyecto.setPlainText(str(proyecto.get("notas", "")))
        # Visor de planta general de la obra: muestra TODOS los ejes y niveles,
        # no el foco de un pórtico particular.
        self.vista_planta.set_portico_activo("")
        self.vista_planta.set_datos(self.ejes, self.niveles, self.estructura)
        self.combo_nivel_planta.blockSignals(True)
        nivel_previo = self.combo_nivel_planta.currentText()
        self.combo_nivel_planta.clear()
        opciones_nivel = ["Todos los niveles"]
        opciones_nivel += [
            str(n.get("nombre", "")) for n in self.niveles if str(n.get("nombre", ""))
        ]
        self.combo_nivel_planta.addItems(opciones_nivel)
        # Conservar el nivel que ya estaba visible: refrescar la planta no
        # debe cambiar el filtro elegido por el usuario.
        if nivel_previo in opciones_nivel:
            self.combo_nivel_planta.setCurrentText(nivel_previo)
        else:
            self.combo_nivel_planta.setCurrentIndex(0)
        self.combo_nivel_planta.blockSignals(False)
        self._cambiar_nivel_planta(self.combo_nivel_planta.currentText())
        cargas_datos = rutas.leer_json(rutas.CARGAS, {}) or {}
        elementos_carga = cargas_datos.get("elementos", {}) or {}
        self.vista_planta.set_losas(elementos_carga)
        # Aplicaciones de cargas.json: flechas de muro perpendicular entre pórticos.
        self.vista_planta.set_aplicaciones(cargas_datos.get("aplicaciones", []))
        estructura_completa = rutas.cargar_estructura()
        estructura = estructura_completa.get(portico, {}) if portico else {}
        referencia = estructura.get("referencia_planta", {}) or {}
        direccion = str(referencia.get("direccion", "x")).upper()
        eje_id = referencia.get("eje_id")
        eje_longitudinal_id = referencia.get("eje_longitudinal_id")
        ejes_por_id = {
            str(item.get("id", f"{familia}:{str(item['nombre']).casefold()}")): item
            for familia, items in self.ejes.items()
            for item in items
        }
        eje = ejes_por_id.get(str(eje_id)) if eje_id else None
        eje_longitudinal = (
            ejes_por_id.get(str(eje_longitudinal_id)) if eje_longitudinal_id else None
        )
        if eje_longitudinal_id and eje_longitudinal is None:
            detalle_primera_columna = " · eje de primera columna no disponible"
        elif eje_longitudinal is not None:
            inicio_longitudinal = float(eje_longitudinal["coordenada_m"])
            detalle_primera_columna = (
                f" · primera columna: eje {eje_longitudinal['nombre']} "
                f"({inicio_longitudinal:.3f} m)"
            )
        else:
            detalle_primera_columna = " · primera columna en origen local"
        if eje_id and eje is None:
            self.referencia_portico_inicio.setText(
                f"Referencia del pórtico: paralelo a {direccion} · referencia a eje "
                "no disponible. Revisá la asociación antes de usar anchos automáticos."
                f"{detalle_primera_columna}."
            )
        elif eje is None:
            self.referencia_portico_inicio.setText(
                f"Referencia del pórtico: paralelo a {direccion} · sin eje asociado "
                "(pórtico legado sin asociación)."
                f"{detalle_primera_columna}."
                if portico else "Seleccioná un pórtico para ver su referencia de planta."
            )
        else:
            self.referencia_portico_inicio.setText(
                f"Referencia del pórtico: paralelo a {direccion} · eje {eje['nombre']} "
                f"({float(eje['coordenada_m']):.3f} m)"
                f"{detalle_primera_columna}."
            )
        elementos = proyecto.get("elementos", {})
        cargas_activas = sum(bool(e.get("activo", True)) for e in elementos.values())
        aplicaciones_activas = sum(
            bool(a.get("activa", True)) for a in proyecto.get("aplicaciones", [])
        )

        def _transferida_a_portico(elemento: dict) -> bool:
            apoyos = elemento.get("apoya_en") or {}
            return (
                all(
                    (apoyos.get(lado) or {}).get("portico")
                    and (apoyos.get(lado) or {}).get("viga")
                    for lado in ("izq", "der")
                )
                and any(
                    (apoyos.get(lado) or {}).get("portico") == portico
                    for lado in ("izq", "der")
                )
            )

        losas_transferidas = sum(
            1
            for elemento in elementos.values()
            if (
                portico
                and elemento.get("tipo") == "losa"
                and elemento.get("activo", True)
                and elemento.get("panel_ejes")
                and _transferida_a_portico(elemento)
            )
        )
        geometria = pipeline.estado_etapa("geometria", portico)
        estado_cargas = pipeline.estado_etapa("cargas", portico)
        estado_motor = pipeline.estado_etapa("portico", portico)
        if portico:
            n_columnas = len(estructura.get("columnas", {}))
            n_tramos = sum(len(v.get("tramos", [])) for v in estructura.get("vigas", {}).values())
            seleccionado = (
                f"Pórtico seleccionado: <b>{html.escape(portico)}</b> "
                f"({n_columnas} columnas, {n_tramos} tramos)"
            )
            proximo = self._siguiente_paso(geometria, estado_cargas, estado_motor)
        else:
            seleccionado = "No hay pórticos cargados en datos/estructura.json."
            proximo = "Cargá la geometría de un pórtico para aplicar cargas y resolverlo."
        self.etiqueta_inicio.setText(
            f"{seleccionado}<br><b>Próximo paso:</b> {html.escape(proximo)}<br>"
            f"{html.escape(resumen_etapas(portico)) if portico else 'Estado: falta geometría'}"
        )

        iconos = pipeline.ICONOS
        nombres_estado = self.NOMBRE_ESTADO
        informes = sorted(
            rutas.listar(rutas.SAL_ANALISIS_CARGAS, "*.txt"),
            key=lambda archivo: archivo.stat().st_mtime,
            reverse=True,
        )
        self.ultimo_analisis = informes[0] if informes else None
        estado_txt = (
            f"{iconos.get(estado_cargas['estado'], '')} "
            f"{nombres_estado.get(estado_cargas['estado'], estado_cargas['estado'])}"
        )
        detalle_txt = estado_cargas["detalle"]
        if informes:
            detalle_txt += f"\nÚltimo TXT: {self._rel(informes[0])}"
        self.lbl_cargas.setText(
            f"{estado_txt} \u00b7 {detalle_txt}\n"
            f"Configuración: {cargas_activas} carga(s) activa(s), "
            f"{aplicaciones_activas} aplicación(es) manual(es) a tramos y "
            f"{losas_transferidas} losa(s) aplicada(s) automáticamente a este pórtico."
        )
        self.boton_ver_cargas.setEnabled(self.ultimo_analisis is not None)

        archivo, datos, filas = datos_solicitaciones(portico)
        diagrama_momentos = (datos or {}).get("diagrama_momentos")
        hay_diagrama = bool(diagrama_momentos and diagrama_momentos.get("barras"))
        self.boton_diagrama_momentos.setEnabled(hay_diagrama)
        estado_motor_txt = (
            f"{iconos.get(estado_motor['estado'], '')} "
            f"{nombres_estado.get(estado_motor['estado'], estado_motor['estado'])} \u00b7 "
            f"{estado_motor['detalle']}"
        )
        if archivo is not None and datos:
            generado = datos.get("generado", "sin fecha")
            self.lbl_motor.setText(
                f"{estado_motor_txt}\nResultado: {self._rel(archivo)} \u00b7 generado {generado} \u00b7 "
                f"{len(datos.get('combinaciones', []))} combinaciones."
            )
            self.boton_ver_motor.setEnabled(True)
            self.boton_ver_motor.setText(
                "Abrir resultado del motor" if estado_motor["estado"] == "ok"
                else "Abrir resultado anterior"
            )
        else:
            self.lbl_motor.setText(
                f"{estado_motor_txt}\nTodavía no hay resultado guardado para este pórtico."
            )
            self.boton_ver_motor.setEnabled(False)
            self.boton_ver_motor.setText("Abrir el JSON del motor")
        # El motor lee cargas.json directamente; el TXT del análisis es un informe,
        # no una entrada necesaria para resolver el pórtico.
        requisitos_motor = geometria["estado"] == "ok"
        motor_al_dia = estado_motor["estado"] == "ok"
        self.boton_resolver_motor.setEnabled(
            bool(portico) and bool(estructura) and requisitos_motor
            and (not motor_al_dia or not hay_diagrama)
        )
        if motor_al_dia and hay_diagrama:
            self.boton_resolver_motor.setText("Solicitaciones al día")
            self.boton_resolver_motor.setToolTip("El resultado vigente está guardado; podés abrirlo con el botón de al lado.")
        elif motor_al_dia:
            self.boton_resolver_motor.setText("Actualizar para ver diagrama")
            self.boton_resolver_motor.setToolTip(
                "Este resultado anterior no incluye el muestreo necesario para dibujar momentos."
            )
        elif not requisitos_motor:
            faltan = []
            if geometria["estado"] != "ok":
                faltan.append("Geometría del pórtico")
            self.boton_resolver_motor.setText("Completar etapas previas")
            self.boton_resolver_motor.setToolTip("Primero completá: " + " y ".join(faltan))
        else:
            self.boton_resolver_motor.setText(
                "Recalcular solicitaciones" if estado_motor["estado"] == "desactualizada"
                else "Resolver solicitaciones del pórtico"
            )
            self.boton_resolver_motor.setToolTip(
                "El motor usa las cargas guardadas en esta obra; no hace falta generar el TXT antes."
            )
        mostrar_envolvente = estado_motor["estado"] in ("ok", "pendiente")
        self._llenar(self.tabla_envolvente, filas if mostrar_envolvente else [])

    def _guardar_ficha_proyecto(self) -> None:
        datos = cargas.datos_cargas()
        datos["obra"] = self.nombre_proyecto.text().strip() or "Obra"
        datos["notas"] = self.notas_proyecto.toPlainText().strip()
        datos["id_proyecto"] = PaginaCargas._id_proyecto(datos)
        rutas.guardar_json(rutas.CARGAS, datos)
        ficha = rutas.leer_json(rutas.CARPETA_OBRA / "obra.json", {}) or {}
        ficha.update({
            "nombre": datos["obra"], "id": datos["id_proyecto"],
            "ubicacion_obra": self.ubicacion_proyecto.text().strip(),
            "comitente": self.comitente_proyecto.text().strip(),
            "responsable": self.responsable_proyecto.text().strip(),
            "ubicacion": datos.get("viento", {}).get("referencia_cirsoc_102_25", {}).get("ubicacion", {}),
        })
        rutas.guardar_json(rutas.CARPETA_OBRA / "obra.json", ficha)
        self.pagina_cargas.recargar()
        self.statusBar().showMessage("Ficha del proyecto guardada en obra.json.", 5000)
        self.refrescar()

    def _ir_a_pestania(self, nombre: str) -> None:
        pestanias = self.findChild(QTabWidget)
        if pestanias is None:
            return
        for indice in range(pestanias.count()):
            if pestanias.tabText(indice) == nombre:
                pestanias.setCurrentIndex(indice)
                return

    def _pedir_referencia_planta_portico(self) -> dict | None:
        niveles = cargar_niveles()
        if len(niveles) < 2:
            QMessageBox.warning(
                self, "Faltan niveles",
                "Definí al menos dos niveles Z antes de crear un pórtico por niveles.",
            )
            return None
        ejes = cargar_ejes()
        if not any(ejes.values()):
            QMessageBox.warning(
                self, "Faltan ejes",
                "Definí los ejes X/Y antes de crear un pórtico por niveles.",
            )
            return None
        dialogo = EditorPorticoPorNiveles(
            ejes, niveles, rutas.cargar_estructura(), self
        )
        if dialogo.exec() != QDialog.DialogCode.Accepted:
            return None
        try:
            return dialogo.geometria()
        except ValueError as exc:
            QMessageBox.warning(self, "Geometría incompleta", str(exc))
            return None

    def _siguiente_numero_portico_disponible(self, estructura: dict | None = None) -> int:
        """Proponé el siguiente número de pórtico contando solo pórticos reales.

        Los muros estructurales (tipo == \"muro\") se excluyen para que su
        numeración no desplace la de los pórticos.
        """
        estructura = estructura or rutas.cargar_estructura()
        usados: set[int] = set()
        for nombre, datos in estructura.items():
            if not isinstance(datos, dict):
                continue
            if str(datos.get("tipo", "")).lower() == "muro":
                continue
            m = re.search(r"(\d+)\s*$", str(nombre))
            if m:
                usados.add(int(m.group(1)))
        return max(usados, default=0) + 1

    def _crear_portico(self) -> None:
        """Crea geometría básica desde Inicio, dentro de la obra actual."""
        estructura = rutas.cargar_estructura()
        siguiente_numero = self._siguiente_numero_portico_disponible(estructura)
        numero_ingresado, ok = QInputDialog.getInt(
            self,
            "Nuevo pórtico",
            "Número acumulativo del pórtico:",
            siguiente_numero,
            1,
            2147483647,
        )
        if not ok:
            return
        nombre = f"Pórtico {numero_ingresado}"
        if nombre in estructura:
            QMessageBox.warning(self, "Nombre existente", f"Ya existe el pórtico {nombre}.")
            return
        try:
            numero_portico = rutas.numero_portico(str(numero_ingresado), estructura)
        except ValueError as exc:
            QMessageBox.warning(self, "Número de pórtico repetido", str(exc))
            return

        ubicacion = self._pedir_referencia_planta_portico()
        if ubicacion is None:
            return
        datos = {
            "vigas": {}, "columnas": {}, "bases": {}, "cargas_puntuales": [],
            "referencia_planta": ubicacion["referencia_planta"],
            "posicion_planta_m": ubicacion["posicion_planta_m"],
            "niveles_geometria": ubicacion["niveles"],
        }
        conteo_transferencias = 0
        ids_columnas_globales: list[str] = []
        ids_bases_globales: list[str] = []
        ref = ubicacion["referencia_planta"]
        direccion = ref.get("direccion")
        desfase_long = float(ref.get("desfase_longitudinal_m", 0.0) or 0.0)
        origen_long = float(ref.get("origen_longitudinal_m", 0.0) or 0.0)
        pos_transversal = float(ubicacion.get("posicion_planta_m", ref.get("desfase_m", 0.0)) or 0.0)
        for piso, tramo_nivel in enumerate(ubicacion["tramos_nivel"]):
            columnas_nivel = tramo_nivel["ejes_columna"]
            coordenadas = [float(eje["x_local_m"]) for eje in columnas_nivel]
            y = float(tramo_nivel["cota_inferior_m"])
            cota_superior = float(tramo_nivel["cota_superior_m"])
            altura = cota_superior - y
            for indice, (eje_columna, x) in enumerate(zip(columnas_nivel, coordenadas)):
                # PRINCIPIO 4: columnas son entidades globales de la obra,
                # identificadas por posición global (X, Y, nivel).
                if direccion == "x":
                    x_global = float(x) + origen_long + desfase_long
                    y_global = pos_transversal
                else:
                    x_global = pos_transversal
                    y_global = float(x) + origen_long + desfase_long
                gid_columna = rutas.obtener_o_crear_columna_global(
                    x=x,
                    nivel=y,
                    altura_m=altura,
                    eje_id=eje_columna["id"],
                    eje_nombre=eje_columna["nombre"],
                    nivel_inicio=tramo_nivel["nivel_inferior"],
                    estructura=estructura,
                    x_global=x_global,
                    y_global=y_global,
                    piso=piso,
                )
                ids_columnas_globales.append(gid_columna)

                if piso == 0:
                    gid_base = rutas.obtener_o_crear_base_global(
                        x=x,
                        nivel=y,
                        eje_id=eje_columna["id"],
                        eje_nombre=eje_columna["nombre"],
                        estructura=estructura,
                        x_global=x_global,
                        y_global=y_global,
                    )
                    ids_bases_globales.append(gid_base)

            nombre_nivel_viga = tramo_nivel["nivel_superior"]
            viga_id = f"V{piso}-{numero_portico}"
            tramos = []
            for indice, (x_inicio, x_fin) in enumerate(
                zip(coordenadas, coordenadas[1:]), 1
            ):
                tramos.append({
                    "id": f"{viga_id} T{indice}",
                    "longitud_m": x_fin - x_inicio,
                    "es_voladizo": False, "x_inicio": x_inicio, "x_fin": x_fin,
                    "cargas_puntuales": [],
                })
            x_final = coordenadas[-1]
            voladizos = {}
            for lado in ("izquierda", "derecha"):
                respuesta = QMessageBox.question(
                    self, f"Voladizo en {nombre_nivel_viga}",
                    f"¿La viga de {nombre_nivel_viga} tiene voladizo a la {lado}?",
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                    QMessageBox.StandardButton.No,
                )
                if respuesta != QMessageBox.StandardButton.Yes:
                    continue
                largo, ok = QInputDialog.getDouble(
                    self, f"Voladizo {lado} · {nombre_nivel_viga}",
                    f"Longitud del voladizo a la {lado} [m]:",
                    1.0, 0.1, 30.0, 2,
                )
                if not ok:
                    return
                voladizos[lado] = largo
            if "izquierda" in voladizos:
                largo = voladizos["izquierda"]
                tramos.append({
                    "id": f"{viga_id} tv_izq", "longitud_m": largo,
                    "es_voladizo": True,
                    "x_inicio": coordenadas[0] - largo, "x_fin": coordenadas[0],
                    "cargas_puntuales": [],
                })
            if "derecha" in voladizos:
                largo = voladizos["derecha"]
                tramos.append({
                    "id": f"{viga_id} tv_der", "longitud_m": largo,
                    "es_voladizo": True,
                    "x_inicio": x_final, "x_fin": x_final + largo,
                    "cargas_puntuales": [],
                })
            datos["vigas"][viga_id] = {
                "tramos": tramos,
                "nivel_nombre": nombre_nivel_viga,
                "cota_m": cota_superior,
                "voladizos_m": dict(voladizos),
            }
            if piso + 1 < len(ubicacion["tramos_nivel"]):
                x_min = min(coordenadas) - voladizos.get("izquierda", 0.0)
                x_max = max(coordenadas) + voladizos.get("derecha", 0.0)
                columnas_superiores = ubicacion["tramos_nivel"][piso + 1]["ejes_columna"]
                fuera_viga = [
                    eje for eje in columnas_superiores
                    if not (x_min - 1e-6 <= float(eje["x_local_m"]) <= x_max + 1e-6)
                ]
                if fuera_viga:
                    nombres = ", ".join(
                        f"{eje['nombre']} ({float(eje['x_local_m']):g} m)"
                        for eje in fuera_viga
                    )
                    QMessageBox.warning(
                        self, "Columna sin apoyo inferior",
                        f"Las columnas de {tramo_nivel['nivel_superior']} en {nombres} "
                        f"quedan fuera de las vigas y voladizos de {nombre_nivel_viga}. "
                        "Extendé la viga/voladizo o cambiá los cruces antes de crear el pórtico.",
                    )
                    return
                for eje in columnas_superiores:
                    if not any(
                        abs(float(eje["x_local_m"]) - x) < 1e-6 for x in coordenadas
                    ):
                        conteo_transferencias += 1

            # Nota: el ingreso manual de cargas puntuales al crear el pórtico se
            # eliminó de este flujo (voladizos inicio/fin se mantienen). Las
            # cargas puntuales se definen ahora desde la pestaña de Cargas.

        datos["columnas"] = ids_columnas_globales
        datos["bases"] = ids_bases_globales
        estructura[nombre] = datos
        rutas.guardar_estructura(estructura)
        self.refrescar()
        self.combo_portico.setCurrentText(nombre)
        self.statusBar().showMessage(f"{nombre} creado en la obra actual.", 5000)
        if conteo_transferencias:
            QMessageBox.information(
                self, "Columnas apoyadas sobre vigas",
                f"El pórtico incluye {conteo_transferencias} columna(s) de planta alta que "
                "descargan en puntos intermedios de vigas o voladizos inferiores. Esos nudos "
                "se conectan directamente en el modelo de análisis.",
            )

    def _eliminar_portico(self) -> None:
        nombre = self._portico()
        estructura = rutas.cargar_estructura() or {}
        if not nombre:
            QMessageBox.information(self, "Sin elemento", "Escogí el elemento que querés eliminar.")
            self.refrescar()
            return

        elemento = estructura.get(nombre)
        es_muro = isinstance(elemento, dict) and str(elemento.get("tipo", "")).lower() == "muro"
        if es_muro:
            respuesta = QMessageBox.question(
                self,
                "Confirmar eliminación",
                f"¿Seguro que querés eliminar el muro {nombre}?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No,
            )
            if respuesta != QMessageBox.StandardButton.Yes:
                return
            if nombre in estructura:
                del estructura[nombre]
                self.estructura = estructura
                rutas.guardar_estructura(self.estructura)
                self.refrescar()
                QMessageBox.information(self, "Muro eliminado", f"Se eliminó el muro {nombre}.")
            return
            QMessageBox.information(
                self, "Sin pórtico", "Elegí el pórtico que querés eliminar."
            )
            self.refrescar()
            return

        datos_cargas = rutas.leer_json(rutas.CARGAS, {}) or {}
        aplicaciones = datos_cargas.get("aplicaciones", [])
        afectadas = sum(
            aplicacion.get("portico") == nombre
            or any(
                reaccion.get("portico") == nombre
                for reaccion in aplicacion.get("reacciones", [])
            )
            for aplicacion in aplicaciones
        )

        losas = rutas.leer_json(rutas.LOSAS, {}) or {}
        nombres_losas_vinculadas = []
        for nombre_losa, elemento in (
            datos_cargas.get("elementos", {}) or {}
        ).items():
            if elemento.get("tipo") != "losa":
                continue
            apoyos = elemento.get("apoya_en", {}) or {}
            if any(
                (apoyo or {}).get("portico") == nombre
                for apoyo in apoyos.values()
            ):
                nombres_losas_vinculadas.append(str(nombre_losa))
        for nombre_losa, datos_losa in (losas.get("losas", {}) or {}).items():
            if any(
                (reaccion.get("apoya_en") or {}).get("portico") == nombre
                for reaccion in datos_losa.get("reacciones", [])
            ):
                nombres_losas_vinculadas.append(str(nombre_losa))
        for archivo in rutas.SAL_LOSAS.rglob("*.json"):
            resultado = rutas.leer_json(archivo, {}) or {}
            if any(
                (reaccion.get("apoya_en") or {}).get("portico") == nombre
                for reaccion in resultado.get("reacciones", [])
            ):
                nombres_losas_vinculadas.append(
                    str(resultado.get("nombre", archivo.stem))
                )
        nombres_losas_vinculadas = sorted(set(nombres_losas_vinculadas))
        if nombres_losas_vinculadas:
            QMessageBox.warning(
                self, "Pórtico usado por losas",
                f"No se puede eliminar {nombre}: está referenciado por "
                f"{', '.join(nombres_losas_vinculadas)}. Cambiá primero el destino de esas losas.",
            )
            return

        respuesta = QMessageBox.question(
            self, "Confirmar eliminación del pórtico",
            f"¿Eliminar {nombre} de la geometría de esta obra?\n\n"
            "Se quitará su geometría y sus aplicaciones de carga "
            f"({afectadas}). El catálogo común de cargas no se borra. Los resultados e informes "
            "ya generados se conservarán como históricos; al crear de nuevo el pórtico habrá "
            "que recalcularlos.",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        if respuesta != QMessageBox.StandardButton.Yes:
            return

        estructura_actualizada = dict(estructura)
        del estructura_actualizada[nombre]
        aplicaciones_actualizadas = []
        aplicaciones_cambiaron = False
        for aplicacion in aplicaciones:
            if aplicacion.get("portico") == nombre:
                aplicaciones_cambiaron = True
                continue
            reacciones = aplicacion.get("reacciones")
            if isinstance(reacciones, list):
                reacciones_restantes = [
                    reaccion for reaccion in reacciones
                    if reaccion.get("portico") != nombre
                ]
                if len(reacciones_restantes) != len(reacciones):
                    aplicaciones_cambiaron = True
                    if not reacciones_restantes:
                        continue
                    aplicacion_actualizada = dict(aplicacion)
                    aplicacion_actualizada["reacciones"] = reacciones_restantes
                    aplicaciones_actualizadas.append(aplicacion_actualizada)
                    continue
            aplicaciones_actualizadas.append(aplicacion)

        try:
            rutas.guardar_estructura(estructura_actualizada)
        except OSError as exc:
            QMessageBox.critical(
                self, "No se pudo eliminar el pórtico",
                f"No se pudo guardar la geometría actualizada: {exc}",
            )
            return
        if aplicaciones_cambiaron:
            datos_cargas_actualizados = dict(datos_cargas)
            datos_cargas_actualizados["aplicaciones"] = aplicaciones_actualizadas
            try:
                rutas.guardar_json(rutas.CARGAS, datos_cargas_actualizados)
            except OSError as exc:
                try:
                    rutas.guardar_estructura(estructura)
                except OSError as error_reversion:
                    QMessageBox.critical(
                        self, "Eliminación incompleta",
                        f"No se pudieron guardar las aplicaciones de carga ({exc}) y tampoco "
                        f"se pudo restaurar la geometría original ({error_reversion}).",
                    )
                    return
                QMessageBox.critical(
                    self, "No se pudo eliminar el pórtico",
                    f"No se pudieron actualizar las aplicaciones de carga; la geometría "
                    f"original fue restaurada. Error: {exc}",
                )
                return

        self.refrescar()
        self.statusBar().showMessage(f"{nombre} eliminado de la obra.", 5000)

    def _ingresar_cargas_puntuales_nivel(
        self, viga_id: str, tramos: list[dict], cota_m: float, nombre_nivel: str,
    ) -> bool:
        numero = 1
        x_min = min(float(tramo["x_inicio"]) for tramo in tramos)
        x_max = max(float(tramo["x_fin"]) for tramo in tramos)
        while QMessageBox.question(
            self, f"Carga puntual · {nombre_nivel}",
            f"¿Agregar una carga puntual a la viga de {nombre_nivel}?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        ) == QMessageBox.StandardButton.Yes:
            x, ok = QInputDialog.getDouble(
                self, "Ubicación de la carga puntual",
                f"Coordenada local X [m] (rango de viga: {x_min:g} a {x_max:g}):",
                max(x_min, min(x_max, 0.5 * (x_min + x_max))),
                x_min, x_max, 3,
            )
            if not ok:
                return False
            valor, ok = QInputDialog.getDouble(
                self, "Valor de la carga puntual",
                "Carga vertical hacia abajo [kN] (usar negativo para hacia arriba):",
                1.0, -100000.0, 100000.0, 3,
            )
            if not ok:
                return False
            tipo, ok = QInputDialog.getText(
                self, "Descripción de la carga puntual",
                "Descripción (por ejemplo, reacción de una viga secundaria o muro):",
                text="Carga puntual",
            )
            if not ok:
                return False
            tramo = next(
                (
                    item for item in tramos
                    if float(item["x_inicio"]) - 1e-6 <= x
                    <= float(item["x_fin"]) + 1e-6
                ),
                None,
            )
            if tramo is None:
                QMessageBox.warning(
                    self, "Carga fuera de la viga",
                    "El punto ingresado no pertenece a una viga o voladizo.",
                )
                continue
            tramo.setdefault("cargas_puntuales", []).append({
                "id": f"{viga_id}-P{numero}",
                "coordenadas": {"x": x, "y": cota_m},
                "valor_kN": valor,
                "tipo": tipo.strip() or "Carga puntual",
            })
            numero += 1
        return True

    @staticmethod
    def _siguiente_paso(geometria: dict, cargas_estado: dict, motor: dict) -> str:
        if geometria["estado"] != "ok":
            return "Completá la geometría: " + geometria["detalle"]
        if motor["estado"] != "ok":
            return "Resolvé las solicitaciones; el motor toma las cargas guardadas en esta obra."
        if cargas_estado["estado"] in ("pendiente", "desactualizada"):
            return "Las solicitaciones están listas; si necesitás el informe TXT, actualizalo desde Cargas."
        return "Revisa las solicitaciones y continúa con el dimensionado de vigas, columnas y bases."

    def _abrir_ultimo_analisis(self) -> None:
        if self.ultimo_analisis is None:
            return
        error = abrir_con_windows(self.ultimo_analisis)
        if error:
            self._aviso("No se pudo abrir", error)

    def _ir_a_definir_cargas(self) -> None:
        pestanias = self.centralWidget().findChild(QTabWidget)
        if pestanias is not None:
            pestanias.setCurrentWidget(self.pagina_cargas)
            self.pagina_cargas.ir_a_elementos()

    def _nueva_losa_desde_inicio(self) -> None:
        pestanias = self.centralWidget().findChild(QTabWidget)
        if pestanias is not None:
            pestanias.setCurrentWidget(self.pagina_cargas)
            self.pagina_cargas.pestanias_cargas.setCurrentIndex(0)
        self.pagina_cargas.nueva_losa_desde_ejes()

    def _ir_a_aplicaciones_cargas(self) -> None:
        pestanias = self.centralWidget().findChild(QTabWidget)
        if pestanias is not None:
            pestanias.setCurrentWidget(self.pagina_cargas)
            self.pagina_cargas.ir_a_aplicaciones()

    def _mostrar_referencias_cargas(self) -> None:
        dialogo = QDialog(self)
        dialogo.setWindowTitle(
            f"Referencias de cargas — {self._portico() or 'sin pórtico seleccionado'}"
        )
        dialogo.resize(620, 360)
        caja = QVBoxLayout(dialogo)
        indicacion = QLabel(
            "Solo se representan aquí las cargas aplicadas al pórtico seleccionado. "
            "Las cargas definidas en el catálogo pero aún no asignadas no aparecen en el esquema."
        )
        indicacion.setWordWrap(True)
        caja.addWidget(indicacion)

        contenido = QLabel(self.vista_portico.leyenda)
        contenido.setTextFormat(Qt.TextFormat.RichText)
        contenido.setWordWrap(True)
        contenido.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        contenido.setStyleSheet(
            "QLabel { color: #334155; background: #f8fafc; padding: 10px; }"
        )
        desplazamiento = QScrollArea()
        desplazamiento.setWidgetResizable(True)
        desplazamiento.setWidget(contenido)
        caja.addWidget(desplazamiento, 1)

        botones = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        botones.rejected.connect(dialogo.reject)
        caja.addWidget(botones)
        dialogo.exec()

    def _abrir_json_motor(self) -> None:
        archivo = buscar_archivo(rutas.SAL_SOLICITACIONES, "{portico}.json", self._portico())
        if archivo is None:
            self._aviso("Sin solicitaciones", "Todavía no hay JSON del motor para este pórtico.")
            return
        error = abrir_con_windows(archivo)
        if error:
            self._aviso("No se pudo abrir", error)

    def _mostrar_diagrama_momentos(self) -> None:
        portico = self._portico()
        _, datos, _ = datos_solicitaciones(portico)
        diagrama = (datos or {}).get("diagrama_momentos")
        if not diagrama or not diagrama.get("barras"):
            self._aviso(
                "Sin diagrama de momentos",
                "Recalculá las solicitaciones para generar los diagramas de momento.",
            )
            return
        estado = pipeline.estado_etapa("portico", portico)
        dialogo = QDialog(self)
        dialogo.setWindowTitle(f"Diagrama de momentos — {portico}")
        dialogo.resize(980, 660)
        caja = QVBoxLayout(dialogo)
        referencia_apoyos = (
            "â–² Apoyo articulado · rayado: empotramiento."
            if diagrama.get("apoyos")
            else "Recalculá las solicitaciones para mostrar los apoyos."
        )
        detalle = QLabel(
            f"Combinación más solicitante: {diagrama['combinacion']}. "
            "Hacé clic en una curva para consultar su valor."
            f"\n{referencia_apoyos} "
            "Escala normalizada; el signo mostrado se invierte respecto de Pynite."
        )
        detalle.setWordWrap(True)
        caja.addWidget(detalle)
        if estado["estado"] != "ok":
            aviso = QLabel(
                "El resultado guardado puede estar desactualizado respecto de la geometría o las cargas."
            )
            aviso.setWordWrap(True)
            caja.addWidget(aviso)
        vista = VistaDiagramaMomentos()
        vista.establecer_diagrama(diagrama)
        selector = QComboBox()
        selector.addItem("Momento", "momento")
        selector.addItem("Corte", "corte")
        hay_corte = all(barra.get("V_kN") for barra in diagrama["barras"])
        if not hay_corte:
            selector.setItemData(
                1,
                "Recalculá las solicitaciones para generar el diagrama de corte.",
                Qt.ItemDataRole.ToolTipRole,
            )
        seleccion = QLabel("Hacé clic sobre una curva para consultar su momento.")
        seleccion.setWordWrap(True)
        seleccion.setStyleSheet(
            "QLabel { color: #1e3a8a; background: #eff6ff; padding: 8px; }"
        )
        vista.solicitacion_seleccionada.connect(seleccion.setText)
        selector.currentIndexChanged.connect(
            lambda indice: (
                vista.establecer_tipo_diagrama(selector.itemData(indice)),
                seleccion.setText(
                    "Hacé clic sobre una curva para consultar su "
                    f"{'momento' if selector.itemData(indice) == 'momento' else 'corte'}."
                ),
            )
        )
        caja.addWidget(selector)
        caja.addWidget(seleccion)
        caja.addWidget(vista, 1)
        botones = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        botones.rejected.connect(dialogo.reject)
        caja.addWidget(botones)
        dialogo.exec()

    def _resolver_motor(self) -> None:
        portico = self._portico()
        if not portico:
            self._aviso("Sin pórtico", "Elegí un pórtico en la barra de arriba.")
            return
        geometria = pipeline.estado_etapa("geometria", portico)
        if geometria["estado"] != "ok":
            self._aviso(
                "Etapas previas incompletas",
                "Antes de resolver el pórtico, completá la geometría: " + geometria["detalle"],
            )
            return
        estado_motor = pipeline.estado_etapa("portico", portico)
        _, resultado_guardado, _ = datos_solicitaciones(portico)
        diagrama_guardado = (resultado_guardado or {}).get("diagrama_momentos") or {}
        tiene_diagrama = bool(
            diagrama_guardado.get("barras")
        )
        if estado_motor["estado"] == "ok" and tiene_diagrama:
            self._aviso(
                "Resultado vigente",
                "Las solicitaciones de este pórtico ya están calculadas. Abrí el resultado guardado; "
                "si cambiaste cargas o geometría, actualizá la ventana y recalculá cuando figure desactualizado.",
            )
            return
        texto = (
            "Se va a resolver el pórtico con el MOTOR (Pynite), por combinaciones:\n\n"
            f'    py -m calc.portico "{portico}" --guardar\n\n'
            f"Escribe salidas/solicitaciones/{rutas.nombre_seguro(portico)}.json"
        )
        if estado_motor["estado"] == "ok":
            texto = (
                "El resultado vigente es anterior al visor de momentos. Se volverá a resolver "
                "el pórtico para guardar los datos del diagrama; no se modifican la geometría "
                "ni las cargas.\n\n" + texto
            )
        if QMessageBox.question(
            self, "Resolver solicitaciones", texto,
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        ) != QMessageBox.StandardButton.Yes:
            return

        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        self.statusBar().showMessage(f"Resolviendo {portico} con el motor…")
        try:
            proceso = subprocess.run(
                [sys.executable, "-m", "calc.portico", portico, "--guardar"],
                cwd=str(rutas.RAIZ),
                env={**os.environ, "PYTHONIOENCODING": "utf-8"},
                capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=1800,
            )
            resultado = {
                "ok": proceso.returncode == 0, "codigo": proceso.returncode,
                "salida": proceso.stdout, "error": proceso.stderr,
            }
        except (subprocess.TimeoutExpired, OSError) as exc:
            resultado = {"ok": False, "codigo": None, "salida": "", "error": str(exc)}
        finally:
            QApplication.restoreOverrideCursor()

        self._mostrar_resultado("Solicitaciones del pórtico (motor)", resultado)
        self.refrescar()

    # ------------------------------------------------------------------
    # Pestaña 1: estado y etapas
    # ------------------------------------------------------------------
    NOMBRE_ESTADO = {
        "ok": "al día",
        "desactualizada": "recalcular",
        "pendiente": "falta calcular",
        "sin_datos": "faltan datos",
        "a_desarrollar": "a desarrollar",
    }

    def _cargar_etapas(self, portico: str) -> None:
        estados = pipeline.semaforo(portico)
        self.tabla_etapas.setRowCount(len(estados))
        for i, dato in enumerate(estados):
            fondo, texto = COLORES_ESTADO.get(dato["estado"], (QColor("#ffffff"), QColor("#101010")))
            celdas = (
                f"{pipeline.ICONOS.get(dato['estado'], '?')}  {self.NOMBRE_ESTADO.get(dato['estado'], dato['estado'])}",
                dato["nombre"],
                dato["detalle"],
                dato["script"] or "—",
            )
            for j, valor in enumerate(celdas):
                item = QTableWidgetItem(str(valor))
                item.setToolTip(str(valor))
                if j == 0:
                    item.setBackground(fondo)
                    item.setForeground(texto)
                self.tabla_etapas.setItem(i, j, item)
            self.tabla_etapas.item(i, 0).setData(Qt.ItemDataRole.UserRole, dato)
        self.tabla_etapas.resizeColumnsToContents()
        self.tabla_etapas.horizontalHeader().setStretchLastSection(True)
        self.detalle_etapa.setPlainText("")
        self.boton_ejecutar.setEnabled(False)
        self.boton_abrir_salida.setEnabled(False)
        self.boton_ejecutar.setText("Calcular etapa")
        self.boton_abrir_salida.setText("Abrir resultado")

    def _dato_etapa(self) -> dict | None:
        fila = self.tabla_etapas.currentRow()
        if fila < 0:
            return None
        item = self.tabla_etapas.item(fila, 0)
        return item.data(Qt.ItemDataRole.UserRole) if item else None

    def _primera_salida(self, dato: dict) -> Path | None:
        for ruta in dato["salidas"]:
            camino = Path(ruta)
            if camino.exists():
                return camino
        return None

    def _texto_etapa(self, dato: dict) -> str:
        lineas = [
            dato["nombre"],
            "=" * 62,
            f"Estado ......... : {self.NOMBRE_ESTADO.get(dato['estado'], dato['estado'])} — {dato['detalle']}",
            f"Qué hace ....... : {dato['descripcion']}",
            f"Programa ....... : {dato['script'] or 'todavía sin programar'}",
            f"Pide datos ..... : {'sí, por consola todavía' if dato['interactiva'] else 'no'}",
            f"Depende de ..... : {', '.join(dato['depende_de']) or '—'}",
            "",
            "Archivos de entrada:",
        ]
        lineas += [f"   - {self._rel(r)}" for r in dato["entradas"]] or ["   (ninguno)"]
        lineas += ["", "Archivos de salida:"]
        lineas += [f"   - {self._rel(r)}" for r in dato["salidas"]] or ["   (ninguno)"]
        if dato.get("nota"):
            lineas += ["", "Nota:", f"   {dato['nota']}"]
        return "\n".join(lineas)

    def _al_elegir_etapa(self) -> None:
        dato = self._dato_etapa()
        if not dato:
            return
        self.detalle_etapa.setPlainText(self._texto_etapa(dato))
        script = dato.get("script")
        estado = dato.get("estado")
        estados = {e["clave"]: e for e in pipeline.semaforo(self._portico())}
        bloqueos = [
            estados[clave]["nombre"]
            for clave in dato.get("depende_de", [])
            if clave in estados and estados[clave]["estado"] != "ok"
        ]
        esta_al_dia = estado == "ok"
        puede = (
            bool(script) and not dato["interactiva"] and not esta_al_dia
            and not bloqueos and (rutas.RAIZ / script).exists()
        )
        self.boton_ejecutar.setEnabled(puede)
        if esta_al_dia:
            self.boton_ejecutar.setText("Etapa al día")
            self.boton_ejecutar.setToolTip("El resultado vigente ya está guardado.")
        elif bloqueos:
            self.boton_ejecutar.setText("Esperando etapas previas")
            self.boton_ejecutar.setToolTip("Completá primero: " + ", ".join(bloqueos))
        elif dato["interactiva"]:
            self.boton_ejecutar.setText("Se ejecuta desde su pestaña")
            self.boton_ejecutar.setToolTip("Esta etapa se inicia desde su pestaña de trabajo.")
        else:
            self.boton_ejecutar.setText(
                "Recalcular etapa" if estado == "desactualizada" else "Calcular etapa"
            )
            self.boton_ejecutar.setToolTip("" if puede else "Esta etapa todavía no se puede ejecutar.")
        salida = self._primera_salida(dato)
        self.boton_abrir_salida.setEnabled(salida is not None)
        self.boton_abrir_salida.setText(
            "Abrir resultado vigente" if esta_al_dia
            else "Abrir resultado anterior" if salida is not None
            else "Resultado todavía no generado"
        )

    # ------------------------------------------------------------------
    # Pestañas 2 y 3: vigas y columnas
    # ------------------------------------------------------------------
    def _cargar_vigas(self, portico: str) -> None:
        archivo, filas = datos_vigas(portico)
        self._llenar(self.tabla_vigas, filas)
        obra = rutas.CARPETA_OBRA.name if rutas.CARPETA_OBRA else "sin obra"
        contexto = f"Obra: {obra}  ·  Pórtico seleccionado: {portico or 'ninguno'}"
        if archivo is None:
            self.etiqueta_vigas.setText(
                f"{contexto}\nSin resultados de vigas: falta calcular la etapa 5."
            )
            self._cargar_planillas_vigas(portico)
            return
        en_rojo = sum(1 for _, malo in filas if malo)
        aviso = (
            f"  ·  {en_rojo} con alguna verificación que NO cumple (en rojo)"
            if en_rojo
            else "  ·  todas las verificaciones cumplen"
        )
        self.etiqueta_vigas.setText(
            f"{contexto}\nFuente: {self._rel(archivo)}  ·  {len(filas)} tramo(s){aviso}"
        )
        if self.advertencias_vigas:
            self.etiqueta_vigas.setText(
                self.etiqueta_vigas.text() + "\nP02: " + " · ".join(self.advertencias_vigas)
            )
        self._cargar_planillas_vigas(portico)

    def _cargar_planillas_vigas(self, portico: str) -> None:
        seleccionada = self.lista_planillas_vigas.currentItem()
        ruta_previa = seleccionada.data(Qt.ItemDataRole.UserRole) if seleccionada else None
        self.lista_planillas_vigas.clear()
        if not portico:
            archivos = []
        else:
            archivos = rutas.listar(rutas.SAL_VIGAS, f"planilla_{portico}_*.txt")
            saneado = rutas.nombre_seguro(portico)
            if saneado != portico:
                archivos += rutas.listar(rutas.SAL_VIGAS, f"planilla_{saneado}_*.txt")
        archivos = sorted(set(archivos), key=lambda p: (p.name.lower(), p.stat().st_mtime))
        for archivo in archivos:
            item = QListWidgetItem(archivo.name)
            item.setData(Qt.ItemDataRole.UserRole, str(archivo))
            item.setToolTip(self._rel(archivo))
            self.lista_planillas_vigas.addItem(item)
        if not archivos:
            self.etiqueta_planillas_vigas.setText(
                "Planillas por tramo: todavía no hay TXT para este pórtico."
            )
            self.texto_planilla_viga.setPlainText("")
            return
        indice = next(
            (i for i in range(self.lista_planillas_vigas.count())
             if self.lista_planillas_vigas.item(i).data(Qt.ItemDataRole.UserRole) == ruta_previa),
            0,
        )
        self.lista_planillas_vigas.setCurrentRow(indice)
        self.etiqueta_planillas_vigas.setText(
            f"{len(archivos)} planilla(s) de armado por tramo · elegí una para verla acá."
        )

    def _al_elegir_planilla_viga(self, actual, _anterior=None) -> None:
        if actual is None:
            self.texto_planilla_viga.setPlainText("")
            return
        ruta = Path(actual.data(Qt.ItemDataRole.UserRole))
        self.texto_planilla_viga.setPlainText(
            rutas.leer_texto(ruta, "No se pudo leer la planilla.")
        )

    def _cargar_planillas_texto(
        self,
        lista: QListWidget,
        etiqueta: QLabel,
        visor: QPlainTextEdit,
        archivos: list[Path],
        tipo: str,
    ) -> None:
        seleccionada = lista.currentItem()
        ruta_previa = seleccionada.data(Qt.ItemDataRole.UserRole) if seleccionada else None
        lista.clear()
        archivos = sorted(set(archivos), key=lambda p: (p.name.lower(), p.stat().st_mtime))
        for archivo in archivos:
            item = QListWidgetItem(archivo.name)
            item.setData(Qt.ItemDataRole.UserRole, str(archivo))
            item.setToolTip(self._rel(archivo))
            lista.addItem(item)
        if not archivos:
            etiqueta.setText(f"{tipo}: todavía no hay TXT para este pórtico.")
            visor.setPlainText("")
            return
        indice = next(
            (i for i in range(lista.count())
             if lista.item(i).data(Qt.ItemDataRole.UserRole) == ruta_previa),
            0,
        )
        lista.setCurrentRow(indice)
        etiqueta.setText(f"{len(archivos)} archivo(s) de {tipo} · elegí uno para verlo acá.")

    def _al_elegir_planilla_columnas(self, actual, _anterior=None) -> None:
        if actual is None:
            self.texto_planilla_columnas.setPlainText("")
            return
        ruta = Path(actual.data(Qt.ItemDataRole.UserRole))
        self.texto_planilla_columnas.setPlainText(
            rutas.leer_texto(ruta, "No se pudo leer la memoria de columnas.")
        )

    def _al_elegir_planilla_bases(self, actual, _anterior=None) -> None:
        if actual is None:
            self.texto_planilla_bases.setPlainText("")
            return
        ruta = Path(actual.data(Qt.ItemDataRole.UserRole))
        self.texto_planilla_bases.setPlainText(
            rutas.leer_texto(ruta, "No se pudo leer la planilla de bases.")
        )

    def _dimensionar_vigas(self) -> None:
        portico = self._portico()
        if not portico:
            QMessageBox.information(self, "Vigas", "Elegí un pórtico primero.")
            return
        resultado = dimensionar_desde_app(self, portico)
        if resultado:
            _salida, self.advertencias_vigas = resultado
            self.refrescar()

    def _dimensionar_columnas(self) -> None:
        portico = self._portico()
        estructura = rutas.cargar_estructura()
        columnas = estructura.get(portico, {}).get("columnas", {})
        if not portico or not columnas:
            QMessageBox.information(
                self, "Columnas", "Elegí un pórtico con columnas cargadas."
            )
            return
        estado_motor = pipeline.estado_etapa("portico", portico)
        if estado_motor["estado"] != "ok":
            QMessageBox.warning(
                self, "Solicitaciones no vigentes",
                "Primero resolvé las solicitaciones del pórtico y verificá que estén al día.",
            )
            return
        respuestas = respuestas_columnas(self, portico, columnas)
        if respuestas is None:
            return
        resultado = self._ejecutar_dimensionado("columnas", portico, respuestas)
        if resultado is not None:
            self._mostrar_resultado("Columnas", resultado)
            self.refrescar()

    def _dimensionar_bases(self) -> None:
        portico = self._portico()
        estructura = rutas.cargar_estructura()
        bases = estructura.get(portico, {}).get("bases", {})
        if not portico or not bases:
            QMessageBox.information(
                self, "Bases", "Elegí un pórtico con bases cargadas."
            )
            return
        estado_motor = pipeline.estado_etapa("portico", portico)
        estado_columnas = pipeline.estado_etapa("columnas", portico)
        if estado_motor["estado"] != "ok" or estado_columnas["estado"] != "ok":
            QMessageBox.warning(
                self, "Faltan etapas previas",
                "Para dimensionar las bases, resolvé primero las solicitaciones "
                "y dimensioná las columnas del pórtico.",
            )
            return
        terreno = rutas.leer_json(rutas.TERRENO, {}) or {}
        if not (terreno.get("profundidad_fundacion_m") and terreno.get("q_adm_kPa")):
            try:
                terreno = pedir_datos_bases(self)
            except (OSError, ValueError) as exc:
                QMessageBox.critical(self, "No se pudieron guardar los datos del terreno", str(exc))
                return
            if terreno is None:
                return
        respuestas = f"{terreno['profundidad_fundacion_m']}\n"
        resultado = self._ejecutar_dimensionado("bases", portico, respuestas)
        if resultado is not None:
            self._mostrar_resultado("Bases", resultado)
            self.refrescar()

    def _ejecutar_dimensionado(self, etapa: str, portico: str, respuestas: str) -> dict | None:
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        self.statusBar().showMessage(f"Dimensionando {etapa} de {portico}…")
        try:
            return pipeline.ejecutar(etapa, portico, respuestas=respuestas)
        finally:
            QApplication.restoreOverrideCursor()

    def _cargar_columnas(self, portico: str) -> None:
        archivo, filas = datos_columnas(portico)
        self._llenar(self.tabla_columnas, filas)
        self.boton_dimensionar_columnas.setEnabled(bool(portico))
        planilla = buscar_archivo(
            rutas.SAL_COLUMNAS, "memoria_{portico}.txt", portico
        ) if portico else None
        self._cargar_planillas_texto(
            self.lista_planillas_columnas,
            self.etiqueta_planillas_columnas,
            self.texto_planilla_columnas,
            [planilla] if planilla else [],
            "memoria de columnas",
        )
        if not filas:
            self.etiqueta_columnas.setText(
                f"El pórtico «{portico}» todavía no tiene columnas dimensionadas. "
                "Si sus solicitaciones están al día, usá «Dimensionar columnas» para generar los resultados."
            )
            return
        en_rojo = sum(1 for _, malo in filas if malo)
        aviso = (
            f"  ·  {en_rojo} con estribo FUERA DE NORMA (en rojo)"
            if en_rojo
            else "  ·  estribos conformes"
        )
        self.etiqueta_columnas.setText(f"Fuente: {self._rel(archivo)}  ·  {len(filas)} columnas{aviso}")

    # ------------------------------------------------------------------
    # Pestaña 4: bases
    # ------------------------------------------------------------------
    def _cargar_bases(self, portico: str) -> None:
        archivo, filas, nota = datos_bases(portico)
        self._llenar(self.tabla_bases, filas)
        self.boton_dimensionar_bases.setEnabled(bool(portico))
        planilla = buscar_archivo(
            rutas.SAL_BASES,
            "bases_completas_{portico}.txt",
            portico,
        ) if portico else None
        self._cargar_planillas_texto(
            self.lista_planillas_bases,
            self.etiqueta_planillas_bases,
            self.texto_planilla_bases,
            [planilla] if planilla else [],
            "planilla de bases",
        )
        if archivo is None:
            estado_columnas = pipeline.estado_etapa("columnas", portico) if portico else None
            if estado_columnas and estado_columnas["estado"] == "ok":
                proximo_paso = "Usá «Dimensionar bases» para generar los resultados."
            else:
                proximo_paso = "Primero dimensioná las columnas del pórtico."
            self.etiqueta_bases.setText(
                f"El pórtico «{portico}» todavía no tiene bases calculadas. {proximo_paso}"
            )
            return
        en_rojo = sum(1 for _, malo in filas if malo)
        aviso = (
            f"  ·  {en_rojo} pide(n) viga de fundación o tienen q mín negativa (en rojo)"
            if en_rojo
            else "  ·  todas las bases cumplen"
        )
        self.etiqueta_bases.setText(
            f"Fuente: {self._rel(archivo)}  ·  {len(filas)} base(s){aviso}"
            + (f"\n{nota}" if nota else "")
        )

    # ------------------------------------------------------------------
    # Pestaña 5: losas
    # ------------------------------------------------------------------
    def _cargar_losas(self) -> None:
        self.lista_losas.clear()
        memorias = memorias_losas()
        for memoria in memorias:
            item = QListWidgetItem(memoria.name)
            item.setData(Qt.ItemDataRole.UserRole, str(memoria))
            item.setToolTip(self._rel(memoria))
            self.lista_losas.addItem(item)
        if memorias:
            self.etiqueta_losas.setText(
                f"{len(memorias)} memoria(s) de losa en salidas/losas — elegí una para leerla acá al lado."
            )
            self.lista_losas.setCurrentRow(0)
        else:
            self.etiqueta_losas.setText("Todavía no hay memorias de losa (etapa 10).")
            self.texto_losa.setPlainText("")

    def _al_elegir_losa(self, actual, _anterior=None) -> None:
        if actual is None:
            return
        ruta = Path(actual.data(Qt.ItemDataRole.UserRole))
        self.texto_losa.setPlainText(rutas.leer_texto(ruta, "(no se pudo leer el archivo)"))

    # ------------------------------------------------------------------
    # Pestaña 6: archivos
    # ------------------------------------------------------------------
    def _cargar_arbol(self) -> None:
        self.arbol.clear()
        for carpeta in (rutas.DATOS_GLOBAL, rutas.DATOS, rutas.SALIDAS):
            raiz = QTreeWidgetItem([carpeta.name + "\\", ""])
            raiz.setData(0, Qt.ItemDataRole.UserRole, str(carpeta))
            self.arbol.addTopLevelItem(raiz)
            self._agregar_hijos(raiz, carpeta)
            raiz.setExpanded(True)
        self.arbol.resizeColumnToContents(0)

    def _agregar_hijos(self, padre: QTreeWidgetItem, carpeta: Path) -> None:
        try:
            contenido = sorted(carpeta.iterdir(), key=lambda p: (p.is_file(), p.name.lower()))
        except OSError:
            return
        for camino in contenido:
            if camino.name.startswith("."):
                continue
            if camino.is_dir():
                nodo = QTreeWidgetItem([camino.name + "\\", ""])
                nodo.setData(0, Qt.ItemDataRole.UserRole, str(camino))
                padre.addChild(nodo)
                self._agregar_hijos(nodo, camino)
            else:
                try:
                    kb = max(1, camino.stat().st_size // 1024)
                except OSError:
                    kb = 0
                nodo = QTreeWidgetItem([camino.name, f"{kb} kB"])
                nodo.setData(0, Qt.ItemDataRole.UserRole, str(camino))
                padre.addChild(nodo)

    # ------------------------------------------------------------------
    # Acciones
    # ------------------------------------------------------------------
    def _aviso(self, titulo: str, texto: str) -> None:
        QMessageBox.warning(self, titulo, texto)

    def _ejecutar_etapa(self) -> None:
        dato = self._dato_etapa()
        if not dato or not dato.get("script"):
            return
        if dato["clave"] == "columnas":
            self._dimensionar_columnas()
            return
        if dato["clave"] == "bases":
            self._dimensionar_bases()
            return
        texto = (
            f"Se va a ejecutar:\n\n    py {dato['script']}\n\n"
            f"Directorio de trabajo: {rutas.RAIZ}"
        )
        advertencia = ADVERTENCIAS.get(dato["clave"])
        if advertencia:
            texto += "\n\nATENCIÓN\n" + advertencia

        respuesta = QMessageBox.question(
            self,
            f"Ejecutar {dato['nombre']}",
            texto,
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        if respuesta != QMessageBox.StandardButton.Yes:
            return

        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        self.statusBar().showMessage(f"Ejecutando {dato['script']}…")
        try:
            resultado = pipeline.ejecutar(dato["clave"], self._portico())
        finally:
            QApplication.restoreOverrideCursor()

        self._mostrar_resultado(dato["nombre"], resultado)
        self.refrescar()

    def _mostrar_resultado(self, titulo: str, resultado: dict) -> None:
        dialogo = QDialog(self)
        dialogo.setWindowTitle(f"Resultado: {titulo}")
        dialogo.resize(900, 540)
        caja = QVBoxLayout(dialogo)
        encabezado = (
            "Terminó bien."
            if resultado.get("ok")
            else f"No terminó bien (código {resultado.get('codigo')})."
        )
        caja.addWidget(QLabel(encabezado))

        salida = QPlainTextEdit()
        salida.setReadOnly(True)
        salida.setPlainText(resultado.get("salida") or "(sin salida en pantalla)")
        if resultado.get("error"):
            salida.appendPlainText("\n--- errores ---\n" + resultado["error"])
        caja.addWidget(salida, 1)

        botones = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        botones.accepted.connect(dialogo.accept)
        botones.rejected.connect(dialogo.reject)
        caja.addWidget(botones)
        dialogo.exec()

    def _abrir_salida_etapa(self) -> None:
        dato = self._dato_etapa()
        if not dato:
            return
        salida = self._primera_salida(dato)
        if salida is None:
            self._aviso("Sin salidas", "Esta etapa todavía no generó ningún archivo.")
            return
        error = abrir_con_windows(salida)
        if error:
            self._aviso("No se pudo abrir", error)

    def _abrir_del_arbol(self, item, _columna: int = 0) -> None:
        if item is None:
            return
        ruta = item.data(0, Qt.ItemDataRole.UserRole)
        if not ruta:
            return
        error = abrir_con_windows(Path(ruta))
        if error:
            self._aviso("No se pudo abrir", error)


# ---------------------------------------------------------------------------
# Arranque
# ---------------------------------------------------------------------------
def main() -> int:
    import traceback
    rutas.inicializar_obras()
    importlib.reload(pipeline)
    rutas.asegurar_directorios()
    aplicacion = QApplication.instance() or QApplication(sys.argv)
    aplicacion.setApplicationName("Calculador")
    try:
        ventana = VentanaPrincipal()
        with open("arranco.txt", "a", encoding="utf-8") as _f:
            _f.write("VentanaPrincipal ok\\n")
    except Exception:
        with open("arranco.txt", "a", encoding="utf-8") as _f:
            _f.write("VentanaPrincipal ERROR:\\n")
            traceback.print_exc(file=_f)
        raise
    try:
        ventana.show()
        with open("arranco.txt", "a", encoding="utf-8") as _f:
            _f.write("show ok\\n")
    except Exception:
        with open("arranco.txt", "a", encoding="utf-8") as _f:
            _f.write("show ERROR:\\n")
            traceback.print_exc(file=_f)
        raise
    return aplicacion.exec()


if __name__ == "__main__":
    raise SystemExit(main())

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
    QDialog,
    QFormLayout,
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

from calc import cargas, pipeline, rutas  # noqa: E402
from app.cargas_proyecto import PaginaCargas  # noqa: E402
from app.diseno_vigas import dimensionar_desde_app  # noqa: E402
from app.visualizacion import VistaPortico  # noqa: E402

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
    "Base", "Lado (m)", "Área (m²)", "Excentricidad (m)", "Espesor (m)",
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
        self.resize(1220, 750)

        # --- barra de arriba ---
        self.combo_portico = QComboBox()
        self.combo_portico.setMinimumWidth(190)
        self.combo_obras = QComboBox()
        self.combo_obras.setMinimumWidth(190)
        self.boton_nueva_obra = QPushButton("Nueva obra")
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
        self.tabla_bases = self._tabla(ENC_BASES)
        self.lista_losas = QListWidget()
        self.texto_losa = QPlainTextEdit()
        self.texto_losa.setReadOnly(True)
        self.arbol = QTreeWidget()
        self.arbol.setHeaderLabels(("Archivo", "Tamaño"))

        # --- pestaña Inicio ---
        self.etiqueta_inicio = QLabel()
        self.etiqueta_inicio.setWordWrap(True)
        self.nombre_proyecto = QLineEdit()
        self.id_proyecto_inicio = QLabel()
        self.id_proyecto_inicio.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self.ubicacion_proyecto = QLabel()
        self.ubicacion_proyecto.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self.notas_proyecto = QPlainTextEdit()
        self.notas_proyecto.setPlaceholderText("Criterios, contexto y notas de la obra")
        self.notas_proyecto.setMaximumHeight(58)
        self.boton_guardar_proyecto = QPushButton("Guardar carátula")
        self.boton_crear_portico = QPushButton("Crear pórtico / estructura")
        self.boton_inicio_cargas = QPushButton("Cargas de la obra")
        self.boton_inicio_losas = QPushButton("Losas alivianadas")
        self.vista_portico = VistaPortico()
        self.leyenda_cargas_visual = QLabel()
        self.leyenda_cargas_visual.setWordWrap(True)
        self.leyenda_cargas_visual.setTextFormat(Qt.TextFormat.RichText)
        self.leyenda_cargas_visual.setStyleSheet(
            "QLabel { color: #0f172a; background: #f8fafc; padding: 6px; }"
        )
        self.scroll_vista_portico = QScrollArea()
        self.scroll_vista_portico.setWidgetResizable(True)
        self.scroll_vista_portico.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.scroll_vista_portico.setMaximumHeight(370)
        self.scroll_vista_portico.setWidget(self.vista_portico)
        self.scroll_leyenda_cargas = QScrollArea()
        self.scroll_leyenda_cargas.setWidgetResizable(True)
        self.scroll_leyenda_cargas.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.scroll_leyenda_cargas.setMaximumHeight(105)
        self.scroll_leyenda_cargas.setWidget(self.leyenda_cargas_visual)
        self.lbl_cargas = QLabel()
        self.lbl_cargas.setWordWrap(True)
        self.lbl_motor = QLabel()
        self.lbl_motor.setWordWrap(True)
        self.boton_ver_cargas = QPushButton("Abrir el último análisis")
        self.boton_ir_cargas = QPushButton("Ver cargas del proyecto")
        self.boton_resolver_motor = QPushButton("Resolver solicitaciones del pórtico")
        self.boton_ver_motor = QPushButton("Abrir el JSON del motor")
        self.tabla_inicio = self._tabla(ENC_ENVOLVENTE)
        self.ultimo_analisis: Path | None = None

        self.pagina_cargas = PaginaCargas(self.refrescar, self._portico)
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

    def _pagina(self, texto_etiqueta: str, tabla: QTableWidget):
        pagina = QWidget()
        caja = QVBoxLayout(pagina)
        etiqueta = QLabel(texto_etiqueta)
        etiqueta.setWordWrap(True)
        caja.addWidget(etiqueta)
        caja.addWidget(tabla)
        return pagina, etiqueta

    def _armar_interfaz(self) -> None:
        contenedor = QWidget()
        principal = QVBoxLayout(contenedor)

        barra = QHBoxLayout()
        barra.addWidget(QLabel("Calculador de estructuras"))
        barra.addWidget(QLabel("Obra:"))
        barra.addWidget(self.combo_obras)
        barra.addWidget(self.boton_nueva_obra)
        barra.addWidget(self.boton_actualizar)
        barra.addWidget(self.etiqueta_resumen, 1)
        principal.addLayout(barra)

        pestanias = QTabWidget()
        pestanias.addTab(self._pagina_inicio(), "Inicio")
        pestanias.addTab(self.pagina_cargas, "Cargas")
        pestanias.addTab(self._pagina_etapas(), "1 · Estado y etapas")

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
        pestanias.addTab(pagina, "2 · Vigas")

        pagina, self.etiqueta_columnas = self._pagina("", self.tabla_columnas)
        pestanias.addTab(pagina, "3 · Columnas")

        pagina, self.etiqueta_bases = self._pagina("", self.tabla_bases)
        pestanias.addTab(pagina, "4 · Bases")

        pestanias.addTab(self._pagina_losas(), "5 · Losas")
        pestanias.addTab(self._pagina_archivos(), "6 · Archivos")

        principal.addWidget(pestanias, 1)
        self.setCentralWidget(contenedor)
        self.statusBar().showMessage(f"Carpeta del proyecto: {rutas.RAIZ}")

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
        self.boton_guardar_proyecto.clicked.connect(self._guardar_ficha_proyecto)
        self.boton_crear_portico.clicked.connect(self._crear_portico)
        self.boton_inicio_cargas.clicked.connect(lambda: self._ir_a_pestania("Cargas"))
        self.boton_inicio_losas.clicked.connect(lambda: self._ir_a_pestania("5 · Losas"))
        self.boton_ver_cargas.clicked.connect(self._abrir_ultimo_analisis)
        self.boton_ir_cargas.clicked.connect(self._ir_a_cargas)
        self.boton_resolver_motor.clicked.connect(self._resolver_motor)
        self.boton_dimensionar_vigas.clicked.connect(self._dimensionar_vigas)
        self.boton_ver_motor.clicked.connect(self._abrir_json_motor)
        self.combo_portico.currentTextChanged.connect(lambda _: self.refrescar())
        self.combo_portico.currentTextChanged.connect(self.pagina_cargas.actualizar_portico)
        self.tabla_etapas.itemSelectionChanged.connect(self._al_elegir_etapa)
        self.boton_ejecutar.clicked.connect(self._ejecutar_etapa)
        self.boton_abrir_salida.clicked.connect(self._abrir_salida_etapa)
        self.lista_losas.currentItemChanged.connect(self._al_elegir_losa)
        self.lista_planillas_vigas.currentItemChanged.connect(self._al_elegir_planilla_viga)
        self.arbol.itemDoubleClicked.connect(self._abrir_del_arbol)

    # ------------------------------------------------------------------
    # Refresco de la información
    # ------------------------------------------------------------------
    def _portico(self) -> str:
        return self.combo_portico.currentText().strip()

    def _rel(self, ruta) -> str:
        """Muestra la ruta recortada respecto de la carpeta del proyecto."""
        try:
            return str(Path(ruta).relative_to(rutas.RAIZ))
        except ValueError:
            return str(ruta)

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

    def refrescar(self) -> None:
        obras = rutas.listar_obras()
        activa = rutas.CARPETA_OBRA
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

        portico = self._portico()
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
        self._cargar_arbol()
        self.statusBar().showMessage(
            f"{len(porticos)} pórtico(s) en {rutas.CARPETA_OBRA or rutas.RAIZ}"
        )

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

    # ------------------------------------------------------------------
    # Pestaña Inicio
    # ------------------------------------------------------------------
    def _pagina_inicio(self) -> QWidget:
        pagina = QWidget()
        caja = QVBoxLayout(pagina)

        grupo_proyecto = QGroupBox("Carátula de la obra")
        form_proyecto = QFormLayout(grupo_proyecto)
        fila_identidad = QHBoxLayout()
        fila_identidad.addWidget(self.nombre_proyecto, 2)
        fila_identidad.addWidget(QLabel("ID:"))
        fila_identidad.addWidget(self.id_proyecto_inicio, 1)
        fila_identidad.addWidget(self.boton_guardar_proyecto)
        form_proyecto.addRow("Nombre:", fila_identidad)
        form_proyecto.addRow("Ubicación:", self.ubicacion_proyecto)
        form_proyecto.addRow("Notas:", self.notas_proyecto)
        caja.addWidget(grupo_proyecto)

        self.etiqueta_inicio.setStyleSheet("font-size: 11pt;")
        caja.addWidget(self.etiqueta_inicio)

        grupo_porticos = QGroupBox("Pórticos de esta obra")
        caja_porticos = QHBoxLayout(grupo_porticos)
        caja_porticos.addWidget(QLabel("Pórtico activo:"))
        caja_porticos.addWidget(self.combo_portico, 1)
        caja_porticos.addWidget(self.boton_crear_portico)
        caja.addWidget(grupo_porticos)

        accesos = QHBoxLayout()
        accesos.addWidget(self.boton_inicio_cargas)
        accesos.addWidget(self.boton_inicio_losas)
        accesos.addStretch(1)
        caja.addLayout(accesos)

        grupo_cargas = QGroupBox("Cargas del proyecto")
        caja_cargas = QVBoxLayout(grupo_cargas)
        caja_cargas.addWidget(self.lbl_cargas)
        fila_cargas = QHBoxLayout()
        fila_cargas.addWidget(self.boton_ir_cargas)
        fila_cargas.addWidget(self.boton_ver_cargas)
        fila_cargas.addStretch(1)
        caja_cargas.addLayout(fila_cargas)
        grupo_motor = QGroupBox("Solicitaciones del pórtico")
        caja_motor = QVBoxLayout(grupo_motor)
        caja_motor.addWidget(self.lbl_motor)
        fila_motor = QHBoxLayout()
        fila_motor.addWidget(self.boton_resolver_motor)
        fila_motor.addWidget(self.boton_ver_motor)
        fila_motor.addStretch(1)
        caja_motor.addLayout(fila_motor)
        estados = QVBoxLayout()
        estados.addWidget(grupo_cargas)
        estados.addWidget(grupo_motor)
        estados.addStretch(1)

        grupo_vista = QGroupBox("Esquema del pórtico y cargas asignadas")
        caja_vista = QVBoxLayout(grupo_vista)
        caja_vista.addWidget(self.scroll_vista_portico)
        caja_vista.addWidget(self.scroll_leyenda_cargas)
        zona = QHBoxLayout()
        zona.addWidget(grupo_vista, 2)
        zona.addLayout(estados, 1)
        caja.addLayout(zona, 2)

        caja.addWidget(QLabel("Envolvente del pórtico elegido (máximos en módulo, según el motor):"))
        caja.addWidget(self.tabla_inicio, 1)
        return pagina

    def _cargar_inicio(self, portico: str) -> None:
        import html

        proyecto = rutas.leer_json(rutas.CARGAS, {}) or {}
        self.nombre_proyecto.setText(str(proyecto.get("obra", "Obra")))
        self.id_proyecto_inicio.setText(PaginaCargas._id_proyecto(proyecto))
        ubicacion = proyecto.get("viento", {}).get("referencia_cirsoc_102_25", {}).get("ubicacion", {})
        self.ubicacion_proyecto.setText(", ".join(
            str(ubicacion.get(k, "")).strip() for k in ("ciudad", "provincia", "pais")
            if str(ubicacion.get(k, "")).strip()
        ) or "Sin ubicación definida")
        self.notas_proyecto.setPlainText(str(proyecto.get("notas", "")))
        self.vista_portico.actualizar(portico)
        self.leyenda_cargas_visual.setText(self.vista_portico.leyenda)
        estructura = rutas.cargar_estructura().get(portico, {}) if portico else {}
        elementos = proyecto.get("elementos", {})
        cargas_activas = sum(bool(e.get("activo", True)) for e in elementos.values())
        aplicaciones_activas = sum(
            bool(a.get("activa", True)) for a in proyecto.get("aplicaciones", [])
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
            f"{aplicaciones_activas} aplicación(es) a tramos."
        )
        self.boton_ver_cargas.setEnabled(self.ultimo_analisis is not None)

        archivo, datos, filas = datos_solicitaciones(portico)
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
        self.boton_resolver_motor.setEnabled(bool(portico) and bool(estructura) and requisitos_motor and not motor_al_dia)
        if motor_al_dia:
            self.boton_resolver_motor.setText("Solicitaciones al día")
            self.boton_resolver_motor.setToolTip("El resultado vigente está guardado; podés abrirlo con el botón de al lado.")
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
        self._llenar(self.tabla_inicio, filas if mostrar_envolvente else [])

    def _guardar_ficha_proyecto(self) -> None:
        datos = cargas.datos_cargas()
        datos["obra"] = self.nombre_proyecto.text().strip() or "Obra"
        datos["notas"] = self.notas_proyecto.toPlainText().strip()
        datos["id_proyecto"] = PaginaCargas._id_proyecto(datos)
        rutas.guardar_json(rutas.CARGAS, datos)
        ficha = rutas.leer_json(rutas.CARPETA_OBRA / "obra.json", {}) or {}
        ficha.update({
            "nombre": datos["obra"], "id": datos["id_proyecto"],
            "ubicacion": datos.get("viento", {}).get("referencia_cirsoc_102_25", {}).get("ubicacion", {}),
        })
        rutas.guardar_json(rutas.CARPETA_OBRA / "obra.json", ficha)
        self.pagina_cargas.recargar()
        self.statusBar().showMessage("Ficha del proyecto guardada en datos/cargas.json.", 5000)
        self.refrescar()

    def _ir_a_pestania(self, nombre: str) -> None:
        pestanias = self.findChild(QTabWidget)
        if pestanias is None:
            return
        for indice in range(pestanias.count()):
            if pestanias.tabText(indice) == nombre:
                pestanias.setCurrentIndex(indice)
                return

    def _crear_portico(self) -> None:
        """Crea geometría básica desde Inicio, dentro de la obra actual."""
        nombre, ok = QInputDialog.getText(self, "Nuevo pórtico", "Nombre del pórtico:")
        if not ok:
            return
        nombre = nombre.strip()
        if not nombre:
            QMessageBox.warning(self, "Nombre requerido", "Ingresá un nombre para el pórtico.")
            return
        estructura = rutas.cargar_estructura()
        if nombre in estructura:
            QMessageBox.warning(self, "Nombre existente", f"Ya existe el pórtico {nombre}.")
            return
        try:
            numero_portico = rutas.numero_portico(nombre, estructura)
        except ValueError as exc:
            QMessageBox.warning(self, "Número de pórtico repetido", str(exc))
            return

        pisos, ok = QInputDialog.getInt(
            self, "Niveles", "Cantidad de niveles sobre planta baja:", 0, 0, 20
        )
        if not ok:
            return
        tramos_texto, ok = QInputDialog.getText(
            self, "Tramos entre columnas",
            "Longitudes de los tramos entre columnas [m] (separadas por punto y coma):",
            text="4; 4",
        )
        if not ok:
            return
        try:
            luces = [float(x.strip().replace(",", ".")) for x in tramos_texto.split(";") if x.strip()]
            if not luces or any(x <= 0 for x in luces):
                raise ValueError
        except ValueError:
            QMessageBox.warning(
                self, "Luces inválidas", "Ingresá longitudes positivas, por ejemplo 4; 5; 3.5."
            )
            return
        voladizos = {}
        for lado in ("izquierda", "derecha"):
            respuesta = QMessageBox.question(
                self, "Voladizo", f"¿El pórtico tiene voladizo a la {lado}?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No,
            )
            if respuesta == QMessageBox.StandardButton.Yes:
                longitud, ok = QInputDialog.getDouble(
                    self, f"Voladizo {lado}", f"Longitud del voladizo a la {lado} [m]:",
                    1.0, 0.1, 30.0, 2,
                )
                if not ok:
                    return
                voladizos[lado] = longitud
        altura, ok = QInputDialog.getDouble(
            self, "Altura", "Altura uniforme entre niveles [m]:", 3.0, 0.1, 30.0, 2
        )
        if not ok:
            return

        datos = {"vigas": {}, "columnas": {}, "bases": {}, "cargas_puntuales": []}
        coordenadas = [0.0]
        for luz in luces:
            coordenadas.append(coordenadas[-1] + luz)
        for nivel in range(pisos + 1):
            y = nivel * altura
            for indice, x in enumerate(coordenadas):
                letra = chr(97 + indice) if indice < 26 else str(indice + 1)
                datos["columnas"][f"C{nivel}-{letra}"] = {
                    "x": x, "altura_m": altura, "nivel": y
                }
                if nivel == 0:
                    datos["bases"][f"B0-{letra}"] = {"x": x, "tipo": "empotramiento"}
            viga_id = f"V{nivel}-{numero_portico}"
            tramos = []
            x = 0.0
            for indice, luz in enumerate(luces, 1):
                tramos.append({
                    "id": f"{viga_id} T{indice}", "longitud_m": luz,
                    "es_voladizo": False, "x_inicio": x, "x_fin": x + luz,
                    "cargas_puntuales": [],
                })
                x += luz
            if "izquierda" in voladizos:
                largo = voladizos["izquierda"]
                tramos.append({
                    "id": f"{viga_id} tv_izq", "longitud_m": largo,
                    "es_voladizo": True, "x_inicio": -largo, "x_fin": 0.0,
                    "cargas_puntuales": [],
                })
            if "derecha" in voladizos:
                largo = voladizos["derecha"]
                tramos.append({
                    "id": f"{viga_id} tv_der", "longitud_m": largo,
                    "es_voladizo": True, "x_inicio": x, "x_fin": x + largo,
                    "cargas_puntuales": [],
                })
            datos["vigas"][viga_id] = {"tramos": tramos}
        estructura[nombre] = datos
        rutas.guardar_estructura(estructura)
        self.refrescar()
        self.combo_portico.setCurrentText(nombre)
        self.statusBar().showMessage(f"{nombre} creado en la obra actual.", 5000)

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

    def _ir_a_cargas(self) -> None:
        pestanias = self.centralWidget().findChild(QTabWidget)
        if pestanias is not None:
            pestanias.setCurrentWidget(self.pagina_cargas)

    def _abrir_json_motor(self) -> None:
        archivo = buscar_archivo(rutas.SAL_SOLICITACIONES, "{portico}.json", self._portico())
        if archivo is None:
            self._aviso("Sin solicitaciones", "Todavía no hay JSON del motor para este pórtico.")
            return
        error = abrir_con_windows(archivo)
        if error:
            self._aviso("No se pudo abrir", error)

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
        if estado_motor["estado"] == "ok":
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

    def _dimensionar_vigas(self) -> None:
        portico = self._portico()
        if not portico:
            QMessageBox.information(self, "Vigas", "Elegí un pórtico primero.")
            return
        resultado = dimensionar_desde_app(self, portico)
        if resultado:
            _salida, self.advertencias_vigas = resultado
            self.refrescar()

    def _cargar_columnas(self, portico: str) -> None:
        archivo, filas = datos_columnas(portico)
        self._llenar(self.tabla_columnas, filas)
        if not filas:
            self.etiqueta_columnas.setText(
                f"Sin columnas de {portico} en salidas/columnas/planilla_columnas.csv: falta calcular la etapa 7."
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
        if archivo is None:
            self.etiqueta_bases.setText(f"Sin bases calculadas para {portico}: falta calcular la etapa 8.")
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
    rutas.inicializar_obras()
    importlib.reload(pipeline)
    rutas.asegurar_directorios()
    aplicacion = QApplication.instance() or QApplication(sys.argv)
    aplicacion.setApplicationName("Calculador")
    ventana = VentanaPrincipal()
    ventana.show()
    return aplicacion.exec()


if __name__ == "__main__":
    raise SystemExit(main())

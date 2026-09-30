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
import os
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
    QDialogButtonBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QMainWindow,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QSplitter,
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from calc import pipeline, rutas  # noqa: E402

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
        cumple_flecha = fle.get("cumple_L360", fle.get("cumple_L180"))
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
                numero(fle.get("lim_L360_mm", fle.get("lim_L180_mm")), 1),
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
# Losas y resumen general
# ---------------------------------------------------------------------------
def memorias_losas() -> list[Path]:
    """Memorias de losa existentes, de la más nueva a la más vieja."""
    memorias = rutas.listar(rutas.SAL_LOSAS, "memoria_losa_*.txt")
    return sorted(memorias, key=lambda p: p.stat().st_mtime, reverse=True)


def resumen_etapas(portico: str) -> str:
    """Una línea con el conteo de etapas por estado."""
    conteo: dict[str, int] = {}
    for dato in pipeline.semaforo(portico):
        conteo[dato["estado"]] = conteo.get(dato["estado"], 0) + 1
    nombres = {
        "ok": "al día",
        "desactualizada": "para recalcular",
        "pendiente": "sin calcular",
        "sin_datos": "sin datos",
        "a_desarrollar": "a desarrollar",
    }
    partes = [f"{cantidad} {nombres.get(estado, estado)}" for estado, cantidad in conteo.items()]
    return "Etapas: " + ", ".join(partes)


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
        self.tabla_columnas = self._tabla(ENC_COLUMNAS)
        self.tabla_bases = self._tabla(ENC_BASES)
        self.lista_losas = QListWidget()
        self.texto_losa = QPlainTextEdit()
        self.texto_losa.setReadOnly(True)
        self.arbol = QTreeWidget()
        self.arbol.setHeaderLabels(("Archivo", "Tamaño"))

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
        barra.addWidget(QLabel("Pórtico:"))
        barra.addWidget(self.combo_portico)
        barra.addWidget(self.boton_actualizar)
        barra.addWidget(self.etiqueta_resumen, 1)
        principal.addLayout(barra)

        pestanias = QTabWidget()
        pestanias.addTab(self._pagina_etapas(), "1 · Estado y etapas")

        pagina, self.etiqueta_vigas = self._pagina("", self.tabla_vigas)
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
        self.combo_portico.currentTextChanged.connect(lambda _: self.refrescar())
        self.tabla_etapas.itemSelectionChanged.connect(self._al_elegir_etapa)
        self.boton_ejecutar.clicked.connect(self._ejecutar_etapa)
        self.boton_abrir_salida.clicked.connect(self._abrir_salida_etapa)
        self.lista_losas.currentItemChanged.connect(self._al_elegir_losa)
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
        porticos = rutas.listar_porticos()
        elegido = self._portico()
        self.combo_portico.blockSignals(True)
        self.combo_portico.clear()
        self.combo_portico.addItems(porticos)
        if elegido in porticos:
            self.combo_portico.setCurrentText(elegido)
        self.combo_portico.blockSignals(False)

        portico = self._portico()
        if porticos:
            self.etiqueta_resumen.setText(resumen_etapas(portico))
        else:
            self.etiqueta_resumen.setText(
                "No hay pórticos cargados en datos/estructura.json (empezá por la etapa 3, geometría)."
            )

        self._cargar_etapas(portico)
        self._cargar_vigas(portico)
        self._cargar_columnas(portico)
        self._cargar_bases(portico)
        self._cargar_losas()
        self._cargar_arbol()
        self.statusBar().showMessage(f"{len(porticos)} pórtico(s) en el proyecto  ·  {rutas.RAIZ}")

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
        puede = bool(script) and not dato["interactiva"] and (rutas.RAIZ / script).exists()
        self.boton_ejecutar.setEnabled(puede)
        self.boton_ejecutar.setToolTip(
            "" if puede else "Esta etapa todavía pide datos por teclado: se corre desde la consola."
        )
        self.boton_abrir_salida.setEnabled(self._primera_salida(dato) is not None)

    # ------------------------------------------------------------------
    # Pestañas 2 y 3: vigas y columnas
    # ------------------------------------------------------------------
    def _cargar_vigas(self, portico: str) -> None:
        archivo, filas = datos_vigas(portico)
        self._llenar(self.tabla_vigas, filas)
        if archivo is None:
            self.etiqueta_vigas.setText(f"Sin resultados de vigas para {portico}: falta calcular la etapa 5.")
            return
        en_rojo = sum(1 for _, malo in filas if malo)
        aviso = (
            f"  ·  {en_rojo} con alguna verificación que NO cumple (en rojo)"
            if en_rojo
            else "  ·  todas las verificaciones cumplen"
        )
        self.etiqueta_vigas.setText(f"Fuente: {self._rel(archivo)}  ·  {len(filas)} tramo(s){aviso}")

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
        for carpeta in (rutas.DATOS, rutas.SALIDAS):
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
    rutas.asegurar_directorios()
    aplicacion = QApplication.instance() or QApplication(sys.argv)
    aplicacion.setApplicationName("Calculador")
    ventana = VentanaPrincipal()
    ventana.show()
    return aplicacion.exec()


if __name__ == "__main__":
    raise SystemExit(main())

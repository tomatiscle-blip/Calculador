"""Edición y visualización de ejes de referencia en planta."""

from __future__ import annotations

import math
import uuid

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QPainter, QPen
from PySide6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPushButton,
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from calc import rutas


def cargar_ejes() -> dict[str, list[dict]]:
    """Lee los ejes de la obra activa; las obras antiguas comienzan sin ejes."""
    datos = rutas.leer_json(rutas.DATOS / "ejes.json", {}) or {}
    resultado = {"x": [], "y": []}
    for direccion in resultado:
        for item in datos.get(direccion, []):
            if not isinstance(item, dict):
                continue
            eje = dict(item)
            eje.setdefault("id", f"{direccion}:{str(eje.get('nombre', '')).casefold()}")
            resultado[direccion].append(eje)
    return resultado


def cargar_niveles() -> list[dict]:
    """Lee niveles de referencia Z de la obra activa."""
    datos = rutas.leer_json(rutas.DATOS / "niveles.json", {}) or {}
    return [dict(item) for item in datos.get("niveles", []) if isinstance(item, dict)]


def resumen_ejes(ejes: dict[str, list[dict]]) -> str:
    """Devuelve un resumen legible de los ejes definidos."""
    partes = []
    for direccion in ("x", "y"):
        etiqueta = direccion.upper()
        items = ejes.get(direccion, [])
        texto = ", ".join(
            f"{item['nombre']} = {float(item['coordenada_m']):g} m"
            for item in items
        )
        partes.append(f"{etiqueta}: {texto or 'sin definir'}")
    return "Ejes de planta — " + " · ".join(partes)


def resumen_niveles(niveles: list[dict]) -> str:
    """Devuelve un resumen legible de las cotas de nivel."""
    items = ", ".join(
        f"{item['nombre']} = {float(item['cota_m']):g} m"
        for item in niveles
    )
    return f"Niveles Z — {items or 'sin definir'}"


class VistaEjes(QWidget):
    """Esquema en planta de las líneas de ejes X/Y definidas."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.ejes: dict[str, list[dict]] = {"x": [], "y": []}
        self.setMinimumSize(300, 220)

    def establecer_ejes(self, ejes: dict[str, list[dict]]) -> None:
        self.ejes = ejes
        self.update()

    def paintEvent(self, event) -> None:
        super().paintEvent(event)
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.fillRect(self.rect(), QColor("#ffffff"))
        painter.setFont(self.font())
        painter.setPen(QPen(QColor("#64748b"), 1))
        painter.drawText(12, 20, "Esquema de referencia (no representa elementos estructurales)")

        x_items = self.ejes.get("x", [])
        y_items = self.ejes.get("y", [])
        left, right = 54.0, max(55.0, self.width() - 36.0)
        top, bottom = 42.0, max(43.0, self.height() - 34.0)
        plot_width, plot_height = right - left, bottom - top
        painter.setPen(QPen(QColor("#cbd5e1"), 1))
        painter.drawRect(int(left), int(top), int(plot_width), int(plot_height))

        posiciones_x = [float(item["coordenada_m"]) for item in x_items]
        posiciones_y = [float(item["coordenada_m"]) for item in y_items]
        min_x, max_x = self._rango(posiciones_x)
        min_y, max_y = self._rango(posiciones_y)

        for item in x_items:
            valor = float(item["coordenada_m"])
            x = left + (valor - min_x) / (max_x - min_x) * plot_width
            painter.setPen(QPen(QColor("#2563eb"), 1, Qt.PenStyle.DashLine))
            painter.drawLine(int(x), int(top), int(x), int(bottom))
            painter.setPen(QColor("#1d4ed8"))
            painter.drawText(int(x - 15), int(bottom + 18), str(item["nombre"]))

        for item in y_items:
            valor = float(item["coordenada_m"])
            y = bottom - (valor - min_y) / (max_y - min_y) * plot_height
            painter.setPen(QPen(QColor("#15803d"), 1, Qt.PenStyle.DashLine))
            painter.drawLine(int(left), int(y), int(right), int(y))
            painter.setPen(QColor("#166534"))
            painter.drawText(8, int(y + 4), str(item["nombre"]))

        if not x_items and not y_items:
            painter.setPen(QColor("#64748b"))
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "Agregá ejes X o Y")

    @staticmethod
    def _rango(posiciones: list[float]) -> tuple[float, float]:
        if not posiciones:
            return 0.0, 1.0
        minimo, maximo = min(posiciones), max(posiciones)
        if minimo == maximo:
            return minimo - 0.5, maximo + 0.5
        return minimo, maximo


class VistaPlantaPorticos(QWidget):
    """Muestra la grilla de referencia y las columnas de cada pórtico."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.ejes: dict[str, list[dict]] = {"x": [], "y": []}
        self.porticos: dict[str, dict] = {}
        self.seleccionado: str | None = None
        self.setMinimumSize(520, 300)

    def establecer_datos(
        self,
        ejes: dict[str, list[dict]],
        estructura: dict,
        asociaciones: dict[str, dict],
        seleccionado: str | None,
    ) -> None:
        self.ejes = ejes
        self.seleccionado = seleccionado
        self.porticos = {}
        for nombre, datos in estructura.items():
            if not isinstance(datos, dict):
                continue
            referencia = asociaciones.get(nombre, {})
            direccion = str(referencia.get("direccion", "x")).lower()
            if direccion not in ("x", "y"):
                direccion = "x"
            eje_transversal = self._buscar_eje(
                ejes, "y" if direccion == "x" else "x",
                referencia.get("eje_id"),
            )
            eje_longitudinal = self._buscar_eje(
                ejes, direccion, referencia.get("eje_longitudinal_id"),
            )
            try:
                posicion = (
                    float(eje_transversal["coordenada_m"])
                    if eje_transversal is not None
                    else float(datos.get("posicion_planta_m", 0.0))
                )
                origen_longitudinal = (
                    float(eje_longitudinal["coordenada_m"])
                    if eje_longitudinal is not None else 0.0
                )
            except (TypeError, ValueError):
                posicion = 0.0
                origen_longitudinal = 0.0
            columnas = datos.get("columnas", {})
            coordenadas_locales = [
                round(float(columna["x"]), 6)
                for columna in columnas.values()
                if isinstance(columna, dict)
                and self._es_numero_finito(columna.get("x"))
            ]
            coordenadas = sorted({
                round(origen_longitudinal + x, 6) for x in coordenadas_locales
            })
            extremos_viga_locales = [
                float(extremo)
                for viga in datos.get("vigas", {}).values()
                if isinstance(viga, dict)
                for tramo in viga.get("tramos", [])
                if isinstance(tramo, dict)
                for extremo in (tramo.get("x_inicio"), tramo.get("x_fin"))
                if self._es_numero_finito(extremo)
            ]
            extremos_viga = (
                (
                    origen_longitudinal + min(extremos_viga_locales),
                    origen_longitudinal + max(extremos_viga_locales),
                )
                if extremos_viga_locales else None
            )
            self.porticos[str(nombre)] = {
                "direccion": direccion,
                "posicion": posicion,
                "columnas": coordenadas,
                "extremos_viga": extremos_viga,
            }
        self.update()

    @staticmethod
    def _buscar_eje(
        ejes: dict[str, list[dict]], familia: str, eje_id: object,
    ) -> dict | None:
        if not eje_id:
            return None
        return next(
            (
                item for item in ejes.get(familia, [])
                if str(item.get(
                    "id", f"{familia}:{str(item.get('nombre', '')).casefold()}"
                )) == str(eje_id)
            ),
            None,
        )

    @staticmethod
    def _es_numero_finito(valor: object) -> bool:
        try:
            return math.isfinite(float(valor))
        except (TypeError, ValueError):
            return False

    def paintEvent(self, event) -> None:
        super().paintEvent(event)
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.fillRect(self.rect(), QColor("#ffffff"))
        painter.setPen(QColor("#334155"))
        painter.drawText(12, 20, "Planta de referencia  ·  X horizontal  ·  Y vertical")

        left, right = 58.0, max(59.0, self.width() - 24.0)
        top, bottom = 36.0, max(37.0, self.height() - 38.0)
        ancho, alto = right - left, bottom - top
        coords_x = [float(item["coordenada_m"]) for item in self.ejes.get("x", [])]
        coords_y = [float(item["coordenada_m"]) for item in self.ejes.get("y", [])]
        for datos in self.porticos.values():
            if datos["direccion"] == "x":
                coords_y.append(datos["posicion"])
                coords_x.extend(datos["columnas"])
                if datos["extremos_viga"]:
                    coords_x.extend(datos["extremos_viga"])
            else:
                coords_x.append(datos["posicion"])
                coords_y.extend(datos["columnas"])
                if datos["extremos_viga"]:
                    coords_y.extend(datos["extremos_viga"])
        min_x, max_x = self._rango_con_margen(coords_x)
        min_y, max_y = self._rango_con_margen(coords_y)

        def punto(x: float, y: float) -> tuple[int, int]:
            px = left + (x - min_x) / (max_x - min_x) * ancho
            py = bottom - (y - min_y) / (max_y - min_y) * alto
            return int(px), int(py)

        painter.setPen(QPen(QColor("#cbd5e1"), 1))
        painter.drawRect(int(left), int(top), int(ancho), int(alto))
        for item in self.ejes.get("x", []):
            x = float(item["coordenada_m"])
            px, _ = punto(x, min_y)
            painter.setPen(QPen(QColor("#2563eb"), 1, Qt.PenStyle.DashLine))
            painter.drawLine(px, int(top), px, int(bottom))
            painter.setPen(QColor("#1d4ed8"))
            painter.drawText(px - 12, int(bottom + 16), str(item.get("nombre", "X")))
        for item in self.ejes.get("y", []):
            y = float(item["coordenada_m"])
            _, py = punto(min_x, y)
            painter.setPen(QPen(QColor("#15803d"), 1, Qt.PenStyle.DashLine))
            painter.drawLine(int(left), py, int(right), py)
            painter.setPen(QColor("#166534"))
            painter.drawText(8, py + 4, str(item.get("nombre", "Y")))

        for nombre, datos in self.porticos.items():
            columnas = datos["columnas"]
            if not columnas:
                continue
            activo = nombre == self.seleccionado
            color = QColor("#dc2626") if activo else QColor("#475569")
            puntos = [
                punto(x, datos["posicion"]) if datos["direccion"] == "x"
                else punto(datos["posicion"], x)
                for x in columnas
            ]
            painter.setPen(QPen(color, 3 if activo else 2))
            if len(puntos) > 1:
                painter.drawLine(*puntos[0], *puntos[-1])
            extremos = datos["extremos_viga"]
            if extremos:
                if datos["direccion"] == "x":
                    p_inicio = punto(extremos[0], datos["posicion"])
                    p_fin = punto(extremos[1], datos["posicion"])
                else:
                    p_inicio = punto(datos["posicion"], extremos[0])
                    p_fin = punto(datos["posicion"], extremos[1])
                if puntos and extremos[0] < columnas[0]:
                    painter.setPen(QPen(color, 2, Qt.PenStyle.DashLine))
                    painter.drawLine(*p_inicio, *puntos[0])
                if puntos and extremos[1] > columnas[-1]:
                    painter.setPen(QPen(color, 2, Qt.PenStyle.DashLine))
                    painter.drawLine(*puntos[-1], *p_fin)
            painter.setBrush(color)
            for px, py in puntos:
                painter.drawEllipse(px - 5, py - 5, 10, 10)
            painter.setPen(color)
            px, py = puntos[0]
            painter.drawText(px + 7, py - 7, nombre)

        painter.setPen(QColor("#64748b"))
        painter.drawText(
            12, self.height() - 10,
            "Ejes punteados  ·  puntos: columnas  ·  voladizos: extensión punteada  ·  "
            "rojo: pórtico seleccionado",
        )

    @staticmethod
    def _rango_con_margen(valores: list[float]) -> tuple[float, float]:
        if not valores:
            return -0.5, 0.5
        minimo, maximo = min(valores), max(valores)
        margen = max((maximo - minimo) * 0.08, 0.5)
        return minimo - margen, maximo + margen


class EditorUbicacionPortico(QDialog):
    """Elige la dirección, el plano y los cruces de ejes que llevan columnas."""

    def __init__(self, ejes: dict[str, list[dict]], parent=None):
        super().__init__(parent)
        self.setWindowTitle("Ubicación y columnas del pórtico")
        self.ejes = ejes
        self.direccion = QComboBox()
        self.direccion.addItem("Paralelo a X", "x")
        self.direccion.addItem("Paralelo a Y", "y")
        self.eje_transversal = QComboBox()

        self.tabla = QTableWidget(0, 3)
        self.tabla.setHorizontalHeaderLabels(
            ("Columna en cruce", "Eje longitudinal", "Coordenada acumulada (m)")
        )
        self.tabla.verticalHeader().setVisible(False)
        self.tabla.horizontalHeader().setStretchLastSection(True)
        self.vista = VistaPlantaPorticos()
        self._actualizar_eje_transversal()
        self._actualizar_tabla_ejes()

        instrucciones = QLabel(
            "Marcá los cruces donde habrá columnas. Las luces se calcularán entre las "
            "columnas seleccionadas; los cruces sin marcar quedan sin apoyo. Si seleccionás "
            "menos de dos ejes, podrás ingresar las luces manualmente. Los voladizos se "
            "definen aparte y pueden pasar más allá de otros ejes."
        )
        instrucciones.setWordWrap(True)
        formulario = QFormLayout()
        formulario.addRow("Dirección del pórtico:", self.direccion)
        formulario.addRow("Eje transversal:", self.eje_transversal)

        botones = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        botones.button(QDialogButtonBox.StandardButton.Ok).setText("Usar esta ubicación")
        botones.accepted.connect(self._aceptar)
        botones.rejected.connect(self.reject)
        self.direccion.currentIndexChanged.connect(self._cambio_direccion)
        self.eje_transversal.currentIndexChanged.connect(self._actualizar_vista)
        self.tabla.itemChanged.connect(self._actualizar_vista)

        caja = QVBoxLayout(self)
        caja.addWidget(instrucciones)
        caja.addLayout(formulario)
        caja.addWidget(self.tabla, 1)
        caja.addWidget(self.vista, 2)
        caja.addWidget(botones)
        self.resize(900, 760)
        self._actualizar_vista()

    def _cambio_direccion(self) -> None:
        self._actualizar_eje_transversal()
        self._actualizar_tabla_ejes()
        self._actualizar_vista()

    def _actualizar_eje_transversal(self) -> None:
        familia = "y" if self.direccion.currentData() == "x" else "x"
        anterior = self.eje_transversal.currentData()
        self.eje_transversal.blockSignals(True)
        self.eje_transversal.clear()
        for item in self.ejes.get(familia, []):
            self.eje_transversal.addItem(
                f"{item['nombre']} ({float(item['coordenada_m']):.3f} m)",
                str(item.get("id", f"{familia}:{str(item['nombre']).casefold()}")),
            )
        indice = self.eje_transversal.findData(anterior)
        if indice < 0:
            indice = 0
        self.eje_transversal.setCurrentIndex(indice)
        self.eje_transversal.blockSignals(False)

    def _actualizar_tabla_ejes(self) -> None:
        familia = str(self.direccion.currentData())
        ejes = sorted(
            self.ejes.get(familia, []),
            key=lambda item: float(item["coordenada_m"]),
        )
        self.tabla.blockSignals(True)
        self.tabla.setRowCount(len(ejes))
        for fila, eje in enumerate(ejes):
            columna = QTableWidgetItem("Columna")
            columna.setFlags(
                columna.flags() | Qt.ItemFlag.ItemIsUserCheckable
            )
            columna.setCheckState(Qt.CheckState.Checked)
            self.tabla.setItem(fila, 0, columna)
            item_eje = QTableWidgetItem(str(eje["nombre"]))
            item_eje.setData(
                Qt.ItemDataRole.UserRole,
                str(eje.get("id", f"{familia}:{str(eje['nombre']).casefold()}")),
            )
            item_eje.setFlags(item_eje.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.tabla.setItem(fila, 1, item_eje)
            item_coordenada = QTableWidgetItem(f"{float(eje['coordenada_m']):g}")
            item_coordenada.setFlags(
                item_coordenada.flags() & ~Qt.ItemFlag.ItemIsEditable
            )
            self.tabla.setItem(fila, 2, item_coordenada)
        self.tabla.blockSignals(False)

    def _ejes_con_columna(self) -> list[dict]:
        seleccionados = []
        for fila in range(self.tabla.rowCount()):
            item_columna = self.tabla.item(fila, 0)
            item_eje = self.tabla.item(fila, 1)
            item_coordenada = self.tabla.item(fila, 2)
            if (
                item_columna is None
                or item_columna.checkState() != Qt.CheckState.Checked
                or item_eje is None
                or item_coordenada is None
            ):
                continue
            seleccionados.append({
                "id": str(item_eje.data(Qt.ItemDataRole.UserRole)),
                "nombre": item_eje.text(),
                "coordenada_m": float(item_coordenada.text()),
            })
        return sorted(seleccionados, key=lambda item: item["coordenada_m"])

    def _actualizar_vista(self) -> None:
        direccion = str(self.direccion.currentData())
        ejes_columna = self._ejes_con_columna()
        eje_id = str(self.eje_transversal.currentData() or "")
        eje_transversal = VistaPlantaPorticos._buscar_eje(
            self.ejes, "y" if direccion == "x" else "x", eje_id
        )
        if eje_transversal is None:
            return
        posicion = float(eje_transversal["coordenada_m"])
        coordenada_inicial = (
            float(ejes_columna[0]["coordenada_m"]) if ejes_columna else 0.0
        )
        columnas = {
            f"C0-{chr(97 + indice)}": {
                "x": float(eje["coordenada_m"]) - coordenada_inicial
            }
            for indice, eje in enumerate(ejes_columna)
        }
        tramos = []
        for indice, (izquierda, derecha) in enumerate(
            zip(ejes_columna, ejes_columna[1:]), 1
        ):
            inicio = float(izquierda["coordenada_m"]) - coordenada_inicial
            fin = float(derecha["coordenada_m"]) - coordenada_inicial
            tramos.append({
                "id": f"vista-{indice}",
                "x_inicio": inicio,
                "x_fin": fin,
            })
        estructura = {
            "Pórtico nuevo": {
                "posicion_planta_m": posicion,
                "columnas": columnas,
                "vigas": {"V": {"tramos": tramos}},
            }
        }
        referencia = {
            "direccion": direccion,
            "eje_id": eje_id or None,
            "desfase_m": 0.0,
            "eje_longitudinal_id": ejes_columna[0]["id"] if ejes_columna else None,
            "desfase_longitudinal_m": 0.0,
        }
        self.vista.establecer_datos(
            self.ejes, estructura, {"Pórtico nuevo": referencia}, "Pórtico nuevo"
        )

    def ubicacion(self) -> dict:
        ejes_columna = self._ejes_con_columna()
        eje_id = str(self.eje_transversal.currentData() or "")
        familia_transversal = (
            "y" if self.direccion.currentData() == "x" else "x"
        )
        eje_transversal = VistaPlantaPorticos._buscar_eje(
            self.ejes, familia_transversal, eje_id
        )
        if eje_transversal is None:
            raise ValueError("Seleccioná un eje transversal para ubicar el pórtico.")
        posicion_transversal = float(eje_transversal["coordenada_m"])
        coordenada_inicial = (
            float(ejes_columna[0]["coordenada_m"]) if ejes_columna else 0.0
        )
        coordenadas = [
            float(eje["coordenada_m"]) - coordenada_inicial
            for eje in ejes_columna
        ]
        luces = (
            [derecha - izquierda for izquierda, derecha in zip(coordenadas, coordenadas[1:])]
            if len(coordenadas) >= 2 else None
        )
        referencia = {
            "direccion": str(self.direccion.currentData()),
            "eje_id": eje_id or None,
            "desfase_m": 0.0,
            "eje_longitudinal_id": ejes_columna[0]["id"] if ejes_columna else None,
            "desfase_longitudinal_m": 0.0,
        }
        return {
            "referencia_planta": referencia,
            "posicion_planta_m": posicion_transversal,
            "coordenadas_columnas_m": coordenadas,
            "luces_m": luces,
        }

    def _aceptar(self) -> None:
        ejes_columna = self._ejes_con_columna()
        if not self.eje_transversal.currentData():
            QMessageBox.warning(
                self, "Eje transversal requerido",
                "Seleccioná un eje transversal para ubicar el pórtico.",
            )
            return
        if not self.eje_transversal.currentData():
            QMessageBox.warning(
                self, "Eje transversal requerido",
                "Seleccioná un eje transversal para ubicar el pórtico.",
            )
            return
        coordenadas = [round(float(eje["coordenada_m"]), 6) for eje in ejes_columna]
        if len(set(coordenadas)) != len(coordenadas):
            QMessageBox.warning(
                self, "Ejes repetidos",
                "Los cruces seleccionados para columnas deben tener coordenadas distintas.",
            )
            return
        self.accept()


class EditorPorticoPorNiveles(QDialog):
    """Define columnas por entrepiso a partir de niveles y ejes de referencia."""

    def __init__(
        self,
        ejes: dict[str, list[dict]],
        niveles: list[dict],
        estructura_existente: dict | None = None,
        parent=None,
    ):
        super().__init__(parent)
        self.setWindowTitle("Geometría del pórtico por niveles")
        self.ejes = ejes
        self.estructura_existente = estructura_existente or {}
        self.niveles = sorted(
            (dict(nivel) for nivel in niveles),
            key=lambda nivel: float(nivel["cota_m"]),
        )
        self.direccion = QComboBox()
        self.direccion.addItem("Paralelo a X", "x")
        self.direccion.addItem("Paralelo a Y", "y")
        self.eje_transversal = QComboBox()

        self.pestanas = QTabWidget()
        self.tablas: list[QTableWidget] = []
        self.vista = VistaPlantaPorticos()
        self._crear_pestanas()

        instrucciones = QLabel(
            "Configurá cada tramo vertical entre niveles. En cada pestaña marcá los ejes "
            "donde nace una columna; la viga queda en el nivel superior de ese tramo. Las "
            "columnas de niveles superiores pueden desplazarse individualmente indicando su "
            "coordenada acumulada. Se conectarán a la viga inferior si caen sobre una luz o "
            "voladizo. Los pórticos existentes aparecen en gris y sus ejes transversales no "
            "se pueden volver a ocupar. El pórtico nuevo debe coincidir con un eje transversal."
        )
        instrucciones.setWordWrap(True)
        formulario = QFormLayout()
        formulario.addRow("Dirección del pórtico:", self.direccion)
        formulario.addRow("Eje transversal del pórtico:", self.eje_transversal)
        botones = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        botones.button(QDialogButtonBox.StandardButton.Ok).setText("Continuar")
        botones.accepted.connect(self._aceptar)
        botones.rejected.connect(self.reject)
        self.direccion.currentIndexChanged.connect(self._cambio_direccion)
        self.eje_transversal.currentIndexChanged.connect(self._actualizar_vista)
        self.pestanas.currentChanged.connect(self._actualizar_vista)

        caja = QVBoxLayout(self)
        caja.addWidget(instrucciones)
        caja.addLayout(formulario)
        caja.addWidget(self.pestanas, 1)
        caja.addWidget(self.vista, 2)
        caja.addWidget(botones)
        self.resize(980, 820)
        self._actualizar_ejes_transversales()
        self._actualizar_vista()

    def _actualizar_ejes_transversales(self) -> None:
        familia = "y" if self.direccion.currentData() == "x" else "x"
        actual = self.eje_transversal.currentData()
        self.eje_transversal.blockSignals(True)
        self.eje_transversal.clear()
        for eje in self.ejes.get(familia, []):
            nombre = str(eje.get("nombre", ""))
            coordenada = float(eje["coordenada_m"])
            ocupantes = self._porticos_en_plano(
                str(self.direccion.currentData()), coordenada
            )
            etiqueta = f"{nombre} ({coordenada:.3f} m)"
            if ocupantes:
                etiqueta += f" — ocupado: {', '.join(ocupantes)}"
            self.eje_transversal.addItem(
                etiqueta,
                str(eje.get("id", f"{familia}:{str(eje['nombre']).casefold()}")),
            )
            modelo = self.eje_transversal.model()
            item = modelo.item(self.eje_transversal.count() - 1) if hasattr(modelo, "item") else None
            if item is not None:
                item.setEnabled(not ocupantes)
        indice = self.eje_transversal.findData(actual)
        if indice < 0 or not self.eje_transversal.model().item(indice).isEnabled():
            indice = next(
                (
                    i for i in range(self.eje_transversal.count())
                    if self.eje_transversal.model().item(i).isEnabled()
                ),
                -1,
            )
        self.eje_transversal.setCurrentIndex(indice)
        self.eje_transversal.blockSignals(False)

    def _porticos_en_plano(self, direccion: str, posicion: float) -> list[str]:
        nombres = []
        for nombre, datos in self.estructura_existente.items():
            referencia = datos.get("referencia_planta", {}) or {}
            direccion_existente = str(referencia.get("direccion", "x")).lower()
            if direccion_existente != direccion:
                continue
            familia_transversal = "y" if direccion == "x" else "x"
            eje = VistaPlantaPorticos._buscar_eje(
                self.ejes, familia_transversal, referencia.get("eje_id")
            )
            posicion_existente = (
                float(eje["coordenada_m"])
                if eje is not None
                else float(datos.get("posicion_planta_m", 0.0))
            )
            if abs(posicion_existente - posicion) <= 1e-4:
                nombres.append(str(nombre))
        return nombres

    def _crear_pestanas(self) -> None:
        seleccion_previa = []
        for tabla in self.tablas:
            seleccion_previa.append({
                str(tabla.item(fila, 1).data(Qt.ItemDataRole.UserRole))
                for fila in range(tabla.rowCount())
                if tabla.item(fila, 0)
                and tabla.item(fila, 0).checkState() == Qt.CheckState.Checked
                and tabla.item(fila, 1)
            })
        while self.pestanas.count():
            self.pestanas.removeTab(0)
        self.tablas = []
        familia = str(self.direccion.currentData())
        ejes = sorted(
            self.ejes.get(familia, []),
            key=lambda eje: float(eje["coordenada_m"]),
        )
        for indice, (inferior, superior) in enumerate(
            zip(self.niveles, self.niveles[1:])
        ):
            nombre_inferior = str(inferior.get("nombre", f"Nivel {indice}"))
            nombre_superior = str(superior.get("nombre", f"Nivel {indice + 1}"))
            dz = float(superior["cota_m"]) - float(inferior["cota_m"])
            tabla = QTableWidget(len(ejes), 3)
            tabla.setHorizontalHeaderLabels(
                ("Columna en este eje", "Eje longitudinal", "Coordenada acumulada (m)")
            )
            tabla.verticalHeader().setVisible(False)
            tabla.horizontalHeader().setStretchLastSection(True)
            ids_ejes = {
                str(eje.get("id", f"{familia}:{str(eje.get('nombre', '')).casefold()}"))
                for eje in ejes
            }
            seleccion = seleccion_previa[indice] if indice < len(seleccion_previa) else set()
            if not seleccion.intersection(ids_ejes):
                seleccion = set()
            for fila, eje in enumerate(ejes):
                item_check = QTableWidgetItem("Columna")
                item_check.setFlags(item_check.flags() | Qt.ItemFlag.ItemIsUserCheckable)
                item_check.setCheckState(
                    Qt.CheckState.Checked
                    if not seleccion or str(eje.get("id")) in seleccion
                    else Qt.CheckState.Unchecked
                )
                tabla.setItem(fila, 0, item_check)
                nombre_eje = QTableWidgetItem(str(eje.get("nombre", "")))
                nombre_eje.setData(
                    Qt.ItemDataRole.UserRole,
                    str(eje.get(
                        "id", f"{familia}:{str(eje.get('nombre', '')).casefold()}"
                    )),
                )
                nombre_eje.setFlags(nombre_eje.flags() & ~Qt.ItemFlag.ItemIsEditable)
                tabla.setItem(fila, 1, nombre_eje)
                cota = QTableWidgetItem(f"{float(eje['coordenada_m']):g}")
                cota.setFlags(cota.flags() & ~Qt.ItemFlag.ItemIsEditable)
                if indice == 0:
                    tabla.setItem(fila, 2, cota)
                else:
                    posicion = QDoubleSpinBox()
                    posicion.setRange(-10000.0, 10000.0)
                    posicion.setDecimals(3)
                    posicion.setSingleStep(0.1)
                    posicion.setValue(float(eje["coordenada_m"]))
                    posicion.valueChanged.connect(self._actualizar_vista)
                    tabla.setCellWidget(fila, 2, posicion)
            tabla.itemChanged.connect(self._actualizar_vista)
            self.tablas.append(tabla)
            self.pestanas.addTab(
                tabla, f"{nombre_inferior} → {nombre_superior}  (Δz = {dz:g} m)"
            )

    def _cambio_direccion(self) -> None:
        self._actualizar_ejes_transversales()
        self._crear_pestanas()
        self._actualizar_vista()

    def _ejes_seleccionados(self, indice: int) -> list[dict]:
        if indice < 0 or indice >= len(self.tablas):
            return []
        tabla = self.tablas[indice]
        familia = str(self.direccion.currentData())
        seleccionados = []
        for fila in range(tabla.rowCount()):
            check = tabla.item(fila, 0)
            item_eje = tabla.item(fila, 1)
            control_cota = tabla.cellWidget(fila, 2)
            item_cota = tabla.item(fila, 2)
            if control_cota is not None:
                coordenada = float(control_cota.value())
            elif item_cota is not None:
                coordenada = float(item_cota.text())
            else:
                continue
            if (
                check is None
                or check.checkState() != Qt.CheckState.Checked
                or item_eje is None
            ):
                continue
            seleccionados.append({
                "id": str(item_eje.data(Qt.ItemDataRole.UserRole)),
                "nombre": item_eje.text(),
                "coordenada_m": coordenada,
                "familia": familia,
            })
        return sorted(seleccionados, key=lambda eje: eje["coordenada_m"])

    def _actualizar_vista(self, *_args) -> None:
        del _args
        indice = self.pestanas.currentIndex()
        if indice < 0 or not self.niveles:
            return
        seleccionados = self._ejes_seleccionados(indice)
        eje_id = str(self.eje_transversal.currentData() or "")
        familia_transversal = "y" if self.direccion.currentData() == "x" else "x"
        eje_transversal = VistaPlantaPorticos._buscar_eje(self.ejes, familia_transversal, eje_id)
        if eje_transversal is None:
            return
        posicion = float(eje_transversal["coordenada_m"])
        coordenadas = [eje["coordenada_m"] for eje in seleccionados]
        ejes_base = self._ejes_seleccionados(0)
        origen = ejes_base[0]["coordenada_m"] if ejes_base else 0.0
        columnas = {
            f"C0-{indice_col}": {"x": x - origen}
            for indice_col, x in enumerate(coordenadas)
        }
        estructura = dict(self.estructura_existente)
        estructura["Pórtico nuevo"] = {
                "posicion_planta_m": posicion,
                "columnas": columnas,
                "vigas": {},
        }
        referencia = {
            "direccion": str(self.direccion.currentData()),
            "eje_id": eje_id,
            "desfase_m": 0.0,
            "eje_longitudinal_id": seleccionados[0]["id"] if seleccionados else None,
            "desfase_longitudinal_m": 0.0,
        }
        asociaciones = {
            nombre: datos.get("referencia_planta", {}) or {}
            for nombre, datos in self.estructura_existente.items()
        }
        asociaciones["Pórtico nuevo"] = referencia
        self.vista.establecer_datos(
            self.ejes, estructura, asociaciones, "Pórtico nuevo"
        )

    def geometria(self) -> dict:
        if len(self.niveles) < 2:
            raise ValueError("Definí al menos dos niveles para crear un pórtico.")
        cotas = [float(nivel["cota_m"]) for nivel in self.niveles]
        if any(cota_fin - cota_inicio <= 1e-6 for cota_inicio, cota_fin in zip(cotas, cotas[1:])):
            raise ValueError("Las cotas de los niveles deben ser estrictamente crecientes.")
        if not self.tablas or any(len(self._ejes_seleccionados(i)) < 2 for i in range(len(self.tablas))):
            raise ValueError(
                "Cada tramo entre niveles necesita columnas en al menos dos ejes distintos."
            )
        for indice in range(len(self.tablas)):
            coordenadas = [
                round(eje["coordenada_m"], 6)
                for eje in self._ejes_seleccionados(indice)
            ]
            if len(coordenadas) != len(set(coordenadas)):
                raise ValueError(
                    f"Los ejes seleccionados entre {self.niveles[indice]['nombre']} y "
                    f"{self.niveles[indice + 1]['nombre']} deben tener coordenadas distintas."
                )
        eje_id = str(self.eje_transversal.currentData() or "")
        familia_transversal = "y" if self.direccion.currentData() == "x" else "x"
        eje_transversal = VistaPlantaPorticos._buscar_eje(self.ejes, familia_transversal, eje_id)
        if eje_transversal is None:
            raise ValueError("Seleccioná un eje transversal disponible para el pórtico.")
        posicion_transversal = float(eje_transversal["coordenada_m"])
        ocupantes = self._porticos_en_plano(
            str(self.direccion.currentData()), posicion_transversal
        )
        if ocupantes:
            raise ValueError(
                f"El eje transversal ya está ocupado por: {', '.join(ocupantes)}."
            )
        niveles = []
        cota_base = float(self.niveles[0]["cota_m"])
        for nivel in self.niveles:
            niveles.append({
                "nombre": str(nivel["nombre"]),
                "cota_m": float(nivel["cota_m"]) - cota_base,
                "cota_origen_m": float(nivel["cota_m"]),
            })
        origen_longitudinal = self._ejes_seleccionados(0)[0]["coordenada_m"]
        tramos_nivel = []
        for indice, (inferior, superior) in enumerate(zip(niveles, niveles[1:])):
            ejes_columnas = self._ejes_seleccionados(indice)
            tramos_nivel.append({
                "nivel_inferior": inferior["nombre"],
                "nivel_superior": superior["nombre"],
                "cota_inferior_m": inferior["cota_m"],
                "cota_superior_m": superior["cota_m"],
                "ejes_columna": [
                    {
                        "id": eje["id"],
                        "nombre": eje["nombre"],
                        "coordenada_m": eje["coordenada_m"],
                        "x_local_m": eje["coordenada_m"] - origen_longitudinal,
                    }
                    for eje in ejes_columnas
                ],
            })
        referencia = {
            "direccion": str(self.direccion.currentData()),
            "eje_id": eje_id,
            "desfase_m": 0.0,
            "eje_longitudinal_id": self._ejes_seleccionados(0)[0]["id"],
            "desfase_longitudinal_m": 0.0,
            "origen_longitudinal_m": origen_longitudinal,
            "niveles": niveles,
            "columnas_por_nivel": [
                {
                    "nivel": tramo["nivel_inferior"],
                    "cota_m": tramo["cota_inferior_m"],
                    "eje_ids": [eje["id"] for eje in tramo["ejes_columna"]],
                }
                for tramo in tramos_nivel
            ],
        }
        return {
            "referencia_planta": referencia,
            "posicion_planta_m": posicion_transversal,
            "niveles": niveles,
            "tramos_nivel": tramos_nivel,
        }

    def _aceptar(self) -> None:
        try:
            self.geometria()
        except ValueError as exc:
            QMessageBox.warning(self, "Faltan columnas por nivel", str(exc))
            return
        self.accept()


class EditorEjes(QDialog):
    """Crea familias de ejes X/Y con coordenadas absolutas no uniformes."""

    DIRECCIONES = {
        "x": ("Ejes X", "X"),
        "y": ("Ejes Y", "Y"),
    }

    def __init__(self, ejes: dict[str, list[dict]], parent=None):
        super().__init__(parent)
        self.setWindowTitle("Ejes de referencia en planta")
        self.ejes_originales = {
            direccion: [
                {
                    **dict(item),
                    "id": item.get(
                        "id", f"{direccion}:{str(item.get('nombre', '')).casefold()}"
                    ),
                }
                for item in ejes.get(direccion, [])
            ]
            for direccion in self.DIRECCIONES
        }
        self.tablas: dict[str, QTableWidget] = {}

        instrucciones = QLabel(
            "Definí las líneas de referencia de la planta. Cada coordenada se mide desde "
            "el mismo origen para su eje; ingresá posiciones acumuladas, no distancias entre "
            "ejes. Las separaciones pueden ser distintas. Los ejes son independientes: no "
            "obligan a que vigas, columnas, muros u otros elementos coincidan con ellos."
        )
        instrucciones.setWordWrap(True)

        pestanias = QTabWidget()
        for direccion, (titulo, prefijo) in self.DIRECCIONES.items():
            pagina = QWidget()
            caja = QVBoxLayout(pagina)
            tabla = QTableWidget(0, 2)
            tabla.setHorizontalHeaderLabels(("Nombre de eje", "Coordenada acumulada (m)"))
            tabla.horizontalHeader().setStretchLastSection(True)
            tabla.verticalHeader().setVisible(False)
            self.tablas[direccion] = tabla
            for eje in self.ejes_originales[direccion]:
                self._agregar_fila(
                    direccion,
                    str(eje.get("nombre", "")),
                    float(eje.get("coordenada_m", 0.0)),
                )
            acciones = QHBoxLayout()
            agregar = QPushButton(f"Agregar eje {prefijo}")
            quitar = QPushButton("Quitar eje seleccionado")
            if direccion == "x":
                agregar.clicked.connect(self._agregar_x)
                quitar.clicked.connect(self._quitar_x)
            else:
                agregar.clicked.connect(self._agregar_y)
                quitar.clicked.connect(self._quitar_y)
            acciones.addWidget(agregar)
            acciones.addWidget(quitar)
            acciones.addStretch(1)
            caja.addLayout(acciones)
            caja.addWidget(tabla)
            pestanias.addTab(pagina, titulo)

        self.vista = VistaEjes()
        self.vista.establecer_ejes(self.ejes_originales)
        self.tablas["x"].itemChanged.connect(self._actualizar_vista)
        self.tablas["y"].itemChanged.connect(self._actualizar_vista)

        botones = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Save | QDialogButtonBox.StandardButton.Cancel
        )
        botones.button(QDialogButtonBox.StandardButton.Save).setText("Guardar ejes")
        botones.accepted.connect(self._guardar)
        botones.rejected.connect(self.reject)

        caja = QVBoxLayout(self)
        caja.addWidget(instrucciones)
        caja.addWidget(pestanias, 1)
        caja.addWidget(self.vista)
        caja.addWidget(botones)
        self.resize(640, 680)

    def _agregar_nuevo(self, direccion: str) -> None:
        tabla = self.tablas[direccion]
        prefijo = self.DIRECCIONES[direccion][1]
        usados = {
            (tabla.item(fila, 0).text().strip().casefold()
             if tabla.item(fila, 0) else "")
            for fila in range(tabla.rowCount())
        }
        indice = 1
        while f"{prefijo}{indice}".casefold() in usados:
            indice += 1
        coordenadas = []
        for fila in range(tabla.rowCount()):
            item = tabla.item(fila, 1)
            try:
                coordenadas.append(float(item.text().replace(",", ".")))
            except (AttributeError, ValueError):
                continue
        coordenada = max(coordenadas) + 5.0 if coordenadas else 0.0
        self._agregar_fila(direccion, f"{prefijo}{indice}", coordenada)
        self._actualizar_vista()

    def _agregar_x(self) -> None:
        self._agregar_nuevo("x")

    def _agregar_y(self) -> None:
        self._agregar_nuevo("y")

    def _agregar_fila(self, direccion: str, nombre: str, coordenada: float) -> None:
        tabla = self.tablas[direccion]
        fila = tabla.rowCount()
        tabla.insertRow(fila)
        item_nombre = QTableWidgetItem(nombre)
        eje_existente = next(
            (
                eje for eje in self.ejes_originales[direccion]
                if eje.get("nombre") == nombre
            ),
            None,
        )
        item_nombre.setData(
            Qt.ItemDataRole.UserRole,
            eje_existente.get("id") if eje_existente else str(uuid.uuid4()),
        )
        tabla.setItem(fila, 0, item_nombre)
        tabla.setItem(fila, 1, QTableWidgetItem(f"{coordenada:g}"))

    def _quitar(self, direccion: str) -> None:
        tabla = self.tablas[direccion]
        fila = tabla.currentRow()
        if fila >= 0:
            tabla.removeRow(fila)
            self._actualizar_vista()

    def _quitar_x(self) -> None:
        self._quitar("x")

    def _quitar_y(self) -> None:
        self._quitar("y")

    def _leer_tablas(self) -> dict[str, list[dict]] | None:
        resultado: dict[str, list[dict]] = {"x": [], "y": []}
        ids: set[str] = set()
        for direccion, tabla in self.tablas.items():
            nombres: set[str] = set()
            coordenadas: set[float] = set()
            for fila in range(tabla.rowCount()):
                item_nombre = tabla.item(fila, 0)
                item_coordenada = tabla.item(fila, 1)
                nombre = item_nombre.text().strip() if item_nombre else ""
                if not nombre:
                    QMessageBox.warning(self, "Nombre requerido", "Cada eje necesita un nombre.")
                    return None
                if nombre.casefold() in nombres:
                    QMessageBox.warning(
                        self, "Nombre repetido",
                        f"El nombre de eje {nombre} está repetido en la familia {direccion.upper()}.",
                    )
                    return None
                try:
                    coordenada = float(item_coordenada.text().strip().replace(",", "."))
                except (AttributeError, ValueError):
                    QMessageBox.warning(
                        self, "Coordenada inválida",
                        f"Ingresá una coordenada numérica para el eje {nombre}.",
                    )
                    return None
                if not math.isfinite(coordenada) or not (-10000.0 <= coordenada <= 10000.0):
                    QMessageBox.warning(
                        self, "Coordenada fuera de rango",
                        f"La coordenada del eje {nombre} debe estar entre -10000 y 10000 m.",
                    )
                    return None
                clave_coordenada = round(coordenada, 6)
                if clave_coordenada in coordenadas:
                    QMessageBox.warning(
                        self, "Coordenada repetida",
                        f"Los ejes de la familia {direccion.upper()} no pueden compartir coordenada.",
                    )
                    return None
                nombres.add(nombre.casefold())
                coordenadas.add(clave_coordenada)
                eje_id = str(item_nombre.data(Qt.ItemDataRole.UserRole) or uuid.uuid4())
                if eje_id in ids:
                    QMessageBox.warning(
                        self, "Identificador de eje repetido",
                        "Cada línea de eje debe conservar una identidad única.",
                    )
                    return None
                ids.add(eje_id)
                resultado[direccion].append({
                    "id": eje_id,
                    "nombre": nombre,
                    "coordenada_m": coordenada,
                })
            resultado[direccion].sort(key=lambda item: item["coordenada_m"])
        return resultado

    def _actualizar_vista(self) -> None:
        ejes = self._leer_tablas_silencioso()
        if ejes is not None:
            self.vista.establecer_ejes(ejes)

    def _leer_tablas_silencioso(self) -> dict[str, list[dict]] | None:
        ejes: dict[str, list[dict]] = {"x": [], "y": []}
        for direccion, tabla in self.tablas.items():
            for fila in range(tabla.rowCount()):
                item_nombre = tabla.item(fila, 0)
                item_coordenada = tabla.item(fila, 1)
                try:
                    nombre = item_nombre.text().strip()
                    coordenada = float(item_coordenada.text().strip().replace(",", "."))
                except (AttributeError, ValueError):
                    return None
                if (
                    not nombre
                    or not math.isfinite(coordenada)
                    or not (-10000.0 <= coordenada <= 10000.0)
                ):
                    return None
                ejes[direccion].append({
                    "id": str(item_nombre.data(Qt.ItemDataRole.UserRole) or uuid.uuid4()),
                    "nombre": nombre, "coordenada_m": coordenada
                })
        return ejes

    def _guardar(self) -> None:
        ejes = self._leer_tablas()
        if ejes is None:
            return
        ejes_por_id = {
            eje["id"]: eje
            for familia in ejes.values()
            for eje in familia
        }
        estructura = rutas.cargar_estructura()
        cargas_proyecto = rutas.leer_json(rutas.CARGAS, {}) or {}
        referencias_losas = []
        for nombre_losa, elemento in (
            cargas_proyecto.get("elementos", {}) or {}
        ).items():
            panel = elemento.get("panel_ejes") or {}
            for familia, lados in (panel.get("ejes", {}) or {}).items():
                for lado, eje_guardado in (lados or {}).items():
                    eje_id = str((eje_guardado or {}).get("id", ""))
                    if eje_id:
                        referencias_losas.append(
                            (str(nombre_losa), str(familia), str(lado), eje_id)
                        )
        ejes_originales_por_id = {
            (familia, str(eje.get("id", ""))): eje
            for familia, lista in self.ejes_originales.items()
            for eje in lista
        }
        for nombre_losa, familia, lado, eje_id in referencias_losas:
            eje_nuevo = next(
                (
                    eje for eje in ejes.get(familia, [])
                    if str(eje.get("id", "")) == eje_id
                ),
                None,
            )
            eje_original = ejes_originales_por_id.get((familia, eje_id))
            if (
                eje_nuevo is None
                or eje_original is None
                or abs(
                    float(eje_nuevo["coordenada_m"])
                    - float(eje_original["coordenada_m"])
                ) > 1e-6
            ):
                QMessageBox.warning(
                    self, "Eje usado por una losa",
                    f"No se puede eliminar, reemplazar o mover el eje {familia.upper()} "
                    f"{(eje_original or {}).get('nombre', eje_id)} porque delimita "
                    f"el paño «{nombre_losa}». "
                    "Reasigná o eliminá primero esa losa.",
                )
                return
        estructura_actualizada = {
            nombre: dict(datos) for nombre, datos in estructura.items()
        }
        for nombre, datos in estructura_actualizada.items():
            referencia = datos.get("referencia_planta", {}) or {}
            eje_transversal_id = referencia.get("eje_id")
            eje_longitudinal_id = referencia.get("eje_longitudinal_id")
            for tipo, eje_id in (
                ("transversal", eje_transversal_id),
                ("longitudinal", eje_longitudinal_id),
            ):
                if eje_id and str(eje_id) not in ejes_por_id:
                    QMessageBox.warning(
                        self, "Eje asociado eliminado",
                        f"El pórtico {nombre} todavía referencia un eje {tipo} que se quitó. "
                        "Reasocialo o dejalo sin eje antes de eliminar esa referencia.",
                    )
                    return
            referencia["desfase_m"] = 0.0
            referencia["desfase_longitudinal_m"] = 0.0
            datos["referencia_planta"] = referencia
            if eje_transversal_id:
                eje = ejes_por_id[str(eje_transversal_id)]
                datos["posicion_planta_m"] = float(eje["coordenada_m"])
        try:
            rutas.guardar_json(rutas.DATOS / "ejes.json", ejes)
            if estructura_actualizada:
                rutas.guardar_estructura(estructura_actualizada)
        except OSError as exc:
            QMessageBox.critical(
                self, "No se pudieron guardar los ejes y actualizar las posiciones", str(exc)
            )
            return
        self.ejes_originales = ejes
        self.accept()


class VistaNiveles(QWidget):
    """Esquema vertical simple de las cotas de nivel."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.niveles: list[dict] = []
        self.setMinimumHeight(150)

    def establecer_niveles(self, niveles: list[dict]) -> None:
        self.niveles = niveles
        self.update()

    def paintEvent(self, event) -> None:
        super().paintEvent(event)
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.fillRect(self.rect(), QColor("#ffffff"))
        left, right = 38, max(39, self.width() - 20)
        top, bottom = 22, max(23, self.height() - 22)
        painter.setPen(QPen(QColor("#cbd5e1"), 1))
        painter.drawLine(left, bottom, right, bottom)
        if not self.niveles:
            painter.setPen(QColor("#64748b"))
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "Agregá niveles Z")
            return
        cotas = [float(nivel["cota_m"]) for nivel in self.niveles]
        minimo, maximo = min(cotas), max(cotas)
        if minimo == maximo:
            minimo -= 0.5
            maximo += 0.5
        for nivel in self.niveles:
            cota = float(nivel["cota_m"])
            y = bottom - (cota - minimo) / (maximo - minimo) * (bottom - top)
            painter.setPen(QPen(QColor("#7c3aed"), 1))
            painter.drawLine(left, int(y), right, int(y))
            painter.setPen(QColor("#5b21b6"))
            painter.drawText(left + 5, int(y) - 3, f"{nivel['nombre']}  {cota:g} m")


class EditorNiveles(QDialog):
    """Define niveles nombrados en cotas absolutas desde un mismo datum."""

    def __init__(self, niveles: list[dict], parent=None):
        super().__init__(parent)
        self.setWindowTitle("Niveles de referencia Z")
        self.tabla = QTableWidget(0, 2)
        self.tabla.setHorizontalHeaderLabels(("Nombre del nivel", "Cota acumulada Z (m)"))
        self.tabla.horizontalHeader().setStretchLastSection(True)
        self.tabla.verticalHeader().setVisible(False)
        for nivel in niveles:
            self._agregar_fila(
                str(nivel.get("nombre", "")),
                float(nivel.get("cota_m", 0.0)),
            )
        instrucciones = QLabel(
            "Definí las cotas de los niveles desde un mismo datum vertical (por ejemplo "
            "PB = 0, PA = 3.0, Azotea = 6.0 m). Son referencias de proyecto: todavía no "
            "cambian las alturas ni la geometría de los pórticos existentes."
        )
        instrucciones.setWordWrap(True)
        botones_fila = QHBoxLayout()
        agregar = QPushButton("Agregar nivel")
        quitar = QPushButton("Quitar nivel seleccionado")
        agregar.clicked.connect(self._agregar_nivel)
        quitar.clicked.connect(self._quitar_nivel)
        botones_fila.addWidget(agregar)
        botones_fila.addWidget(quitar)
        botones_fila.addStretch(1)
        self.vista = VistaNiveles()
        self.vista.establecer_niveles(self._leer_silencioso())
        self.tabla.itemChanged.connect(self._actualizar_vista)
        botones = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Save | QDialogButtonBox.StandardButton.Cancel
        )
        botones.button(QDialogButtonBox.StandardButton.Save).setText("Guardar niveles")
        botones.accepted.connect(self._guardar)
        botones.rejected.connect(self.reject)
        caja = QVBoxLayout(self)
        caja.addWidget(instrucciones)
        caja.addLayout(botones_fila)
        caja.addWidget(self.tabla, 1)
        caja.addWidget(self.vista)
        caja.addWidget(botones)
        self.resize(560, 580)

    def _agregar_fila(self, nombre: str, cota: float) -> None:
        fila = self.tabla.rowCount()
        self.tabla.insertRow(fila)
        self.tabla.setItem(fila, 0, QTableWidgetItem(nombre))
        self.tabla.setItem(fila, 1, QTableWidgetItem(f"{cota:g}"))

    def _agregar_nivel(self) -> None:
        usados = {
            self.tabla.item(fila, 0).text().strip().casefold()
            for fila in range(self.tabla.rowCount())
            if self.tabla.item(fila, 0)
        }
        if not self.tabla.rowCount() and "pb" not in usados:
            nombre, cota = "PB", 0.0
        else:
            indice = 1
            while f"N{indice}".casefold() in usados:
                indice += 1
            cotas = []
            for fila in range(self.tabla.rowCount()):
                try:
                    cotas.append(float(self.tabla.item(fila, 1).text().replace(",", ".")))
                except (AttributeError, ValueError):
                    continue
            nombre = f"N{indice}"
            cota = max(cotas) + 3.0 if cotas else 0.0
        self._agregar_fila(nombre, cota)
        self._actualizar_vista()

    def _quitar_nivel(self) -> None:
        fila = self.tabla.currentRow()
        if fila >= 0:
            self.tabla.removeRow(fila)
            self._actualizar_vista()

    def _leer_silencioso(self) -> list[dict]:
        niveles = []
        for fila in range(self.tabla.rowCount()):
            try:
                nombre = self.tabla.item(fila, 0).text().strip()
                cota = float(self.tabla.item(fila, 1).text().strip().replace(",", "."))
            except (AttributeError, ValueError):
                return []
            if nombre and math.isfinite(cota):
                niveles.append({"nombre": nombre, "cota_m": cota})
        return sorted(niveles, key=lambda nivel: nivel["cota_m"])

    def _leer_tabla(self) -> list[dict] | None:
        niveles = []
        nombres: set[str] = set()
        cotas: set[float] = set()
        for fila in range(self.tabla.rowCount()):
            nombre_item, cota_item = self.tabla.item(fila, 0), self.tabla.item(fila, 1)
            nombre = nombre_item.text().strip() if nombre_item else ""
            if not nombre:
                QMessageBox.warning(self, "Nombre requerido", "Cada nivel necesita un nombre.")
                return None
            try:
                cota = float(cota_item.text().strip().replace(",", ".")) if cota_item else math.nan
            except ValueError:
                cota = math.nan
            if not math.isfinite(cota) or not (-1000.0 <= cota <= 10000.0):
                QMessageBox.warning(
                    self, "Cota inválida",
                    f"Ingresá una cota válida entre -1000 y 10000 m para {nombre}.",
                )
                return None
            if nombre.casefold() in nombres:
                QMessageBox.warning(self, "Nombre repetido", f"El nivel {nombre} está repetido.")
                return None
            if round(cota, 6) in cotas:
                QMessageBox.warning(
                    self, "Cota repetida", "Dos niveles no pueden tener la misma cota."
                )
                return None
            nombres.add(nombre.casefold())
            cotas.add(round(cota, 6))
            niveles.append({"nombre": nombre, "cota_m": cota})
        return sorted(niveles, key=lambda nivel: nivel["cota_m"])

    def _actualizar_vista(self) -> None:
        self.vista.establecer_niveles(self._leer_silencioso())

    def _guardar(self) -> None:
        niveles = self._leer_tabla()
        if niveles is None:
            return
        try:
            rutas.guardar_json(rutas.DATOS / "niveles.json", {"niveles": niveles})
        except OSError as exc:
            QMessageBox.critical(self, "No se pudieron guardar los niveles", str(exc))
            return
        self.accept()


class EditorAsociacionPorticos(QDialog):
    """Asocia cada pórtico directamente a ejes de planta, sin desfases."""

    def __init__(self, estructura: dict, ejes: dict[str, list[dict]], parent=None):
        super().__init__(parent)
        self.setWindowTitle("Asociar pórticos a ejes")
        self.estructura = estructura
        self.ejes = ejes
        self.controles: dict[str, tuple[QComboBox, QComboBox, QComboBox]] = {}
        self.nombres_filas: list[str] = []
        instrucciones = QLabel(
            "La dirección indica el desarrollo del pórtico. El eje transversal ubica su plano; "
            "el eje longitudinal ubica la primera columna. Las posiciones coinciden con los "
            "ejes seleccionados: no se admiten desfases. Las luces se miden desde la primera "
            "columna; los voladizos pueden extenderse más allá de otros ejes."
        )
        instrucciones.setWordWrap(True)
        self.tabla = QTableWidget(len(estructura), 4)
        self.tabla.setHorizontalHeaderLabels(
            (
                "Pórtico", "Dirección", "Eje transversal", "Eje de primera columna",
            )
        )
        self.tabla.verticalHeader().setVisible(False)
        self.tabla.horizontalHeader().setStretchLastSection(True)
        for fila, (nombre, datos) in enumerate(estructura.items()):
            self.nombres_filas.append(nombre)
            item = QTableWidgetItem(nombre)
            item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.tabla.setItem(fila, 0, item)
            referencia = datos.get("referencia_planta", {}) or {}
            direccion = QComboBox()
            direccion.addItem("Paralelo a X", "x")
            direccion.addItem("Paralelo a Y", "y")
            actual = str(referencia.get("direccion", "x"))
            direccion.setCurrentIndex(max(0, direccion.findData(actual)))
            eje = QComboBox()
            eje_longitudinal = QComboBox()
            self.controles[nombre] = (direccion, eje, eje_longitudinal)
            self.tabla.setCellWidget(fila, 1, direccion)
            self.tabla.setCellWidget(fila, 2, eje)
            self.tabla.setCellWidget(fila, 3, eje_longitudinal)
            direccion.currentIndexChanged.connect(
                lambda indice, n=nombre: self._direccion_portico_cambiada(n, indice)
            )
            self._actualizar_ejes_portico(
                nombre,
                str(referencia.get("eje_id", "")),
                str(referencia.get("eje_longitudinal_id", "")),
            )
            eje.currentIndexChanged.connect(self._eje_transversal_cambiado)
            eje_longitudinal.currentIndexChanged.connect(self._actualizar_vista)

        self.vista = VistaPlantaPorticos()
        self.tabla.itemSelectionChanged.connect(self._actualizar_vista)
        self._actualizar_ejes_transversales_ocupados()
        if self.tabla.rowCount():
            self.tabla.setCurrentCell(0, 0)
        self._actualizar_vista()
        botones = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Save | QDialogButtonBox.StandardButton.Cancel
        )
        botones.button(QDialogButtonBox.StandardButton.Save).setText("Guardar asociaciones")
        botones.accepted.connect(self._guardar)
        botones.rejected.connect(self.reject)
        caja = QVBoxLayout(self)
        caja.addWidget(instrucciones)
        caja.addWidget(self.tabla, 1)
        caja.addWidget(self.vista, 2)
        caja.addWidget(botones)
        self.resize(920, 720)

    def _direccion_portico_cambiada(self, nombre: str, indice: int) -> None:
        if indice >= 0:
            selector, selector_longitudinal = self.controles[nombre][1:]
            eje_transversal = str(selector.currentData() or "")
            eje_longitudinal = str(selector_longitudinal.currentData() or "")
            self._actualizar_ejes_portico(
                nombre, eje_longitudinal, eje_transversal
            )
            self._eje_transversal_cambiado()

    def _eje_transversal_cambiado(self, *_args) -> None:
        del _args
        self._actualizar_ejes_transversales_ocupados()
        self._actualizar_vista()

    def _actualizar_ejes_transversales_ocupados(self) -> None:
        ocupaciones = []
        for nombre, controles in self.controles.items():
            direccion, selector = controles[:2]
            eje_id = str(selector.currentData() or "")
            familia = "y" if direccion.currentData() == "x" else "x"
            eje = VistaPlantaPorticos._buscar_eje(self.ejes, familia, eje_id)
            if eje is not None:
                ocupaciones.append((
                    nombre,
                    str(direccion.currentData()),
                    float(eje["coordenada_m"]),
                ))

        for nombre, controles in self.controles.items():
            direccion, selector = controles[:2]
            familia = "y" if direccion.currentData() == "x" else "x"
            modelo = selector.model()
            for indice in range(selector.count()):
                item = modelo.item(indice) if hasattr(modelo, "item") else None
                if item is None:
                    continue
                eje = VistaPlantaPorticos._buscar_eje(
                    self.ejes, familia, str(selector.itemData(indice) or "")
                )
                ocupado = eje is not None and any(
                    otro_nombre != nombre
                    and otra_direccion == str(direccion.currentData())
                    and abs(float(eje["coordenada_m"]) - otra_posicion) <= 1e-4
                    for otro_nombre, otra_direccion, otra_posicion in ocupaciones
                )
                item.setEnabled(not ocupado)

    def _actualizar_ejes_portico(
        self, nombre: str, eje_id: str = "", eje_longitudinal_id: str = ""
    ) -> None:
        direccion, selector, selector_longitudinal = self.controles[nombre]
        familia_transversal = "y" if direccion.currentData() == "x" else "x"
        familias = (
            (selector, familia_transversal, eje_id),
            (selector_longitudinal, str(direccion.currentData()), eje_longitudinal_id),
        )
        for combo, familia, eje_actual in familias:
            combo.blockSignals(True)
            combo.clear()
            combo.addItem("Seleccionar eje…", "")
            for item in self.ejes.get(familia, []):
                combo.addItem(
                    f"{item['nombre']} ({float(item['coordenada_m']):.3f} m)",
                    str(item.get("id", f"{familia}:{str(item['nombre']).casefold()}")),
                )
            indice = combo.findData(eje_actual)
            combo.setCurrentIndex(indice if indice >= 0 else 0)
            combo.blockSignals(False)

    def _actualizar_vista(self) -> None:
        asociaciones = {}
        for nombre, (direccion, selector, selector_longitudinal) in self.controles.items():
            asociaciones[nombre] = {
                "direccion": direccion.currentData(),
                "eje_id": selector.currentData(),
                "desfase_m": 0.0,
                "eje_longitudinal_id": selector_longitudinal.currentData(),
                "desfase_longitudinal_m": 0.0,
            }
        fila = self.tabla.currentRow()
        seleccionado = self.nombres_filas[fila] if 0 <= fila < len(self.nombres_filas) else None
        self.vista.establecer_datos(
            self.ejes, self.estructura, asociaciones, seleccionado
        )

    def _guardar(self) -> None:
        ejes_por_id = {
            str(item.get("id", f"{direccion}:{str(item['nombre']).casefold()}")): item
            for direccion, items in self.ejes.items()
            for item in items
        }
        nuevos_datos = {nombre: dict(datos) for nombre, datos in self.estructura.items()}
        planos_ocupados = []
        for nombre, (direccion, selector, selector_longitudinal) in self.controles.items():
            eje_id = str(selector.currentData() or "")
            eje_longitudinal_id = str(selector_longitudinal.currentData() or "")
            if not eje_id or not eje_longitudinal_id:
                QMessageBox.warning(
                    self, "Faltan ejes",
                    f"Seleccioná el eje transversal y el eje de primera columna para {nombre}.",
                )
                return
            referencia = {
                "direccion": str(direccion.currentData()),
                "eje_id": eje_id,
                "desfase_m": 0.0,
                "eje_longitudinal_id": eje_longitudinal_id,
                "desfase_longitudinal_m": 0.0,
            }
            nuevos_datos[nombre]["referencia_planta"] = referencia
            if eje_id:
                eje = ejes_por_id.get(eje_id)
                if eje is None:
                    QMessageBox.warning(
                        self, "Eje no disponible",
                        f"No se encontró el eje seleccionado para {nombre}.",
                    )
                    return
                posicion = float(eje["coordenada_m"])
                if any(
                    otra_direccion == str(direccion.currentData())
                    and abs(otra_posicion - posicion) <= 1e-4
                    for otra_direccion, otra_posicion in planos_ocupados
                ):
                    QMessageBox.warning(
                        self, "Eje transversal ocupado",
                        f"El eje transversal seleccionado para {nombre} ya está ocupado "
                        "por otro pórtico paralelo.",
                    )
                    return
                planos_ocupados.append((str(direccion.currentData()), posicion))
                nuevos_datos[nombre]["posicion_planta_m"] = (
                    posicion
                )
            if eje_longitudinal_id and eje_longitudinal_id not in ejes_por_id:
                QMessageBox.warning(
                    self, "Eje no disponible",
                    f"No se encontró el eje longitudinal seleccionado para {nombre}.",
                )
                return
        try:
            rutas.guardar_estructura(nuevos_datos)
        except OSError as exc:
            QMessageBox.critical(
                self, "No se pudieron guardar las asociaciones", str(exc)
            )
            return
        self.estructura = nuevos_datos
        self.accept()

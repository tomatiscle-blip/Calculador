"""Alta geométrica de paños de losa vinculados a pórticos, muros y ejes."""

from __future__ import annotations

import math

from PySide6.QtCore import QPoint, Qt
from PySide6.QtGui import QColor, QPainter, QPen, QPolygon
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QVBoxLayout,
    QWidget,
)

from app.ejes import VistaPlantaPorticos


class VistaPanelLosa(QWidget):
    """Vista de planta para elegir cuatro ejes y revisar apoyos unidireccionales."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.ejes: dict[str, list[dict]] = {"x": [], "y": []}
        self.porticos: dict[str, dict] = {}
        self.muros: dict[str, dict] = {}
        self.ejes_seleccionados: dict[str, list[str]] = {"x": [], "y": []}
        self.direccion_luz = "x"
        self.coordenadas_paño: tuple[float, float, float, float] | None = None
        self.apoyos: dict[str, str] = {}
        self.muros_apoyo: set[str] = set()
        self.setMinimumSize(560, 320)

    def establecer_datos(
        self,
        ejes: dict[str, list[dict]],
        estructura: dict,
        ejes_seleccionados: dict[str, list[str]],
        direccion_luz: str,
        coordenadas_paño: tuple[float, float, float, float] | None,
        apoyos: dict[str, dict],
        cota_apoyo: float | None = None,
    ) -> None:
        self.ejes = ejes
        self.ejes_seleccionados = {
            familia: list(ids) for familia, ids in ejes_seleccionados.items()
        }
        self.direccion_luz = direccion_luz
        self.coordenadas_paño = coordenadas_paño
        self.cota_apoyo = cota_apoyo
        self.porticos = {}
        self.muros = {}
        self.muros_apoyo = set()
        self.apoyos = {}
        for lado, apoyo in (apoyos or {}).items():
            if not isinstance(apoyo, dict):
                continue
            if str(apoyo.get("tipo", "")) == "portico":
                self.apoyos[str(apoyo.get("portico", ""))] = str(apoyo.get("viga", ""))
            elif str(apoyo.get("tipo", "")) == "muro":
                self.muros_apoyo.add(str(apoyo.get("nombre", "")))
        for nombre, datos in estructura.items():
            referencia = datos.get("referencia_planta", {}) or {}
            direccion = str(referencia.get("direccion", "x")).lower()
            if direccion not in ("x", "y"):
                continue
            familia = "y" if direccion == "x" else "x"
            eje = VistaPlantaPorticos._buscar_eje(
                ejes, familia, referencia.get("eje_id")
            )
            if eje is None:
                continue
            longitudinal = VistaPlantaPorticos._buscar_eje(
                ejes, direccion, referencia.get("eje_longitudinal_id")
            )
            origen = float(
                referencia.get(
                    "origen_longitudinal_m",
                    longitudinal["coordenada_m"] if longitudinal else 0.0,
                )
            )
            if str(datos.get("tipo", "")).lower() == "muro":
                # Los muros son flexibles: se muestran independientemente del nivel de apoyo.
                try:
                    fin_muro = float(datos.get("x_fin", origen))
                except (TypeError, ValueError):
                    fin_muro = origen
                self.muros[str(nombre)] = {
                    "direccion": direccion,
                    "posicion": float(eje["coordenada_m"]),
                    "origen": origen,
                    "fin": fin_muro,
                    "espesor_m": float(datos.get("espesor_m", 0.0)),
                    "altura_m": float(datos.get("altura_m", 0.0)),
                }
                continue
            # Pórtico: solo mostrar vigas cuya cota coincida con el nivel de apoyo seleccionado.
            for viga_id, viga in (datos.get("vigas", {}) or {}).items():
                try:
                    cota_viga = float(viga.get("cota_m"))
                except (TypeError, ValueError):
                    continue
                if cota_apoyo is not None and abs(cota_viga - cota_apoyo) > 1e-4:
                    continue
                tramos = [
                    tramo
                    for tramo in viga.get("tramos", [])
                    if isinstance(tramo, dict)
                    and math.isfinite(float(tramo.get("x_inicio", 0.0)))
                    and math.isfinite(float(tramo.get("x_fin", 0.0)))
                ]
                if not tramos:
                    continue
                inicio = origen + min(
                    min(float(t["x_inicio"]), float(t["x_fin"])) for t in tramos
                )
                fin = origen + max(
                    max(float(t["x_inicio"]), float(t["x_fin"])) for t in tramos
                )
                self.porticos.setdefault(str(nombre), {
                    "direccion": direccion,
                    "posicion": float(eje["coordenada_m"]),
                    "origen": origen,
                    "vigas": {},
                })["vigas"][str(viga_id)] = {
                    "nivel_nombre": str(viga.get("nivel_nombre", "")),
                    "cota_m": viga.get("cota_m"),
                    "inicio": inicio,
                    "fin": fin,
                }
                tramos = [
                    tramo for tramo in viga.get("tramos", [])
                    if isinstance(tramo, dict)
                    and math.isfinite(float(tramo.get("x_inicio", 0.0)))
                    and math.isfinite(float(tramo.get("x_fin", 0.0)))
                ]
                if not tramos:
                    continue
                inicio = origen + min(
                    min(float(t["x_inicio"]), float(t["x_fin"])) for t in tramos
                )
                fin = origen + max(
                    max(float(t["x_inicio"]), float(t["x_fin"])) for t in tramos
                )
                self.porticos.setdefault(str(nombre), {
                    "direccion": direccion,
                    "posicion": float(eje["coordenada_m"]),
                    "origen": origen,
                    "vigas": {},
                })["vigas"][str(viga_id)] = {
                    "nivel_nombre": str(viga.get("nivel_nombre", "")),
                    "cota_m": viga.get("cota_m"),
                    "inicio": inicio,
                    "fin": fin,
                }
        self.update()

    def mousePressEvent(self, event) -> None:
        if event.button() != Qt.MouseButton.LeftButton:
            return
        izquierda, derecha, arriba, abajo = self._limites_pantalla()
        x_items = self.ejes.get("x", [])
        y_items = self.ejes.get("y", [])
        coords_x, coords_y = self._coordenadas_referencia()
        min_x, max_x = VistaPlantaPorticos._rango_con_margen(coords_x)
        min_y, max_y = VistaPlantaPorticos._rango_con_margen(coords_y)
        px, py = float(event.position().x()), float(event.position().y())
        candidatos: list[tuple[float, str, str]] = []
        if max_x > min_x:
            for eje in x_items:
                posicion = izquierda + (
                    (float(eje["coordenada_m"]) - min_x) / (max_x - min_x)
                ) * (derecha - izquierda)
                candidatos.append((
                    abs(px - posicion), "x",
                    str(eje.get("id", f"x:{str(eje.get('nombre', '')).casefold()}")),
                ))
        if max_y > min_y:
            for eje in y_items:
                posicion = abajo - (
                    (float(eje["coordenada_m"]) - min_y) / (max_y - min_y)
                ) * (abajo - arriba)
                candidatos.append((
                    abs(py - posicion), "y",
                    str(eje.get("id", f"y:{str(eje.get('nombre', '')).casefold()}")),
                ))
        if not candidatos:
            return
        distancia, familia, eje_id = min(candidatos)
        if distancia > 14:
            return
        seleccionados = self.ejes_seleccionados[familia]
        if eje_id in seleccionados:
            seleccionados.remove(eje_id)
        elif len(seleccionados) < 2:
            seleccionados.append(eje_id)
        self.update()
        callback = getattr(self, "al_cambiar_ejes", None)
        if callback:
            callback()

    def _coordenadas_referencia(self) -> tuple[list[float], list[float]]:
        coords_x = [float(e["coordenada_m"]) for e in self.ejes.get("x", [])]
        coords_y = [float(e["coordenada_m"]) for e in self.ejes.get("y", [])]
        for datos in self.porticos.values():
            familia_longitudinal = "x" if datos["direccion"] == "x" else "y"
            extremos = [
                extremo
                for viga in datos["vigas"].values()
                for extremo in (viga["inicio"], viga["fin"])
            ]
            if familia_longitudinal == "x":
                coords_y.append(datos["posicion"])
                coords_x.extend(extremos)
            else:
                coords_x.append(datos["posicion"])
                coords_y.extend(extremos)
        for datos in self.muros.values():
            familia_longitudinal = "x" if datos["direccion"] == "x" else "y"
            if familia_longitudinal == "x":
                coords_y.append(datos["posicion"])
            else:
                coords_x.append(datos["posicion"])
        return coords_x, coords_y

    def _limites_pantalla(self) -> tuple[float, float, float, float]:
        return 66.0, max(67.0, self.width() - 24.0), 38.0, max(39.0, self.height() - 45.0)

    def paintEvent(self, event) -> None:
        super().paintEvent(event)
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.fillRect(self.rect(), QColor("#ffffff"))
        painter.setPen(QColor("#334155"))
        painter.drawText(12, 20, "Planta Â· clic en dos ejes X y dos ejes Y para delimitar el paño")

        izquierda, derecha, arriba, abajo = self._limites_pantalla()
        coords_x, coords_y = self._coordenadas_referencia()
        min_x, max_x = VistaPlantaPorticos._rango_con_margen(coords_x)
        min_y, max_y = VistaPlantaPorticos._rango_con_margen(coords_y)

        def punto(x: float, y: float) -> tuple[int, int]:
            px = izquierda + (x - min_x) / (max_x - min_x) * (derecha - izquierda)
            py = abajo - (y - min_y) / (max_y - min_y) * (abajo - arriba)
            return int(px), int(py)

        painter.setPen(QPen(QColor("#cbd5e1"), 1))
        painter.drawRect(
            int(izquierda), int(arriba), int(derecha - izquierda), int(abajo - arriba)
        )
        for familia, color, etiqueta, vertical in (
            ("x", QColor("#2563eb"), "#1d4ed8", True),
            ("y", QColor("#15803d"), "#166534", False),
        ):
            for eje in self.ejes.get(familia, []):
                eje_id = str(eje.get("id", f"{familia}:{str(eje.get('nombre', '')).casefold()}"))
                seleccionado = eje_id in self.ejes_seleccionados[familia]
                lapiz = QPen(color, 2 if seleccionado else 1, Qt.PenStyle.DashLine)
                painter.setPen(lapiz)
                coord = float(eje["coordenada_m"])
                if vertical:
                    px, _ = punto(coord, min_y)
                    painter.drawLine(px, int(arriba), px, int(abajo))
                    painter.setPen(QColor(etiqueta))
                    painter.drawText(px - 12, int(abajo + 16), str(eje["nombre"]))
                else:
                    _, py = punto(min_x, coord)
                    painter.drawLine(int(izquierda), py, int(derecha), py)
                    painter.setPen(QColor(etiqueta))
                    painter.drawText(8, py + 4, str(eje["nombre"]))

        if self.coordenadas_paño is not None:
            x0, x1, y0, y1 = self.coordenadas_paño
            p0 = punto(x0, y1)
            p1 = punto(x1, y0)
            painter.setPen(QPen(QColor("#2563eb"), 2))
            painter.setBrush(QColor(59, 130, 246, 48))
            painter.drawRect(p0[0], p0[1], p1[0] - p0[0], p1[1] - p0[1])
            centro = punto((x0 + x1) / 2, (y0 + y1) / 2)
            fin = punto(
                (x0 + x1) / 2 + (x1 - x0) * 0.18
                if self.direccion_luz == "x" else (x0 + x1) / 2,
                (y0 + y1) / 2
                if self.direccion_luz == "x"
                else (y0 + y1) / 2 + (y1 - y0) * 0.18,
            )
            painter.setPen(QPen(QColor("#dc2626"), 3))
            painter.drawLine(*centro, *fin)
            painter.drawText(centro[0] + 6, centro[1] - 6, f"Luz {self.direccion_luz.upper()}")

        for nombre, datos in self.porticos.items():
            for viga_id, viga in datos["vigas"].items():
                if datos["direccion"] == "x":
                    p_inicio = punto(viga["inicio"], datos["posicion"])
                    p_fin = punto(viga["fin"], datos["posicion"])
                else:
                    p_inicio = punto(datos["posicion"], viga["inicio"])
                    p_fin = punto(datos["posicion"], viga["fin"])
                es_apoyo = self.apoyos.get(nombre) == viga_id
                painter.setPen(QPen(
                    QColor("#15803d") if es_apoyo else QColor("#64748b"),
                    4 if es_apoyo else 1,
                ))
                painter.drawLine(*p_inicio, *p_fin)
                if es_apoyo:
                    centro_x = (p_inicio[0] + p_fin[0]) // 2
                    centro_y = (p_inicio[1] + p_fin[1]) // 2
                    painter.setBrush(QColor("#15803d"))
                    painter.setPen(QPen(QColor("#166534"), 1))
                    painter.drawPolygon(QPolygon([
                        QPoint(centro_x - 7, centro_y - 9),
                        QPoint(centro_x + 7, centro_y - 9),
                        QPoint(centro_x, centro_y),
                    ]))
                    painter.drawText(centro_x + 8, centro_y - 8, f"{nombre} Â· {viga_id}")

        for nombre, datos in self.muros.items():
            direccion = datos["direccion"]
            espesor = datos["espesor_m"]
            origen = float(datos.get("origen", 0.0))
            fin = float(datos.get("fin", origen))
            if direccion == "x":
                p_inicio = punto(origen, datos["posicion"])
                p_fin = punto(fin, datos["posicion"])
            else:
                p_inicio = punto(datos["posicion"], origen)
                p_fin = punto(datos["posicion"], fin)
            es_apoyo_muro = nombre in self.muros_apoyo
            color_linea = QColor("#2563eb") if es_apoyo_muro else QColor("#94a3b8")
            pintor_linea = QPen(color_linea, 5 if es_apoyo_muro else 2)
            painter.setPen(pintor_linea)
            painter.drawLine(*p_inicio, *p_fin)
            if es_apoyo_muro:
                centro_x = (p_inicio[0] + p_fin[0]) // 2
                centro_y = (p_inicio[1] + p_fin[1]) // 2
                painter.setBrush(QColor("#2563eb"))
                painter.setPen(QPen(QColor("#1d4ed8"), 1))
                painter.drawPolygon(QPolygon([
                    QPoint(centro_x - 7, centro_y - 9),
                    QPoint(centro_x + 7, centro_y - 9),
                    QPoint(centro_x, centro_y),
                ]))
                espesor_txt = f"{espesor:.2f} m" if espesor > 0 else ""
                painter.drawText(
                    centro_x + 8, centro_y - 8,
                    f"Muro {nombre}{' Â· ' + espesor_txt if espesor_txt else ''}"
                )

        pintor_leyenda = QColor("#64748b")
        partes_leyenda = [
            "Verde: viga de pórtico seleccionada como apoyo",
            "Azul: muro estructural seleccionado como apoyo",
            "Gris: elemento disponible no seleccionado",
            "flecha roja: dirección de trabajo",
        ]
        painter.setPen(pintor_leyenda)
        painter.drawText(
            12, self.height() - 12,
            "  Â·  ".join(partes_leyenda),
        )

class EditorLosaDesdeEjes(QDialog):
    """Delimita un paño, vincula apoyos (pórtico+viga o muro) y devuelve su geometría."""

    def __init__(
        self,
        ejes: dict[str, list[dict]],
        niveles: list[dict],
        estructura: dict,
        parent=None,
        elemento: dict | None = None,
        nombre: str = "",
    ):
        super().__init__(parent)
        self.setWindowTitle("Geometría y apoyos de losa")
        self.ejes = ejes
        self.niveles = niveles
        self.estructura = estructura
        self.elemento = dict(elemento or {})
        self.nombre = QLineEdit(nombre)
        self.nombre.setPlaceholderText("Nombre del paño")
        self.tipologia = QComboBox()
        self.tipologia.addItem("Alivianada unidireccional", "alivianada")
        self.tipologia.addItem("Maciza", "maciza")
        self.tipologia.addItem("Casetonada (cálculo pendiente)", "casetonada")
        tipo_actual = str(self.elemento.get("tipologia", "alivianada"))
        indice_tipo = self.tipologia.findData(tipo_actual)
        self.tipologia.setCurrentIndex(max(indice_tipo, 0))
        self.direccion_luz = QComboBox()
        self.direccion_luz.addItem("Luz en X Â· apoyos en ejes X", "x")
        self.direccion_luz.addItem("Luz en Y Â· apoyos en ejes Y", "y")
        direccion_previa = str(
            (self.elemento.get("panel_ejes") or {}).get("direccion_luz", "x")
        )
        self.direccion_luz.setCurrentIndex(
            max(self.direccion_luz.findData(direccion_previa), 0)
        )
        self.nivel_apoyo = QComboBox()
        for cota, nivel in self._niveles_de_vigas():
            self.nivel_apoyo.addItem(f"{nivel} Â· {cota:g} m", cota)
        nivel_previo = (self.elemento.get("panel_ejes") or {}).get("cota_apoyo_m")
        if nivel_previo is not None:
            indice = self.nivel_apoyo.findData(float(nivel_previo))
            if indice >= 0:
                self.nivel_apoyo.setCurrentIndex(indice)

        self.desfase_inicio = self._spin()
        self.desfase_fin = self._spin()
        panel_previo = self.elemento.get("panel_ejes") or {}
        self.desfase_inicio.setValue(float(panel_previo.get("desfase_inicio_m", 0.0)))
        self.desfase_fin.setValue(float(panel_previo.get("desfase_fin_m", 0.0)))
        self.selector_vigas = {
            "izq": QComboBox(),
            "der": QComboBox(),
        }
        self.etiqueta_apoyos = QLabel()
        self.etiqueta_apoyos.setWordWrap(True)
        self.estado = QLabel()
        self.estado.setWordWrap(True)
        self.vista = VistaPanelLosa()
        self.vista.al_cambiar_ejes = self._actualizar
        ejes_iniciales = panel_previo.get("ejes", {}) or {}
        self.vista.ejes_seleccionados = {
            familia: [
                str(ejes_iniciales.get(familia, {}).get(lado, {}).get("id"))
                for lado in ("inicio", "fin")
                if ejes_iniciales.get(familia, {}).get(lado, {}).get("id")
            ]
            for familia in ("x", "y")
        }
        self.formulario = QFormLayout()
        self.formulario.addRow("Nombre del paño:", self.nombre)
        self.formulario.addRow("Tipología:", self.tipologia)
        self.formulario.addRow("Dirección de trabajo:", self.direccion_luz)
        self.formulario.addRow("Nivel de las vigas de apoyo:", self.nivel_apoyo)
        self.formulario.addRow("Apoyo en el primer borde:", self.selector_vigas["izq"])
        self.formulario.addRow("Apoyo en el segundo borde:", self.selector_vigas["der"])
        self.formulario.addRow("Desfase del borde libre inicial (m):", self.desfase_inicio)
        self.formulario.addRow("Desfase del borde libre final (m):", self.desfase_fin)
        self.formulario.addRow("Apoyos encontrados:", self.etiqueta_apoyos)
        self.formulario.addRow(self.estado)
        self.direccion_luz.currentIndexChanged.connect(self._actualizar)
        self.nivel_apoyo.currentIndexChanged.connect(self._actualizar)
        for selector in self.selector_vigas.values():
            selector.currentIndexChanged.connect(self._actualizar)
        self.desfase_inicio.valueChanged.connect(self._actualizar)
        self.desfase_fin.valueChanged.connect(self._actualizar)
        self.nombre.textChanged.connect(self._actualizar)
        botones = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        botones.accepted.connect(self._aceptar)
        botones.rejected.connect(self.reject)
        caja = QVBoxLayout(self)
        instrucciones = QLabel(
            "Hacé clic en dos ejes X y dos ejes Y para dibujar el paño. Elegí la dirección "
            "de trabajo: los bordes de apoyo (izq/der) se detectan automáticamente según "
            "lo que haya en esos ejes (muro estructural o pórtico con viga). Los bordes "
            "libres pueden desplazarse con los desfases."
        )
        instrucciones.setWordWrap(True)
        caja.addLayout(self.formulario)
        caja.addWidget(instrucciones)
        caja.addWidget(self.vista, 1)
        caja.addWidget(botones)
        self.resize(1040, 760)
        self._actualizar()

    @staticmethod
    def _spin() -> QDoubleSpinBox:
        control = QDoubleSpinBox()
        control.setRange(-10000.0, 10000.0)
        control.setDecimals(3)
        control.setSingleStep(0.1)
        control.setSuffix(" m")
        return control

    def _niveles_de_vigas(self) -> list[tuple[float, str]]:
        niveles: dict[float, str] = {}
        for datos in self.estructura.values():
            for viga in (datos.get("vigas", {}) or {}).values():
                try:
                    cota = float(viga["cota_m"])
                except (KeyError, TypeError, ValueError):
                    continue
                nombre = str(viga.get("nivel_nombre", f"Nivel {cota:g}"))
                niveles.setdefault(cota, nombre)
        return sorted((cota, nombre) for cota, nombre in niveles.items())

    def _elemento_en_eje(self, eje: dict, direccion_apoyo: str) -> dict | None:
        """Detecta qué hay apoyado en un eje: muro estructural o pórtico.

        Busca en self.estructura un elemento cuya referencia_planta coincida
        con el eje (mismo eje_id y coordenada) y cuya dirección de planta
        sea paralela al borde de apoyo.
        """
        eje_id = str(eje.get("id", ""))
        for nombre, datos in self.estructura.items():
            if not isinstance(datos, dict):
                continue
            ref = datos.get("referencia_planta", {}) or {}
            ref_eje_id = str(ref.get("eje_id", ""))
            if ref_eje_id != eje_id:
                continue
            # El eje transversal ya queda identificado por eje_id; no hay que
            # comparar la coordenada longitudinal del origen porque pertenece a
            # otro sistema de ejes (familia de columnas, no de apoyos).
            ref_direccion = str(ref.get("direccion", "x")).lower()
            if ref_direccion == direccion_apoyo:
                return {"nombre": str(nombre), "datos": datos, "tipo": str(datos.get("tipo", "portico"))}
        return None

    def _poblar_vigas(self, lado: str, eje: dict, direccion: str, cota: float) -> str | None:
        selector = self.selector_vigas[lado]
        previo = str(selector.currentData() or "")
        selector.blockSignals(True)
        selector.clear()
        elemento = self._elemento_en_eje(eje, direccion)
        if elemento and str(elemento["tipo"]).lower() == "portico":
            portico_nombre = str(elemento["nombre"])
            datos = self.estructura.get(portico_nombre, {})
            for viga_id, viga in (datos.get("vigas", {}) or {}).items():
                try:
                    cota_viga = float(viga.get("cota_m"))
                except (TypeError, ValueError):
                    continue
                if abs(cota_viga - cota) > 1e-4:
                    continue
                extremos = self._viga_extremos(portico_nombre, str(viga_id))
                if extremos is None:
                    continue
                selector.addItem(
                    f"{viga_id} Â· {extremos[0]:g}â€“{extremos[1]:g} m",
                    f"portico:{portico_nombre}:{viga_id}",
                )
            indice = selector.findData(previo)
            selector.setCurrentIndex(indice if indice >= 0 else (0 if selector.count() else -1))
        selector.blockSignals(False)
        return portico_nombre if (elemento and str(elemento["tipo"]).lower() == "portico") else None

    def _poblar_muro(self, lado: str, eje: dict, direccion: str) -> str | None:
        """Si en el eje hay un muro estructural, lo agrega al selector como única opción."""
        selector = self.selector_vigas[lado]
        previo = str(selector.currentData() or "")
        selector.blockSignals(True)
        selector.clear()
        elemento = self._elemento_en_eje(eje, direccion)
        if elemento and str(elemento["tipo"]).lower() == "muro":
            datos = elemento["datos"]
            ref = datos.get("referencia_planta", {}) or {}
            ref_direccion = str(ref.get("direccion", "x")).lower()
            espesor = float(datos.get("espesor_m", 0.0))
            altura = float(datos.get("altura_m", 0.0))
            selector.addItem(
                f"Muro {str(elemento['nombre'])} ({ref_direccion}, {espesor:.2f}Ã—{altura:.2f} m)",
                f"muro:{str(elemento['nombre'])}",
            )
            selector.setCurrentIndex(0 if selector.count() else -1)
        selector.blockSignals(False)
        return str(elemento["nombre"]) if (elemento and str(elemento["tipo"]).lower() == "muro") else None

    def _eje_por_id(self, familia: str, eje_id: str) -> dict | None:
        return VistaPlantaPorticos._buscar_eje(self.ejes, familia, eje_id)

    def _coordenadas_ejes(self) -> tuple[dict, dict] | None:
        elegidos = {}
        for familia in ("x", "y"):
            ejes = [
                self._eje_por_id(familia, eje_id)
                for eje_id in self.vista.ejes_seleccionados[familia]
            ]
            ejes = [eje for eje in ejes if eje is not None]
            if len(ejes) != 2:
                return None
            elegidos[familia] = sorted(
                ejes, key=lambda eje: float(eje["coordenada_m"])
            )
            if abs(
                float(elegidos[familia][1]["coordenada_m"])
                - float(elegidos[familia][0]["coordenada_m"])
            ) <= 1e-6:
                return None
        return elegidos["x"], elegidos["y"]

    def _portico_en_eje(
        self, direccion_portico: str, familia_eje: str, eje: dict
    ) -> str | None:
        ocupantes = []
        for nombre, datos in self.estructura.items():
            referencia = datos.get("referencia_planta", {}) or {}
            if str(referencia.get("direccion", "x")) != direccion_portico:
                continue
            if float(referencia.get("desfase_m", 0.0) or 0.0) != 0.0:
                continue
            eje_actual = self._eje_por_id(
                familia_eje, str(referencia.get("eje_id", ""))
            )
            if eje_actual and abs(
                float(eje_actual["coordenada_m"]) - float(eje["coordenada_m"])
            ) <= 1e-4:
                ocupantes.append(str(nombre))
        return ocupantes[0] if len(ocupantes) == 1 else None

    def _viga_extremos(self, portico: str, viga_id: str) -> tuple[float, float] | None:
        datos = self.estructura.get(portico, {})
        referencia = datos.get("referencia_planta", {}) or {}
        direccion = str(referencia.get("direccion", "x"))
        familia = direccion
        eje = self._eje_por_id(familia, str(referencia.get("eje_longitudinal_id", "")))
        origen = float(
            referencia.get(
                "origen_longitudinal_m",
                eje["coordenada_m"] if eje else 0.0,
            )
        )
        viga = (datos.get("vigas", {}) or {}).get(viga_id, {})
        tramos = viga.get("tramos", [])
        if not tramos:
            return None
        intervalos = [
            (
                min(float(t["x_inicio"]), float(t["x_fin"])),
                max(float(t["x_inicio"]), float(t["x_fin"])),
            )
            for t in tramos
            if isinstance(t, dict)
            and math.isfinite(float(t.get("x_inicio", 0.0)))
            and math.isfinite(float(t.get("x_fin", 0.0)))
        ]
        if not intervalos:
            return None
        return origen + min(i[0] for i in intervalos), origen + max(i[1] for i in intervalos)

    def _viga_cubre(
        self, portico: str, viga_id: str, inicio: float, fin: float
    ) -> bool:
        extremos = self._viga_extremos(portico, viga_id)
        if extremos is None:
            return False
        vi, vf = extremos
        cubierto = vi
        for ti, tf in sorted(((vi, vf),), key=lambda p: p[0]):
            if ti > cubierto + 1e-6:
                break
            cubierto = max(cubierto, tf)
        return cubierto >= fin - 1e-6

    def _apoyo_texto(self, lado: str) -> str:
        dato = str(self.selector_vigas[lado].currentData() or "")
        if dato.startswith("muro:"):
            return f"{lado}: muro {dato.split(':', 1)[1]}"
        if dato.startswith("portico:"):
            partes = dato.split(":")
            return f"{lado}: {partes[1]} Â· viga {partes[2]}"
        return f"{lado}: sin apoyo"

    def _apoyo_dato(self, lado: str) -> tuple[str, str | None, str | None]:
        dato = str(self.selector_vigas[lado].currentData() or "")
        if dato.startswith("muro:"):
            return "muro", dato.split(":", 1)[1], None
        if dato.startswith("portico:"):
            partes = dato.split(":")
            return "portico", partes[1], partes[2]
        return "ninguno", None, None

    def _actualizar(self) -> None:
        direccion = str(self.direccion_luz.currentData() or "x")
        coordenadas = self._coordenadas_ejes()
        cota_apoyo = self.nivel_apoyo.currentData()
        nivel_data = float(cota_apoyo) if cota_apoyo is not None else None
        apoyos: dict[str, dict] = {}
        ejes_soporte = []
        # La dirección del elemento de apoyo es perpendicular a la luz:
        # si la luz corre en X, los bordes de apoyo están en ejes X y el
        # muro/pórtico que los soporta corre en Y (y viceversa).
        direccion_apoyo = "y" if direccion == "x" else "x"
        if coordenadas and nivel_data is not None:
            x_ejes, y_ejes = coordenadas
            ejes_soporte = x_ejes if direccion == "x" else y_ejes
            for lado, eje in zip(("izq", "der"), ejes_soporte):
                elemento = self._elemento_en_eje(eje, direccion_apoyo)
                tipo_el = str((elemento or {}).get("tipo", "")).lower()
                if tipo_el == "portico":
                    portico = self._poblar_vigas(lado, eje, direccion_apoyo, float(nivel_data))
                    dato = str(self.selector_vigas[lado].currentData() or "")
                    if dato.startswith("portico:"):
                        partes = dato.split(":")
                        apoyos[lado] = {"tipo": "portico", "portico": partes[1], "viga": partes[2]}
                elif tipo_el == "muro":
                    nombre_muro = self._poblar_muro(lado, eje, direccion_apoyo)
                    dato = str(self.selector_vigas[lado].currentData() or "")
                    if dato.startswith("muro:"):
                        apoyos[lado] = {"tipo": "muro", "nombre": dato.split(":", 1)[1]}

        intervalo_comun = None
        if coordenadas and nivel_data is not None:
            x_ejes, y_ejes = coordenadas
            ejes_libres = y_ejes if direccion == "x" else x_ejes
            coord_inicio_ejes = float(ejes_libres[0]["coordenada_m"])
            coord_fin_ejes = float(ejes_libres[1]["coordenada_m"])

            rangos = []
            hay_apoyo_valido = all(apoyos.get(lado) for lado in ("izq", "der"))
            if hay_apoyo_valido:
                for lado in ("izq", "der"):
                    if apoyos[lado].get("tipo") == "portico":
                        portico_nombre = str(apoyos[lado]["portico"])
                        viga_id = str(apoyos[lado]["viga"])
                        extremos = self._viga_extremos(portico_nombre, viga_id)
                        if extremos:
                            rangos.append(extremos)
                    # muro: no restringe el intervalo (usa el eje completo)

            if rangos:
                # Al menos un lado es pórtico: intersección de los rangos de viga
                inicio = max(r[0] for r in rangos)
                fin = min(r[1] for r in rangos)
                if fin > inicio:
                    intervalo_comun = (inicio, fin)
            elif hay_apoyo_valido:
                # Ambos lados son muro: el paño ocupa todo el tramo entre ejes
                intervalo_comun = (coord_inicio_ejes, coord_fin_ejes)

        self.desfase_inicio.blockSignals(True)
        self.desfase_fin.blockSignals(True)
        if intervalo_comun and coordenadas:
            coord_inicio = coord_inicio_ejes
            coord_fin = coord_fin_ejes
            inicio_min = max(intervalo_comun[0], coord_inicio - 10000)
            inicio_max = min(intervalo_comun[1], coord_fin - 1e-6)
            fin_min = max(intervalo_comun[0], inicio_min + 1e-6)
            fin_max = min(intervalo_comun[1], coord_fin + 10000)
            self.desfase_inicio.setRange(
                coord_inicio - inicio_max, coord_inicio - inicio_min
            )
            self.desfase_fin.setRange(
                fin_min - coord_fin, fin_max - coord_fin
            )
            panel_previo = self.elemento.get("panel_ejes") or {}
            if not getattr(self, "_rangos_inicializados", False):
                self.desfase_inicio.setValue(
                    float(panel_previo.get("desfase_inicio_m", 0.0))
                )
                self.desfase_fin.setValue(
                    float(panel_previo.get("desfase_fin_m", 0.0))
                )
                self._rangos_inicializados = True
            inicio_paño = coord_inicio - float(self.desfase_inicio.value())
            fin_paño = coord_fin + float(self.desfase_fin.value())
            if (
                inicio_paño < intervalo_comun[0] - 1e-6
                or fin_paño > intervalo_comun[1] + 1e-6
                or inicio_paño >= fin_paño
            ):
                intervalo_comun = None
        else:
            self.desfase_inicio.setRange(0.0, 0.0)
            self.desfase_fin.setRange(0.0, 0.0)
            inicio_paño = fin_paño = 0.0
        self.desfase_inicio.blockSignals(False)
        self.desfase_fin.blockSignals(False)

        coordenadas_paño = None
        nombres_ejes = "Seleccioná dos ejes X y dos ejes Y."
        if coordenadas:
            x_ejes, y_ejes = coordenadas
            x0 = float(x_ejes[0]["coordenada_m"])
            x1 = float(x_ejes[1]["coordenada_m"])
            y0 = float(y_ejes[0]["coordenada_m"])
            y1 = float(y_ejes[1]["coordenada_m"])
            coordenadas_paño = (
                x0 if direccion == "x" else inicio_paño,
                x1 if direccion == "x" else fin_paño,
                inicio_paño if direccion == "x" else y0,
                fin_paño if direccion == "x" else y1,
            )
            nombres_ejes = (
                f"Paño X: {x_ejes[0]['nombre']}â€“{x_ejes[1]['nombre']} Â· "
                f"Y: {y_ejes[0]['nombre']}â€“{y_ejes[1]['nombre']}"
            )

        self.etiqueta_apoyos.setText(" Â· ".join(
            self._apoyo_texto(lado) for lado in ("izq", "der")
        ))
        if not coordenadas:
            estado = nombres_ejes
        elif not all(apoyos.get(lado) for lado in ("izq", "der")):
            faltan = [lado for lado in ("izq", "der") if not apoyos.get(lado)]
            estado = f"Falta apoyo en el borde {', '.join(faltan)} (no hay muro ni pórtico en ese eje)."
        elif intervalo_comun is None:
            ambos_muros = all(
                apoyos.get(lado, {}).get("tipo") == "muro" for lado in ("izq", "der")
            )
            if not ambos_muros:
                estado = "Las vigas de apoyo no cubren el ancho completo del paño."
            else:
                estado = (
                    f"{nombres_ejes} Â· apoyos en muro estructural"
                )
        else:
            estado = (
                f"{nombres_ejes} Â· luz {direccion.upper()} "
                f"{abs(float(ejes_soporte[1]['coordenada_m']) - float(ejes_soporte[0]['coordenada_m'])):.3f} m"
            )
        self.estado.setText(estado)
        self.vista.establecer_datos(
            self.ejes, self.estructura, self.vista.ejes_seleccionados,
            direccion, coordenadas_paño, apoyos, cota_apoyo,
        )

    def resultado(self) -> tuple[str, dict]:
        nombre = self.nombre.text().strip()
        if not nombre:
            raise ValueError("Ingresá un nombre para el paño.")
        coordenadas = self._coordenadas_ejes()
        if coordenadas is None:
            raise ValueError("Seleccioná exactamente dos ejes X y dos ejes Y.")
        x_ejes, y_ejes = coordenadas
        direccion = str(self.direccion_luz.currentData())
        ejes_soporte = x_ejes if direccion == "x" else y_ejes
        apoyos: dict[str, dict] = {}
        for lado, eje in zip(("izq", "der"), ejes_soporte):
            tipo, nombre_el, viga_id = self._apoyo_dato(lado)
            if tipo == "muro":
                if not nombre_el:
                    raise ValueError(
                        f"Seleccioná un muro estructural para el apoyo {lado}."
                    )
                apoyos[lado] = {"tipo": "muro", "nombre": nombre_el}
            elif tipo == "portico":
                if not nombre_el or not viga_id:
                    raise ValueError(
                        f"Falta un pórtico y una viga de apoyo en el eje {eje['nombre']}."
                    )
                apoyos[lado] = {"tipo": "portico", "portico": nombre_el, "viga": viga_id}
            else:
                raise ValueError(
                    f"No hay apoyo válido en el borde {lado} (eje {eje['nombre']})."
                )
        span = abs(
            float(ejes_soporte[1]["coordenada_m"])
            - float(ejes_soporte[0]["coordenada_m"])
        )
        ejes_libres = y_ejes if direccion == "x" else x_ejes
        eje_inicio, eje_fin = ejes_libres
        coord_inicio = float(eje_inicio["coordenada_m"]) - self.desfase_inicio.value()
        coord_fin = float(eje_fin["coordenada_m"]) + self.desfase_fin.value()
        if coord_fin <= coord_inicio:
            raise ValueError("Los desfases hacen que el ancho del paño sea nulo o negativo.")
        es_portico = lambda lado: apoyos.get(lado, {}).get("tipo") == "portico"
        if es_portico("izq") and es_portico("der"):
            if not all(
                self._viga_cubre(
                    str(apoyos[lado]["portico"]),
                    str(apoyos[lado]["viga"]),
                    coord_inicio,
                    coord_fin,
                )
                for lado in ("izq", "der")
            ):
                raise ValueError(
                    "Las vigas de apoyo no cubren el ancho completo del paño y sus desfases."
                )
        ancho = coord_fin - coord_inicio
        if not math.isfinite(span) or not math.isfinite(ancho) or span <= 0 or ancho <= 0:
            raise ValueError("La luz y el ancho del paño deben ser positivos.")
        cota = float(self.nivel_apoyo.currentData())
        nivel_nombre = self.nivel_apoyo.currentText().split(" Â· ", 1)[0]
        geometria = {
            "tipologia": str(self.tipologia.currentData()),
            "luz_transversal_m": span,
            "ancho_losa_m": ancho,
            "apoya_en": apoyos,
            "panel_ejes": {
                "direccion_luz": direccion,
                "nivel_apoyo_nombre": nivel_nombre,
                "cota_apoyo_m": cota,
                "ejes": {
                    "x": {
                        "inicio": {
                            "id": str(x_ejes[0].get("id", "")),
                            "nombre": str(x_ejes[0]["nombre"]),
                            "coordenada_m": float(x_ejes[0]["coordenada_m"]),
                        },
                        "fin": {
                            "id": str(x_ejes[1].get("id", "")),
                            "nombre": str(x_ejes[1]["nombre"]),
                            "coordenada_m": float(x_ejes[1]["coordenada_m"]),
                        },
                    },
                    "y": {
                        "inicio": {
                            "id": str(y_ejes[0].get("id", "")),
                            "nombre": str(y_ejes[0]["nombre"]),
                            "coordenada_m": float(y_ejes[0]["coordenada_m"]),
                        },
                        "fin": {
                            "id": str(y_ejes[1].get("id", "")),
                            "nombre": str(y_ejes[1]["nombre"]),
                            "coordenada_m": float(y_ejes[1]["coordenada_m"]),
                        },
                    },
                },
                "desfase_inicio_m": float(self.desfase_inicio.value()),
                "desfase_fin_m": float(self.desfase_fin.value()),
                "coordenada_libre_inicio_m": coord_inicio,
                "coordenada_libre_fin_m": coord_fin,
            },
        }
        return nombre, geometria

    def _aceptar(self) -> None:
        try:
            self.resultado()
        except (TypeError, ValueError) as exc:
            QMessageBox.warning(self, "Geometría de losa incompleta", str(exc))
            return
        self.accept()

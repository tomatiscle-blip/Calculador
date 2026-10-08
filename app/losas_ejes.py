"""Alta geométrica de paños de losa vinculados a pórticos y ejes."""

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
        self.ejes_seleccionados: dict[str, list[str]] = {"x": [], "y": []}
        self.direccion_luz = "x"
        self.coordenadas_paño: tuple[float, float, float, float] | None = None
        self.apoyos: dict[str, str] = {}
        self.setMinimumSize(560, 320)

    def establecer_datos(
        self,
        ejes: dict[str, list[dict]],
        estructura: dict,
        ejes_seleccionados: dict[str, list[str]],
        direccion_luz: str,
        coordenadas_paño: tuple[float, float, float, float] | None,
        apoyos: dict[str, str],
    ) -> None:
        self.ejes = ejes
        self.ejes_seleccionados = {
            familia: list(ids) for familia, ids in ejes_seleccionados.items()
        }
        self.direccion_luz = direccion_luz
        self.coordenadas_paño = coordenadas_paño
        self.apoyos = apoyos
        self.porticos = {}
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
            for viga_id, viga in (datos.get("vigas", {}) or {}).items():
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
        return coords_x, coords_y

    def _limites_pantalla(self) -> tuple[float, float, float, float]:
        return 66.0, max(67.0, self.width() - 24.0), 38.0, max(39.0, self.height() - 45.0)

    def paintEvent(self, event) -> None:
        super().paintEvent(event)
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.fillRect(self.rect(), QColor("#ffffff"))
        painter.setPen(QColor("#334155"))
        painter.drawText(12, 20, "Planta · clic en dos ejes X y dos ejes Y para delimitar el paño")

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
                    painter.drawText(centro_x + 8, centro_y - 8, f"{nombre} · {viga_id}")

        painter.setPen(QColor("#64748b"))
        painter.drawText(
            12, self.height() - 12,
            "Verde: vigas de apoyo seleccionadas  ·  flecha roja: dirección de trabajo",
        )


class EditorLosaDesdeEjes(QDialog):
    """Delimita un paño, vincula dos vigas de pórtico y devuelve su geometría."""

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
        self.direccion_luz.addItem("Luz en X · apoyos en X mínimo/máximo", "x")
        self.direccion_luz.addItem("Luz en Y · apoyos en Y mínimo/máximo", "y")
        direccion_previa = str(
            (self.elemento.get("panel_ejes") or {}).get("direccion_luz", "x")
        )
        self.direccion_luz.setCurrentIndex(
            max(self.direccion_luz.findData(direccion_previa), 0)
        )
        self.nivel_apoyo = QComboBox()
        for cota, nivel in self._niveles_de_vigas():
            self.nivel_apoyo.addItem(f"{nivel} · {cota:g} m", cota)
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
        self.formulario.addRow("Viga de apoyo en el primer borde:", self.selector_vigas["izq"])
        self.formulario.addRow("Viga de apoyo en el segundo borde:", self.selector_vigas["der"])
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
            "de trabajo: los apoyos deben ser dos pórticos paralelos al borde y ubicados "
            "exactamente sobre sus ejes. Los bordes libres pueden desplazarse solo dentro "
            "del tramo de viga disponible en ambos apoyos."
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
        return (
            origen + min(
                min(float(t["x_inicio"]), float(t["x_fin"])) for t in tramos
            ),
            origen + max(
                max(float(t["x_inicio"]), float(t["x_fin"])) for t in tramos
            ),
        )

    def _viga_cubre(
        self, portico: str, viga_id: str, inicio: float, fin: float
    ) -> bool:
        datos = self.estructura.get(portico, {})
        referencia = datos.get("referencia_planta", {}) or {}
        direccion = str(referencia.get("direccion", "x"))
        eje = self._eje_por_id(
            direccion, str(referencia.get("eje_longitudinal_id", ""))
        )
        origen = float(
            referencia.get(
                "origen_longitudinal_m",
                eje["coordenada_m"] if eje else 0.0,
            )
        )
        viga = (datos.get("vigas", {}) or {}).get(viga_id, {})
        intervalos = sorted(
            (
                origen + min(float(tramo["x_inicio"]), float(tramo["x_fin"])),
                origen + max(float(tramo["x_inicio"]), float(tramo["x_fin"])),
            )
            for tramo in viga.get("tramos", [])
        )
        cubierto_hasta = inicio
        for tramo_inicio, tramo_fin in intervalos:
            if tramo_inicio > cubierto_hasta + 1e-6:
                break
            cubierto_hasta = max(cubierto_hasta, tramo_fin)
        return cubierto_hasta >= fin - 1e-6

    def _poblar_vigas(self, lado: str, eje: dict, direccion: str, cota: float) -> str | None:
        selector = self.selector_vigas[lado]
        previo = str(selector.currentData() or "")
        selector.blockSignals(True)
        selector.clear()
        portico = self._portico_en_eje(
            "y" if direccion == "x" else "x",
            "x" if direccion == "x" else "y",
            eje,
        )
        if portico:
            datos = self.estructura[portico]
            for viga_id, viga in (datos.get("vigas", {}) or {}).items():
                try:
                    cota_viga = float(viga.get("cota_m"))
                except (TypeError, ValueError):
                    continue
                if abs(cota_viga - cota) > 1e-4:
                    continue
                extremos = self._viga_extremos(portico, str(viga_id))
                if extremos is None:
                    continue
                selector.addItem(
                    f"{viga_id} · {extremos[0]:g}–{extremos[1]:g} m",
                    str(viga_id),
                )
            indice = selector.findData(previo)
            selector.setCurrentIndex(indice if indice >= 0 else (0 if selector.count() else -1))
        selector.blockSignals(False)
        return portico

    def _actualizar(self) -> None:
        direccion = str(self.direccion_luz.currentData() or "x")
        coordenadas = self._coordenadas_ejes()
        nivel_data = self.nivel_apoyo.currentData()
        apoyos: dict[str, str] = {}
        porticos_borde: dict[str, str | None] = {}
        ejes_soporte = []
        if coordenadas and nivel_data is not None:
            x_ejes, y_ejes = coordenadas
            ejes_soporte = x_ejes if direccion == "x" else y_ejes
            for lado, eje in zip(("izq", "der"), ejes_soporte):
                portico = self._poblar_vigas(
                    lado, eje, direccion, float(nivel_data)
                )
                porticos_borde[lado] = portico
                viga_id = str(self.selector_vigas[lado].currentData() or "")
                if portico and viga_id:
                    apoyos[portico] = viga_id

        intervalos = []
        for lado in ("izq", "der"):
            portico = porticos_borde.get(lado)
            viga_id = str(self.selector_vigas[lado].currentData() or "")
            intervalos.append(
                self._viga_extremos(portico, viga_id)
                if portico and viga_id else None
            )
        intervalo_comun = None
        if all(intervalos):
            inicio = max(intervalo[0] for intervalo in intervalos if intervalo)
            fin = min(intervalo[1] for intervalo in intervalos if intervalo)
            if fin > inicio:
                intervalo_comun = (inicio, fin)

        self.desfase_inicio.blockSignals(True)
        self.desfase_fin.blockSignals(True)
        if intervalo_comun and coordenadas:
            x_ejes, y_ejes = coordenadas
            ejes_libres = y_ejes if direccion == "x" else x_ejes
            coord_inicio = float(ejes_libres[0]["coordenada_m"])
            coord_fin = float(ejes_libres[1]["coordenada_m"])
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
                or not all(
                    self._viga_cubre(
                        str(porticos_borde.get(lado) or ""),
                        str(self.selector_vigas[lado].currentData() or ""),
                        inicio_paño,
                        fin_paño,
                    )
                    for lado in ("izq", "der")
                )
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
                f"Paño X: {x_ejes[0]['nombre']}–{x_ejes[1]['nombre']} · "
                f"Y: {y_ejes[0]['nombre']}–{y_ejes[1]['nombre']}"
            )

        self.etiqueta_apoyos.setText(
            " · ".join(
                f"{lado}: {porticos_borde.get(lado) or 'sin pórtico compatible'}"
                + (
                    f" / {self.selector_vigas[lado].currentData()}"
                    if self.selector_vigas[lado].currentData() else ""
                )
                for lado in ("izq", "der")
            )
        )
        if not coordenadas:
            estado = nombres_ejes
        elif not all(porticos_borde.get(lado) for lado in ("izq", "der")):
            estado = "Falta un pórtico sobre uno de los ejes de apoyo. Asociá los pórticos a esos ejes."
        elif not all(self.selector_vigas[lado].currentData() for lado in ("izq", "der")):
            estado = "Los pórticos seleccionados no tienen una viga en el nivel elegido."
        elif intervalo_comun is None:
            estado = "Las vigas de apoyo no cubren el ancho completo del paño."
        else:
            estado = (
                f"{nombres_ejes} · luz {direccion.upper()} "
                f"{abs(float(ejes_soporte[1]['coordenada_m']) - float(ejes_soporte[0]['coordenada_m'])):.3f} m"
            )
        self.estado.setText(estado)
        self.vista.establecer_datos(
            self.ejes, self.estructura, self.vista.ejes_seleccionados,
            direccion, coordenadas_paño, apoyos,
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
        apoyos = {}
        for lado, eje in zip(("izq", "der"), ejes_soporte):
            portico = self._portico_en_eje(
                "y" if direccion == "x" else "x",
                "x" if direccion == "x" else "y",
                eje,
            )
            viga_id = str(self.selector_vigas[lado].currentData() or "")
            if not portico or not viga_id:
                raise ValueError(
                    f"Falta un pórtico y una viga de apoyo en el eje {eje['nombre']}."
                )
            apoyos[lado] = {"portico": portico, "viga": viga_id}
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
        nivel_nombre = self.nivel_apoyo.currentText().split(" · ", 1)[0]
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

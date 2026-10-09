"""Esquema de geometría y cargas asignadas para la pantalla Inicio."""

from __future__ import annotations

from PySide6.QtCore import QPointF, QRectF, Qt, Signal
from PySide6.QtGui import QColor, QFont, QFontMetricsF, QPainter, QPen, QPolygonF
from PySide6.QtWidgets import QWidget

from calc import cargas, materiales, rutas


class VistaDiagramaMomentos(QWidget):
    """Dibuja la geometría del pórtico y sus diagramas de momento y corte."""

    solicitacion_seleccionada = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumSize(620, 360)
        self.diagrama: dict | None = None
        self.tipo_diagrama = "momento"
        self.setMouseTracking(True)
        self.setToolTip("Hacé clic sobre una curva para consultar el valor en esa estación.")
        self._puntos_seleccionables: list[tuple[QPointF, str, float, float]] = []
        self._seleccion: tuple[str, float, float] | None = None

    def establecer_diagrama(self, diagrama: dict | None) -> None:
        self.diagrama = diagrama
        self.update()

    def establecer_tipo_diagrama(self, tipo: str) -> None:
        if tipo not in ("momento", "corte"):
            raise ValueError(f"Tipo de diagrama no válido: {tipo}")
        self.tipo_diagrama = tipo
        self._seleccion = None
        self.update()

    def paintEvent(self, event) -> None:
        super().paintEvent(event)
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.fillRect(self.rect(), QColor("#f8fafc"))
        painter.setPen(QPen(QColor("#dbe4f0"), 1))
        painter.setBrush(Qt.BrushStyle.NoBrush)
        painter.drawRoundedRect(self.rect().adjusted(1, 1, -2, -2), 10, 10)
        self._puntos_seleccionables = []
        barras = (self.diagrama or {}).get("barras", [])
        if not barras:
            painter.setPen(QColor("#64748b"))
            painter.drawText(
                self.rect(),
                Qt.AlignmentFlag.AlignCenter,
                "No hay datos de momentos para representar.",
            )
            return

        puntos_extremos = [
            extremo
            for barra in barras
            for extremo in (barra["inicio"], barra["fin"])
        ]
        xs = [float(punto["x"]) for punto in puntos_extremos]
        ys = [float(punto["y"]) for punto in puntos_extremos]
        xmin, xmax = min(xs), max(xs)
        ymin, ymax = min(ys), max(ys)
        rango_x, rango_y = max(xmax - xmin, 1.0), max(ymax - ymin, 1.0)
        margen_x, margen_y = 64.0, 42.0
        escala = min(
            (self.width() - 2 * margen_x) / rango_x,
            (self.height() - 2 * margen_y) / rango_y,
        )
        ancho_dibujo, alto_dibujo = rango_x * escala, rango_y * escala
        origen_x = (self.width() - ancho_dibujo) / 2
        origen_y = self.height() - (self.height() - alto_dibujo) / 2

        def punto(x: float, y: float) -> QPointF:
            return QPointF(
                origen_x + (x - xmin) * escala,
                origen_y - (y - ymin) * escala,
            )

        clave_valores = "M_kNm" if self.tipo_diagrama == "momento" else "V_kN"
        clave_maximo = (
            "max_abs_kNm" if self.tipo_diagrama == "momento" else "max_abs_kN"
        )
        maximo = float((self.diagrama or {}).get(clave_maximo, 0.0))
        if maximo <= 0:
            maximo = max((
                abs(float(valor))
                for barra in barras
                for valor in barra.get(clave_valores, [])
            ), default=0.0)
        amplitud_px = min(self.width() - 2 * margen_x, self.height() - 2 * margen_y) * 0.16
        etiquetas = []
        painter.setPen(QPen(QColor("#475569"), 2))
        painter.setBrush(Qt.BrushStyle.NoBrush)

        for barra in barras:
            inicio, fin = barra["inicio"], barra["fin"]
            a = punto(float(inicio["x"]), float(inicio["y"]))
            b = punto(float(fin["x"]), float(fin["y"]))
            estaciones = barra.get("x_m", [])
            valores = barra.get(clave_valores, [])
            longitud = float(estaciones[-1]) if estaciones else 0.0
            if not estaciones or len(estaciones) != len(valores) or longitud <= 0:
                continue

            dx, dy = b.x() - a.x(), b.y() - a.y()
            largo_pantalla = (dx * dx + dy * dy) ** 0.5
            if largo_pantalla <= 0:
                continue
            normal_x, normal_y = -dy / largo_pantalla, dx / largo_pantalla
            curva = []
            for x_local, valor in zip(estaciones, valores):
                fraccion = min(1.0, max(0.0, float(x_local) / longitud))
                base_x = a.x() + dx * fraccion
                base_y = a.y() + dy * fraccion
                desfase = (
                    -float(valor) / maximo * amplitud_px if maximo > 0 else 0.0
                )
                punto_curva = QPointF(
                    base_x + normal_x * desfase,
                    base_y + normal_y * desfase,
                )
                curva.append(punto_curva)
                self._puntos_seleccionables.append((
                    punto_curva,
                    str(barra["id"]),
                    float(x_local),
                    -float(valor),
                ))
            area = QPolygonF([a, *curva, b])
            painter.setPen(Qt.PenStyle.NoPen)
            color_area = QColor("#60a5fa")
            color_area.setAlpha(72)
            painter.setBrush(color_area)
            painter.drawPolygon(area)
            painter.setPen(QPen(QColor("#1d4ed8"), 2))
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawPolyline(QPolygonF(curva))
            painter.setPen(QPen(QColor("#475569"), 2))
            painter.drawLine(a, b)

            punto_medio = QPointF((a.x() + b.x()) / 2, (a.y() + b.y()) / 2)
            etiquetas.append((punto_medio, str(barra["id"]), QColor("#1e3a8a")))

        if self.tipo_diagrama == "corte" and not any(
            barra.get("V_kN") for barra in barras
        ):
            painter.setPen(QColor("#64748b"))
            painter.drawText(
                self.rect(),
                Qt.AlignmentFlag.AlignCenter,
                "No hay datos de corte. Recalculá las solicitaciones para generarlos.",
            )
            return

        for apoyo in (self.diagrama or {}).get("apoyos", []):
            base = punto(float(apoyo["x"]), float(apoyo.get("y", 0.0)))
            painter.setPen(QPen(QColor("#334155"), 1.5))
            painter.setBrush(QColor("#f8fafc"))
            if apoyo.get("tipo") == "articulado":
                painter.drawPolygon(QPolygonF([
                    QPointF(base.x(), base.y() + 1),
                    QPointF(base.x() - 8, base.y() + 12),
                    QPointF(base.x() + 8, base.y() + 12),
                ]))
                painter.drawLine(
                    QPointF(base.x() - 11, base.y() + 13),
                    QPointF(base.x() + 11, base.y() + 13),
                )
            else:
                painter.drawLine(
                    QPointF(base.x() - 10, base.y() + 5),
                    QPointF(base.x() + 10, base.y() + 5),
                )
                for desplazamiento in (-7, 0, 7):
                    painter.drawLine(
                        QPointF(base.x() + desplazamiento, base.y() + 5),
                        QPointF(base.x() + desplazamiento - 5, base.y() + 11),
                    )

        painter.setFont(QFont(self.font().family(), 9, QFont.Weight.DemiBold))
        metricas = QFontMetricsF(painter.font())
        area_segura = QRectF(self.rect()).adjusted(10, 10, -10, -10)
        ocupadas: list[QRectF] = []
        for ancla, texto, color in etiquetas:
            ancho = metricas.horizontalAdvance(texto) + 18
            alto = metricas.height() + 8
            posiciones = []
            for radio in (10, 34, 58, 82, 106):
                posiciones.extend((
                    (ancla.x() + radio, ancla.y() - alto / 2),
                    (ancla.x() - ancho - radio, ancla.y() - alto / 2),
                    (ancla.x() - ancho / 2, ancla.y() - alto - radio),
                    (ancla.x() - ancho / 2, ancla.y() + radio),
                    (ancla.x() + radio, ancla.y() - alto - radio),
                    (ancla.x() + radio, ancla.y() + radio),
                    (ancla.x() - ancho - radio, ancla.y() - alto - radio),
                    (ancla.x() - ancho - radio, ancla.y() + radio),
                ))
            rectangulo = next((
                QRectF(x, y, ancho, alto)
                for x, y in posiciones
                if area_segura.contains(QRectF(x, y, ancho, alto))
                and not any(
                    QRectF(x, y, ancho, alto).adjusted(-3, -2, 3, 2).intersects(otra)
                    for otra in ocupadas
                )
            ), None)
            if rectangulo is None:
                candidatos = [
                    QRectF(x, y, ancho, alto)
                    for y in range(
                        int(area_segura.top()),
                        int(area_segura.bottom() - alto),
                        int(alto + 5),
                    )
                    for x in range(
                        int(area_segura.left()),
                        int(area_segura.right() - ancho),
                        int(ancho + 5),
                    )
                ]
                rectangulo = next((
                    rect for rect in sorted(
                        candidatos,
                        key=lambda r: (r.center().x() - ancla.x()) ** 2
                        + (r.center().y() - ancla.y()) ** 2,
                    )
                    if not any(
                        rect.adjusted(-3, -2, 3, 2).intersects(otra)
                        for otra in ocupadas
                    )
                ), None)
            if rectangulo is None:
                continue

            objetivo = QPointF(
                min(max(ancla.x(), rectangulo.left()), rectangulo.right()),
                min(max(ancla.y(), rectangulo.top()), rectangulo.bottom()),
            )
            painter.setPen(QPen(QColor(color.red(), color.green(), color.blue(), 170), 1))
            painter.drawLine(ancla, objetivo)
            fondo = QColor("#ffffff")
            fondo.setAlpha(242)
            painter.setBrush(fondo)
            painter.setPen(QPen(color, 1.3))
            painter.drawRoundedRect(rectangulo, 7, 7)
            painter.setPen(QColor("#172554"))
            painter.drawText(rectangulo, Qt.AlignmentFlag.AlignCenter, texto)
            ocupadas.append(rectangulo.adjusted(-3, -2, 3, 2))

        if self._seleccion:
            barra_seleccionada, estacion_seleccionada, _ = self._seleccion
            candidatos = [
                (punto, estacion, momento)
                for punto, barra_id, estacion, momento in self._puntos_seleccionables
                if barra_id == barra_seleccionada
            ]
            if candidatos:
                punto, _, _ = min(
                    candidatos,
                    key=lambda candidato: abs(candidato[1] - estacion_seleccionada),
                )
                painter.setPen(QPen(QColor("#ffffff"), 2))
                painter.setBrush(QColor("#dc2626"))
                painter.drawEllipse(punto, 6, 6)

    def mousePressEvent(self, event) -> None:
        if event.button() == Qt.MouseButton.LeftButton and self._puntos_seleccionables:
            posicion = event.position()
            punto, barra_id, estacion, momento = min(
                self._puntos_seleccionables,
                key=lambda muestra: (
                    (muestra[0].x() - posicion.x()) ** 2
                    + (muestra[0].y() - posicion.y()) ** 2
                ),
            )
            distancia = (
                (punto.x() - posicion.x()) ** 2
                + (punto.y() - posicion.y()) ** 2
            ) ** 0.5
            if distancia <= 24:
                self._seleccion = (barra_id, estacion, momento)
                texto = (
                    f"{barra_id} · x = {estacion:.2f} m · "
                    f"{'M' if self.tipo_diagrama == 'momento' else 'V'} = "
                    f"{momento:+.2f} "
                    f"{'kN·m' if self.tipo_diagrama == 'momento' else 'kN'}"
                )
                self.solicitacion_seleccionada.emit(texto)
                self.update()
        super().mousePressEvent(event)

class VistaPortico(QWidget):
    """Dibuja el pórtico seleccionado y los intervalos de carga del proyecto."""

    PALETA = ("#2563eb", "#15803d", "#d97706", "#7c3aed", "#0891b2", "#c2410c")

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumSize(380, 250)
        self.setToolTip("Esquema de tramos y cargas asignadas; no está a escala de dibujo.")
        self.nombre = ""
        self.columnas: dict = {}
        self.vigas: dict = {}
        self.tramos: list[dict] = []
        self.cargas: dict[str, list[dict]] = {}
        self.aviso = ""
        self.leyenda = "Sin cargas asignadas a este pórtico."
        self.ancho_viento_m = 0.0
        self.viento_lineal_kN_m = 0.0
        self.viento_nodos: list[dict] = []
        self.capas_carga = 0

    def actualizar(self, nombre: str) -> None:
        self.nombre = nombre
        estructura = rutas.cargar_estructura().get(nombre, {}) if nombre else {}
        self.columnas = estructura.get("columnas", {})
        self.vigas = estructura.get("vigas", {})
        self.tramos = self._geometria_tramos()
        self.cargas = {}
        self.capas_carga = 0
        self.aviso = ""
        self.ancho_viento_m = 0.0
        self.viento_lineal_kN_m = 0.0
        self.viento_nodos = []
        if nombre:
            try:
                from calc.portico import _cargas_asignadas

                self.cargas, avisos = _cargas_asignadas(nombre)
                self.aviso = "; ".join(avisos)
            except Exception as exc:
                self.aviso = f"No se pudieron leer las cargas: {exc}"
        self.capas_carga = max((
            len({c.get("aplicacion_id", "") for c in cargas_tramo})
            for cargas_tramo in self.cargas.values()
        ), default=0)
        self.setMinimumHeight(max(250, 190 + self.capas_carga * 24))
        proyecto = rutas.leer_json(rutas.CARGAS, {}) or {}
        configuracion_viento = proyecto.get("viento", {}) or {}
        if configuracion_viento.get("activo"):
            try:
                items_viento = cargas.items_viento_general(proyecto, materiales.cargar())
                self.viento_lineal_kN_m = sum(float(item["valor"]) for item in items_viento)
                self.ancho_viento_m = float(configuracion_viento.get("ancho_tributario_m", 0.0))
                columnas_por_piso = {}
                for cid, columna in self.columnas.items():
                    prefijo = str(cid).split("-")[0]
                    if prefijo.startswith("C") and prefijo[1:].isdigit():
                        columnas_por_piso.setdefault(prefijo, []).append(columna)
                for columnas in columnas_por_piso.values():
                    if not columnas:
                        continue
                    altura = float(columnas[0].get("altura_m", 0.0))
                    fuerza_nudo = self.viento_lineal_kN_m * altura / len(columnas)
                    for columna in columnas:
                        self.viento_nodos.append({
                            "x": float(columna["x"]),
                            "y": float(columna.get("nivel", 0.0)) + float(columna.get("altura_m", 0.0)),
                            "valor_kN": fuerza_nudo,
                        })
            except Exception as exc:
                self.aviso = f"No se pudo representar el viento: {exc}"
        ids = sorted({c.get("carga_id", "") for grupo in self.cargas.values() for c in grupo})
        self.colores = {cid: QColor(self.PALETA[i % len(self.PALETA)]) for i, cid in enumerate(ids)}
        lineas = []
        for cid in ids:
            carga = next(c for grupo in self.cargas.values() for c in grupo if c.get("carga_id", "") == cid)
            tipos = ", ".join(dict.fromkeys(c["tipo"] for grupo in self.cargas.values() for c in grupo if c.get("carga_id") == cid))
            formas = ", ".join(dict.fromkeys(
                "puntual" if c.get("tipo_aplicacion") == "puntual" else "lineal"
                for grupo in self.cargas.values() for c in grupo if c.get("carga_id") == cid
            ))
            lineas.append(
                f'<span style="color:{self.colores[cid].name()};font-weight:bold">━━</span> '
                f'{carga.get("descripcion", cid)} · {formas} · {tipos}'
            )
        if self.viento_lineal_kN_m:
            altura_total = max(
                (float(c.get("nivel", 0.0)) + float(c.get("altura_m", 0.0)) for c in self.columnas.values()),
                default=0.0,
            ) - min((float(c.get("nivel", 0.0)) for c in self.columnas.values()), default=0.0)
            fuerza_total = sum(nudo["valor_kN"] for nudo in self.viento_nodos)
            nudos = len(self.viento_nodos)
            lineas.append(
                f'<span style="color:#d97706;font-weight:bold">━━ W →</span> '
                f'Viento: {self.viento_lineal_kN_m:.2f} kN/m × H={altura_total:.2f} m '
                f'= {fuerza_total:.2f} kN en el pórtico; '
                f'b tributario={self.ancho_viento_m:.2f} m; '
                f'el motor aplica fuerzas W en los {nudos} nudos superiores.'
            )
        self.leyenda = "<br>".join(lineas) if lineas else "Sin cargas asignadas a este pórtico."
        self.update()

    def _geometria_tramos(self) -> list[dict]:
        tramos = []
        for viga_id, viga in self.vigas.items():
            cabeza = str(viga_id).split("-")[0]
            piso = cabeza[1:] if cabeza[:1].upper() == "V" else ""
            if not piso.isdigit():
                continue
            columnas = sorted(
                (c for cid, c in self.columnas.items() if str(cid).startswith(f"C{piso}-")),
                key=lambda c: float(c.get("x", 0.0)),
            )
            if not columnas:
                continue
            y = float(columnas[0].get("nivel", 0.0)) + float(columnas[0].get("altura_m", 0.0))
            for indice, tramo in enumerate(viga.get("tramos", [])):
                longitud = float(tramo.get("longitud_m", 0.0))
                tramo_id = str(tramo.get("id", f"{viga_id} tramo {indice + 1}"))
                es_voladizo = bool(tramo.get("es_voladizo"))
                if es_voladizo:
                    if "izq" in tramo_id.lower():
                        x0, x1 = float(columnas[0]["x"]) - longitud, float(columnas[0]["x"])
                        invertir_carga = False
                    else:
                        x0, x1 = float(columnas[-1]["x"]), float(columnas[-1]["x"]) + longitud
                        invertir_carga = True
                    apoyado = True
                elif indice + 1 < len(columnas):
                    x0, x1 = float(columnas[indice]["x"]), float(columnas[indice + 1]["x"])
                    invertir_carga = False
                    apoyado = True
                else:
                    x0 = float(columnas[-1]["x"])
                    x1 = x0 + longitud
                    invertir_carga = False
                    apoyado = False
                tramos.append({
                    "id": tramo_id, "x0": x0, "x1": x1, "y": y,
                    "apoyado": apoyado, "invertir_carga": invertir_carga,
                })
        return tramos

    def paintEvent(self, _event) -> None:
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.fillRect(self.rect(), QColor("#ffffff"))
        if not self.columnas:
            painter.setPen(QColor("#64748b"))
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "Seleccioná un pórtico con geometría")
            return

        xs = [float(c.get("x", 0.0)) for c in self.columnas.values()]
        ys = [float(c.get("nivel", 0.0)) for c in self.columnas.values()]
        ys.extend(float(c.get("nivel", 0.0)) + float(c.get("altura_m", 0.0)) for c in self.columnas.values())
        for tramo in self.tramos:
            xs.extend((tramo["x0"], tramo["x1"]))
            ys.append(tramo["y"])
        if not xs or not ys:
            return
        ancho, alto = self.width(), self.height()
        margen_x = 56.0
        margen_y_superior = max(44.0, 30.0 + self.capas_carga * 24.0)
        margen_y_inferior = 38.0
        xmin, xmax = min(xs), max(xs)
        ymin, ymax = min(ys), max(ys)
        rango_x, rango_y = max(xmax - xmin, 1.0), max(ymax - ymin, 1.0)
        alto_util = alto - margen_y_superior - margen_y_inferior
        escala = min((ancho - 2 * margen_x) / rango_x, alto_util / rango_y)
        dibujar_ancho = rango_x * escala
        dibujar_alto = rango_y * escala
        origen_x = (ancho - dibujar_ancho) / 2
        origen_y = alto - margen_y_inferior - (alto_util - dibujar_alto) / 2

        def punto(x, y):
            return QPointF(origen_x + (x - xmin) * escala, origen_y - (y - ymin) * escala)

        painter.setPen(QPen(QColor("#334155"), 3))
        for columna in self.columnas.values():
            x = float(columna.get("x", 0.0))
            y0 = float(columna.get("nivel", 0.0))
            y1 = y0 + float(columna.get("altura_m", 0.0))
            painter.drawLine(punto(x, y0), punto(x, y1))
            painter.setBrush(QColor("#334155"))
            painter.drawEllipse(punto(x, y0), 3.5, 3.5)

        font = QFont()
        font.setPointSize(8)
        painter.setFont(font)
        for tramo in self.tramos:
            a, b = punto(tramo["x0"], tramo["y"]), punto(tramo["x1"], tramo["y"])
            if tramo["apoyado"]:
                painter.setPen(QPen(QColor("#334155"), 4))
            else:
                painter.setPen(QPen(QColor("#dc2626"), 3, Qt.PenStyle.DashLine))
            painter.drawLine(a, b)
            color_texto = QColor("#334155") if tramo["apoyado"] else QColor("#b91c1c")
            painter.setPen(color_texto)
            texto_tramo = f"{tramo['id']} · {abs(tramo['x1'] - tramo['x0']):.2f} m"
            ancho_texto_tramo = painter.fontMetrics().horizontalAdvance(texto_tramo)
            x_texto_tramo = (a.x() + b.x() - ancho_texto_tramo) / 2
            painter.drawText(QPointF(x_texto_tramo, a.y() + 17), texto_tramo)

        grupos = {}
        por_id = {t["id"]: t for t in self.tramos}
        for tramo_id, cargas_tramo in self.cargas.items():
            geometria = por_id.get(tramo_id)
            if not geometria:
                continue
            for carga in cargas_tramo:
                clave = (tramo_id, carga.get("aplicacion_id", ""))
                grupos.setdefault(clave, []).append(carga)

        niveles_tramo = {}
        for (tramo_id, aplicacion_id), cargas_aplicacion in grupos.items():
            geometria = por_id[tramo_id]
            carga_ref = cargas_aplicacion[0]
            color_linea = self.colores.get(carga_ref.get("carga_id", ""), QColor("#7c3aed"))
            if carga_ref.get("tipo_aplicacion") == "puntual":
                x_local = float(carga_ref["x_m"])
                x_aplicado = (
                    geometria["x1"] - x_local
                    if geometria["invertir_carga"]
                    else geometria["x0"] + x_local
                )
                x_aplicado = max(geometria["x0"], min(geometria["x1"], x_aplicado))
                posicion = punto(x_aplicado, geometria["y"])
                nivel = niveles_tramo.get(tramo_id, 0)
                niveles_tramo[tramo_id] = nivel + 1
                y_inicio = posicion.y() - 16 - nivel * 24
                valor = float(carga_ref["valor_kN"])
                sentido = float(carga_ref.get("signo", -1))
                y_fin = y_inicio + (10 if sentido < 0 else -10)
                painter.setPen(QPen(color_linea, 1.7))
                painter.drawLine(QPointF(posicion.x(), y_inicio), QPointF(posicion.x(), y_fin))
                delta = 4 if sentido < 0 else -4
                painter.setBrush(color_linea)
                painter.drawPolygon(QPolygonF([
                    QPointF(posicion.x(), y_fin),
                    QPointF(posicion.x() - 3, y_fin - delta),
                    QPointF(posicion.x() + 3, y_fin - delta),
                ]))
                influencia = carga_ref.get("influencia_m")
                detalle_influencia = (
                    f" · L/2={float(influencia):.2f} m"
                    if influencia is not None else ""
                )
                if abs(geometria["x1"] - geometria["x0"]) < 2.0:
                    etiqueta = f"{carga_ref['tipo']} {valor:.1f} kN"
                else:
                    etiqueta = f"{carga_ref['tipo']} {valor:.2f} kN{detalle_influencia}"
                painter.setPen(QColor("#0f172a"))
                texto_x = min(posicion.x() + 5, ancho - painter.fontMetrics().horizontalAdvance(etiqueta) - 4)
                painter.drawText(QPointF(max(4.0, texto_x), y_inicio - 2), etiqueta)
                continue

            x0 = float(carga_ref["x_inicio_m"])
            x1 = float(carga_ref["x_fin_m"])
            if geometria["invertir_carga"]:
                inicio, fin = geometria["x1"] - x1, geometria["x1"] - x0
            else:
                inicio, fin = geometria["x0"] + x0, geometria["x0"] + x1
            inicio = max(geometria["x0"], min(geometria["x1"], inicio))
            fin = max(geometria["x0"], min(geometria["x1"], fin))
            if fin <= inicio:
                continue
            nivel = niveles_tramo.get(tramo_id, 0)
            niveles_tramo[tramo_id] = nivel + 1
            base = punto(inicio, geometria["y"])
            final = punto(fin, geometria["y"])
            y_linea = base.y() - 14 - nivel * 24
            tipos = {}
            for carga in cargas_aplicacion:
                tipo = carga["tipo"]
                tipos[tipo] = tipos.get(tipo, 0.0) + float(carga["valor_kN_m"])
            painter.setPen(QPen(color_linea, 1.5))
            painter.drawLine(QPointF(base.x(), y_linea), QPointF(final.x(), y_linea))

            categorias = list(tipos)
            for indice, tipo in enumerate(categorias):
                fraccion = (indice + 1) / (len(categorias) + 1)
                x_arrow = base.x() + (final.x() - base.x()) * fraccion
                sentido = next(
                    (c.get("signo", -1) for c in cargas_aplicacion if c["tipo"] == tipo), -1
                )
                y_fin = y_linea + (8 if sentido < 0 else -8)
                color_flecha = color_linea
                painter.setPen(QPen(color_flecha, 1.5))
                painter.drawLine(QPointF(x_arrow, y_linea), QPointF(x_arrow, y_fin))
                delta = 4 if sentido < 0 else -4
                painter.setBrush(color_flecha)
                painter.drawPolygon(QPolygonF([
                    QPointF(x_arrow, y_fin),
                    QPointF(x_arrow - 3, y_fin - delta),
                    QPointF(x_arrow + 3, y_fin - delta),
                ]))

            painter.setPen(QColor("#0f172a") if geometria["apoyado"] else QColor("#b91c1c"))
            longitud_tramo = abs(geometria["x1"] - geometria["x0"])
            decimales = 1 if longitud_tramo < 2.0 else 2
            valores = " / ".join(
                f"{t} {v:.{decimales}f}" for t, v in tipos.items()
            )
            sin_apoyo = (
                " · voladizo"
                if not geometria["apoyado"] and longitud_tramo >= 2.0
                else ""
            )
            etiqueta = f"{valores} kN/m{sin_apoyo}"
            ancho_texto = painter.fontMetrics().horizontalAdvance(etiqueta)
            x_texto = max(4.0, min((base.x() + final.x() - ancho_texto) / 2, ancho - ancho_texto - 4))
            painter.drawText(QPointF(x_texto, y_linea - 4), etiqueta)

        # El viento se indica una vez fuera del marco; el motor conserva su
        # distribución interna en los nudos superiores.
        if self.viento_nodos:
            x_marco = min(float(c.get("x", 0.0)) for c in self.columnas.values())
            y_base = min(float(c.get("nivel", 0.0)) for c in self.columnas.values())
            y_cima = max(float(c.get("nivel", 0.0)) + float(c.get("altura_m", 0.0)) for c in self.columnas.values())
            centro = punto(x_marco, (y_base + y_cima) / 2)
            y_flecha = centro.y()
            x1 = centro.x() - 8
            x0 = x1 - 32
            painter.setPen(QPen(QColor("#d97706"), 2.5))
            painter.drawLine(QPointF(x0, y_flecha), QPointF(x1, y_flecha))
            painter.setBrush(QColor("#d97706"))
            painter.drawPolygon(QPolygonF([
                QPointF(x1, y_flecha), QPointF(x1 - 6, y_flecha - 4), QPointF(x1 - 6, y_flecha + 4)
            ]))
            painter.setPen(QColor("#9a3412"))
            total = sum(nudo["valor_kN"] for nudo in self.viento_nodos)
            painter.drawText(QPointF(x0 - 8, y_flecha - 7), f"W = {total:.2f} kN")

        if self.aviso:
            painter.setPen(QColor("#b91c1c"))
            painter.drawText(8, 17, self.aviso)

"""Esquema de geometría y cargas asignadas para la pantalla Inicio."""

from __future__ import annotations

from PySide6.QtCore import QPointF, QRect, Qt
from PySide6.QtGui import QColor, QFont, QPainter, QPen, QPolygonF
from PySide6.QtWidgets import QWidget

from calc import cargas, materiales, rutas


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

    def actualizar(self, nombre: str) -> None:
        self.nombre = nombre
        estructura = rutas.cargar_estructura().get(nombre, {}) if nombre else {}
        self.columnas = estructura.get("columnas", {})
        self.vigas = estructura.get("vigas", {})
        self.tramos = self._geometria_tramos()
        self.cargas = {}
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
            lineas.append(f'<span style="color:{self.colores[cid].name()};font-weight:bold">━━</span> {carga.get("descripcion", cid)} · {tipos}')
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
        margen_x, margen_y_superior, margen_y_inferior = 42.0, 92.0, 38.0
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

        if self.columnas:
            y_base = min(float(c.get("nivel", 0.0)) for c in self.columnas.values())
            y_cima = max(float(c.get("nivel", 0.0)) + float(c.get("altura_m", 0.0)) for c in self.columnas.values())
            extremo_izquierdo = min(float(c.get("x", 0.0)) for c in self.columnas.values())
            base_marco, cima_marco = punto(extremo_izquierdo, y_base), punto(extremo_izquierdo, y_cima)
            p_base = QPointF(base_marco.x() - 16, base_marco.y())
            p_cima = QPointF(cima_marco.x() - 16, cima_marco.y())
            painter.setPen(QPen(QColor("#64748b"), 1))
            painter.drawLine(p_base, p_cima)
            painter.drawLine(QPointF(p_base.x() - 4, p_base.y()), QPointF(p_base.x() + 5, p_base.y()))
            painter.drawLine(QPointF(p_cima.x() - 4, p_cima.y()), QPointF(p_cima.x() + 5, p_cima.y()))
            painter.drawText(QPointF(p_cima.x() - 4, p_cima.y() - 5), f"H={y_cima - y_base:.2f} m")

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
            texto_tramo = f"{tramo['id']} ({abs(tramo['x1'] - tramo['x0']):.2f} m)"
            painter.drawText(QPointF((a.x() + b.x()) / 2 - 35, a.y() + 18), texto_tramo)

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
            y_linea = base.y() - 16 - nivel * 30
            tipos = {}
            for carga in cargas_aplicacion:
                tipo = carga["tipo"]
                tipos[tipo] = tipos.get(tipo, 0.0) + float(carga["valor_kN_m"])
            color_linea = self.colores.get(carga_ref.get("carga_id", ""), QColor("#7c3aed"))
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
            valores = " · ".join(f"{t}: {v:.2f} kN/m" for t, v in tipos.items())
            sin_apoyo = " · SIN APOYO" if not geometria["apoyado"] else ""
            etiqueta = f"{valores}{sin_apoyo}"
            ancho_texto = painter.fontMetrics().horizontalAdvance(etiqueta)
            x_texto = max(4.0, min((base.x() + final.x() - ancho_texto) / 2, ancho - ancho_texto - 4))
            painter.drawText(QPointF(x_texto, y_linea - 4), etiqueta)

        # El motor actual concentra el resultante de viento en los nudos superiores.
        for nudo in self.viento_nodos:
            p = punto(nudo["x"], nudo["y"])
            y_flecha = p.y() - 8
            x0, x1 = p.x() - 11, p.x() + 14
            painter.setPen(QPen(QColor("#d97706"), 2.5))
            painter.drawLine(QPointF(x0, y_flecha), QPointF(x1, y_flecha))
            painter.setBrush(QColor("#d97706"))
            painter.drawPolygon(QPolygonF([
                QPointF(x1, y_flecha), QPointF(x1 - 6, y_flecha - 4), QPointF(x1 - 6, y_flecha + 4)
            ]))
            painter.setPen(QColor("#9a3412"))
            painter.drawText(QPointF(x1 + 3, y_flecha - 3), f"W {nudo['valor_kN']:.2f} kN")

        if self.aviso:
            painter.setPen(QColor("#b91c1c"))
            painter.drawText(8, 17, self.aviso)

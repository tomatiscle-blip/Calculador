"""Esquema de geometría y cargas asignadas para la pantalla Inicio."""

from __future__ import annotations

import math

from PySide6.QtCore import QPointF, QRectF, Qt, Signal
from PySide6.QtGui import QColor, QFont, QFontMetricsF, QPainter, QPen, QPolygonF
from PySide6.QtWidgets import QWidget

from calc import cargas, materiales, rutas

#-------MOMENTOS Y CORTES-----
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
            painter.end()
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
            painter.end()
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
                painter.end()

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

#--VISTA DE PORTICOS----
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
        self.muros: dict = {}
        self.tipo = "pórtico"
        self.muro_nombre = ""
        self.aviso = ""
        self.leyenda = "Sin cargas asignadas a este pórtico."
        self.ancho_viento_m = 0.0
        self.viento_lineal_kN_m = 0.0
        self.viento_nodos: list[dict] = []
        self.capas_carga = 0
        self.niveles: list[dict] = []
        self.colores: dict = {}

    def establecer_niveles(self, niveles: list[dict]) -> None:
        """Recibe los niveles de la obra activa para mostrarlos como referencia en el visor."""
        self.niveles = list(niveles or [])
        self.update()

    # ------------------------------------------------------------------
    # Helpers de referencia cruzada muro ↔ pórtico (estructura mixta)
    # ------------------------------------------------------------------
    def _datos_portico_activo(self) -> dict:
        if not self.nombre:
            return {}
        estructura = rutas.cargar_estructura() or {}
        return estructura.get(self.nombre, {})

    def _muros_referencia_para_portico(self) -> list[dict]:
        portico = self._datos_portico_activo()
        if not portico:
            return []
        ref_p = portico.get("referencia_planta", {}) or {}
        dir_p = str(ref_p.get("direccion", "x")).lower()
        pos_p = float(portico.get("posicion_planta_m", 0.0))
        resultado: list[dict] = []
        for muro in self.muros.values():
            ref_m = muro.get("referencia_planta", {}) or {}
            dir_m = str(ref_m.get("direccion", "x")).lower()
            pos_m = float(muro.get("posicion_planta_m", 0.0))
            if dir_m == dir_p and abs(pos_m - pos_p) <= 1e-4:
                resultado.append(muro)
        return resultado

    def _portico_referencia_para_muro(self) -> dict | None:
        if self.tipo != "muro" or not self.muro_nombre:
            return None
        muro = self.muros.get(self.muro_nombre)
        if not muro:
            return None
        ref_m = muro.get("referencia_planta", {}) or {}
        dir_m = str(ref_m.get("direccion", "x")).lower()
        pos_m = float(muro.get("posicion_planta_m", 0.0))
        estructura = rutas.cargar_estructura() or {}
        for nombre, datos in estructura.items():
            # Un pórtico en planta se reconoce por tener referencia_planta +
            # posicion_planta_m (no siempre tiene campo "tipo"). Se excluyen
            # los muros para no auto-referenciarse.
            if str(datos.get("tipo", "")).lower() == "muro":
                continue
            ref_p = datos.get("referencia_planta", {}) or {}
            if not ref_p or "direccion" not in ref_p:
                continue
            pos_p = float(datos.get("posicion_planta_m", 0.0))
            dir_p = str(ref_p.get("direccion", "x")).lower()
            if dir_p == dir_m and abs(pos_p - pos_m) <= 1e-4:
                datos = dict(datos)
                datos["nombre"] = nombre
                return datos
        return None

    def actualizar(self, nombre: str) -> None:
        self.nombre = nombre
        estructura = rutas.cargar_estructura().get(nombre, {}) if nombre else {}
        es_muro = isinstance(estructura, dict) and str(estructura.get("tipo", "")).lower() == "muro"
        if nombre and es_muro:
            self.tipo = "muro"
            self.muro_nombre = nombre
            self.columnas = {}
            self.vigas = {}
            # Solo este muro en el visor de inicio
            self.muros = {nombre: estructura}
            self.tramos = []
            self.cargas = {}
        else:
            self.tipo = "pórtico"
            self.muro_nombre = ""
            self.vigas = estructura.get("vigas", {})
            # Contexto: se conservan los muros existentes del proyecto para
            # dibujarlos tenues detrás del pórtico (ver paintEvent). Esto no
            # altera la escala del pórtico.
            estructura_completa = rutas.cargar_estructura() or {}
            # --- Construir self.columnas como diccionario {gid: datos} ---
            # datos["columnas"] puede ser:
            #  a) dict {cid: datos} (formato antiguo de pórticos creados antes del PRINCIPIO 4)
            #  b) lista de IDs globales ["COL-001", "COL-002"] (formato nuevo, PRINCIPIO 4)
            columnas_raw = estructura.get("columnas", {}) or {}
            columnas_globales = estructura_completa.get("columnas", {}) or {}
            if isinstance(columnas_raw, dict):
                self.columnas = columnas_raw
            elif isinstance(columnas_raw, list):
                self.columnas = {
                    gid: dict(columnas_globales.get(gid, {}))
                    for gid in columnas_raw
                    if gid in columnas_globales
                }
            else:
                self.columnas = {}
            # Coordenada local de columnas compartidas: recalcular a partir de
            # su posición global (X/Y) y la geometría del pórtico activo, para
            # que pórticos transversales no hereden la coordenada del creador.
            if self.tipo != "muro" and self.columnas:
                ref = estructura.get("referencia_planta", {}) or {}
                direccion = str(ref.get("direccion", "x")).lower()
                if direccion not in ("x", "y"):
                    direccion = "x"
                origen = VistaPlantaPorticos._numero(ref.get("origen_longitudinal_m"))
                if origen is None:
                    eje_longitudinal = VistaPlantaPorticos._buscar_eje(
                        self.ejes, direccion, ref.get("eje_longitudinal_id")
                    )
                    origen = (
                        VistaPlantaPorticos._numero(eje_longitudinal.get("coordenada_m"))
                        if eje_longitudinal is not None
                        else 0.0
                    )
                origen = (origen or 0.0) + (VistaPlantaPorticos._numero(ref.get("desfase_longitudinal_m")) or 0.0)
                for _cid, _col in self.columnas.items():
                    if not isinstance(_col, dict):
                        continue
                    if "X" not in _col and "Y" not in _col:
                        continue
                    try:
                        if direccion == "x":
                            _global = _col.get("X", _col.get("x"))
                        else:
                            _global = _col.get("Y", _col.get("x"))
                        if _global is not None:
                            _col["x"] = float(_global) - float(origen)
                    except (TypeError, ValueError):
                        pass
            self.muros = {
                n: d for n, d in estructura_completa.items()
                if isinstance(d, dict) and str(d.get("tipo", "")).lower() == "muro"
            }
            self.tramos = self._geometria_tramos()
        self.cargas = {}
        self.capas_carga = 0
        self.aviso = ""
        self.ancho_viento_m = 0.0
        self.viento_lineal_kN_m = 0.0
        self.viento_nodos = []
        if nombre and self.tipo != "muro":
            try:
                from calc.portico import _cargas_asignadas

                self.cargas, avisos = _cargas_asignadas(nombre)
                self.aviso = "; ".join(avisos)
            except Exception as exc:
                self.aviso = f"No se pudieron leer las cargas: {exc}"
        elif self.tipo == "muro":
            # El visor de Inicio solo dibuja pórticos. Los muros se visualizan
            # en la pestaña «Muros» con su propio visor y cargas aplicadas.
            self.cargas = {}
            self.aviso = "Muro seleccionado — ver pestaña «Muros» para su esquema y cargas."
        self.capas_carga = max((
            len({c.get("aplicacion_id", "") for c in cargas_tramo})
            for cargas_tramo in self.cargas.values()
        ), default=0)
        self.setMinimumHeight(max(250, 190 + self.capas_carga * 24))
        proyecto = rutas.leer_json(rutas.CARGAS, {}) or {}
        configuracion_viento = proyecto.get("viento", {}) or {}
        self.direccion_viento = str(configuracion_viento.get("direccion_viento", "x")).strip().lower()
        if self.direccion_viento not in ("x", "y"):
            self.direccion_viento = "x"
        self.anchos_viento_por_portico = configuracion_viento.get("anchos_por_portico", {}) or {}
        porticos_activos = configuracion_viento.get("porticos_activos")
        activo = configuracion_viento.get("activo")
        if porticos_activos is None:
            activo = bool(activo)
        else:
            activo = bool(porticos_activos) and (self.nombre in (porticos_activos or []))
        if activo:
            try:
                items_viento = cargas.items_viento_general(proyecto, materiales.cargar())
                self.viento_lineal_kN_m = sum(float(item["valor"]) for item in items_viento)
                if self.nombre in self.anchos_viento_por_portico and self.anchos_viento_por_portico[self.nombre]:
                    self.ancho_viento_m = float(self.anchos_viento_por_portico[self.nombre])
                else:
                    self.ancho_viento_m = float(configuracion_viento.get("ancho_tributario_m", 0.0))
                columnas_por_piso = {}
                for cid, columna in self.columnas.items():
                    # PRINCIPIO 4: las columnas ahora son globales (COL-XXX), ya no
                    # hay prefijo "C{piso}-". Se agrupan por su cota de base (nivel)
                    # para repartir el viento piso por piso.
                    clave = round(float(columna.get("nivel", 0.0)), 4)
                    columnas_por_piso.setdefault(clave, []).append(columna)
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
            detalle_por_portico = ""
            afectados = self.anchos_viento_por_portico
            if afectados:
                resumen = ", ".join(f"{n}={b:.2f}m" for n, b in afectados.items())
                detalle_por_portico = f"; b por pórtico afectado: {resumen}"
            lineas.append(
                f'<span style="color:#d97706;font-weight:bold">━━ W →</span> '
                f'Viento dirección {self.direccion_viento.upper()}: '
                f'{self.viento_lineal_kN_m:.2f} kN/m × H={altura_total:.2f} m '
                f'= {fuerza_total:.2f} kN en toda la estructura; '
                f'b tributario={self.ancho_viento_m:.2f} m{detalle_por_portico}; '
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
            cota_viga = viga.get("cota_m")
            if cota_viga is not None:
                columnas = sorted(
                    (c for c in self.columnas.values()
                     if abs(float(c.get("nivel", 0.0)) + float(c.get("altura_m", 0.0)) - float(cota_viga)) < 1e-6),
                    key=lambda c: float(c.get("x", 0.0)),
                )
            else:
                columnas = sorted(
                    (c for cid, c in self.columnas.items() if str(cid).startswith(f"C{piso}-")),
                    key=lambda c: float(c.get("x", 0.0)),
                )
                # PRINCIPIO 4: fallback para pórticos con IDs globales (COL-XXX)
                # que no tienen prefijo "C{piso}-": usar todas las columnas.
                if not columnas and self.columnas:
                    columnas = sorted(
                        self.columnas.values(),
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
        if self.tipo == "muro":
            xs = []
            ys = []
            for muro in self.muros.values():
                xs.append(float(muro.get("x_inicio", 0.0)))
                xs.append(float(muro.get("x_fin", 0.0)))
                ys.append(float(muro.get("nivel_inicio", 0.0)))
                ys.append(float(muro.get("nivel_fin", 0.0)))
            # Pórtico de referencia: incluir sus columnas en el escalado
            portico_ref = self._portico_referencia_para_muro()
            if portico_ref:
                cols_ref = portico_ref.get("columnas", {}) or {}
                if isinstance(cols_ref, (list, tuple)):
                    cat_col = self.estructura.get("columnas", {}) if isinstance(self.estructura, dict) else {}
                    cols_ref = {gid: (cat_col.get(gid, {}) if isinstance(cat_col, dict) else {}) for gid in cols_ref}
                for col in cols_ref.values():
                    cx = float(col.get("x", 0.0))
                    cy0 = float(col.get("nivel", 0.0))
                    cy1 = cy0 + float(col.get("altura_m", 0.0))
                    xs.append(cx)
                    ys.extend((cy0, cy1))
            if not xs or not ys:
                painter.setPen(QColor("#64748b"))
                painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "No hay geometría de muro seleccionada")
                painter.end()
                return
        elif not self.columnas:
            painter.setPen(QColor("#64748b"))
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "Seleccioná un pórtico con geometría")
            painter.end()
            return

        if self.tipo != "muro":
            xs = [float(c.get("x", 0.0)) for c in self.columnas.values()]
            ys = [float(c.get("nivel", 0.0)) for c in self.columnas.values()]
            ys.extend(float(c.get("nivel", 0.0)) + float(c.get("altura_m", 0.0)) for c in self.columnas.values())
            for tramo in self.tramos:
                xs.extend((tramo["x0"], tramo["x1"]))
                ys.append(tramo["y"])
            # Muros que apoyan en este pórtico: se incluyen en el escalado
            # para que se vean como referencia dentro del visor.
            for muro_ref in self._muros_referencia_para_portico():
                xs.extend((float(muro_ref.get("x_inicio", 0.0)), float(muro_ref.get("x_fin", 0.0))))
                ys.extend((float(muro_ref.get("nivel_inicio", 0.0)), float(muro_ref.get("nivel_fin", 0.0))))
            if not xs or not ys:
                painter.end()
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

        # ------------------------------------------------------------------
        # Líneas de nivel como referencia visual (el elemento no queda flotando)
        # ------------------------------------------------------------------
        if self.niveles:
            cota_minima = min(float(nn["cota_m"]) for nn in self.niveles)
            for n in self.niveles:
                cota = float(n["cota_m"])
                y_pantalla = punto(0.0, cota).y()
                nombre = n.get("nombre", "?")
                es_pb = abs(cota - cota_minima) < 1e-6
                color = QColor("#7c3aed") if es_pb else QColor("#15803d")
                grosor = 2 if es_pb else 1
                painter.setPen(QPen(color, grosor, Qt.PenStyle.DashLine))
                painter.drawLine(
                    int(margen_x), int(y_pantalla), int(ancho - margen_x), int(y_pantalla)
                )
                etiqueta = (
                    f"PLANTA BAJA · {nombre}  {cota:g} m" if es_pb
                    else f"{nombre}  {cota:g} m"
                )
                painter.setPen(color)
                painter.drawText(
                    int(margen_x + 4), int(y_pantalla) - 4, etiqueta
                )

        # ------------------------------------------------------------------
        # Dibujar elemento de referencia (estructura mixta) detrás del activo
        # ------------------------------------------------------------------
        if self.tipo != "muro":
            # Muros que apoyan en este pórtico: se dibujan como referencia tenue
            for muro_ref in self._muros_referencia_para_portico():
                x0 = float(muro_ref.get("x_inicio", 0.0))
                x1 = float(muro_ref.get("x_fin", 0.0))
                y0 = float(muro_ref.get("nivel_inicio", 0.0))
                y1 = float(muro_ref.get("nivel_fin", 0.0))
                espesor = float(muro_ref.get("espesor_m", 0.20))
                offset = float(muro_ref.get("offset", 0.0))
                nombre_muro = muro_ref.get("nombre", "muro")
                a = punto(x0 + offset, y0)
                b = punto(x1 + offset, y0)
                c = punto(x1 + offset, y1)
                d = punto(x0 + offset, y1)
                painter.setBrush(QColor(13, 148, 136, 90))  # teal muy tenue
                painter.setPen(QPen(QColor(13, 148, 136, 160), 2))
                painter.drawPolygon(QPolygonF([a, b, c, d]))
                painter.setBrush(Qt.BrushStyle.NoBrush)
                painter.setPen(QColor(13, 148, 136, 160))
                cx = (a.x() + c.x()) / 2
                painter.drawText(
                    QPointF(cx, min(a.y(), b.y()) - 8),
                    f"{nombre_muro} · muro de ref.",
                )
        else:
            # Muro activo: dibujar el pórtico que lo soporta como referencia
            portico_ref = self._portico_referencia_para_muro()
            if portico_ref:
                columnas_ref = portico_ref.get("columnas", {}) or {}
                color_ref = QColor(120, 130, 145, 150)
                painter.setPen(QPen(color_ref, 2))
                painter.setBrush(color_ref)
                for col in columnas_ref.values():
                    cx = float(col.get("x", 0.0))
                    cy0 = float(col.get("nivel", 0.0))
                    cy1 = cy0 + float(col.get("altura_m", 0.0))
                    painter.drawLine(punto(cx, cy0), punto(cx, cy1))
                    painter.drawEllipse(punto(cx, cy0), 3, 3)
                nombre_portico = portico_ref.get("nombre", "pórtico")
                if columnas_ref:
                    primera_col = next(iter(columnas_ref.values()))
                    cx = float(primera_col.get("x", 0.0))
                    cy0 = float(primera_col.get("nivel", 0.0))
                    cy1 = cy0 + float(primera_col.get("altura_m", 0.0))
                    painter.setPen(color_ref)
                    painter.drawText(
                        QPointF(punto(cx, cy1).x(), punto(cx, cy1).y() + 10),
                        f"{nombre_portico} · pórtico de ref.",
                    )

        if self.tipo == "muro":
            for indice, (nombre, muro) in enumerate(self.muros.items()):
                x0 = float(muro.get("x_inicio", 0.0))
                x1 = float(muro.get("x_fin", 0.0))
                y0 = float(muro.get("nivel_inicio", 0.0))
                y1 = float(muro.get("nivel_fin", 0.0))
                espesor = float(muro.get("espesor_m", 0.20))
                offset = float(muro.get("offset", 0.0))
                longitud = abs(x1 - x0)

                # convertir a puntos en pantalla
                a = punto(x0 + offset, y0)
                b = punto(x1 + offset, y0)
                c = punto(x1 + offset, y1)
                d = punto(x0 + offset, y1)

                # dibujar rectángulo del muro (relleno suave + borde)
                painter.setBrush(QColor("#0f766e"))
                painter.setPen(QPen(QColor("#0f766e"), 3))
                painter.drawPolygon(QPolygonF([a, b, c, d]))
                painter.setBrush(Qt.BrushStyle.NoBrush)

                # etiqueta rica; se apila en niveles para que no se superponga
                etiqueta = (
                    f"{nombre} · L={longitud:.2f} m · e={espesor:.2f} m · "
                    f"niveles {y0:.2f}-{y1:.2f} m"
                )
                painter.setPen(QColor("#0f172a"))
                cx = (a.x() + c.x()) / 2
                painter.drawText(
                    QPointF(cx, min(a.y(), b.y()) - 10 - indice * 16),
                    etiqueta,
                )

            # Un muro se dibuja por sí solo. Si tiene cargas aplicadas, se
            # representan como flechas sobre el paño (exclusivo de muros).
            if self.cargas:
                self._dibujar_cargas_muro(painter, punto)
            painter.end()
            return

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
        # distribución interna en los nudos superiores. La flecha apunta en el
        # sentido de la dirección del viento (X = horizontal, Y = vertical).
        if self.viento_nodos:
            x_marco = min(float(c.get("x", 0.0)) for c in self.columnas.values())
            y_base = min(float(c.get("nivel", 0.0)) for c in self.columnas.values())
            y_cima = max(float(c.get("nivel", 0.0)) + float(c.get("altura_m", 0.0)) for c in self.columnas.values())
            centro = punto(x_marco, (y_base + y_cima) / 2)
            color_viento = QColor("#d97706")
            if self.direccion_viento == "y":
                # Flecha vertical: viento sobre ejes Y (pórticos planos X al frente).
                x_flecha = centro.x()
                y1 = centro.y() + 8
                y0 = y1 + 32
                painter.setPen(QPen(color_viento, 2.5))
                painter.drawLine(QPointF(x_flecha, y0), QPointF(x_flecha, y1))
                painter.setBrush(color_viento)
                painter.drawPolygon(QPolygonF([
                    QPointF(x_flecha, y1),
                    QPointF(x_flecha - 4, y1 - 6),
                    QPointF(x_flecha + 4, y1 - 6),
                ]))
                painter.setPen(QColor("#9a3412"))
                total = sum(nudo["valor_kN"] for nudo in self.viento_nodos)
                painter.drawText(QPointF(x_flecha - 26, y0 - 7), f"W = {total:.2f} kN · dir. Y · toda la estructura")
            else:
                # Flecha horizontal: viento sobre ejes X (pórticos planos Y al frente).
                y_flecha = centro.y()
                x1 = centro.x() - 8
                x0 = x1 - 32
                painter.setPen(QPen(color_viento, 2.5))
                painter.drawLine(QPointF(x0, y_flecha), QPointF(x1, y_flecha))
                painter.setBrush(color_viento)
                painter.drawPolygon(QPolygonF([
                    QPointF(x1, y_flecha), QPointF(x1 - 6, y_flecha - 4), QPointF(x1 - 6, y_flecha + 4)
                ]))
                painter.setPen(QColor("#9a3412"))
                total = sum(nudo["valor_kN"] for nudo in self.viento_nodos)
                painter.drawText(QPointF(x0 - 12, y_flecha - 7), f"W = {total:.2f} kN · dir. X · toda la estructura")

        if self.aviso:
            painter.setPen(QColor("#b91c1c"))
            painter.drawText(8, 17, self.aviso)
        painter.end()

    # ---------------------------------------------------------------------------
    # Dibujo de cargas sobre muro (exclusivo del visor de muro)
    # ---------------------------------------------------------------------------
    def _dibujar_cargas_muro(self, painter, punto_fn) -> None:
        """Dibuja las cargas lineales aplicadas sobre el muro activo."""
        if not self.muros:
            return
        # Tomamos el primer (y único) muro del visor
        nombre, muro = next(iter(self.muros.items()))
        x0 = float(muro.get("x_inicio", 0.0))
        x1 = float(muro.get("x_fin", 0.0))
        y0 = float(muro.get("nivel_inicio", 0.0))
        y1 = float(muro.get("nivel_fin", 0.0))
        longitud = abs(x1 - x0)
        if longitud <= 0:
            return

        # Escala vertical: dibujar las cargas en la parte superior del muro
        altura_pantalla = abs(punto_fn(x0, y1).y() - punto_fn(x0, y0).y())
        if altura_pantalla <= 0:
            return

        for tramo_id, cargas_tramo in self.cargas.items():
            for carga in cargas_tramo:
                cx0 = float(carga.get("x_inicio_m", 0.0))
                cx1 = float(carga.get("x_fin_m", longitud))
                valor = float(carga.get("valor_kN_m", 0.0))
                tipo = str(carga.get("tipo", "D"))
                # Clamp al rango del muro
                cx0 = max(0.0, min(cx0, longitud))
                cx1 = max(0.0, min(cx1, longitud))
                if cx1 <= cx0:
                    continue
                # Convertir coordenadas del muro a pantalla
                p_inicio = punto_fn(x0 + cx0, y1)
                p_fin = punto_fn(x0 + cx1, y1)
                # Flecha hacia abajo (carga vertical): en coordenadas de pantalla
                # de Qt, y crece hacia abajo, así que y_flecha debe ser mayor.
                y_flecha = p_inicio.y() + 12
                painter.setPen(QPen(QColor("#b91c1c"), 2))
                painter.drawLine(QPointF((p_inicio.x() + p_fin.x()) / 2, p_inicio.y()),
                                 QPointF((p_inicio.x() + p_fin.x()) / 2, y_flecha))
                # Punta de flecha: apunta hacia abajo (y_flecha es la punta más abajo).
                painter.setBrush(QColor("#b91c1c"))
                painter.drawPolygon(QPolygonF([
                    QPointF((p_inicio.x() + p_fin.x()) / 2, y_flecha),
                    QPointF((p_inicio.x() + p_fin.x()) / 2 - 4, y_flecha - 6),
                    QPointF((p_inicio.x() + p_fin.x()) / 2 + 4, y_flecha - 6),
                ]))
                # Etiqueta con valor (debajo de la punta de la flecha).
                etiqueta = f"{valor:.2f} kN/m ({tipo})"
                ancho_txt = painter.fontMetrics().horizontalAdvance(etiqueta)
                x_txt = (p_inicio.x() + p_fin.x()) / 2 - ancho_txt / 2
                painter.setPen(QColor("#b91c1c"))
                painter.drawText(QPointF(x_txt, y_flecha + 14), etiqueta)

# ---------------------------------------------------------------------------
# Vista de planta: ejes X/Y con muros, pórticos y losas posicionados
# ---------------------------------------------------------------------------
class VistaPlantaPorticos(QWidget):
    """Planta legible por nivel: todos los elementos existentes en el nivel
    visible se dibujan en sólido; lo que pertenece a otro nivel va en
    punteado como proyección. Los rótulos se desplazan para no solaparse y
    cada losa muestra su tipo (alivianada, macisa, casetonada, cruzada...).
    """

    TOL = 0.01  # tolerancia (m) al comparar cotas de nivel

    # --- Paleta visual de la planta ---------------------------------------
    # Cada familia de elementos tiene color y símbolo propios para que se
    # distingan de un vistazo (y coincidan con la leyenda).
    COLOR_EJE = "#cbd5e1"            # grilla de ejes: referencia, punteada
    COLOR_EJE_BURBUJA = "#ffffff"    # fondo de la burbuja de identificación
    COLOR_EJE_BORDE = "#475569"      # borde de la burbuja de eje
    COLOR_PORTICO = "#64748b"        # pórtico inactivo
    COLOR_PORTICO_ACTIVO = "#1e3a8a"  # pórtico activo
    COLOR_VIGA = "#b45309"           # viga del nivel visible
    COLOR_VIGA_BORDE = "#7c2d12"     # contorno de la viga
    COLOR_MURO = "#dc2626"           # muro: banda maciza
    COLOR_MURO_BORDE = "#7f1d1d"     # contorno del muro
    COLOR_COLUMNA = "#2563eb"        # columna activa: nodo de planta
    COLOR_COLUMNA_INACTIVA = "#93c5fd"  # columna de pórtico inactivo
    COLOR_BASE = "#7c3aed"           # base: triángulo
    COLOR_APOYO = "#db2777"          # apoyo de losa: círculo hueco

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumSize(500, 400)
        self.setToolTip("Planta de ejes con muros, pórticos y losas (vista superior).")
        self.ejes = {}
        self.niveles = []
        self.estructura = {}
        self.losas = {}
        self.aplicaciones: list[dict] = []  # flechas de muro perpendicular (cargas.json)
        self.muros_seleccionados = set()
        self.portico_activo = ""
        self.nivel_activo = None  # cota Z del nivel filtrado; None = todos
        self.zoom_factor = 1.0  # multiplicador de escala para zoom con rueda
        # Paneo por arrastre: desplazamiento de la vista en píxeles.
        self.desplazamiento_x = 0.0
        self.desplazamiento_y = 0.0
        self._arrastrando = False
        self._posicion_arrastre = None
        # Datos de viento para dibujar la flecha direccional en planta
        self.viento_nodos: list[dict] = []
        self.direccion_viento: str = "x"
        self.viento_lineal_kN_m: float = 0.0
        self.ancho_viento_m: float = 0.0

    # ------------------------------------------------------------------
    # Vista: zoom, paneo y reencuadre
    # ------------------------------------------------------------------
    def reiniciar_vista(self) -> None:
        """Vuelve al encuadre original (sin zoom ni desplazamiento)."""
        self.zoom_factor = 1.0
        self.desplazamiento_x = 0.0
        self.desplazamiento_y = 0.0
        self.update()

    def mousePressEvent(self, event):
        """Arrastre con el botón izquierdo: mueve la planta (paneo)."""
        if event.button() != Qt.MouseButton.LeftButton:
            super().mousePressEvent(event)
            return
        self._arrastrando = True
        self._posicion_arrastre = event.position()
        self.setCursor(Qt.CursorShape.ClosedHandCursor)

    def mouseMoveEvent(self, event):
        """Sigue el arrastre: desplaza el encuadre en píxeles."""
        if not self._arrastrando or self._posicion_arrastre is None:
            super().mouseMoveEvent(event)
            return
        posicion = event.position()
        delta = posicion - self._posicion_arrastre
        self.desplazamiento_x += delta.x()
        self.desplazamiento_y += delta.y()
        self._posicion_arrastre = posicion
        self.update()

    def mouseReleaseEvent(self, event):
        """Termina el arrastre."""
        if event.button() != Qt.MouseButton.LeftButton:
            super().mouseReleaseEvent(event)
            return
        self._arrastrando = False
        self._posicion_arrastre = None
        self.unsetCursor()
        self.update()

    # ------------------------------------------------------------------
    # Datos de entrada
    # ------------------------------------------------------------------
    def set_datos(self, ejes, niveles, estructura):
        self.ejes = ejes or {}
        self.niveles = niveles or []
        self.estructura = estructura or {}
        self.update()

    def set_losas(self, elementos_carga):
        """Recibe los elementos de cargas.json; filtra los tipo losa/cubierta."""
        self.losas = {
            k: v for k, v in (elementos_carga or {}).items()
            if isinstance(v, dict) and str(v.get("tipo", "")).lower() in ("losa", "cubierta")
        }
        self.update()

    def set_aplicaciones(self, aplicaciones):
        """Recibe las aplicaciones de cargas.json (flechas de muro perpendicular).

        Se conservan, entre todo lo recibido, las aplicaciones de tipo
        ``muro_perpendicular`` que están activas. Cada una vincula dos
        pórticos mediante sus ``reacciones`` (porticó + tramo_id + x_m) y se
        resuelve al vuelo en ``_escena`` para dibujar la flecha entre los
        puntos de apoyo.
        """
        self.aplicaciones = [
            a for a in (aplicaciones or [])
            if isinstance(a, dict)
            and a.get("tipo_aplicacion") == "muro_perpendicular"
            and a.get("activa", True)
        ]
        self.update()

    def set_nivel_activo(self, cota):
        """Filtra la vista al nivel con esa cota; None muestra todos."""
        self.nivel_activo = cota
        self.update()

    def set_portico_activo(self, nombre):
        self.portico_activo = (nombre or "").strip()
        self.update()

    def set_datos_viento(self, viento_nodos: list[dict], direccion: str,
                         lineal_kN_m: float, ancho_m: float) -> None:
        """Recibe los datos de viento del proyecto para dibujar la flecha en planta."""
        self.viento_nodos = viento_nodos or []
        self.direccion_viento = str(direccion or "x").strip().lower()
        if self.direccion_viento not in ("x", "y"):
            self.direccion_viento = "x"
        self.viento_lineal_kN_m = float(lineal_kN_m or 0.0)
        self.ancho_viento_m = float(ancho_m or 0.0)
        self.update()

    # ------------------------------------------------------------------
    # Utilidades de planta
    # ------------------------------------------------------------------
    @staticmethod
    def _buscar_eje(ejes, familia, eje_id):
        """Devuelve el eje pedido en la familia (admite ids legados familia:nombre)."""
        if not eje_id:
            return None
        return next(
            (
                item for item in (ejes or {}).get(familia, [])
                if str(item.get("id", f"{familia}:{str(item.get('nombre', '')).casefold()}"))
                == str(eje_id)
            ),
            None,
        )

    @staticmethod
    def _numero(valor):
        try:
            numero = float(valor)
        except (TypeError, ValueError):
            return None
        return numero if math.isfinite(numero) else None

    def _geometria(self, datos):
        """(dirección, posición transversal, origen longitudinal) del elemento.

        La dirección es la familia de ejes por la que avanza el elemento; el eje
        de aplicación (posicion_planta_m) queda siempre en la familia opuesta.
        """
        ref = datos.get("referencia_planta", {}) or {}
        direccion = str(ref.get("direccion", "x")).lower()
        if direccion not in ("x", "y"):
            direccion = "x"
        eje_transversal = self._buscar_eje(
            self.ejes, "y" if direccion == "x" else "x", ref.get("eje_id")
        )
        posicion = (
            self._numero(eje_transversal.get("coordenada_m"))
            if eje_transversal is not None
            else None
        )
        if posicion is None:
            posicion = self._numero(datos.get("posicion_planta_m")) or 0.0
        if str(datos.get("tipo", "")).lower() == "muro":
            # Un muro guarda sus extremos en coordenadas absolutas de su eje.
            origen = 0.0
        else:
            origen = self._numero(ref.get("origen_longitudinal_m"))
            if origen is None:
                eje_longitudinal = self._buscar_eje(
                    self.ejes, direccion, ref.get("eje_longitudinal_id")
                )
                origen = (
                    self._numero(eje_longitudinal.get("coordenada_m"))
                    if eje_longitudinal is not None
                    else 0.0
                )
            origen += self._numero(ref.get("desfase_longitudinal_m")) or 0.0
        return direccion, posicion, origen

    @staticmethod
    def _punto(geom, longitud):
        """Coordenada de planta (x, y) para una medida sobre el elemento."""
        direccion, posicion, origen = geom
        valor = origen + longitud
        if direccion == "y":
            return posicion, valor
        return valor, posicion

    def _en_nivel(self, cota, nombre_nivel=None):
        """True si el elemento pertenece al nivel visible (o si no hay filtro)."""
        if self.nivel_activo is None:
            return True
        valor = self._numero(cota)
        if valor is None and nombre_nivel is not None:
            for nivel in self.niveles:
                if str(nivel.get("nombre", "")).strip() == str(nombre_nivel).strip():
                    valor = self._numero(nivel.get("cota_m"))
                    break
        if valor is None:
            # Sin cota conocida preferimos mostrarlo antes que esconderlo.
            return True
        return abs(valor - self.nivel_activo) <= self.TOL

    def _en_rango_nivel(self, cota_ini, cota_fin):
        """True si el rango de cotas de un muro toca el nivel visible."""
        if self.nivel_activo is None:
            return True
        ini = self._numero(cota_ini)
        fin = self._numero(cota_fin)
        if ini is None or fin is None:
            return True
        bajo, alto = sorted((ini, fin))
        return (bajo - self.TOL) <= self.nivel_activo <= (alto + self.TOL)

    def _cota_base(self):
        """Cota de las bases: el nivel más bajo de la obra."""
        cotas = [
            cota for cota in (self._numero(n.get("cota_m")) for n in self.niveles)
            if cota is not None
        ]
        return min(cotas) if cotas else 0.0

    def _rect_panel(self, panel):
        """Rectángulo (x0, y0, x1, y1) de una losa delimitada por ejes."""
        lados = panel.get("ejes", {}) or {}
        rect = {}
        for familia in ("x", "y"):
            lado = lados.get(familia) or {}
            inicio = (lado.get("inicio") or {}).get("coordenada_m")
            fin = (lado.get("fin") or {}).get("coordenada_m")
            a, b = self._numero(inicio), self._numero(fin)
            if a is None or b is None:
                return None
            rect[familia] = (min(a, b), max(a, b))
        (x0, x1), (y0, y1) = rect["x"], rect["y"]
        if x1 - x0 <= 0 or y1 - y0 <= 0:
            return None
        return x0, y0, x1, y1

    def _rect_legado(self, datos):
        """Banda de una losa sin panel de ejes, ubicada por sus pórticos de apoyo."""
        refs = datos.get("porticos", datos.get("aplicar_a", []))
        if isinstance(refs, str):
            refs = [refs]
        geoms, puntos = [], []
        for ref in refs:
            datos_p = self.estructura.get(str(ref))
            if not isinstance(datos_p, dict) or str(datos_p.get("tipo", "")).lower() == "muro":
                continue
            geom = self._geometria(datos_p)
            geoms.append(geom)
            cols = datos_p.get("columnas", {}) or {}
            if isinstance(cols, (list, tuple)):
                cat_col = self.estructura.get("columnas", {}) if isinstance(self.estructura, dict) else {}
                cols = {gid: (cat_col.get(gid, {}) if isinstance(cat_col, dict) else {}) for gid in cols}
            for columna in cols.values():
                if not isinstance(columna, dict):
                    continue
                x = self._numero(columna.get("x"))
                if x is not None:
                    puntos.append(self._punto(geom, x))
        if not geoms or not puntos:
            return None
        try:
            ancho = self._numero(datos.get("ancho_losa_m")) or 1.0
        except (TypeError, ValueError):
            ancho = 1.0
        xs = [p[0] for p in puntos]
        ys = [p[1] for p in puntos]
        if all(geom[0] == "x" for geom in geoms):
            y = sum(geom[1] for geom in geoms) / len(geoms)
            return min(xs), y - ancho / 2, max(xs), y + ancho / 2
        if all(geom[0] == "y" for geom in geoms):
            x = sum(geom[1] for geom in geoms) / len(geoms)
            return x - ancho / 2, min(ys), x + ancho / 2, max(ys)
        return None

    def _rotulo(self, x, y, texto, direccion="x", color="#334155", prioridad=5,
                desplazamiento=0.0, tamano=None, flotante=True):
        """Rótulo de planta orientado según el eje del elemento.

        ``direccion`` es la familia de ejes por la que corre el elemento: con
        "x" el texto queda horizontal (lee en el sentido X) y con "y" queda
        vertical, para que el rótulo siga la dirección de la viga, el muro o el
        pórtico en lugar de cruzarse con él.

        ``flotante`` marca los rótulos técnicos dispersos (nombres de vigas,
        muros, columnas, bases, losas...). Se conservan en la escena para
        poder consultarlos, pero no se dibujan sobre el plano: el dibujo se
        mantiene limpio y la información queda en la leyenda.
        """
        return {
            "x": float(x),
            "y": float(y),
            "texto": str(texto),
            "giro": 0.0 if str(direccion).lower() == "x" else 90.0,
            "color": color,
            "prioridad": int(prioridad),
            "desplazamiento": float(desplazamiento),
            "tamano": tamano,
            "flotante": bool(flotante),
        }

    def _colocar_rotulos(self, rotulos, px, py):
        """Calcula la posición final de cada rótulo evitando superposiciones.

        Devuelve una lista de tuplas ``(rotulo, rectangulo)`` con los rótulos
        que se pudieron colocar. Los que no tienen lugar libre se descartan:
        es preferible omitir un texto antes que tapar otro.
        """
        if not rotulos:
            return []
        fuente_base = QFont(self.font())
        fuente_base.setPointSize(max(9, fuente_base.pointSize() or 9))
        metricas = QFontMetricsF(fuente_base)
        colocados: list[QRectF] = []
        resultado: list[tuple[dict, QRectF]] = []
        for rot in sorted(rotulos, key=lambda item: item["prioridad"]):
            texto = rot["texto"]
            if not texto:
                continue
            fuente = QFont(fuente_base)
            if rot.get("tamano"):
                fuente.setPointSize(int(rot["tamano"]))
            fuente.setBold(int(rot["prioridad"]) <= 2)
            metricas = QFontMetricsF(fuente)
            ancho = metricas.horizontalAdvance(texto)
            alto = metricas.height()
            if rot["giro"]:
                ancho, alto = alto, ancho
            cx = px(rot["x"])
            cy = py(rot["y"])
            rect = QRectF(cx - ancho / 2, cy - alto / 2, ancho, alto)
            paso = rot["desplazamiento"] or 14.0
            if rot["desplazamiento"]:
                # El rótulo se aparta de la línea del elemento (perpendicular).
                if rot["giro"]:
                    rect.translate(paso, 0)
                else:
                    rect.translate(0, -paso)

            limite = 2.0
            ancho_widget = float(self.width())
            alto_widget = float(self.height())

            def _ajustar(caja):
                """Mueve la caja para que no se salga del widget."""
                dx = dy = 0.0
                if caja.left() < limite:
                    dx = limite - caja.left()
                elif caja.right() > ancho_widget - limite:
                    dx = ancho_widget - limite - caja.right()
                if caja.top() < limite:
                    dy = limite - caja.top()
                elif caja.bottom() > alto_widget - limite:
                    dy = alto_widget - limite - caja.bottom()
                if not dx and not dy:
                    return caja
                return QRectF(caja).translated(dx, dy)

            # Búsqueda de hueco: se recorre una espiral de posiciones alrededor
            # del rótulo (empezando por la perpendicular al texto) hasta
            # encontrar un lugar libre. Con muchos rótulos cerca (columnas en
            # grilla, pórticos que se cruzan) hace falta explorar en las dos
            # dimensiones, no solo a lo largo de un eje.
            candidato = rect
            intentos = 0
            max_intentos = 48
            offsets = (
                (0, -1), (0, 1), (-1, 0), (1, 0),
                (-1, -1), (1, 1), (-1, 1), (1, -1),
            )
            while True:
                candidato = _ajustar(candidato)
                if not any(candidato.intersects(otro) for otro in colocados):
                    break
                intentos += 1
                if intentos >= max_intentos:
                    candidato = None
                    break
                fila, columna = offsets[(intentos - 1) % len(offsets)]
                anillo = intentos // len(offsets) + 1
                salto = paso * anillo
                if rot["giro"]:
                    # Texto vertical: la perpendicular al texto es la horizontal.
                    dx = salto * fila
                    dy = salto * columna * 0.6
                else:
                    # Texto horizontal: la perpendicular al texto es la vertical.
                    dy = salto * fila
                    dx = salto * columna * 0.6
                candidato = QRectF(rect)
                candidato.translate(QPointF(dx, dy))
            if candidato is None:
                continue  # no hay lugar: mejor omitirlo que tapar otro rótulo
            colocados.append(candidato)
            resultado.append((rot, candidato))
        return resultado

    def _pintar_rotulos(self, painter, rotulos, px, py):
        """Dibuja los rótulos de planta centrados y sin superponerse.

        Cada rótulo se dibuja centrado en su elemento y con la orientación del
        eje al que pertenece. Si dos rótulos se superponen, el de menor
        prioridad (el más reciente) se corre lo necesario para separarse; si
        no hay lugar, se omite para no ensuciar el dibujo.

        Los rótulos marcados como ``flotante`` (rótulos técnicos dispersos) se
        conservan en la escena pero no se dibujan: el plano queda limpio y la
        información sigue disponible en la leyenda y en las consultas.
        """
        if not rotulos:
            return
        visibles = [rot for rot in rotulos if not rot.get("flotante", False)]
        for rot, rect in self._colocar_rotulos(visibles, px, py):
            texto = rot["texto"]
            fuente = QFont(self.font())
            fuente.setPointSize(max(9, fuente.pointSize() or 9))
            if rot.get("tamano"):
                fuente.setPointSize(int(rot["tamano"]))
            fuente.setBold(int(rot["prioridad"]) <= 2)
            metricas = QFontMetricsF(fuente)
            ancho = metricas.horizontalAdvance(texto)
            alto = metricas.height()
            if rot["giro"]:
                ancho, alto = alto, ancho
            painter.save()
            painter.setFont(fuente)
            painter.setPen(QPen(QColor(rot["color"]), 1))
            painter.translate(rect.center().x(), rect.center().y())
            if rot["giro"]:
                painter.rotate(rot["giro"])
            painter.drawText(
                QRectF(-ancho / 2, -alto / 2, ancho, alto),
                Qt.AlignmentFlag.AlignCenter,
                texto,
            )
            painter.restore()

    def _apoyo_losa(self, datos_losa, panel, rect, lado):
        """Punto de planta (x, y) y texto de un apoyo de losa (izq/der).

        La losa tiene dos bordes de apoyo, dados por la dirección de la luz:
        con luz en X los apoyos están sobre los bordes X del paño y con luz en
        Y sobre los bordes Y. El punto devuelto está en el medio del borde, que
        es donde se apoya la losa. Se devuelve ``None`` si el apoyo no está
        declarado o si no se puede verificar (pórtico inexistente, viga que no
        cubre el borde, muro sin posición).
        """
        apoyos = datos_losa.get("apoya_en") or {}
        destino = apoyos.get(lado) or {}
        if not isinstance(destino, dict) or not destino:
            return None
        x0, y0, x1, y1 = rect
        direccion = str((panel.get("direccion_luz", "x"))).lower()
        if direccion not in ("x", "y"):
            direccion = "x"
        centro_x, centro_y = (x0 + x1) / 2, (y0 + y1) / 2
        if direccion == "x":
            borde = x0 if lado == "izq" else x1
            punto = (borde, centro_y)
        else:
            borde = y0 if lado == "izq" else y1
            punto = (centro_x, borde)

        tipo = str(destino.get("tipo", "")).lower()
        if tipo == "muro":
            nombre_m = str(destino.get("nombre", ""))
            datos_m = self.estructura.get(nombre_m)
            if not isinstance(datos_m, dict):
                return None
            geom = self._geometria(datos_m)
            if geom is None:
                return None
            # El muro debe contener el borde del paño (dentro de tolerancia).
            x_ini = self._numero(datos_m.get("x_inicio"))
            x_fin = self._numero(datos_m.get("x_fin"))
            if x_ini is None or x_fin is None:
                return None
            muro_direccion, muro_posicion, _ = geom
            if muro_direccion == "x":
                if not (min(x_ini, x_fin) - 0.5 <= borde <= max(x_ini, x_fin) + 0.5):
                    return None
                return muro_posicion, punto[1], nombre_m
            if not (min(x_ini, x_fin) - 0.5 <= borde <= max(x_ini, x_fin) + 0.5):
                return None
            return punto[0], muro_posicion, nombre_m
        if tipo == "portico":
            nombre_p = str(destino.get("portico", ""))
            viga_id = str(destino.get("viga", ""))
            datos_p = self.estructura.get(nombre_p)
            if not isinstance(datos_p, dict):
                return None
            viga = (datos_p.get("vigas", {}) or {}).get(viga_id) or {}
            if not isinstance(viga, dict):
                return None
            geom = self._geometria(datos_p)
            if geom is None:
                return None
            # La viga de apoyo corre perpendicular a la luz de la losa.
            esperado = "y" if direccion == "x" else "x"
            if geom[0] != esperado:
                return None
            # La viga debe tener un tramo que cubra el borde del paño.
            cubre = False
            for tramo in (viga.get("tramos") or []):
                if not isinstance(tramo, dict):
                    continue
                a = self._numero(tramo.get("x_inicio"))
                b = self._numero(tramo.get("x_fin"))
                if a is None or b is None:
                    continue
                if min(a, b) - 0.5 <= borde <= max(a, b) + 0.5:
                    cubre = True
                    break
            if not cubre:
                return None
            texto = f"{nombre_p} · {viga_id}" if viga_id else nombre_p
            return punto[0], punto[1], texto
        return None


    def _escena(self):
        """Arma la planta del nivel visible: todos los elementos existentes
        en sólido, las proyecciones de otros niveles en punteado y los rótulos
        con evasión de solapamiento. No hay enfoque único: se dibuja todo lo
        que pertenece al nivel elegido (o todo, si no hay filtro).
        """
        porticos, muros, losas, etiquetas = [], [], [], []
        cargas_muro: list[dict] = []
        rotulos: list[dict] = []
        apoyos_losa: list[dict] = []
        proyecciones: list[dict] = []
        puntos: list[tuple[float, float]] = []
        conteo = {"columnas": 0, "vigas": 0, "bases": 0, "muros": 0, "losas": 0}
        # Registro global de columnas por coordenada de planta (x, y).
        # Las columnas en cruces de ejes pertenecen a dos pórticos (uno en X,
        # otro en Y) pero son la misma columna física: se reutiliza el nodo.
        columnas_globales: dict[tuple[float, float], list[dict]] = {}
        # Registro global de bases por coordenada de planta: dos bases que
        # comparten el mismo punto global son la misma base física (aunque
        # tengan base_id distintos). Se unifican por identidad geométrica.
        bases_globales: dict[tuple[float, float], list[dict]] = {}
        cota_base = self._cota_base()
        # Niveles ordenados por cota: el nivel superior al visible tiene su
        # viga como techo punteado de este nivel.
        niveles_ordenados = sorted(
            (n for n in (self.niveles or []) if isinstance(n, dict)),
            key=lambda n: self._numero(n.get("cota_m")) or 0.0,
        )
        cotas = [self._numero(n.get("cota_m")) for n in niveles_ordenados]
        cotas = [c for c in cotas if c is not None]
        cota_superior = None
        if self.nivel_activo is not None and cotas:
            mayores = [c for c in cotas if c > self.nivel_activo + self.TOL]
            cota_superior = min(mayores) if mayores else None

        def _nombre_eje_aplicacion(datos: dict, direccion: str) -> str:
            """Nombre del eje transversal donde corre el elemento (Y3, X2...)."""
            ref = datos.get("referencia_planta", {}) or {}
            eje = self._buscar_eje(
                self.ejes, "y" if direccion == "x" else "x", ref.get("eje_id")
            )
            if eje is not None:
                return str(eje.get("nombre", ""))
            posicion = self._numero(datos.get("posicion_planta_m"))
            if posicion is None:
                return ""
            return f"{posicion:.2f} m"

        # Tipo de losa por nombre, para pintar y para la leyenda.
        tipos_losa: dict[str, str] = {}

        for nombre, datos in self.estructura.items():
            if not isinstance(datos, dict):
                continue
            geom = self._geometria(datos)
            if geom is None:
                continue
            direccion, posicion, _ = geom
            es_muro = str(datos.get("tipo", "")).lower() == "muro"

            if es_muro:
                # ---------------- muros ----------------
                eje_nombre = _nombre_eje_aplicacion(datos, direccion)
                x0 = self._numero(datos.get("x_inicio"))
                x1 = self._numero(datos.get("x_fin"))
                if x0 is None or x1 is None:
                    continue
                puntos.extend([self._punto(geom, x0), self._punto(geom, x1)])
                ni = self._numero(datos.get("nivel_inicio"))
                en_nivel = self._en_rango_nivel(
                    datos.get("nivel_inicio"), datos.get("nivel_fin")
                )
                # SÓLIDO = nace en este nivel. Si solo lo atraviesa (nació
                # abajo), va punteado con aviso para no confundir.
                nace_aca = (
                    self.nivel_activo is not None
                    and ni is not None
                    and abs(ni - self.nivel_activo) <= self.TOL
                )
                longitud = abs(max(x0, x1) - min(x0, x1))
                espesor = self._numero(datos.get("espesor_m")) or 0.15
                if en_nivel and (self.nivel_activo is None or nace_aca):
                    muros.append({
                        "nombre": str(nombre),
                        "direccion": direccion,
                        "posicion": posicion,
                        "x0": min(x0, x1),
                        "x1": max(x0, x1),
                        "espesor": espesor,
                        "activo": str(nombre) == str(self.portico_activo),
                    })
                    rotulos.append(self._rotulo(
                        x=(min(x0, x1) + max(x0, x1)) / 2,
                        y=posicion,
                        texto=f"{nombre} · eje {eje_nombre}",
                        direccion=direccion,
                        color="#b91c1c",
                        prioridad=2,
                        desplazamiento=12.0,
                        # Solo se rotula el muro activo.
                        flotante=str(nombre) != str(self.portico_activo),
                    ))
                    conteo["muros"] += 1
                elif en_nivel:
                    # Atraviesa el nivel pero nació abajo: proyección.
                    proyecciones.append({
                        "tipo": "muro",
                        "nombre": str(nombre),
                        "direccion": direccion,
                        "posicion": posicion,
                        "x0": min(x0, x1),
                        "x1": max(x0, x1),
                        "nota": "nace abajo",
                    })
                    rotulos.append(self._rotulo(
                        x=(min(x0, x1) + max(x0, x1)) / 2,
                        y=posicion,
                        texto=f"{nombre} (viene de abajo)",
                        direccion=direccion,
                        color="#b91c1c",
                        prioridad=4,
                        desplazamiento=12.0,
                    ))
                continue  # los muros no son pórticos: no seguir en este bloque

            # ---------------- pórticos ----------------
            columnas_datos = datos.get("columnas", {}) or {}
            vigas_datos = datos.get("vigas", {}) or {}
            bases_datos = datos.get("bases", {}) or {}
            # PRINCIPIO 4: los pórticos nuevos almacenan listas de IDs globales
            # (COL-XXX / BASE-XXX) en lugar de dicts locales. Resolverlas contra
            # el catálogo global de la estructura.
            if isinstance(columnas_datos, (list, tuple)):
                cat_col = self.estructura.get("columnas", {}) if isinstance(self.estructura, dict) else {}
                columnas_datos = {
                    gid: (cat_col.get(gid, {}) if isinstance(cat_col, dict) else {})
                    for gid in columnas_datos
                }
            if isinstance(bases_datos, (list, tuple)):
                cat_base = self.estructura.get("bases", {}) if isinstance(self.estructura, dict) else {}
                bases_datos = {
                    bid: (cat_base.get(bid, {}) if isinstance(cat_base, dict) else {})
                    for bid in bases_datos
                }
            eje_nombre = _nombre_eje_aplicacion(datos, direccion)
            # Nombre del eje longitudinal inicial (X1, Y2...): dice desde
            # dónde se miden las x locales de columnas/vigas.
            ref = datos.get("referencia_planta", {}) or {}
            eje_long = self._buscar_eje(
                self.ejes, direccion, ref.get("eje_longitudinal_id")
            )
            eje_long_nombre = (
                str(eje_long.get("nombre", "")) if eje_long is not None else ""
            )

            # Eje del pórtico: siempre se dibuja, con todas sus columnas de
            # referencia (incluso si en este nivel no hay ninguna).
            largos = [
                self._numero(col.get("x"))
                for col in columnas_datos.values()
                if isinstance(col, dict)
            ]
            largos = [x for x in largos if x is not None]
            eje_ini = min(largos) if largos else 0.0
            eje_fin = max(largos) if largos else 0.0
            eje_p0, eje_p1 = self._punto(geom, eje_ini), self._punto(geom, eje_fin)
            puntos.extend([eje_p0, eje_p1])

            columnas = []
            for cid, col in columnas_datos.items():
                if not isinstance(col, dict):
                    continue
                x = self._numero(col.get("x"))
                if x is None or not self._en_nivel(col.get("nivel"), col.get("nivel_inicio")):
                    continue
                punto = self._punto(geom, x)
                puntos.append(punto)
                # Eje de la columna (X2, Y1...): dice en qué cruce está parada.
                eje_col = ""
                col_eje = col.get("eje_nombre") or ""
                if col_eje:
                    eje_col = str(col_eje)
                else:
                    buscado = self._buscar_eje(self.ejes, direccion, col.get("eje_id"))
                    if buscado is not None:
                        eje_col = str(buscado.get("nombre", ""))
                columnas.append((str(cid), punto, eje_col))
                conteo["columnas"] += 1
                # Registro global por coordenada de planta: deduplica columnas
                # que están en el mismo punto pero definidas en dos pórticos
                # diferentes (cruce de ejes X y Y).
                clave = (round(punto[0], 3), round(punto[1], 3))
                columnas_globales.setdefault(clave, []).append({
                    "nombre": str(cid),
                    "punto": punto,
                    "eje": eje_col,
                    "portico": str(nombre),
                })

            vigas = []
            for vid, viga in vigas_datos.items():
                if not isinstance(viga, dict):
                    continue
                extremos = [
                    (self._numero(t.get("x_inicio")), self._numero(t.get("x_fin")))
                    for t in (viga.get("tramos") or [])
                    if isinstance(t, dict)
                ]
                extremos = [(a, b) for a, b in extremos if a is not None and b is not None]
                if not extremos:
                    continue
                inicio = min(min(a, b) for a, b in extremos)
                fin = max(max(a, b) for a, b in extremos)
                p0, p1 = self._punto(geom, inicio), self._punto(geom, fin)
                puntos.extend([p0, p1])
                if self._en_nivel(viga.get("cota_m"), viga.get("nivel_nombre")):
                    vigas.append((str(vid), p0, p1))
                    rotulos.append(self._rotulo(
                        x=(p0[0] + p1[0]) / 2,
                        y=(p0[1] + p1[1]) / 2,
                        texto=f"{vid} · eje {eje_nombre}",
                        direccion=direccion,
                        color="#b45309",
                        prioridad=3,
                        desplazamiento=10.0,
                        # Rótulo técnico: no se dibuja sobre el plano.
                        flotante=True,
                    ))
                    conteo["vigas"] += 1
                elif (
                    cota_superior is not None
                    and self._numero(viga.get("cota_m")) is not None
                    and abs(
                        (self._numero(viga.get("cota_m")) or 0.0) - cota_superior
                    ) <= self.TOL
                ):
                    # Viga del nivel SUPERIOR = techo de este nivel: se dibuja
                    # punteada para ubicar, sin contar como viga del nivel.
                    proyecciones.append({
                        "tipo": "viga",
                        "nombre": str(vid),
                        "direccion": direccion,
                        "p0": p0,
                        "p1": p1,
                        "nota": "techo del nivel",
                    })
                    rotulos.append(self._rotulo(
                        x=(p0[0] + p1[0]) / 2,
                        y=(p0[1] + p1[1]) / 2,
                        texto=f"{vid} (techo)",
                        direccion=direccion,
                        color="#b45309",
                        prioridad=4,
                        desplazamiento=10.0,
                        # Rótulo técnico: no se dibuja sobre el plano.
                        flotante=True,
                    ))

            bases = []
            if self._en_nivel(cota_base):
                for bid, base in bases_datos.items():
                    if not isinstance(base, dict):
                        continue
                    x = self._numero(base.get("x"))
                    if x is None:
                        continue
                    punto = self._punto(geom, x)
                    puntos.append(punto)
                    bases.append((str(bid), punto))
                    # Identidad global única de la base: punto global redondeado.
                    clave = (round(punto[0], 3), round(punto[1], 3))
                    bases_globales.setdefault(clave, []).append({
                        "bid": str(bid),
                        "punto": punto,
                        "portico": str(nombre),
                    })

            porticos.append({
                "nombre": str(nombre),
                "activo": str(nombre) == str(self.portico_activo),
                "direccion": direccion,
                "eje_p0": eje_p0,
                "eje_p1": eje_p1,
                "columnas": columnas,
                "vigas": vigas,
                "bases": bases,
            })
            # Rótulo con el eje de inserción: dice sobre qué eje corre.
            detalle_eje = f" · {eje_nombre}" if eje_nombre else ""
            rotulos.append(self._rotulo(
                x=(eje_p0[0] + eje_p1[0]) / 2,
                y=(eje_p0[1] + eje_p1[1]) / 2,
                texto=f"{nombre}{detalle_eje}",
                direccion=direccion,
                color="#1d4ed8" if str(nombre) == str(self.portico_activo) else "#475569",
                prioridad=0,
                desplazamiento=14.0,
                # Solo se rotula el pórtico activo: los demás nombres son
                # rótulos técnicos que ensucian el plano.
                flotante=str(nombre) != str(self.portico_activo),
            ))

        # Columnas únicas por coordenada global: las columnas en cruces de ejes
        # pertenecen a dos pórticos pero son un solo nodo. Se dibuja una sola vez.
        columnas_unicas = []
        conteo_columnas_unicas = 0
        for clave, entradas in columnas_globales.items():
            conteo_columnas_unicas += 1
            punto = entradas[0]["punto"]
            ejes = sorted({e["eje"] for e in entradas if e["eje"]})
            nombre_principal = entradas[0]["nombre"]
            if len(entradas) > 1:
                nombre_principal = f"{nombre_principal} (+{len(entradas)-1})"
            columnas_unicas.append({
                "nombre": nombre_principal,
                "punto": punto,
                "ejes": ejes,
                "porticos": [e["portico"] for e in entradas],
            })
        conteo["columnas"] = conteo_columnas_unicas

        # Rótulos de columnas: van al sistema de rótulos (con evasión de
        # solapamiento) para que dos columnas vecinas no se tapen.
        # Para no saturar la planta se rotula solo lo relevante: si hay un
        # pórtico activo se rotulan sus columnas; si no hay ninguno (vista
        # general de la obra) se rotulan únicamente cuando son pocas.
        if self.portico_activo:
            columnas_rotuladas = [
                c for c in columnas_unicas
                if any(str(pn) == str(self.portico_activo) for pn in c["porticos"])
            ]
        elif len(columnas_unicas) <= 10:
            columnas_rotuladas = list(columnas_unicas)
        else:
            columnas_rotuladas = []

        # Rótulos de columnas: van al sistema de rótulos (con evasión de
        # solapamiento) para que dos columnas vecinas no se tapen.
        for col in columnas_rotuladas:
            punto = col["punto"]
            ejes_txt = ", ".join(col["ejes"]) if col["ejes"] else ""
            texto = f"{col['nombre']} ({ejes_txt})" if ejes_txt else col["nombre"]
            es_activa = any(
                str(pn) == str(self.portico_activo) for pn in col["porticos"]
            )
            rotulos.append(self._rotulo(
                x=punto[0],
                y=punto[1],
                texto=texto,
                direccion="x",
                color="#1d4ed8" if es_activa else "#475569",
                prioridad=2,
                desplazamiento=13.0,
                # Solo se rotulan las columnas del pórtico activo.
                flotante=not es_activa,
            ))
        # Bases únicas por coordenada global: bases que comparten el mismo
        # punto físico (aunque tengan base_id distintos) se dibujan una vez.
        bases_unicas = []
        for clave, entradas in bases_globales.items():
            punto = entradas[0]["punto"]
            bids = [e["bid"] for e in entradas]
            nombre_principal = entradas[0]["bid"]
            if len(entradas) > 1:
                nombre_principal = f"{nombre_principal} (+{len(entradas)-1})"
            bases_unicas.append({
                "nombre": nombre_principal,
                "punto": punto,
                "porticos": [e["portico"] for e in entradas],
                "bases": bids,
            })
        conteo["bases"] = len(bases_unicas)

        # Rótulos de bases: también con evasión de solapamiento. Cuando la base
        # está bajo una columna se omite su rótulo: la columna ya identifica
        # ese punto y los dos textos se amontonarían.
        puntos_columna = {
            (round(c["punto"][0], 3), round(c["punto"][1], 3))
            for c in columnas_unicas
        }
        for base in bases_unicas:
            punto = base["punto"]
            clave = (round(punto[0], 3), round(punto[1], 3))
            if clave in puntos_columna:
                continue
            if not self.portico_activo and len(bases_unicas) > 10:
                # Vista general con muchas bases: no se rotula ninguna para
                # no tapar los nombres de pórticos y muros.
                break
            rotulos.append(self._rotulo(
                x=punto[0],
                y=punto[1],
                texto=base["nombre"],
                direccion="x",
                color="#5b21b6",
                prioridad=3,
                desplazamiento=20.0,
                # Rótulo técnico: no se dibuja sobre el plano.
                flotante=True,
            ))

        # ---------------- flechas de muro perpendicular entre pórticos ----------------
        # Cada aplicación muro_perpendicular víncula dos pórticos mediante sus
        # reacciones (porticó + tramo_id + x_m sobre la viga). Se resuelve la
        # geometría de cada pórtico al momento y se dibuja la flecha doble entre
        # los dos puntos de apoyo, en la dirección de los pórticos. El tramo_id
        # no altera la geometría (x_m ya es la posición longitudinal sobre la
        # viga), pero se conserva por si se quiere validar la cobertura.
        for aplicacion in self.aplicaciones:
            reacciones = aplicacion.get("reacciones") or []
            if len(reacciones) < 2:
                continue
            puntos_apoyo: list[tuple[float, float]] = []
            direccion = None
            for reaccion in reacciones:
                nombre_p = str(reaccion.get("portico", "")).strip()
                datos_p = self.estructura.get(nombre_p)
                if not isinstance(datos_p, dict):
                    puntos_apoyo = []
                    break
                geom = self._geometria(datos_p)
                if geom is None:
                    puntos_apoyo = []
                    break
                if direccion is None:
                    direccion = geom[0]
                x_m = self._numero(reaccion.get("x_m"))
                if x_m is None:
                    puntos_apoyo = []
                    break
                puntos_apoyo.append(self._punto(geom, x_m))
            if len(puntos_apoyo) < 2:
                continue
            cargas_muro.append({
                "p1": puntos_apoyo[0],
                "p2": puntos_apoyo[1],
                "direccion": direccion or "x",
                "influencia_m": self._numero(
                    aplicacion.get("separacion_m")
                ) or 0.0,
            })
            puntos.extend(puntos_apoyo)

        # ---------------- losas ----------------
        for nombre, datos in self.losas.items():
            if not isinstance(datos, dict):
                continue
            # Tipo de losa: alivianada, macisa, casetonada, cruzada...
            tipo_losa = str(datos.get("tipologia", "alivianada")).lower()
            tipos_losa[str(nombre)] = tipo_losa
            panel = datos.get("panel_ejes") or {}
            rect = None
            if panel:
                # Losa delimitada por ejes: pertenece a un solo nivel de apoyo.
                if not self._en_nivel(panel.get("cota_apoyo_m"), panel.get("nivel_apoyo_nombre")):
                    continue
                rect = self._rect_panel(panel)
            else:
                rect = self._rect_legado(datos)
            if rect is None:
                # Sin datos de ejes no se dibuja, pero se avisa para que no
                # desaparezca en silencio.
                etiquetas.append(str(nombre))
                continue
            x0, y0, x1, y1 = rect
            puntos.extend([(x0, y0), (x1, y1)])
            panel_datos = datos.get("panel_ejes") or {}
            direccion_luz = str(panel_datos.get("direccion_luz", "x")).lower()
            if direccion_luz not in ("x", "y"):
                direccion_luz = "x"
            losas.append((str(nombre), rect, tipo_losa))
            # Nombre de la losa al centro del paño, con tipo y apoyos.
            apoyos_txt = []
            for lado in ("izq", "der"):
                apoyo = self._apoyo_losa(datos, panel_datos, rect, lado)
                if apoyo is None:
                    continue
                ax, ay, texto = apoyo
                apoyos_txt.append(texto)
                apoyos_losa.append({
                    "x": ax,
                    "y": ay,
                    "texto": texto,
                    "lado": lado,
                    "direccion_luz": direccion_luz,
                })
            detalle_apoyos = (
                f" → {' + '.join(apoyos_txt)}" if apoyos_txt else ""
            )
            rotulos.append(self._rotulo(
                x=(x0 + x1) / 2,
                y=(y0 + y1) / 2,
                texto=f"{nombre} ({tipo_losa}) carga:{direccion_luz.upper()}{detalle_apoyos}",
                direccion=direccion_luz,
                color="#047857",
                prioridad=4,
                # Rótulo técnico: no se dibuja sobre el plano.
                flotante=True,
            ))
            # Marcas de apoyo en los bordes del paño.
            conteo["losas"] += 1

        return {
            "porticos": porticos,
            "muros": muros,
            "cargas_muro": cargas_muro,
            "losas": losas,
            "etiquetas": etiquetas,
            "rotulos": rotulos,
            "apoyos_losa": apoyos_losa,
            "proyecciones": proyecciones,
            "puntos": puntos,
            "conteo": conteo,
            "tipos_losa": tipos_losa,
            "columnas_unicas": columnas_unicas,
            "bases_unicas": bases_unicas,
        }

    def _encuadre(self, escena):
        """Devuelve el encuadre de la planta: (px, py, xmin, xmax, ...).

        ``px`` y ``py`` son las funciones que convierten coordenadas de planta
        (en metros) a píxeles del widget. Se separan del ``paintEvent`` para
        poder reutilizar la misma transformación al colocar los rótulos.
        """
        xs = [p[0] for p in escena["puntos"]]
        ys = [p[1] for p in escena["puntos"]]
        if not xs or not ys:
            for familia, items in self.ejes.items():
                for item in items:
                    c = self._numero(item.get("coordenada_m"))
                    if c is None:
                        continue
                    (xs if familia == "x" else ys).append(c)
        if not xs or not ys:
            xs = [0.0, 10.0]
            ys = [0.0, 10.0]

        xmin, xmax = min(xs), max(xs)
        ymin, ymax = min(ys), max(ys)
        aire = max(xmax - xmin, ymax - ymin, 5.0) * 0.05 + 0.5
        xmin -= aire
        xmax += aire
        ymin -= aire
        ymax += aire
        # Margenes para las burbujas de los ejes (van fuera del área de obra)
        # y para la leyenda y el pie de texto.
        margen_x, margen_y = 78.0, 62.0
        rango_x = max(xmax - xmin, 3.0)
        rango_y = max(ymax - ymin, 3.0)
        # Zoom acotado para que la planta ocupe mejor el área disponible.
        # Factor reducido para que la estructura respire y no quede amontanada.
        escala = min(
            max(self.width() - 2 * margen_x, 1.0) / rango_x,
            max(self.height() - 2 * margen_y, 1.0) / rango_y,
        ) * 0.92 * self.zoom_factor
        # El desplazamiento (paneo por arrastre) se aplica al origen: mover la
        # vista no cambia la escala ni el encuadre de la obra.
        ox = (self.width() - rango_x * escala) / 2 + self.desplazamiento_x
        oy = (self.height() - rango_y * escala) / 2 + self.desplazamiento_y

        def px(x):
            return ox + (x - xmin) * escala

        def py(y):
            return oy + (ymax - y) * escala

        return px, py, xmin, xmax, ymin, ymax, escala

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.fillRect(self.rect(), QColor("#ffffff"))

        escena = self._escena()

        # --- encuadre: todos los elementos del nivel visible ---
        px, py, xmin, xmax, ymin, ymax, _escala = self._encuadre(escena)

        # --- base / terreno: recuadro semitransparente que abarca la planta ---
        # Se pinta primero, detrás de losas, muros y columnas, para dar
        # contexto de superficie enterrada sin ocultar los elementos.
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(QColor(200, 200, 210, 16))
        painter.drawRect(QRectF(QPointF(px(xmin), py(ymax)), QPointF(px(xmax), py(ymin))))

        margen = 10.0

        # --- grilla de ejes: fina, punteada y muy clara (solo referencia) ---
        painter.setPen(QPen(QColor(self.COLOR_EJE), 1, Qt.PenStyle.DashLine))
        fuente_eje = QFont(painter.font().family(), 10, QFont.Weight.Bold)
        painter.setFont(fuente_eje)
        metricas_eje = QFontMetricsF(fuente_eje)

        def _burbuja(ancla, texto, offset=QPointF(0.0, 0.0)):
            """Dibuja la burbuja de un eje con radio según el texto.

            ``ancla`` es el punto del eje y ``offset`` separa la burbuja del
            perímetro exterior de la obra (siempre hacia afuera, nunca sobre
            el plano). Se traza una guía fina que une la burbuja con su eje.
            """
            ancho = metricas_eje.horizontalAdvance(texto) + 10.0
            radio = max(13.0, min(ancho / 2 + 2.0, 30.0))
            centro = QPointF(ancla.x() + offset.x(), ancla.y() + offset.y())
            # Guía del eje hacia la burbuja.
            painter.setPen(QPen(QColor(self.COLOR_EJE_BORDE), 1, Qt.PenStyle.DotLine))
            painter.drawLine(ancla, centro)
            painter.setPen(QPen(QColor(self.COLOR_EJE_BORDE), 1.5))
            painter.setBrush(QColor(self.COLOR_EJE_BURBUJA))
            painter.drawEllipse(centro, radio, radio)
            painter.setPen(QColor("#1e293b"))
            painter.drawText(
                QRectF(centro.x() - radio, centro.y() - radio, radio * 2, radio * 2),
                Qt.AlignmentFlag.AlignCenter, texto,
            )
            return radio

        # Burbujas de ejes X: por debajo del perímetro exterior. Si dos ejes
        # quedan demasiado juntos, la segunda baja otra fila (más afuera).
        ultimo_x = None
        fila_x = 0
        for item in self.ejes.get("x", []):
            c = self._numero(item.get("coordenada_m"))
            if c is None:
                continue
            painter.setPen(QPen(QColor(self.COLOR_EJE), 1, Qt.PenStyle.DashLine))
            painter.drawLine(QPointF(px(c), py(ymin)), QPointF(px(c), py(ymax)))
            nombre = str(item.get("nombre", item.get("id", "?")))
            x_px = px(c)
            if ultimo_x is not None and x_px - ultimo_x < 34.0:
                fila_x = 1 - fila_x
            else:
                fila_x = 0
            ultimo_x = x_px
            _burbuja(
                QPointF(x_px, py(ymin)),
                nombre,
                QPointF(0.0, 26.0 + 26.0 * fila_x),
            )

        # Burbujas de ejes Y: a la izquierda del perímetro exterior, con el
        # mismo escalonado hacia afuera.
        ultimo_y = None
        fila_y = 0
        for item in self.ejes.get("y", []):
            c = self._numero(item.get("coordenada_m"))
            if c is None:
                continue
            painter.setPen(QPen(QColor(self.COLOR_EJE), 1, Qt.PenStyle.DashLine))
            painter.drawLine(QPointF(px(xmin), py(c)), QPointF(px(xmax), py(c)))
            nombre = str(item.get("nombre", item.get("id", "?")))
            y_px = py(c)
            if ultimo_y is not None and abs(y_px - ultimo_y) < 34.0:
                fila_y = 1 - fila_y
            else:
                fila_y = 0
            ultimo_y = y_px
            _burbuja(
                QPointF(px(xmin), y_px),
                nombre,
                QPointF(-30.0 - 26.0 * fila_y, 0.0),
            )

        # --- losas primero (fondo): color según tipo de losa ---
        color_por_tipo = {
            "alivianada": QColor("#059669"),
            "maciza": QColor("#0f766e"),
            "casetonada": QColor("#0891b2"),
            "cruzada": QColor("#0d9488"),
        }
        for nombre, (x0, y0, x1, y1), tipo_losa in escena["losas"]:
            color = color_por_tipo.get(tipo_losa, QColor("#059669"))
            painter.setPen(QPen(color, 1.5, Qt.PenStyle.SolidLine))
            painter.setBrush(QColor(color.red(), color.green(), color.blue(), 28))
            painter.drawRect(QRectF(QPointF(px(x0), py(y1)), QPointF(px(x1), py(y0))))
            # --- flecha doble tipo cota en el centro: sentido de carga de la losa ---
            datos_losa = self.losas.get(nombre, {})
            panel_datos = datos_losa.get("panel_ejes") or {}
            direccion_luz = str(panel_datos.get("direccion_luz", "x")).lower()
            if direccion_luz not in ("x", "y"):
                direccion_luz = "x"
            cx = (px(x0) + px(x1)) / 2
            cy = (py(y0) + py(y1)) / 2
            largo = 12
            punt = 4
            painter.setPen(QPen(QColor("#047857"), 1.5))
            painter.setBrush(QColor("#047857"))
            if direccion_luz == "x":
                # Flecha horizontal doble (carga en sentido X)
                painter.drawLine(QPointF(cx - largo, cy), QPointF(cx + largo, cy))
                painter.drawPolygon(QPolygonF([
                    QPointF(cx - largo, cy),
                    QPointF(cx - largo + punt, cy - punt / 2),
                    QPointF(cx - largo + punt, cy + punt / 2),
                ]))
                painter.drawPolygon(QPolygonF([
                    QPointF(cx + largo, cy),
                    QPointF(cx + largo - punt, cy - punt / 2),
                    QPointF(cx + largo - punt, cy + punt / 2),
                ]))
            else:
                # Flecha vertical doble (carga en sentido Y)
                painter.drawLine(QPointF(cx, cy - largo), QPointF(cx, cy + largo))
                painter.drawPolygon(QPolygonF([
                    QPointF(cx, cy - largo),
                    QPointF(cx - punt / 2, cy - largo + punt),
                    QPointF(cx + punt / 2, cy - largo + punt),
                ]))
                painter.drawPolygon(QPolygonF([
                    QPointF(cx, cy + largo),
                    QPointF(cx - punt / 2, cy + largo - punt),
                    QPointF(cx + punt / 2, cy + largo - punt),
                ]))

        # --- muros: banda roja con borde oscuro (se ve sobre la losa) ---
        for muro in escena["muros"]:
            espesor = muro["espesor"]
            if muro["direccion"] == "x":
                p1 = QPointF(px(muro["x0"]), py(muro["posicion"] - espesor / 2))
                p2 = QPointF(px(muro["x1"]), py(muro["posicion"] - espesor / 2))
                p3 = QPointF(px(muro["x1"]), py(muro["posicion"] + espesor / 2))
                p4 = QPointF(px(muro["x0"]), py(muro["posicion"] + espesor / 2))
            else:
                p1 = QPointF(px(muro["posicion"] - espesor / 2), py(muro["x0"]))
                p2 = QPointF(px(muro["posicion"] + espesor / 2), py(muro["x0"]))
                p3 = QPointF(px(muro["posicion"] + espesor / 2), py(muro["x1"]))
                p4 = QPointF(px(muro["posicion"] - espesor / 2), py(muro["x1"]))
            painter.setPen(QPen(QColor(self.COLOR_MURO_BORDE), 2))
            color_muro = QColor(self.COLOR_MURO)
            color_muro.setAlpha(120)
            painter.setBrush(color_muro)
            painter.drawPolygon(QPolygonF([p1, p2, p3, p4]))

        # --- flechas de muro perpendicular entre pórticos: dobles tipo cota ---
        # Cada flecha conecta los dos puntos de apoyo (un pórtico cada uno) y
        # lleva puntas hacia afuera, indicando la dirección transversal del
        # muro perpendicular que los carga.
        for carga in escena["cargas_muro"]:
            p1, p2 = carga["p1"], carga["p2"]
            x_inicio, y_inicio = px(p1[0]), py(p1[1])
            x_fin, y_fin = px(p2[0]), py(p2[1])
            dx = x_fin - x_inicio
            dy = y_fin - y_inicio
            distancia = (dx * dx + dy * dy) ** 0.5
            if distancia < 4:
                continue
            ux = dx / distancia
            uy = dy / distancia
            # Reforma los extremos para que las puntas no tapen los ejes.
            margen = 6.0
            bx0 = x_inicio + ux * margen
            by0 = y_inicio + uy * margen
            bx1 = x_fin - ux * margen
            by1 = y_fin - uy * margen
            color = QColor("#0ea5e9")
            color.setAlpha(220)
            painter.setPen(QPen(color, 2.5))
            painter.setBrush(color)
            # Barra central.
            painter.drawLine(QPointF(bx0, by0), QPointF(bx1, by1))
            # Punta de flecha en cada extremo, apuntando hacia afuera.
            punt = 6.0
            for cx, cy, vx, vy in (
                (bx1, by1, ux, uy),
                (bx0, by0, -ux, -uy),
            ):
                base_x = cx - punt * vx
                base_y = cy - punt * vy
                px_n = -vy * (punt / 2)
                py_n = vx * (punt / 2)
                painter.drawPolygon(QPolygonF([
                    QPointF(cx, cy),
                    QPointF(base_x + px_n, base_y + py_n),
                    QPointF(base_x - px_n, base_y - py_n),
                ]))

        # --- pórticos: eje fino + elementos gruesos del nivel ---
        fuente_el = QFont(painter.font().family(), 9)
        painter.setFont(fuente_el)
        for portico in escena["porticos"]:
            activo = portico["activo"]
            color = QColor(self.COLOR_PORTICO_ACTIVO if activo else self.COLOR_PORTICO)
            p0 = portico["eje_p0"]
            p1 = portico["eje_p1"]

            # Eje del pórtico: línea fina punteada (referencia, no elemento).
            painter.setPen(QPen(color, 1, Qt.PenStyle.DotLine))
            painter.drawLine(QPointF(px(p0[0]), py(p0[1])), QPointF(px(p1[0]), py(p1[1])))

            # Vigas del nivel (banda gruesa, color propio).
            for vid, a, b in portico["vigas"]:
                painter.setPen(QPen(QColor(self.COLOR_VIGA), 6 if activo else 5))
                painter.drawLine(QPointF(px(a[0]), py(a[1])), QPointF(px(b[0]), py(b[1])))
                painter.setPen(QPen(QColor(self.COLOR_VIGA_BORDE), 1, Qt.PenStyle.SolidLine))
                painter.drawLine(QPointF(px(a[0]), py(a[1])), QPointF(px(b[0]), py(b[1])))

            # Bases: se dibujan una sola vez por coordenada global (bases_unicas),
            # más abajo, para no duplicar las que comparten el mismo punto físico.


        # --- columnas: una sola vez por coordenada global (nodos compartidos en cruces) ---
        # Color: si alguno de los pórticos que contienen la columna está activo, se usa el color activo.
        # El nombre de la columna NO se dibuja aquí: lo aporta el sistema de
        # rótulos, que evita superposiciones entre textos.
        porticos_por_nombre = {p["nombre"]: p for p in escena["porticos"]}
        for col in escena.get("columnas_unicas", []):
            punto = col["punto"]
            cx, cy = px(punto[0]), py(punto[1])
            es_activa = any(
                porticos_por_nombre.get(pn, {}).get("activo", False)
                for pn in col.get("porticos", [])
            )
            col_color = QColor(self.COLOR_COLUMNA if es_activa else self.COLOR_COLUMNA_INACTIVA)
            painter.setPen(QPen(QColor("#0f172a"), 1))
            painter.setBrush(col_color)
            painter.drawRect(cx - 6, cy - 6, 12, 12)

        # --- bases: una sola vez por coordenada global (el nombre lo dibuja
        # el sistema de rótulos) ---
        for base in escena.get("bases_unicas", []):
            punto = base["punto"]
            cx, cy = px(punto[0]), py(punto[1])
            painter.setPen(QPen(QColor("#4c1d95"), 1.5))
            painter.setBrush(QColor(self.COLOR_BASE))
            painter.drawPolygon(QPolygonF([
                QPointF(cx, cy + 6),
                QPointF(cx - 7, cy + 18),
                QPointF(cx + 7, cy + 18),
            ]))

        # --- proyecciones: lo que es de otro nivel va punteado y apagado ---
        trazo_proy = QPen(QColor("#94a3b8"), 1, Qt.PenStyle.DashLine)
        for proy in escena.get("proyecciones", []):
            if proy.get("tipo") == "viga":
                a, b = proy["p0"], proy["p1"]
                painter.setPen(QPen(QColor("#d6a05c"), 2, Qt.PenStyle.DashLine))
                painter.drawLine(QPointF(px(a[0]), py(a[1])), QPointF(px(b[0]), py(b[1])))
            else:
                espesor = 0.15
                if proy["direccion"] == "x":
                    p1 = QPointF(px(proy["x0"]), py(proy["posicion"] - espesor / 2))
                    p2 = QPointF(px(proy["x1"]), py(proy["posicion"] - espesor / 2))
                    p3 = QPointF(px(proy["x1"]), py(proy["posicion"] + espesor / 2))
                    p4 = QPointF(px(proy["x0"]), py(proy["posicion"] + espesor / 2))
                else:
                    p1 = QPointF(px(proy["posicion"] - espesor / 2), py(proy["x0"]))
                    p2 = QPointF(px(proy["posicion"] + espesor / 2), py(proy["x0"]))
                    p3 = QPointF(px(proy["posicion"] + espesor / 2), py(proy["x1"]))
                    p4 = QPointF(px(proy["posicion"] - espesor / 2), py(proy["x1"]))
                painter.setPen(trazo_proy)
                painter.setBrush(Qt.BrushStyle.NoBrush)
                painter.drawPolygon(QPolygonF([p1, p2, p3, p4]))


        # --- rótulos orientados según el eje, sin superposición ---
        self._pintar_rotulos(painter, escena["rotulos"], px, py)

        # --- textos: nivel, leyenda y conteo de elementos del nivel ---
        conteo = escena["conteo"]
        if self.nivel_activo is not None:
            nombre_nivel = "?"
            for nivel in self.niveles:
                cota = self._numero(nivel.get("cota_m"))
                if cota is not None and abs(cota - self.nivel_activo) <= self.TOL:
                    nombre_nivel = str(nivel.get("nombre", "?"))
                    break
            titulo = f"Nivel visible: {nombre_nivel} (cota {self.nivel_activo:.2f} m)"
        else:
            titulo = "Todos los niveles"
        if self.portico_activo:
            titulo += f"   ·   Pórtico activo: {self.portico_activo}"
        painter.setPen(QPen(QColor("#7c3aed"), 1))
        painter.drawText(QPointF(margen, 22), titulo)

        # --- leyenda visual: muestra de color + símbolo por tipo de elemento ---
        tipos_presentes = sorted({t for t in escena["tipos_losa"].values()})
        leyenda_tipos = " · ".join(t.title() for t in tipos_presentes) if tipos_presentes else "sin losas"

        # La leyenda se dibuja como un panel en la esquina inferior derecha:
        # las burbujas de ejes ocupan los márgenes superior e izquierdo, y el
        # pie de texto ocupa la esquina inferior izquierda.
        entradas_leyenda = [
            (self.COLOR_COLUMNA, "cuadrado", "Columna activa (nodo de planta)"),
            (self.COLOR_COLUMNA_INACTIVA, "cuadrado", "Columna de pórtico inactivo"),
            (self.COLOR_VIGA, "linea", "Viga del nivel visible"),
            (self.COLOR_BASE, "triangulo", "Base / cimiento"),
            (self.COLOR_BASE, "triangulo", "Base compartida: bid (+N) = N bases/columnas"),
            ("#059669", "cuadrado", f"Losa ({leyenda_tipos})"),
            (self.COLOR_MURO, "linea", "Muro"),
            (self.COLOR_APOYO, "circulo", "Apoyo de losa"),
            ("#94a3b8", "linea", "Proyección de otro nivel (punteado)"),
            (self.COLOR_PORTICO, "linea", "Eje de pórtico inactivo (punteado fino)"),
            (self.COLOR_PORTICO_ACTIVO, "linea", "Eje de pórtico activo"),
            (self.COLOR_EJE, "linea", "Ejes de la obra (punteado claro)"),
        ]

        _tam_leyenda = 9

        def _muestra_leyenda(y, color, simbolo, texto, x_texto=margen + 20):
            painter.save()
            painter.setPen(QPen(QColor(color), 2))
            painter.setBrush(QColor(QColor(color).red(), QColor(color).green(), QColor(color).blue(), 40))
            if simbolo == "cuadrado":
                painter.drawRect(QRectF(QPointF(x_texto - 20, y - 8), QPointF(x_texto - 6, y + 6)))
            elif simbolo == "linea":
                painter.drawLine(QPointF(x_texto - 18, y - 1), QPointF(x_texto - 8, y - 1))
            elif simbolo == "triangulo":
                painter.drawPolygon(QPolygonF([
                    QPointF(x_texto - 13, y - 6),
                    QPointF(x_texto - 18, y + 5),
                    QPointF(x_texto - 8, y + 5),
                ]))
            elif simbolo == "circulo":
                painter.setBrush(Qt.BrushStyle.NoBrush)
                painter.drawEllipse(QPointF(x_texto - 13, y - 1), 5, 5)
                painter.setBrush(QColor(QColor(color).red(), QColor(color).green(), QColor(color).blue(), 40))
            painter.setPen(QColor("#334155"))
            painter.setFont(QFont(painter.font().family(), _tam_leyenda))
            painter.drawText(QPointF(x_texto, y + 3), texto)
            painter.restore()

        # --- panel de leyenda (esquina inferior derecha) ---
        # Se busca el tamaño de fuente más grande con el que la leyenda entre
        # en el ancho disponible; si no entra de ninguna forma, se dibuja
        # arriba a la derecha (zona libre de burbujas de ejes).
        _ancho_disponible = self.width() - 2 * margen - 24.0
        _fuente_leyenda = QFont(painter.font().family(), 9)
        _metricas_leyenda = QFontMetricsF(_fuente_leyenda)
        _ancho_leyenda = 0
        for _color, _sim, _texto in entradas_leyenda:
            _ancho_leyenda = max(_ancho_leyenda, _metricas_leyenda.horizontalAdvance(_texto))
        for _tam in (9, 8, 7):
            _fuente_leyenda = QFont(painter.font().family(), _tam)
            _metricas_leyenda = QFontMetricsF(_fuente_leyenda)
            _ancho_leyenda = 0
            for _color, _sim, _texto in entradas_leyenda:
                _ancho_leyenda = max(_ancho_leyenda, _metricas_leyenda.horizontalAdvance(_texto))
            _tam_leyenda = _tam
            if _ancho_leyenda + 46 <= _ancho_disponible:
                break
        fila = 15 if _tam_leyenda <= 7 else 16
        _ancho_panel = _ancho_leyenda + 46
        _alto_panel = fila * len(entradas_leyenda) + 26
        if _ancho_panel <= _ancho_disponible:
            # Abajo a la derecha, con margen respecto de los bordes.
            _x_panel = self.width() - _ancho_panel - 12.0
            _y_panel = max(self.height() - _alto_panel - 10.0, 40.0)
        else:
            # No entra en ancho: se dibuja arriba a la derecha, bajo el título.
            _x_panel = max(self.width() - _ancho_panel - 12.0, margen)
            _y_panel = 40.0

        painter.save()
        painter.setPen(QPen(QColor("#cbd5e1"), 1))
        painter.setBrush(QColor(255, 255, 255, 235))
        painter.drawRoundedRect(
            QRectF(QPointF(_x_panel, _y_panel),
                   QPointF(_x_panel + _ancho_panel, _y_panel + _alto_panel)),
            6, 6,
        )
        painter.setPen(QPen(QColor("#475569"), 1))
        painter.setFont(QFont(painter.font().family(), _tam, QFont.Weight.Bold))
        painter.drawText(QPointF(_x_panel + 12, _y_panel + 14), "Leyenda:")
        _y = _y_panel + 14 + fila
        for _color, _sim, _texto in entradas_leyenda:
            _muestra_leyenda(_y, _color, _sim, _texto, x_texto=_x_panel + 32)
            _y += fila
        painter.restore()

        painter.setPen(QPen(QColor("#334155"), 1))
        painter.drawText(
            QPointF(margen, self.height() - 30),
            f"Columnas: {conteo['columnas']}  ·  Vigas: {conteo['vigas']}  ·  "
            f"Bases: {conteo['bases']}  ·  Losas: {conteo['losas']}  ·  "
            f"Muros: {conteo['muros']}",
        )
        nota = "Solo se dibujan los elementos del nivel elegido"
        if escena["etiquetas"]:
            nota += "  ·  Sin posición en planta: " + ", ".join(escena["etiquetas"][:4])
            if len(escena["etiquetas"]) > 4:
                nota += f" (+{len(escena['etiquetas']) - 4})"
        painter.setPen(QPen(QColor("#94a3b8"), 1))
        painter.drawText(QPointF(margen, self.height() - 12), nota)

        # --- flecha de dirección del viento en planta (esquina superior derecha) ---
        if self.viento_lineal_kN_m > 0 or self.viento_nodos:
            color_viento = QColor("#d97706")
            painter.setPen(QPen(color_viento, 2.5))
            painter.setBrush(color_viento)
            x_fin = self.width() - 40.0
            y_fin = 44.0
            largo = 34.0
            ancho = 14.0
            if self.direccion_viento == "y":
                # Flecha vertical: viento sobre ejes Y (pórticos planos X al frente).
                painter.drawPolygon(QPolygonF([
                    QPointF(x_fin, y_fin),
                    QPointF(x_fin - ancho / 2, y_fin + largo),
                    QPointF(x_fin + ancho / 2, y_fin + largo),
                ]))
                painter.setPen(QColor("#9a3412"))
                painter.drawText(
                    QPointF(x_fin - 78, y_fin + largo + 12),
                    f"Viento → dir. Y  ({self.viento_lineal_kN_m:.2f} kN/m · W total: {sum(n['valor_kN'] for n in self.viento_nodos):.2f} kN)",
                )
            else:
                # Flecha horizontal: viento sobre ejes X (pórticos planos Y al frente).
                painter.drawPolygon(QPolygonF([
                    QPointF(x_fin - largo, y_fin),
                    QPointF(x_fin, y_fin - ancho / 2),
                    QPointF(x_fin, y_fin + ancho / 2),
                ]))
                painter.setPen(QColor("#9a3412"))
                painter.drawText(
                    QPointF(x_fin + 10, y_fin - 4),
                    f"Viento → dir. X  ({self.viento_lineal_kN_m:.2f} kN/m · W total: {sum(n['valor_kN'] for n in self.viento_nodos):.2f} kN)",
                )
        else:
            # Indicador explícito cuando el viento NO aplica en este pórtico/nivel.
            painter.setPen(QPen(QColor("#94a3b8"), 1.5))
            painter.drawText(
                QPointF(self.width() - 150, 44),
                "Viento: no aplica",
            )

        painter.end()

    def wheelEvent(self, event):
        """Zoom con la rueda del mouse: adelante acerca, atrás aleja."""
        delta = event.angleDelta().y()
        if delta > 0:
            self.zoom_factor = min(self.zoom_factor * 1.15, 4.0)
        else:
            self.zoom_factor = max(self.zoom_factor / 1.15, 0.25)
        self.update()









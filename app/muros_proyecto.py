"""Pantalla «Muros»: cálculo de los muros ESTRUCTURALES (estructura mixta).

Solo lista los muros creados en Inicio sobre un pórtico (los que en
``estructura.json`` tienen ``tipo == "muro"``). Los muros del catálogo de Cargas
NO son estructurales y no aparecen acá.

Permite componer cada muro en capas (revoque, hoja de ladrillo, cámara/aislación,
revoque interior) y calcula su peso propio con los pesos específicos de
``materiales.json``. La resistencia térmica, el diseño estructural y el cimiento
quedan marcados como pendientes hasta cargar conductividades λ y la
reglamentación de mampostería.
"""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
    QAbstractSpinBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QMessageBox,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from calc import cargas, mamposteria, muros, rutas
from app.visualizacion import VistaPortico


class EditorCapasMuro(QDialog):
    """Compone un muro en capas y calcula su peso propio en vivo."""

    def __init__(self, nombre: str, capas: list[dict], altura_m: float,
                 longitud_m: float, parent=None, altura_de_niveles: bool = False):
        super().__init__(parent)
        self.setWindowTitle(f"Composición del muro — {nombre}")
        self.nombre = nombre
        self._grupos = muros.grupos_disponibles()  # {grupo: {clave: {...}}}
        self._ultimo: dict = {}

        self.tabla = QTableWidget(0, 3)
        self.tabla.setHorizontalHeaderLabels(("Capa (material)", "Espesor (cm)", "γ (kN/m³)"))
        self.tabla.verticalHeader().setVisible(False)
        self.tabla.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)

        self.boton_agregar = QPushButton("Agregar capa")
        self.boton_quitar = QPushButton("Quitar capa")
        self.boton_agregar.clicked.connect(lambda: self._agregar_capa())
        self.boton_quitar.clicked.connect(self._quitar_capa)

        self.altura = QDoubleSpinBox()
        self.altura.setSuffix(" m")
        self.altura.setDecimals(3)
        self.altura.setRange(0.1, 30.0)
        self.altura.setValue(max(0.1, altura_m))
        # La altura sale de los niveles (nivel de arranque → nivel final); se
        # muestra en solo-lectura para que no se pise a mano por error.
        self.altura_de_niveles = altura_de_niveles
        self.altura.setReadOnly(altura_de_niveles)
        self.altura.setButtonSymbols(
            QAbstractSpinBox.ButtonSymbols.NoButtons if altura_de_niveles
            else QAbstractSpinBox.ButtonSymbols.UpDownArrows
        )
        etiqueta_altura = "Altura del paño (de los niveles):" if altura_de_niveles else "Altura del paño:"
        self.longitud = QDoubleSpinBox()
        self.longitud.setSuffix(" m")
        self.longitud.setDecimals(3)
        self.longitud.setRange(0.1, 200.0)
        self.longitud.setValue(max(0.1, longitud_m))

        self.resultado = QTextEdit()
        self.resultado.setReadOnly(True)
        self.resultado.setMinimumHeight(120)

        botones = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        botones.button(QDialogButtonBox.StandardButton.Ok).setText("Guardar composición")
        botones.accepted.connect(self.accept)
        botones.rejected.connect(self.reject)

        fila = QHBoxLayout()
        fila.addWidget(self.boton_agregar)
        fila.addWidget(self.boton_quitar)
        fila.addStretch(1)

        forma = QFormLayout()
        forma.addRow(etiqueta_altura, self.altura)
        forma.addRow("Longitud del paño:", self.longitud)

        caja = QVBoxLayout(self)
        caja.addWidget(QLabel("Componé el muro de afuera hacia adentro (revoque → hoja → aislación → hoja → revoque)."))
        caja.addWidget(self.tabla, 1)
        caja.addLayout(fila)
        caja.addLayout(forma)
        caja.addWidget(QLabel("Resultado del peso propio:"))
        caja.addWidget(self.resultado)
        caja.addWidget(botones)
        self.resize(560, 620)

        for capa in capas:
            self._agregar_capa(capa)
        self.altura.valueChanged.connect(self._recalcular)
        self.longitud.valueChanged.connect(self._recalcular)
        self.tabla.itemChanged.connect(self._recalcular)
        self._recalcular()

    def _agregar_capa(self, capa: dict | None = None) -> None:
        self.tabla.blockSignals(True)
        fila = self.tabla.rowCount()
        self.tabla.insertRow(fila)

        combo = QComboBox()
        for grupo, tabla in self._grupos.items():
            for clave, info in tabla.items():
                grupo_visible = "Morteros / revoques" if grupo == "Morteros" else grupo
                combo.addItem(f"{info['nombre']} ({grupo_visible})", (grupo, clave))
        if capa:
            # findData con una tupla no matchea de forma confiable en PySide6,
            # así que buscamos el índice comparando el dato guardado a mano.
            objetivo = (str(capa.get("grupo")), str(capa.get("clave")))
            for i in range(combo.count()):
                if combo.itemData(i) == objetivo:
                    combo.setCurrentIndex(i)
                    break
        combo.currentIndexChanged.connect(self._recalcular)
        self.tabla.setCellWidget(fila, 0, combo)

        # Espesor: se carga en CENTÍMETROS (cómodo para mampostería) y se
        # convierte a metros al guardar. Rango 0,5 cm a 100 cm.
        espesor = QDoubleSpinBox()
        espesor.setSuffix(" cm")
        espesor.setDecimals(1)
        espesor.setRange(0.5, 100.0)
        espesor.setSingleStep(0.5)
        espesor_cm = float(capa.get("espesor_m", 0.02)) * 100.0 if capa else 2.0
        espesor.setValue(espesor_cm)
        espesor.valueChanged.connect(self._recalcular)
        self.tabla.setCellWidget(fila, 1, espesor)

        self.tabla.setItem(fila, 2, QTableWidgetItem(""))
        self.tabla.blockSignals(False)
        self._recalcular()

    def _quitar_capa(self) -> None:
        fila = self.tabla.currentRow()
        if fila >= 0:
            self.tabla.removeRow(fila)
            self._recalcular()

    def _capas_actuales(self) -> list[dict]:
        capas: list[dict] = []
        for fila in range(self.tabla.rowCount()):
            combo = self.tabla.cellWidget(fila, 0)
            espesor = self.tabla.cellWidget(fila, 1)
            if not isinstance(combo, QComboBox) or not isinstance(espesor, QDoubleSpinBox):
                continue
            datos = combo.currentData()
            if not datos:
                continue
            grupo, clave = datos
            # El espesor se muestra/carga en cm; el cálculo trabaja en metros.
            espesor_m = espesor.value() / 100.0
            capas.append({"grupo": grupo, "clave": clave, "espesor_m": espesor_m})
        return capas

    def _recalcular(self, *_args) -> None:
        # Bloqueamos las señales de la tabla mientras escribimos la columna de
        # γ: si no, cada setItem emite itemChanged, que vuelve a llamar acá y
        # desborda la pila (la app se cierra de golpe al tocar un material).
        self.tabla.blockSignals(True)
        try:
            for fila in range(self.tabla.rowCount()):
                combo = self.tabla.cellWidget(fila, 0)
                datos = combo.currentData() if isinstance(combo, QComboBox) else None
                texto = ""
                if datos:
                    try:
                        texto = f"{muros._gamma_de(datos[1], datos[0]):.2f}"
                    except ValueError:
                        texto = "—"
                self.tabla.setItem(fila, 2, QTableWidgetItem(texto))
        finally:
            self.tabla.blockSignals(False)

        capas = self._capas_actuales()
        try:
            resultado = muros.calcular_muro_compuesto(
                self.nombre, capas, self.altura.value(), self.longitud.value(),
            )
        except ValueError as exc:
            self.resultado.setPlainText(f"Revisá las capas: {exc}")
            self._ultimo = {}
            return
        self._ultimo = resultado
        lineas = [
            f"Espesor total: {resultado['espesor_total_m']:.3f} m",
            f"Peso propio: {resultado['peso_propio_kN_m2']:.3f} kN/m²",
            f"Carga lineal (por m de alto): {resultado['peso_lineal_kN_m']:.3f} kN/m",
            f"Peso total del paño ({resultado['longitud_m']:.2f} m × {resultado['altura_m']:.2f} m): "
            f"{resultado['peso_total_kN']:.2f} kN",
        ]
        termica = resultado.get("resistencia_termica")
        if termica:
            lineas += [
                "",
                "Resistencia térmica (ISO 6946, muro vertical; valor informativo):",
                f"  R de las capas: {termica['R_capas_m2K_W']:.3f} m²·K/W",
                f"  R total (con Rsi={termica['R_si_m2K_W']:.2f}, Rse={termica['R_se_m2K_W']:.2f}): "
                f"{termica['R_total_m2K_W']:.3f} m²·K/W",
                f"  Transmitancia térmica U = 1/Rtotal: {termica['U_W_m2K']:.3f} W/m²·K",
                "  Sin un límite U objetivo configurado, no corresponde marcarla como cumple/no cumple.",
            ]
        else:
            lineas += [
                "",
                "Resistencia térmica: pendiente (faltan conductividades λ de "
                "alguna capa en materiales.json).",
            ]
        # Verificación estructural (CIRSOC 501) con la altura del paño.
        try:
            carga_losas = muros.cargas_losas_en_muro(self.nombre)
            combinaciones_axiales = muros.combinaciones_axiales_basicas(
                resultado["peso_lineal_kN_m"] + carga_losas["D_kN_m"],
                carga_losas["L_kN_m"],
            )
            axial_servicio = combinaciones_axiales["D_kN_m"] + combinaciones_axiales["L_kN_m"]
            est = mamposteria.verificar_muro(
                capas, altura_m=self.altura.value(), mortero="mortero_normal_n",
                carga_axial_kN_m=combinaciones_axiales["Pu_kN_m"],
            )
        except (FileNotFoundError, KeyError):
            est = {"ok": False,
                   "pendiente": "Falta datos/mamposteria.json (biblioteca CIRSOC 501)."}
        if est.get("ok"):
            comp = est["compresion"]
            cort = est["corte"]
            zapata = muros.verificar_zapata_corrida(axial_servicio)
            lineas += [
                "",
                "Verificación estructural (CIRSOC 501) por metro de muro:",
                "  Se toma como portante la hoja de mampostería más gruesa; "
                "revoques y aislantes no aportan resistencia estructural.",
                "  El estado PASA/NO PASA mostrado abajo es solo por compresión axial "
                "con peso propio y losas; aún no incluye viento/sismo ni excentricidad.",
                f"  Hoja portante: {est['hoja_portante']['clave']} "
                f"(fm = {est['hoja_portante']['fm_mpa']:.2f} MPa, "
                f"e = {est['hoja_portante']['espesor_m']:.2f} m)",
                f"  Esbeltez h/r = {est['esbeltez_hr']:.1f}"
                f"{'  ⚠ ESBELTO' if est['es_esbelta'] else ''}",
                f"  Carga axial de servicio: propio={resultado['peso_lineal_kN_m']:.2f} + "
                f"losa D={carga_losas['D_kN_m']:.2f} + L={carga_losas['L_kN_m']:.2f} "
                f"= {axial_servicio:.2f} kN/m",
                f"  Combinaciones gravitatorias: "
                + ", ".join(f"{k}={v:.2f}" for k, v in combinaciones_axiales["combinaciones_kN_m"].items())
                + f" kN/m; gobierna {combinaciones_axiales['combinacion_critica']} "
                + f"→ Pu={combinaciones_axiales['Pu_kN_m']:.2f} kN/m",
                f"  Compresión: φ·Pn = {comp['phi_Pn_kN']:.1f} kN "
                f"(φ={comp['phi']:.2f}, {comp['rama']}) — "
                f"{'PASA' if comp['pasa'] else 'NO PASA'}",
                f"  Corte:      φ·Vn = {cort['phi_Vn_kN']:.1f} kN (φ={cort['phi']:.2f})",
                (f"  Zapata corrida preliminar: B ≥ {zapata['ancho_minimo_m']:.2f} m, "
                 f"q_adm={zapata['q_adm_kPa']:.1f} kPa, profundidad "
                 f"{zapata['profundidad_m']} m. {zapata['pendiente']}"
                 if zapata.get("ok") else f"  Zapata corrida: {zapata.get('pendiente')}"),
                "  Recomendaciones:",
            ]
            lineas += [f"    · {r}" for r in est["recomendaciones"]]
        else:
            capas_mamp = [c for c in capas if str(c.get("grupo", "")).lower() == "mamposteria"]
            if capas_mamp:
                portante = max(capas_mamp, key=lambda c: float(c.get("espesor_m", 0) or 0))
                lineas += [
                    "", f"Hoja considerada portante: {portante.get('clave')} "
                    f"({float(portante.get('espesor_m', 0)):.3f} m).",
                    "Revoques y aislaciones se consideran capas no portantes.",
                ]
            lineas += [
                "",
                "Verificación estructural: " + str(est.get("pendiente", "no disponible.")),
            ]
        self._ultimo["cargas_losas"] = carga_losas if "carga_losas" in locals() else None
        self._ultimo["combinaciones_axiales"] = (
            combinaciones_axiales if "combinaciones_axiales" in locals() else None
        )
        self._ultimo["diseno_estructural"] = est
        self._ultimo["cimiento"] = (
            muros.verificar_zapata_corrida(axial_servicio)
            if "axial_servicio" in locals() else None
        )
        if 'carga_losas' in locals() and carga_losas.get("losas"):
            lineas += ["", "Losas activas que descargan en este muro:"]
            lineas += [f"  {l['nombre']}: D={l['D_kN_m']:.2f}, L={l['L_kN_m']:.2f} kN/m"
                       for l in carga_losas["losas"]]
        if 'carga_losas' in locals() and carga_losas.get("avisos"):
            lineas += ["", "Avisos de reacciones de losas:"]
            lineas += [f"  {aviso}" for aviso in carga_losas["avisos"]]
        self.resultado.setPlainText("\n".join(lineas))

    def datos_guardar(self) -> tuple[list[dict], float, float, dict]:
        return (
            self._capas_actuales(),
            self.altura.value(),
            self.longitud.value(),
            self._ultimo,
        )

class PaginaMuros(QWidget):
    """Pestaña «Muros»: lista los muros estructurales y permite componerlos."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.intro = QLabel(
            "Muros ESTRUCTURALES de esta obra (los creados sobre un pórtico en "
            "Inicio, estructura mixta). Los muros del catálogo de Cargas no son "
            "estructurales y no aparecen acá. Componé cada muro en capas para "
            "calcular peso propio, reacciones de losas apoyadas, verificación axial "
            "CIRSOC 501 y predimensionado inicial de zapata corrida."
        )
        self.intro.setWordWrap(True)

        # --- selector de muro activo ---
        self.combo_muro_activo = QComboBox()
        self.combo_muro_activo.setMinimumWidth(180)
        self.combo_muro_activo.currentIndexChanged.connect(self._actualizar_visor_muro)
        self.lbl_muro_activo = QLabel("Muro activo:")

        # --- visor de muro (VistaPortico dibuja el paño y sus cargas) ---
        self.visor_muro = VistaPortico()
        self.visor_muro.setMinimumHeight(300)
        self.visor_muro.setMinimumWidth(400)

        self.resumen_calculo = QTextEdit()
        self.resumen_calculo.setReadOnly(True)
        self.resumen_calculo.setMaximumHeight(155)
        self.resumen_calculo.setPlaceholderText("Resumen guardado del cálculo del muro")

        self.tabla = QTableWidget(0, 11)
        self.tabla.setHorizontalHeaderLabels((
            "Muro", "Longitud (m)", "Espesor (m)", "Altura (m)",
            "Peso propio (kN/m)", "Composición", "U (W/m²·K)", "R (m²·K/W)",
            "Compresión axial", "Zapata corrida sugerida", "Refuerzos / notas",
        ))
        self.tabla.verticalHeader().setVisible(False)
        self.tabla.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tabla.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.tabla.setWordWrap(True)
        self.tabla.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        self.tabla.horizontalHeader().setSectionResizeMode(5, QHeaderView.ResizeMode.Stretch)
        self.tabla.horizontalHeader().setSectionResizeMode(10, QHeaderView.ResizeMode.Stretch)
        self.tabla.setColumnWidth(8, 165)
        self.tabla.setColumnWidth(9, 210)
        self.tabla.setColumnWidth(10, 270)
        self.tabla.itemSelectionChanged.connect(self._seleccionar_muro_desde_tabla)

        self.boton_componer = QPushButton("Componer muro / calcular peso propio…")
        self.boton_componer.clicked.connect(self._componer)
        self.boton_actualizar = QPushButton("Actualizar lista")
        self.boton_actualizar.clicked.connect(self.refrescar)
        self.boton_eliminar_muro = QPushButton("Eliminar muro…")
        self.boton_eliminar_muro.clicked.connect(self._eliminar_muro)

        fila = QHBoxLayout()
        fila.addWidget(self.lbl_muro_activo)
        fila.addWidget(self.combo_muro_activo, 1)
        fila.addStretch(1)
        fila.addWidget(self.boton_componer)
        fila.addWidget(self.boton_actualizar)
        fila.addWidget(self.boton_eliminar_muro)

        caja = QVBoxLayout(self)
        caja.addWidget(self.intro)
        caja.addLayout(fila)
        caja.addWidget(self.visor_muro, 1)
        caja.addWidget(QLabel("Resumen del último cálculo guardado:"))
        caja.addWidget(self.resumen_calculo)
        caja.addWidget(self.tabla, 1)
        self.refrescar()

    # ------------------------------------------------------------------
    def _eliminar_muro(self) -> None:
        nombre = self.combo_muro_activo.currentText().strip()
        if not nombre:
            QMessageBox.information(self, "Sin muro", "No hay muro seleccionado.")
            return
        respuesta = QMessageBox.question(
            self, "Eliminar muro", f"¿Eliminar «{nombre}» de la obra?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
        )
        if respuesta != QMessageBox.StandardButton.Yes:
            return
        estructura = rutas.cargar_estructura() or {}
        if nombre in estructura:
            del estructura[nombre]
            try:
                rutas.guardar_estructura(estructura)
            except OSError as exc:
                QMessageBox.critical(self, "No se pudo guardar", str(exc))
                return
        self.refrescar()
        parent = self.parent()
        if hasattr(parent, "refrescar"):
            try:
                parent.refrescar()
            except Exception:
                pass

    def _muros_estructurales(self) -> dict[str, dict]:
        estructura = rutas.cargar_estructura() or {}
        return {
            nombre: datos
            for nombre, datos in estructura.items()
            if isinstance(datos, dict) and str(datos.get("tipo", "")).lower() == "muro"
        }

    def refrescar(self) -> None:
        muro_activo_anterior = self.combo_muro_activo.currentText().strip()
        muros_est = self._muros_estructurales()
        self.tabla.setRowCount(len(muros_est))
        for fila, (nombre, datos) in enumerate(muros_est.items()):
            longitud = float(datos.get("longitud_m", 0.0))
            espesor = float(datos.get("espesor_m", 0.0))
            capas = datos.get("capas", []) or []
            peso_lineal = None
            altura = None
            if capas:
                try:
                    peso_m2 = muros.peso_propio_capas(capas)["peso_propio_kN_m2"]
                    altura = float(datos.get("altura_m", 0.0))
                    peso_lineal = peso_m2 * altura if altura else None
                except ValueError:
                    peso_lineal = None
            guardado = datos.get("resultado_calculo") or {}
            termica = guardado.get("resistencia_termica") or {}
            estructural = guardado.get("diseno_estructural") or {}
            comp = estructural.get("compresion") or {}
            if estructural.get("ok"):
                demanda = comp.get("carga_axial_kN_m")
                capacidad = comp.get("phi_Pn_kN")
                estado = "PASA" if comp.get("pasa") else "NO PASA"
                axial = f"{estado}\n{float(demanda or 0):.1f} / {float(capacidad or 0):.1f} kN/m"
            else:
                axial = "Sin verificar"
            cimiento = guardado.get("cimiento") or {}
            if cimiento.get("ok"):
                base = (
                    f"B ≥ {float(cimiento.get('ancho_minimo_m', 0)):.2f} m\n"
                    f"qadm {float(cimiento.get('q_adm_kPa', 0)):.1f} kPa · "
                    f"Df {cimiento.get('profundidad_m')} m"
                )
            else:
                base = "Sin predimensionar"
            capas_txt = "\n".join(
                f"{c.get('clave', '?')} · {float(c.get('espesor_m', 0)) * 100:.1f} cm"
                for c in capas
            ) or "Sin composición"
            recomendaciones = estructural.get("recomendaciones") or []
            notas_refuerzo = (
                "Orientaciones CIRSOC 501; detalle de armadura pendiente"
                if recomendaciones else "Sin recomendación guardada"
            )
            valores = (
                nombre, f"{longitud:.2f}", f"{espesor:.3f}",
                f"{altura:.2f}" if altura else "—",
                f"{peso_lineal:.3f}" if peso_lineal is not None else "—",
                capas_txt,
                f"{float(termica['U_W_m2K']):.3f}" if termica.get("U_W_m2K") is not None else "—",
                f"{float(termica['R_total_m2K_W']):.3f}" if termica.get("R_total_m2K_W") is not None else "—",
                axial, base, notas_refuerzo,
            )
            for col, valor in enumerate(valores):
                item = QTableWidgetItem(valor)
                if col == 10 and recomendaciones:
                    item.setToolTip("\n".join(recomendaciones))
                elif col == 9 and cimiento.get("ok"):
                    item.setToolTip(str(cimiento.get("pendiente", "")))
                elif col == 8 and estructural.get("ok"):
                    item.setToolTip(
                        "Verificación por compresión axial solamente. No incluye "
                        "viento/sismo ni excentricidad."
                    )
                elif col == 5:
                    item.setToolTip(capas_txt)
                self.tabla.setItem(fila, col, item)
        if not muros_est:
            self.tabla.setRowCount(1)
            self.tabla.setItem(
                0, 0,
                QTableWidgetItem("Todavía no hay muros estructurales en esta obra."),
            )
        self.tabla.resizeRowsToContents()
        # Llenar combo de muro activo
        self.combo_muro_activo.blockSignals(True)
        self.combo_muro_activo.clear()
        self.combo_muro_activo.addItems(list(muros_est.keys()))
        indice_activo = self.combo_muro_activo.findText(muro_activo_anterior)
        if indice_activo >= 0:
            self.combo_muro_activo.setCurrentIndex(indice_activo)
        self.combo_muro_activo.blockSignals(False)
        self._actualizar_visor_muro()

    def _actualizar_visor_muro(self) -> None:
        """Configura el visor para mostrar el muro activo con sus cargas."""
        nombre = self.combo_muro_activo.currentText().strip()
        if not nombre:
            return
        fila_activa = next(
            (
                fila for fila in range(self.tabla.rowCount())
                if self.tabla.item(fila, 0)
                and self.tabla.item(fila, 0).text() == nombre
            ),
            -1,
        )
        if fila_activa >= 0 and self.tabla.currentRow() != fila_activa:
            self.tabla.blockSignals(True)
            self.tabla.selectRow(fila_activa)
            self.tabla.blockSignals(False)
        muros_est = self._muros_estructurales()
        datos = muros_est.get(nombre)
        if not datos:
            return
        self.resumen_calculo.setPlainText(self._texto_resumen_calculo(nombre, datos))

        # Configurar el visor para el muro
        self.visor_muro.tipo = "muro"
        self.visor_muro.muro_nombre = nombre
        self.visor_muro.muros = {nombre: datos}
        self.visor_muro.columnas = {}
        self.visor_muro.vigas = {}
        self.visor_muro.tramos = []
        # Mantener las referencias de planta baja y niveles superiores en el
        # alzado del muro, igual que en el visor de pórticos.
        from app.ejes import cargar_niveles
        self.visor_muro.niveles = cargar_niveles()

        # Cargar cargas aplicadas al muro desde cargas.json
        self.visor_muro.cargas = self._cargar_cargas_muro(nombre)
        # Las reacciones de losas se transfieren por el vínculo geométrico
        # apoya_en y no necesariamente tienen una aplicación manual asociada.
        try:
            reacciones = muros.cargas_losas_en_muro(nombre)
            for indice, losa in enumerate(reacciones.get("losas", [])):
                for tipo, campo in (("D", "D_kN_m"), ("L", "L_kN_m")):
                    valor = float(losa.get(campo, 0.0))
                    if valor <= 0:
                        continue
                    self.visor_muro.cargas.setdefault("MURO", []).append({
                        "aplicacion_id": f"reaccion_losa_{indice}_{tipo}",
                        "carga_id": f"reaccion_losa_{indice}_{tipo}",
                        "nombre": f"{losa['nombre']} · reacción de losa",
                        "x_inicio_m": 0.0,
                        "x_fin_m": float(datos.get("longitud_m", 0.0)),
                        "valor_kN_m": valor,
                        "tipo": tipo,
                        "cota_m": losa.get("cota_m"),
                    })
        except (FileNotFoundError, KeyError, ValueError):
            pass
        self.visor_muro.capas_carga = max(
            (len({c.get("aplicacion_id", "") for c in cargas_tramo})
             for cargas_tramo in self.visor_muro.cargas.values()),
            default=0,
        )
        self.visor_muro.setMinimumHeight(max(250, 190 + self.visor_muro.capas_carga * 24))
        self.visor_muro.update()

    def _seleccionar_muro_desde_tabla(self) -> None:
        fila = self.tabla.currentRow()
        item = self.tabla.item(fila, 0) if fila >= 0 else None
        if item is None:
            return
        indice = self.combo_muro_activo.findText(item.text())
        if indice >= 0 and indice != self.combo_muro_activo.currentIndex():
            self.combo_muro_activo.blockSignals(True)
            self.combo_muro_activo.setCurrentIndex(indice)
            self.combo_muro_activo.blockSignals(False)
        self._actualizar_visor_muro()

    @staticmethod
    def _texto_resumen_calculo(nombre: str, datos: dict) -> str:
        """Resume en la pestaña el último cálculo que se guardó con el muro."""
        resultado = datos.get("resultado_calculo") or {}
        if not resultado:
            return (
                f"{nombre}: todavía no hay un resultado guardado. "
                "Abrí ‘Componer muro / calcular peso propio…’ y guardá la composición."
            )
        lineas = []
        termica = resultado.get("resistencia_termica") or {}
        if termica:
            lineas.append(
                f"Térmica: U = {float(termica.get('U_W_m2K', 0)):.3f} W/m²·K "
                "(informativo; sin límite U para dictamen)."
            )
        losas = resultado.get("cargas_losas") or {}
        lineas.append(
            f"Cargas de losas: D = {float(losas.get('D_kN_m', 0)):.2f}, "
            f"L = {float(losas.get('L_kN_m', 0)):.2f} kN/m."
        )
        estructural = resultado.get("diseno_estructural") or {}
        if estructural.get("ok"):
            comp = estructural.get("compresion", {})
            combos = resultado.get("combinaciones_axiales") or {}
            combo_txt = (
                f"; gobierna {combos.get('combinacion_critica')} "
                f"(Pu={float(combos.get('Pu_kN_m') or 0):.2f} kN/m)"
                if combos else ""
            )
            lineas.append(
                f"Estructural, compresión axial: {'PASA' if comp.get('pasa') else 'NO PASA'}; "
                f"demanda {float(comp.get('carga_axial_kN_m') or 0):.2f} / "
                f"capacidad φPn {float(comp.get('phi_Pn_kN') or 0):.2f} kN/m{combo_txt}. "
                "No incluye viento/sismo ni excentricidad."
            )
        else:
            lineas.append("Estructural: " + str(estructural.get("pendiente", "sin verificación disponible.")))
        cimiento = resultado.get("cimiento") or {}
        if cimiento.get("ok"):
            lineas.append(
                f"Zapata corrida preliminar: B ≥ {float(cimiento.get('ancho_minimo_m', 0)):.2f} m, "
                f"q_adm = {float(cimiento.get('q_adm_kPa', 0)):.1f} kPa, "
                f"profundidad {cimiento.get('profundidad_m')} m. "
                "Faltan excentricidad, flexión, corte, asentamiento y acciones accidentales."
            )
        else:
            lineas.append("Cimiento: " + str(cimiento.get("pendiente", "sin predimensionado guardado.")))
        return "\n".join(lineas)

    def _cargar_cargas_muro(self, nombre_muro: str) -> dict:
        """Lee las cargas activas destinadas a un muro desde datos/cargas.json."""
        try:
            datos_cargas = rutas.leer_json(rutas.CARGAS, {}) or {}
            aplicaciones = datos_cargas.get("aplicaciones", [])
            elementos = datos_cargas.get("elementos", {})
            por_id = {
                elemento.get("id"): (nombre, elemento)
                for nombre, elemento in elementos.items()
                if elemento.get("id")
            }
            resultado: dict[str, list[dict]] = {}
            for aplicacion in aplicaciones:
                if not aplicacion.get("activa", True):
                    continue
                if aplicacion.get("portico") != nombre_muro:
                    continue
                encontrado = por_id.get(aplicacion.get("carga_id"))
                if not encontrado:
                    continue
                nom_carga, fuente = encontrado
                if not fuente.get("activo", True):
                    continue
                tramo_id = aplicacion.get("tramo_id", "MURO")
                # Calcular valor lineal del elemento fuente
                valor = self._valor_lineal_elemento(fuente)
                resultado.setdefault(tramo_id, []).append({
                    "aplicacion_id": aplicacion.get("id", ""),
                    "carga_id": aplicacion.get("carga_id", ""),
                    "nombre": nom_carga,
                    "x_inicio_m": float(aplicacion.get("x_inicio_m", 0.0)),
                    "x_fin_m": float(aplicacion.get("x_fin_m", 0.0)),
                    "valor_kN_m": valor,
                    "tipo": "D",
                })
            return resultado
        except Exception:
            return {}

    def _valor_lineal_elemento(self, elemento: dict) -> float:
        """Calcula el valor lineal (kN/m) de un elemento de carga."""
        biblioteca = rutas.leer_json(rutas.MATERIALES, {}) or {}
        try:
            items = cargas.items_de_elemento(
                elemento.get("id", "?"), elemento, biblioteca
            )
            return sum(float(it["valor"]) for it in items if it.get("tipo") == "D")
        except Exception:
            return 0.0

    def _nombre_muro_a_componer(self, muros_est: dict[str, dict]) -> str | None:
        fila = self.tabla.currentRow()
        item_nombre = self.tabla.item(fila, 0) if fila >= 0 else None
        nombre = (
            item_nombre.text().strip()
            if item_nombre is not None
            else self.combo_muro_activo.currentText().strip()
        )
        return nombre if nombre in muros_est else None

    def _componer(self) -> None:
        muros_est = self._muros_estructurales()
        if not muros_est:
            QMessageBox.information(
                self, "Sin muros",
                "Creá un muro sobre un pórtico desde Inicio (estructura mixta) "
                "para poder calcularlo acá.",
            )
            return
        nombre = self._nombre_muro_a_componer(muros_est)
        if nombre is None:
            QMessageBox.information(
                self, "Seleccioná un muro",
                "Elegí un muro de la tabla o del selector de muro activo.",
            )
            return
        datos = muros_est[nombre]
        capas = datos.get("capas", []) or []
        # La altura sale de los niveles definidos en el editor de muros
        # (nivel de arranque → nivel final). Si no hay niveles, se cae al
        # valor manual guardado en altura_m.
        nivel_inicio = float(datos.get("nivel_inicio", 0.0) or 0.0)
        nivel_fin = float(datos.get("nivel_fin", 0.0) or 0.0)
        altura_de_niveles = nivel_fin > nivel_inicio
        altura = (nivel_fin - nivel_inicio) if altura_de_niveles else float(datos.get("altura_m", 0.0))
        longitud = float(datos.get("longitud_m", 1.0))
        dialogo = EditorCapasMuro(
            nombre, capas, altura, longitud, self,
            altura_de_niveles=altura_de_niveles,
        )
        if dialogo.exec() != QDialog.DialogCode.Accepted:
            return
        nuevas_capas, nueva_altura, nueva_longitud, resultado = dialogo.datos_guardar()
        estructura = rutas.cargar_estructura() or {}
        muro = estructura.get(nombre, {})
        muro["capas"] = nuevas_capas
        muro["altura_m"] = nueva_altura
        muro["espesor_m"] = resultado.get("espesor_total_m", muro.get("espesor_m", 0.0))
        muro["peso_propio_kN_m2"] = resultado.get("peso_propio_kN_m2")
        muro["peso_lineal_kN_m"] = resultado.get("peso_lineal_kN_m")
        # Guardar la memoria resumida junto con la composición en estructura.json
        # para que el último resultado no dependa de mantener abierto el diálogo.
        muro["resultado_calculo"] = {
            "resistencia_termica": resultado.get("resistencia_termica"),
            "cargas_losas": resultado.get("cargas_losas"),
            "combinaciones_axiales": resultado.get("combinaciones_axiales"),
            "diseno_estructural": resultado.get("diseno_estructural"),
            "cimiento": resultado.get("cimiento"),
        }
        # La longitud en planta se conserva (la define la estructura mixta).
        estructura[nombre] = muro
        try:
            rutas.guardar_estructura(estructura)
        except OSError as exc:
            QMessageBox.critical(self, "No se pudo guardar", str(exc))
            return
        indice_activo = self.combo_muro_activo.findText(nombre)
        if indice_activo >= 0:
            self.combo_muro_activo.setCurrentIndex(indice_activo)
        self.refrescar()
        for fila in range(self.tabla.rowCount()):
            item = self.tabla.item(fila, 0)
            if item is not None and item.text() == nombre:
                self.tabla.selectRow(fila)
                break
        QMessageBox.information(
            self, "Muro actualizado",
            f"Se guardó la composición de «{nombre}». "
            f"Espesor {muro['espesor_m']:.3f} m · peso propio "
            f"{(muro.get('peso_propio_kN_m2') or 0):.3f} kN/m².",
        )

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

from calc import mamposteria, muros, rutas


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
                combo.addItem(f"{info['nombre']} ({grupo})", (grupo, clave))
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
                "Resistencia térmica (ISO 6946, muro vertical):",
                f"  R de las capas: {termica['R_capas_m2K_W']:.3f} m²·K/W",
                f"  R total (con Rsi={termica['R_si_m2K_W']:.2f}, Rse={termica['R_se_m2K_W']:.2f}): "
                f"{termica['R_total_m2K_W']:.3f} m²·K/W",
                f"  Transmitancia térmica U = 1/Rtotal: {termica['U_W_m2K']:.3f} W/m²·K",
            ]
        else:
            lineas += [
                "",
                "Resistencia térmica: pendiente (faltan conductividades λ de "
                "alguna capa en materiales.json).",
            ]
        lineas += [
            "",
            "Pendiente de implementar: cimiento (requiere reglamentación / datos "
            "geotécnicos).",
        ]

        # Verificación estructural (CIRSOC 501) con la altura del paño.
        try:
            est = mamposteria.verificar_muro(
                capas, altura_m=self.altura.value(), mortero="mortero_normal_n",
            )
        except (FileNotFoundError, KeyError):
            est = {"ok": False,
                   "pendiente": "Falta datos/mamposteria.json (biblioteca CIRSOC 501)."}
        if est.get("ok"):
            comp = est["compresion"]
            cort = est["corte"]
            lineas += [
                "",
                "Verificación estructural (CIRSOC 501) por metro de muro:",
                f"  Hoja portante: {est['hoja_portante']['clave']} "
                f"(fm = {est['hoja_portante']['fm_mpa']:.2f} MPa, "
                f"e = {est['hoja_portante']['espesor_m']:.2f} m)",
                f"  Esbeltez h/r = {est['esbeltez_hr']:.1f}"
                f"{'  ⚠ ESBELTO' if est['es_esbelta'] else ''}",
                f"  Compresión: φ·Pn = {comp['phi_Pn_kN']:.1f} kN "
                f"(φ={comp['phi']:.2f}, {comp['rama']})",
                f"  Corte:      φ·Vn = {cort['phi_Vn_kN']:.1f} kN (φ={cort['phi']:.2f})",
                "  Recomendaciones:",
            ]
            lineas += [f"    · {r}" for r in est["recomendaciones"]]
        else:
            lineas += [
                "",
                "Verificación estructural: " + str(est.get("pendiente", "no disponible.")),
            ]
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
            "calcular su peso propio; la resistencia térmica, el diseño "
            "estructural y el cimiento se habilitan al cargar las "
            "conductividades λ y la reglamentación de mampostería."
        )
        self.intro.setWordWrap(True)

        self.tabla = QTableWidget(0, 5)
        self.tabla.setHorizontalHeaderLabels((
            "Muro", "Longitud (m)", "Espesor (m)", "Altura (m)", "Peso propio (kN/m)",
        ))
        self.tabla.verticalHeader().setVisible(False)
        self.tabla.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tabla.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.tabla.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)

        self.boton_componer = QPushButton("Componer muro / calcular peso propio…")
        self.boton_componer.clicked.connect(self._componer)
        self.boton_actualizar = QPushButton("Actualizar lista")
        self.boton_actualizar.clicked.connect(self.refrescar)

        fila = QHBoxLayout()
        fila.addWidget(self.boton_componer)
        fila.addWidget(self.boton_actualizar)
        fila.addStretch(1)

        caja = QVBoxLayout(self)
        caja.addWidget(self.intro)
        caja.addWidget(self.tabla, 1)
        caja.addLayout(fila)
        self.refrescar()

    # ------------------------------------------------------------------
    def _muros_estructurales(self) -> dict[str, dict]:
        estructura = rutas.cargar_estructura() or {}
        return {
            nombre: datos
            for nombre, datos in estructura.items()
            if isinstance(datos, dict) and str(datos.get("tipo", "")).lower() == "muro"
        }

    def refrescar(self) -> None:
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
            valores = (
                nombre,
                f"{longitud:.2f}",
                f"{espesor:.3f}",
                f"{altura:.2f}" if altura else "—",
                f"{peso_lineal:.3f}" if peso_lineal is not None else "—",
            )
            for col, valor in enumerate(valores):
                self.tabla.setItem(fila, col, QTableWidgetItem(valor))
        if not muros_est:
            self.tabla.setRowCount(1)
            self.tabla.setItem(
                0, 0,
                QTableWidgetItem("Todavía no hay muros estructurales en esta obra."),
            )

    def _componer(self) -> None:
        muros_est = self._muros_estructurales()
        if not muros_est:
            QMessageBox.information(
                self, "Sin muros",
                "Creá un muro sobre un pórtico desde Inicio (estructura mixta) "
                "para poder calcularlo acá.",
            )
            return
        nombre = self.tabla.currentItem().text() if self.tabla.currentItem() else None
        if nombre not in muros_est:
            nombre = next(iter(muros_est))
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
        # La longitud en planta se conserva (la define la estructura mixta).
        estructura[nombre] = muro
        try:
            rutas.guardar_estructura(estructura)
        except OSError as exc:
            QMessageBox.critical(self, "No se pudo guardar", str(exc))
            return
        self.refrescar()
        QMessageBox.information(
            self, "Muro actualizado",
            f"Se guardó la composición de «{nombre}». "
            f"Espesor {muro['espesor_m']:.3f} m · peso propio "
            f"{(muro.get('peso_propio_kN_m2') or 0):.3f} kN/m².",
        )



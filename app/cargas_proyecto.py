"""Pantalla para revisar y vincular las cargas definidas para la obra."""

from __future__ import annotations

import re
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QComboBox,
    QCheckBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QHBoxLayout,
    QInputDialog,
    QLabel,
    QMessageBox,
    QPushButton,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QFormLayout,
    QLineEdit,
    QVBoxLayout,
    QWidget,
)

from calc import cargas, rutas


class EditorElemento(QDialog):
    """Edición de capas de losas/cubiertas y dimensiones de muros/encadenados."""

    def __init__(self, nombre: str, elemento: dict, biblioteca: dict, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Componer carga — {nombre}")
        self.elemento = dict(elemento)
        self.biblioteca = biblioteca
        self.nombre = QLineEdit(nombre)
        self.formulario = QFormLayout()
        self.formulario.addRow("Descripción:", self.nombre)
        self.caja_datos = QVBoxLayout()
        self.formulario.addRow(self.caja_datos)
        self.setLayout(self.formulario)

        self.componentes = QTableWidget(0, 3)
        self.componentes.setHorizontalHeaderLabels(("Material", "Espesor (m)", "Tipo"))
        self.componentes.horizontalHeader().setStretchLastSection(True)
        self.grupo = QComboBox()
        self.material = QComboBox()
        for grupo, entradas in biblioteca.items():
            if isinstance(entradas, list) and any(e.get("clave") for e in entradas if isinstance(e, dict)):
                self.grupo.addItem(grupo)
        self.grupo.currentTextChanged.connect(self._cargar_materiales)
        self._cargar_materiales(self.grupo.currentText())
        self.espesor = self._numero(0.01, 5.0, 0.01, 2)
        self.boton_agregar = QPushButton("Agregar capa")
        self.boton_quitar = QPushButton("Quitar capa seleccionada")
        self.boton_agregar.clicked.connect(self._agregar_capa)
        self.boton_quitar.clicked.connect(lambda: self.componentes.removeRow(self.componentes.currentRow())
                                          if self.componentes.currentRow() >= 0 else None)

        self.tipo = self.elemento.get("tipo")
        if self.tipo in ("losa", "cubierta"):
            if self.tipo == "losa":
                self.tipologia = QComboBox()
                self.tipologia.addItems(("alivianada", "maciza", "casetonada"))
                self.tipologia.setCurrentText(self.elemento.get("tipologia", "alivianada"))
                self.caja_datos.addWidget(self._fila("Tipología:", self.tipologia))
            self.luz_completa = self._numero(0.0, 100.0, 0.1, 2)
            self.luz_completa.setValue(float(self.elemento.get("luz_transversal_m", 0.0)))
            self.caja_datos.addWidget(self._fila("Luz completa de la losa (m):", self.luz_completa))
            self.caja_datos.addWidget(QLabel(
                "El ancho tributario se elige en cada aplicación a una viga o tramo."
            ))
            self.sobrecarga = QComboBox()
            for uso in biblioteca.get("Sobrecargas", []):
                if uso.get("clave"):
                    self.sobrecarga.addItem(uso.get("uso", uso["clave"]), uso["clave"])
            idx = self.sobrecarga.findData(self.elemento.get("sobrecarga", "vivienda"))
            if idx >= 0:
                self.sobrecarga.setCurrentIndex(idx)
            self.caja_datos.addWidget(self._fila("Sobrecarga de uso:", self.sobrecarga))
            if self.tipo == "cubierta":
                self.succion_activa = QCheckBox("Incluir succión de viento en esta cubierta")
                self.succion_activa.setChecked(bool(self.elemento.get("viento_activo", 0)))
                self.pendiente = self._numero(0.0, 89.0, 1.0, 1)
                self.pendiente.setValue(float(self.elemento.get("pendiente_grados", 0.0)))
                self.altura_manual_activa = QCheckBox("Definir altura vertical manualmente")
                altura = self.elemento.get("altura_vertical_m")
                self.altura_manual_activa.setChecked(altura is not None)
                self.altura_manual = self._numero(0.0, 100.0, 0.1, 2)
                self.altura_manual.setValue(float(altura or 0.0))
                self.caja_datos.addWidget(self.succion_activa)
                self.caja_datos.addWidget(self._fila("Pendiente (grados):", self.pendiente))
                self.caja_datos.addWidget(self.altura_manual_activa)
                self.caja_datos.addWidget(self._fila("Altura vertical (m):", self.altura_manual))
            fila_capas = QHBoxLayout()
            fila_capas.addWidget(self.grupo)
            fila_capas.addWidget(self.material, 1)
            fila_capas.addWidget(QLabel("Espesor (m):"))
            fila_capas.addWidget(self.espesor)
            fila_capas.addWidget(self.boton_agregar)
            self.caja_datos.addLayout(fila_capas)
            self.caja_datos.addWidget(self.componentes)
            self.caja_datos.addWidget(self.boton_quitar)
            self._cargar_componentes()
        elif self.tipo == "muro":
            self.material_muro = QComboBox()
            self.claves_muro: list[str] = []
            for e in biblioteca.get("Mamposteria", []):
                if e.get("clave"):
                    self.material_muro.addItem(e["nombre"], e["clave"])
                    self.claves_muro.append(e["clave"])
            idx = self.material_muro.findData(self.elemento.get("clave"))
            if idx >= 0:
                self.material_muro.setCurrentIndex(idx)
            self.espesor_muro = self._numero(0.01, 2.0, 0.01, 2)
            self.altura_muro = self._numero(0.1, 30.0, 0.1, 2)
            self.espesor_muro.setValue(float(self.elemento.get("espesor_m", 0.18)))
            self.altura_muro.setValue(float(self.elemento.get("altura_m", 2.8)))
            self.caja_datos.addWidget(self._fila("Mampostería:", self.material_muro))
            self.caja_datos.addWidget(self._fila("Espesor (m):", self.espesor_muro))
            self.caja_datos.addWidget(self._fila("Altura (m):", self.altura_muro))
        elif self.tipo == "encadenado":
            self.base = self._numero(0.05, 2.0, 0.05, 2)
            self.altura = self._numero(0.05, 3.0, 0.05, 2)
            self.cantidad = QSpinBox()
            self.cantidad.setRange(1, 100)
            self.base.setValue(float(self.elemento.get("base_m", 0.2)))
            self.altura.setValue(float(self.elemento.get("altura_m", 0.2)))
            self.cantidad.setValue(int(self.elemento.get("cantidad", 1)))
            self.caja_datos.addWidget(self._fila("Base (m):", self.base))
            self.caja_datos.addWidget(self._fila("Altura (m):", self.altura))
            self.caja_datos.addWidget(self._fila("Cantidad:", self.cantidad))

        botones = QDialogButtonBox(QDialogButtonBox.StandardButton.Save | QDialogButtonBox.StandardButton.Cancel)
        botones.accepted.connect(self.accept)
        botones.rejected.connect(self.reject)
        self.caja_datos.addWidget(botones)
        self.resize(760, 520)

    @staticmethod
    def _numero(minimo: float, maximo: float, paso: float, decimales: int) -> QDoubleSpinBox:
        control = QDoubleSpinBox()
        control.setRange(minimo, maximo)
        control.setSingleStep(paso)
        control.setDecimals(decimales)
        return control

    @staticmethod
    def _fila(etiqueta: str, control) -> QWidget:
        fila = QWidget()
        layout = QHBoxLayout(fila)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(QLabel(etiqueta))
        layout.addWidget(control)
        return fila

    def _cargar_materiales(self, grupo: str) -> None:
        self.material.clear()
        for material in self.biblioteca.get(grupo, []):
            if material.get("clave"):
                self.material.addItem(material["nombre"], material["clave"])

    def _agregar_capa(self) -> None:
        grupo = self.grupo.currentText()
        clave = self.material.currentData()
        if not grupo or not clave:
            return
        material = next(m for m in self.biblioteca[grupo] if m.get("clave") == clave)
        componente = {"grupo": grupo, "clave": clave}
        if material.get("tipo") == "volumetrico":
            componente["espesor_m"] = self.espesor.value()
        fila = self.componentes.rowCount()
        self.componentes.insertRow(fila)
        nombre = material["nombre"]
        self.componentes.setItem(fila, 0, QTableWidgetItem(nombre))
        self.componentes.setItem(fila, 1, QTableWidgetItem(str(componente.get("espesor_m", "—"))))
        self.componentes.setItem(fila, 2, QTableWidgetItem(material.get("tipo", "")))
        self.componentes.item(fila, 0).setData(Qt.ItemDataRole.UserRole, componente)

    def _cargar_componentes(self) -> None:
        for componente in self.elemento.get("componentes", []):
            grupo, clave = componente.get("grupo"), componente.get("clave")
            material = next((m for m in self.biblioteca.get(grupo, []) if m.get("clave") == clave), None)
            if not material:
                continue
            fila = self.componentes.rowCount()
            self.componentes.insertRow(fila)
            self.componentes.setItem(fila, 0, QTableWidgetItem(componente.get("nombre", material["nombre"])))
            self.componentes.setItem(fila, 1, QTableWidgetItem(str(componente.get("espesor_m", "—"))))
            self.componentes.setItem(fila, 2, QTableWidgetItem(material.get("tipo", "")))
            self.componentes.item(fila, 0).setData(Qt.ItemDataRole.UserRole, dict(componente))

    def resultado(self) -> tuple[str, dict]:
        elemento = dict(self.elemento)
        nombre = self.nombre.text().strip()
        if self.tipo in ("losa", "cubierta"):
            elemento["componentes"] = [
                self.componentes.item(f, 0).data(Qt.ItemDataRole.UserRole)
                for f in range(self.componentes.rowCount())
            ]
            elemento["luz_transversal_m"] = self.luz_completa.value()
            elemento["sobrecarga"] = self.sobrecarga.currentData()
            if self.tipo == "losa":
                elemento["tipologia"] = self.tipologia.currentText()
            elif self.tipo == "cubierta":
                elemento["viento_activo"] = int(self.succion_activa.isChecked())
                elemento["pendiente_grados"] = self.pendiente.value()
                elemento["altura_vertical_m"] = (
                    self.altura_manual.value() if self.altura_manual_activa.isChecked() else None
                )
        elif self.tipo == "muro":
            elemento.update({
                "grupo": "Mamposteria", "clave": self.material_muro.currentData(),
                "espesor_m": self.espesor_muro.value(), "altura_m": self.altura_muro.value(),
            })
        elif self.tipo == "encadenado":
            elemento.update({"base_m": self.base.value(), "altura_m": self.altura.value(),
                             "cantidad": self.cantidad.value()})
        return nombre, elemento


class EditorAplicaciones(QDialog):
    """Permite aplicar una misma carga varias veces, con distinto tramo y ancho."""

    def __init__(self, carga_id: str, elemento: dict, estructura: dict, aplicaciones: list[dict], parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Aplicaciones de {carga_id}")
        self.carga_id = carga_id
        self.elemento = elemento
        self.aplicaciones = [dict(a) for a in aplicaciones if a.get("carga_id") == carga_id]
        self.todos = [dict(a) for a in aplicaciones if a.get("carga_id") != carga_id]
        self.estructura = estructura
        self.es_superficial = elemento.get("tipo") in ("losa", "cubierta")

        self.tabla = QTableWidget(0, 5)
        self.tabla.setHorizontalHeaderLabels(("Pórtico", "Tramo", "Desde (m)", "Hasta (m)", "Ancho trib. (m)"))
        self.tabla.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.tabla.setSelectionMode(QTableWidget.SelectionMode.SingleSelection)
        self.tabla.horizontalHeader().setStretchLastSection(True)
        self.portico = QComboBox()
        self.tramo = QComboBox()
        self.desde = QDoubleSpinBox()
        self.hasta = QDoubleSpinBox()
        for control in (self.desde, self.hasta):
            control.setRange(0.0, 100.0)
            control.setDecimals(2)
            control.setSingleStep(0.1)
        self.ancho_modo = QComboBox()
        self.ancho_modo.addItem("Media luz de la losa", "media_luz")
        self.ancho_modo.addItem("Luz completa de la losa", "luz_completa")
        self.ancho_modo.addItem("Ancho personalizado", "manual")
        self.ancho_manual = QDoubleSpinBox()
        self.ancho_manual.setRange(0.01, 100.0)
        self.ancho_manual.setDecimals(2)
        self.ancho_manual.setSingleStep(0.1)
        self.modo_legacy = QComboBox()
        self.modo_legacy.addItem("Reemplazar D/L distribuidas antiguas de P00", "reemplazar")
        self.modo_legacy.addItem("Sumar a D/L distribuidas antiguas de P00", "sumar")
        modo_existente = next((a.get("modo_cargas_previas") for a in self.aplicaciones
                               if a.get("modo_cargas_previas")), "reemplazar")
        self.modo_legacy.setCurrentIndex(max(0, self.modo_legacy.findData(modo_existente)))
        self.boton_agregar = QPushButton("Agregar aplicación")
        self.boton_quitar = QPushButton("Quitar aplicación seleccionada")
        self.boton_agregar.clicked.connect(self._agregar)
        self.boton_quitar.clicked.connect(self._quitar)
        self.portico.currentIndexChanged.connect(self._cargar_tramos)
        self.tramo.currentIndexChanged.connect(self._actualizar_rango)

        for nombre in estructura:
            self.portico.addItem(nombre, nombre)
        self._cargar_tramos()
        self.ancho_modo.setEnabled(self.es_superficial)
        self.ancho_manual.setEnabled(self.es_superficial)
        self.ancho_modo.currentIndexChanged.connect(self._actualizar_ancho)
        self._actualizar_ancho()

        caja = QVBoxLayout(self)
        caja.addWidget(QLabel(
            "Cada fila es una aplicación independiente. Las posiciones se miden desde el inicio local del tramo."
        ))
        caja.addWidget(self.tabla, 1)
        formulario = QFormLayout()
        formulario.addRow("Pórtico:", self.portico)
        formulario.addRow("Viga / tramo:", self.tramo)
        formulario.addRow("Posición inicial local (m):", self.desde)
        formulario.addRow("Posición final local (m):", self.hasta)
        if self.es_superficial:
            formulario.addRow("Ancho tributario:", self.ancho_modo)
            formulario.addRow("Ancho personalizado (m):", self.ancho_manual)
        formulario.addRow("Cargas D/L que ya estaban guardadas por P00:", self.modo_legacy)
        caja.addLayout(formulario)
        fila = QHBoxLayout()
        fila.addWidget(self.boton_agregar)
        fila.addWidget(self.boton_quitar)
        caja.addLayout(fila)
        botones = QDialogButtonBox(QDialogButtonBox.StandardButton.Save | QDialogButtonBox.StandardButton.Cancel)
        botones.accepted.connect(self.accept)
        botones.rejected.connect(self.reject)
        caja.addWidget(botones)
        self._refrescar_tabla()
        self.resize(850, 480)

    def _cargar_tramos(self, *_args) -> None:
        nombre = self.portico.currentData()
        self.tramo.blockSignals(True)
        self.tramo.clear()
        p = self.estructura.get(nombre, {})
        for viga_id, viga in p.get("vigas", {}).items():
            for tramo in viga.get("tramos", []):
                tid = tramo.get("id")
                if tid:
                    self.tramo.addItem(f"{tid} ({float(tramo.get('longitud_m', 0)):.2f} m)",
                                       (tid, float(tramo.get("longitud_m", 0))))
        self.tramo.blockSignals(False)
        self._actualizar_rango()

    def _actualizar_rango(self, *_args) -> None:
        tramo = self.tramo.currentData()
        longitud = float(tramo[1]) if tramo else 0.0
        self.desde.setMaximum(longitud)
        self.hasta.setRange(0.0, longitud)
        self.hasta.setValue(longitud)

    def _actualizar_ancho(self, *_args) -> None:
        self.ancho_manual.setEnabled(self.es_superficial and self.ancho_modo.currentData() == "manual")

    def _agregar(self) -> None:
        portico = self.portico.currentData()
        tramo = self.tramo.currentData()
        if not portico or not tramo:
            QMessageBox.warning(self, "Sin tramo", "Elegí un tramo del pórtico.")
            return
        x0, x1 = self.desde.value(), self.hasta.value()
        if x1 <= x0:
            QMessageBox.warning(self, "Intervalo inválido", "La posición final debe ser mayor que la inicial.")
            return
        luz = float(self.elemento.get("luz_transversal_m", 0.0))
        modo = self.ancho_modo.currentData() if self.es_superficial else None
        if self.es_superficial:
            if modo != "manual" and luz <= 0:
                QMessageBox.warning(self, "Falta la luz de la losa", "Editá la composición e ingresá la luz completa de la losa.")
                return
            ancho = self.ancho_manual.value() if modo == "manual" else (luz / 2 if modo == "media_luz" else luz)
        else:
            ancho = None
        modos_mismo_tramo = {
            a.get("modo_cargas_previas", "reemplazar")
            for a in self.todos + self.aplicaciones
            if a.get("portico") == portico and a.get("tramo_id") == tramo[0]
        }
        if len(modos_mismo_tramo) > 1:
            QMessageBox.warning(
                self, "Cargas anteriores", 
                "Este tramo ya tiene aplicaciones con criterios distintos para las cargas de P00. "
                "Quitá y volvé a crear esas aplicaciones con un mismo criterio antes de agregar otra."
            )
            return
        if modos_mismo_tramo and self.modo_legacy.currentData() not in modos_mismo_tramo:
            QMessageBox.warning(
                self, "Cargas anteriores",
                "Todas las aplicaciones del mismo tramo deben usar el mismo criterio para las cargas de P00."
            )
            return
        numero = len(self.todos) + len(self.aplicaciones) + 1
        while any(a.get("id") == f"Aplicacion 0-{numero}" for a in self.todos + self.aplicaciones):
            numero += 1
        self.aplicaciones.append({
            "id": f"Aplicacion 0-{numero}", "carga_id": self.carga_id,
            "portico": portico, "tramo_id": tramo[0],
            "x_inicio_m": x0, "x_fin_m": x1,
            "ancho_modo": modo, "ancho_tributario_m": ancho,
            "modo_cargas_previas": self.modo_legacy.currentData(),
            "activa": True,
        })
        self._refrescar_tabla()

    def _quitar(self) -> None:
        fila = self.tabla.currentRow()
        if fila >= 0:
            self.aplicaciones.pop(fila)
            self._refrescar_tabla()

    def _refrescar_tabla(self) -> None:
        self.tabla.setRowCount(len(self.aplicaciones))
        for i, a in enumerate(self.aplicaciones):
            valores = (a.get("portico", ""), a.get("tramo_id", ""), a.get("x_inicio_m", 0),
                       a.get("x_fin_m", 0), a.get("ancho_tributario_m", "—"))
            for j, valor in enumerate(valores):
                self.tabla.setItem(i, j, QTableWidgetItem(str(valor)))
        self.tabla.resizeColumnsToContents()

    def resultado(self) -> list[dict]:
        return self.todos + self.aplicaciones


class PaginaCargas(QWidget):
    """Catálogo del proyecto y asignación de cargas a intervalos de tramos."""

    def __init__(self, al_guardar=None):
        super().__init__()
        self.al_guardar = al_guardar
        self.etiqueta = QLabel()
        self.etiqueta.setWordWrap(True)
        self.viento_activo = QCheckBox("Aplicar viento horizontal general al pórtico")
        self.ancho_viento = QDoubleSpinBox()
        self.ancho_viento.setRange(0.0, 100.0)
        self.ancho_viento.setDecimals(2)
        self.ancho_viento.setSingleStep(0.1)
        catalogo_viento = rutas.leer_json(rutas.DATOS / "viento_cirsoc_102_25.json", {}) or {}
        self.catalogo_viento = catalogo_viento.get("ciudades", {})
        self.ciudad_viento = QComboBox()
        for ciudad in self.catalogo_viento:
            self.ciudad_viento.addItem(ciudad, ciudad)
        self.categoria_riesgo_viento = QComboBox()
        self.categoria_riesgo_viento.addItem("Elegir categoría…", None)
        for categoria in catalogo_viento.get("categorias_riesgo", []):
            self.categoria_riesgo_viento.addItem(f"Categoría {categoria}", categoria)
        self.velocidad_viento_ref = QLabel()
        self.velocidad_viento_ref.setMinimumWidth(150)
        self.nota_viento_ref = QLabel(
            "Referencia guardada en el proyecto; el motor todavía calcula el viento con la configuración anterior."
        )
        self.nota_viento_ref.setWordWrap(True)
        self.ciudad_viento.currentIndexChanged.connect(self._actualizar_velocidad_referencia)
        self.categoria_riesgo_viento.currentIndexChanged.connect(self._actualizar_velocidad_referencia)
        self.tabla = QTableWidget(0, 6)
        self.tabla.setHorizontalHeaderLabels(
            ("ID", "Descripción", "Aplicación", "Categoría", "Activa", "N.º aplicaciones")
        )
        self.tabla.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.tabla.setSelectionMode(QTableWidget.SelectionMode.SingleSelection)
        self.tabla.setAlternatingRowColors(True)
        self.tabla.verticalHeader().setVisible(False)
        self.tabla.horizontalHeader().setStretchLastSection(True)

        self.boton_nueva = QPushButton("Nueva carga")
        self.boton_guardar = QPushButton("Guardar cargas")
        self.boton_editar = QPushButton("Editar composición")
        self.boton_aplicaciones = QPushButton("Aplicar a tramos…")
        self.boton_informe = QPushButton("Generar informe TXT")
        self.boton_abrir = QPushButton("Abrir último informe")
        self.boton_abrir.setEnabled(False)
        self.ultimo_informe: Path | None = None

        caja = QVBoxLayout(self)
        caja.addWidget(self.etiqueta)
        fila_viento = QHBoxLayout()
        fila_viento.addWidget(self.viento_activo)
        fila_viento.addWidget(QLabel("Ancho tributario del viento (m):"))
        fila_viento.addWidget(self.ancho_viento)
        fila_viento.addStretch(1)
        caja.addLayout(fila_viento)
        fila_referencia_viento = QHBoxLayout()
        fila_referencia_viento.addWidget(QLabel("Ubicación CIRSOC 102-25:"))
        fila_referencia_viento.addWidget(self.ciudad_viento)
        fila_referencia_viento.addWidget(QLabel("Categoría de riesgo:"))
        fila_referencia_viento.addWidget(self.categoria_riesgo_viento)
        fila_referencia_viento.addWidget(self.velocidad_viento_ref)
        fila_referencia_viento.addStretch(1)
        caja.addLayout(fila_referencia_viento)
        caja.addWidget(self.nota_viento_ref)
        caja.addWidget(self.tabla, 1)
        fila = QHBoxLayout()
        fila.addWidget(self.boton_aplicaciones)
        fila.addWidget(self.boton_nueva)
        fila.addWidget(self.boton_editar)
        fila.addWidget(self.boton_guardar)
        caja.addLayout(fila)
        acciones = QHBoxLayout()
        acciones.addWidget(self.boton_informe)
        acciones.addWidget(self.boton_abrir)
        acciones.addStretch(1)
        caja.addLayout(acciones)

        self.tabla.itemChanged.connect(self._actualizar_estado_activo)
        self.boton_guardar.clicked.connect(self.guardar)
        self.boton_nueva.clicked.connect(self.nueva)
        self.boton_editar.clicked.connect(self.editar)
        self.boton_aplicaciones.clicked.connect(self.editar_aplicaciones)
        self.boton_informe.clicked.connect(self.generar_informe)
        self.boton_abrir.clicked.connect(self.abrir_informe)
        self.recargar()

    def _actualizar_estado_activo(self, item: QTableWidgetItem) -> None:
        if item.column() == 4:
            estado = "Sí" if item.checkState() == Qt.CheckState.Checked else "No"
            if item.text() != estado:
                item.setText(estado)

    def _actualizar_velocidad_referencia(self, *_args) -> None:
        ciudad = self.ciudad_viento.currentData()
        categoria = self.categoria_riesgo_viento.currentData()
        valor = (self.catalogo_viento.get(ciudad, {}).get("velocidades", {}) or {}).get(categoria)
        self.velocidad_viento_ref.setText(
            f"V básica: {valor:.1f} m/s" if valor is not None else "V básica: elegir categoría"
        )

    @staticmethod
    def _id_nuevo(elementos: dict) -> str:
        usados = {str(e.get("id", "")) for e in elementos.values()}
        n = 1
        while f"Carga 0-{n}" in usados:
            n += 1
        return f"Carga 0-{n}"

    @staticmethod
    def _id_proyecto(datos: dict) -> str:
        actual = str(datos.get("id_proyecto", "")).strip()
        if actual:
            return actual
        obra = str(datos.get("obra", "obra")).strip().lower()
        return re.sub(r"[^a-z0-9]+", "_", obra).strip("_") or "obra"

    def recargar(self) -> None:
        datos = cargas.datos_cargas()
        id_proyecto = self._id_proyecto(datos)
        viento = datos.get("viento", {})
        self.viento_activo.setChecked(bool(viento.get("activo", 0)))
        self.ancho_viento.setValue(float(viento.get("ancho_tributario_m", 0.0)))
        referencia = viento.get("referencia_cirsoc_102_25", {})
        ciudad = referencia.get("ubicacion", {}).get("ciudad", "Santa Fe")
        indice_ciudad = self.ciudad_viento.findData(ciudad)
        self.ciudad_viento.setCurrentIndex(max(indice_ciudad, 0))
        indice_categoria = self.categoria_riesgo_viento.findData(referencia.get("categoria_riesgo"))
        self.categoria_riesgo_viento.setCurrentIndex(max(indice_categoria, 0))
        self._actualizar_velocidad_referencia()
        elementos = datos.get("elementos", {})
        aplicaciones = datos.get("aplicaciones", [])
        self.tabla.setRowCount(len(elementos))
        self.tabla.setProperty("carga_names", list(elementos))
        for fila, (nombre, elemento) in enumerate(elementos.items()):
            if not elemento.get("id"):
                elemento["id"] = self._id_nuevo(elementos)
            tipo = str(elemento.get("tipo", "?"))
            forma = {
                "losa": "superficial → lineal",
                "cubierta": "superficial → lineal",
                "muro": "lineal",
                "encadenado": "lineal",
            }.get(tipo, tipo)
            categorias = "D + L" if tipo in ("losa", "cubierta") else "D"
            if elemento.get("viento_activo"):
                categorias += " + W"
            n_aplicaciones = sum(1 for a in aplicaciones if a.get("carga_id") == elemento.get("id"))
            valores = (
                elemento["id"], nombre, forma, categorias,
                "Sí" if elemento.get("activo") else "No",
                str(n_aplicaciones),
            )
            for columna, valor in enumerate(valores):
                item = QTableWidgetItem(str(valor))
                item.setData(Qt.ItemDataRole.UserRole, nombre)
                if columna == 4:
                    item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
                    item.setCheckState(
                        Qt.CheckState.Checked if elemento.get("activo") else Qt.CheckState.Unchecked
                    )
                self.tabla.setItem(fila, columna, item)
        self.etiqueta.setText(
            f"{len(elementos)} elementos de carga. "
            "Losas y cubiertas convierten kN/m² a kN/m con el ancho tributario de cada aplicación; "
            "muros y encadenados generan carga lineal. Marcá las cargas que entran "
            "al informe. Los cambios se guardan con los botones de abajo; los TXT "
            "anteriores quedan como estaban. Una carga puede tener varias aplicaciones "
            "en distintos tramos; las puntuales siguen ingresándose desde P00."
        )
        self.tabla.resizeColumnsToContents()

    def _fila_nombre(self) -> str | None:
        fila = self.tabla.currentRow()
        item = self.tabla.item(fila, 0) if fila >= 0 else None
        return item.data(Qt.ItemDataRole.UserRole) if item else None

    def editar(self) -> None:
        nombre = self._fila_nombre()
        if not nombre:
            QMessageBox.information(self, "Cargas", "Elegí primero una carga de la lista.")
            return
        self._guardar_datos()
        datos = cargas.datos_cargas()
        elemento = datos.get("elementos", {}).get(nombre, {})
        editor = EditorElemento(nombre, elemento, rutas.leer_json(rutas.MATERIALES, {}) or {}, self)
        if editor.exec() != QDialog.DialogCode.Accepted:
            return
        nuevo_nombre, actualizado = editor.resultado()
        if not nuevo_nombre:
            QMessageBox.warning(self, "Descripción requerida", "La carga necesita una descripción.")
            return
        elementos = datos["elementos"]
        if nuevo_nombre != nombre and nuevo_nombre in elementos:
            QMessageBox.warning(self, "Descripción existente", "Ya existe una carga con ese nombre.")
            return
        elementos.pop(nombre)
        elementos[nuevo_nombre] = actualizado
        rutas.guardar_json(rutas.CARGAS, datos)
        self.recargar()
        if self.al_guardar:
            self.al_guardar()

    def editar_aplicaciones(self) -> None:
        nombre = self._fila_nombre()
        if not nombre:
            QMessageBox.information(self, "Cargas", "Elegí una carga de la lista.")
            return
        self._guardar_datos()
        datos = cargas.datos_cargas()
        elemento = datos.get("elementos", {}).get(nombre, {})
        carga_id = elemento.get("id")
        if not carga_id:
            QMessageBox.warning(self, "Falta ID", "No se pudo identificar la carga del proyecto.")
            return
        editor = EditorAplicaciones(
            carga_id, elemento, rutas.cargar_estructura(), datos.get("aplicaciones", []), self
        )
        if editor.exec() != QDialog.DialogCode.Accepted:
            return
        datos["aplicaciones"] = editor.resultado()
        rutas.guardar_json(rutas.CARGAS, datos)
        self.recargar()
        if self.al_guardar:
            self.al_guardar()

    def nueva(self) -> None:
        tipos = ("losa", "cubierta", "muro", "encadenado")
        tipo, ok = QInputDialog.getItem(self, "Nueva carga", "Tipo de elemento:", tipos, 0, False)
        if not ok:
            return
        datos = cargas.datos_cargas()
        biblioteca = rutas.leer_json(rutas.MATERIALES, {}) or {}
        if tipo in ("losa", "cubierta"):
            elemento = {
                "tipo": tipo, "activo": 1, "ancho_tributario_m": 1.0,
                "componentes": [],
                "sobrecarga": "vivienda", "viento_activo": 0,
            }
            if tipo == "losa":
                elemento["tipologia"] = "alivianada"
            else:
                elemento.update({"pendiente_grados": 0, "altura_vertical_m": None})
        elif tipo == "muro":
            primero = next((m for m in biblioteca.get("Mamposteria", []) if m.get("clave")), {})
            elemento = {"tipo": tipo, "activo": 1, "grupo": "Mamposteria",
                        "clave": primero.get("clave", "ladrillo_hueco"),
                        "espesor_m": 0.18, "altura_m": 2.8}
        else:
            elemento = {"tipo": tipo, "activo": 1, "base_m": 0.2,
                        "altura_m": 0.2, "cantidad": 1}
        editor = EditorElemento("Carga nueva", elemento, biblioteca, self)
        if editor.exec() != QDialog.DialogCode.Accepted:
            return
        nombre, elemento = editor.resultado()
        if not nombre:
            QMessageBox.warning(self, "Descripción requerida", "La carga necesita una descripción.")
            return
        elementos = datos.setdefault("elementos", {})
        if nombre in elementos:
            QMessageBox.warning(self, "Descripción existente", "Ya existe una carga con ese nombre.")
            return
        elemento["id"] = self._id_nuevo(elementos)
        elementos[nombre] = elemento
        datos["id_proyecto"] = self._id_proyecto(datos)
        rutas.guardar_json(rutas.CARGAS, datos)
        self.recargar()
        if self.al_guardar:
            self.al_guardar()

    def _guardar_datos(self) -> None:
        datos = cargas.datos_cargas()
        datos["id_proyecto"] = self._id_proyecto(datos)
        viento = datos.setdefault("viento", {})
        viento["activo"] = int(self.viento_activo.isChecked())
        viento["ancho_tributario_m"] = self.ancho_viento.value()
        ciudad = self.ciudad_viento.currentData()
        categoria = self.categoria_riesgo_viento.currentData()
        datos_ciudad = self.catalogo_viento.get(ciudad, {})
        velocidad = (datos_ciudad.get("velocidades", {}) or {}).get(categoria)
        viento["referencia_cirsoc_102_25"] = {
            "ubicacion": {
                "pais": "Argentina",
                "provincia": datos_ciudad.get("provincia", ""),
                "ciudad": ciudad,
            },
            "categoria_riesgo": categoria,
            "V_basica_ms": velocidad,
            "fuente": "CIRSOC 102-25, Figura 1.5-1D",
        }
        elementos = datos.get("elementos", {})
        usados: set[str] = set()
        siguiente = 1
        for elemento in elementos.values():
            identificador = str(elemento.get("id", "")).strip()
            if not identificador or identificador in usados:
                while f"Carga 0-{siguiente}" in usados:
                    siguiente += 1
                identificador = f"Carga 0-{siguiente}"
                siguiente += 1
                elemento["id"] = identificador
            usados.add(identificador)
        nombre = self._fila_nombre()
        for fila in range(self.tabla.rowCount()):
            item = self.tabla.item(fila, 0)
            activo = self.tabla.item(fila, 4)
            if item and activo:
                elementos[item.data(Qt.ItemDataRole.UserRole)]["activo"] = (
                    activo.checkState() == Qt.CheckState.Checked
                )
        rutas.guardar_json(rutas.CARGAS, datos)

    def guardar(self) -> None:
        self._guardar_datos()
        self.recargar()
        if self.al_guardar:
            self.al_guardar()
        QMessageBox.information(self, "Cargas", "Asignación guardada en datos/cargas.json.")

    def generar_informe(self) -> None:
        try:
            self._guardar_datos()
            datos = cargas.datos_cargas()
            aplicaciones = datos.get("aplicaciones", [])
            if any(a.get("activa", True) for a in aplicaciones):
                analisis = cargas.AnalisisCargas(datos.get("obra", "Obra"))
                analisis.id_proyecto = self._id_proyecto(datos)
                texto = cargas.informe_aplicaciones(datos)
            else:
                analisis = cargas.cargas_del_conjunto()
                texto = cargas.informe_completo(analisis)
            self.ultimo_informe = cargas.guardar_informe(analisis, texto)
        except (OSError, KeyError, ValueError) as exc:
            QMessageBox.critical(self, "No se pudo generar el informe", str(exc))
            return
        self.boton_abrir.setEnabled(True)
        QMessageBox.information(self, "Informe generado", f"Se guardó:\n{self.ultimo_informe}")
        if self.al_guardar:
            self.al_guardar()

    def abrir_informe(self) -> None:
        if self.ultimo_informe and self.ultimo_informe.exists():
            from app.principal import abrir_con_windows

            error = abrir_con_windows(self.ultimo_informe)
            if error:
                QMessageBox.warning(self, "No se pudo abrir", error)

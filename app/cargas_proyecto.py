"""Pantalla para revisar y vincular las cargas definidas para la obra."""

from __future__ import annotations

import re
import math
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
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QFormLayout,
    QLineEdit,
    QVBoxLayout,
    QWidget,
)

from calc import cargas, losas, losas_macizas, rutas
from app.ejes import cargar_ejes, cargar_niveles
from app.losas_ejes import EditorLosaDesdeEjes


class EditorElemento(QDialog):
    """Edición de capas de losas/cubiertas y dimensiones de muros/encadenados."""

    def __init__(self, nombre: str, elemento: dict, biblioteca: dict, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Componer carga — {nombre}")
        self.elemento = dict(elemento)
        self.tipo = self.elemento.get("tipo")
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
        grupos = list(biblioteca.items())
        preferidos = (
            ("Forjados", "Hormigon") if self.tipo == "losa"
            else ("Cubiertas", "EstructuraCubierta") if self.tipo == "cubierta"
            else ()
        )
        grupos.sort(key=lambda par: preferidos.index(par[0]) if par[0] in preferidos else len(preferidos))
        nombres_grupo = {
            "Forjados": "Forjado estructural",
            "Hormigon": "Hormigón (materiales y capas)",
            "Morteros_Revoques": "Morteros y revoques",
            "Mamposteria": "Mampostería",
            "Viguetas_Bovedillas": "Viguetas y bovedillas",
        }
        for grupo, entradas in grupos:
            if isinstance(entradas, list) and any(e.get("clave") for e in entradas if isinstance(e, dict)):
                self.grupo.addItem(nombres_grupo.get(grupo, grupo), grupo)
        self.grupo.currentIndexChanged.connect(self._cargar_materiales)
        self._cargar_materiales()
        self.espesor = self._numero(0.01, 5.0, 0.01, 2)
        self.boton_agregar = QPushButton("Agregar capa")
        self.boton_quitar = QPushButton("Quitar capa seleccionada")
        self.boton_agregar.clicked.connect(self._agregar_capa)
        self.boton_quitar.clicked.connect(lambda: self.componentes.removeRow(self.componentes.currentRow())
                                          if self.componentes.currentRow() >= 0 else None)

        if self.tipo in ("losa", "cubierta"):
            if self.tipo == "losa":
                self.tipologia = QComboBox()
                self.tipologia.addItems(("alivianada", "maciza", "casetonada"))
                self.tipologia.setCurrentText(self.elemento.get("tipologia", "alivianada"))
                self.caja_datos.addWidget(self._fila("Tipología:", self.tipologia))
            self.luz_completa = self._numero(0.0, 100.0, 0.1, 2)
            self.luz_completa.setValue(float(self.elemento.get("luz_transversal_m", 0.0)))
            if self.elemento.get("panel_ejes"):
                self.luz_completa.setReadOnly(True)
                self.luz_completa.setToolTip(
                    "Dimensión derivada de los ejes y apoyos; editá el paño desde el plano."
                )
            etiqueta_luz = (
                (
                    "Luz del paño según ejes (m):"
                    if self.elemento.get("panel_ejes")
                    else "Luz libre entre apoyos (m):"
                ) if self.tipo == "losa"
                else "Luz transversal de cubierta (m):"
            )
            self.caja_datos.addWidget(self._fila(etiqueta_luz, self.luz_completa))
            if self.tipo == "losa":
                self.ancho_losa = self._numero(0.0, 100.0, 0.1, 2)
                self.ancho_losa.setValue(float(self.elemento.get("ancho_losa_m", 0.0)))
                if self.elemento.get("panel_ejes"):
                    self.ancho_losa.setReadOnly(True)
                    self.ancho_losa.setToolTip(
                        "Dimensión derivada de los ejes y desfases de bordes libres; "
                        "editá el paño desde el plano."
                    )
                self.caja_datos.addWidget(self._fila(
                    (
                        "Ancho del paño según ejes (m):"
                        if self.elemento.get("panel_ejes")
                        else "Ancho del paño paralelo a los apoyos (m):"
                    ),
                    self.ancho_losa,
                ))
                componentes_previos = self.elemento.get("componentes", [])
                comp_maciza = next((
                    c for c in componentes_previos
                    if (c.get("grupo"), c.get("clave")) in {
                        ("Forjados", "Losa_maciza"), ("Hormigon", "armado")
                    }
                ), {})
                self.espesor_maciza = self._numero(0.05, 2.0, 0.01, 2)
                espesor_guardado = float(comp_maciza.get("espesor_m", 0.0) or 0.0)
                if espesor_guardado > 0:
                    self.espesor_maciza.setValue(espesor_guardado)
                else:
                    self.espesor_maciza.setValue(max(0.05, self._espesor_cirsoc()))
                self.espesor_automatico = QCheckBox(
                    "Predimensionar automáticamente según CIRSOC 201-25 (luz/20)"
                )
                self.espesor_automatico.setChecked(bool(
                    self.elemento.get("espesor_automatico", True)
                ))
                self.campo_espesor_maciza = self._fila(
                    "Espesor estructural de losa maciza (m):", self.espesor_maciza
                )
                self.caja_datos.addWidget(self.espesor_automatico)
                self.caja_datos.addWidget(self.campo_espesor_maciza)
                base_alivianada = next((
                    c for c in componentes_previos
                    if c.get("grupo") == "Forjados" and c.get("clave") == "Losa_alivianada"
                ), {})
                catalogo_alivianada = next((
                    m for m in self.biblioteca.get("Forjados", [])
                    if m.get("clave") == "Losa_alivianada"
                ), {})
                peso_catalogo = float(catalogo_alivianada.get("valor", 1.81))
                peso_guardado = float(base_alivianada.get("peso_propio_override_kNm2", peso_catalogo))
                self.peso_propio_alivianada = self._numero(0.01, 20.0, 0.05, 2)
                self.peso_propio_alivianada.setValue(peso_guardado)
                self.campo_peso_alivianada = self._fila(
                    "Peso propio base de losa alivianada (kN/m²):", self.peso_propio_alivianada
                )
                self.caja_datos.addWidget(self.campo_peso_alivianada)
                self.caja_datos.addWidget(QLabel(
                    "La alivianada agrega el forjado base de catálogo y la maciza calcula "
                    "su peso con el espesor. La casetonada no agrega un forjado automáticamente: "
                    "incluí sus componentes permanentes en esta composición; su motor de cálculo "
                    "todavía está pendiente. El peso de catálogo de la alivianada es editable "
                    "y aún no se recalibra según la vigueta elegida."
                ))
                self.luz_completa.valueChanged.connect(self._actualizar_espesor_cirsoc)
                self.espesor_automatico.toggled.connect(self._actualizar_espesor_cirsoc)
                self.tipologia.currentTextChanged.connect(self._actualizar_campos_tipologia)
                self._actualizar_campos_tipologia()
            self.caja_datos.addWidget(QLabel(
                "En paños creados desde ejes, las cargas D/L se transfieren automáticamente "
                "a las vigas de apoyo; las losas legadas usan el ancho tributario de su aplicación."
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
            fila_capas.addWidget(QLabel("Grupo:"))
            fila_capas.addWidget(self.grupo)
            fila_capas.addWidget(QLabel("Material:"))
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

    def _cargar_materiales(self, *_args) -> None:
        grupo = self.grupo.currentData()
        self.material.clear()
        for material in self.biblioteca.get(grupo, []):
            if material.get("clave"):
                self.material.addItem(material["nombre"], material["clave"])

    def _agregar_capa(self) -> None:
        grupo = self.grupo.currentData()
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

    def _espesor_cirsoc(self) -> float:
        """Predimensionado conservador para una losa maciza unidireccional simple: h = L/20."""
        luz = float(self.luz_completa.value()) if hasattr(self, "luz_completa") else 0.0
        return math.ceil(luz / 20.0 * 100.0 - 1e-9) / 100.0 if luz > 0 else 0.05

    def _actualizar_espesor_cirsoc(self, *_args) -> None:
        if self.espesor_automatico.isChecked():
            self.espesor_maciza.setValue(self._espesor_cirsoc())
        self.espesor_maciza.setEnabled(not self.espesor_automatico.isChecked())

    def _actualizar_campos_tipologia(self, *_args) -> None:
        visible = self.tipologia.currentText() == "maciza"
        alivianada = self.tipologia.currentText() == "alivianada"
        self.espesor_automatico.setVisible(visible)
        self.campo_espesor_maciza.setVisible(visible)
        self.campo_peso_alivianada.setVisible(alivianada)
        if visible:
            self._actualizar_espesor_cirsoc()

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
                elemento["ancho_losa_m"] = self.ancho_losa.value()
                tipologia = self.tipologia.currentText()
                elemento["tipologia"] = tipologia
                sistemas_forjado = {"Losa_alivianada", "Losa_maciza", "Losa_casetonada"}
                componentes = [
                    dict(c) for c in elemento["componentes"]
                    if not (c.get("grupo") == "Forjados" and c.get("clave") in sistemas_forjado)
                ]
                if tipologia == "alivianada":
                    catalogo = next((
                        m for m in self.biblioteca.get("Forjados", [])
                        if m.get("clave") == "Losa_alivianada"
                    ), {})
                    peso_catalogo = float(catalogo.get("valor", 1.81))
                    componente_forjado = {"grupo": "Forjados", "clave": "Losa_alivianada"}
                    peso_propio = self.peso_propio_alivianada.value()
                    if abs(peso_propio - peso_catalogo) > 1e-9:
                        componente_forjado["peso_propio_override_kNm2"] = peso_propio
                    componentes.append(componente_forjado)
                elif tipologia == "maciza":
                    # Migra el espesor estructural legado a Forjados y evita contar dos veces el hormigón.
                    componentes = [
                        c for c in componentes
                        if not (c.get("grupo") == "Hormigon" and c.get("clave") == "armado")
                    ]
                    componentes.append({
                        "grupo": "Forjados",
                        "clave": "Losa_maciza",
                        "espesor_m": self.espesor_maciza.value(),
                    })
                    elemento["espesor_automatico"] = self.espesor_automatico.isChecked()
                elemento["componentes"] = componentes
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
    """Asigna cargas lineales, reacciones de muros y anchos tributarios."""

    def __init__(self, carga_id: str, elemento: dict, estructura: dict, aplicaciones: list[dict], parent=None,
                 portico_inicial: str = ""):
        super().__init__(parent)
        self.setWindowTitle(f"Aplicaciones de {carga_id}")
        self.carga_id = carga_id
        self.elemento = elemento
        self.aplicaciones = [dict(a) for a in aplicaciones if a.get("carga_id") == carga_id]
        self.todos = [dict(a) for a in aplicaciones if a.get("carga_id") != carga_id]
        self.estructura = estructura
        self.es_superficial = elemento.get("tipo") in ("losa", "cubierta")
        self.es_muro = elemento.get("tipo") == "muro"
        self.pares_porticos = self._pares_porticos(portico_inicial)

        self.tabla = QTableWidget(0, 5)
        self.tabla.setHorizontalHeaderLabels(("Pórtico(s)", "Aplicación", "Tramo(s)", "Posición / intervalo (m)", "Ancho trib. (m)"))
        self.tabla.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.tabla.setSelectionMode(QTableWidget.SelectionMode.SingleSelection)
        self.tabla.horizontalHeader().setStretchLastSection(True)

        self.modo_aplicacion = QComboBox()
        self.modo_aplicacion.addItem("Muro lineal sobre una viga", "lineal")
        self.modo_aplicacion.addItem("Muro perpendicular entre pórticos", "perpendicular")
        self.modo_aplicacion.setVisible(self.es_muro)

        self.portico = QComboBox()
        self.tramo = QComboBox()
        self.desde = self._numero(0.0, 100.0)
        self.hasta = self._numero(0.0, 100.0)
        self.panel_lineal = QWidget()
        formulario_lineal = QFormLayout(self.panel_lineal)
        formulario_lineal.setContentsMargins(0, 0, 0, 0)
        formulario_lineal.addRow("Pórtico:", self.portico)
        formulario_lineal.addRow("Viga / tramo:", self.tramo)
        formulario_lineal.addRow("Posición inicial local (m):", self.desde)
        formulario_lineal.addRow("Posición final local (m):", self.hasta)

        self.par_porticos = QComboBox()
        self.tramo_a = QComboBox()
        self.tramo_b = QComboBox()
        self.posicion_a = self._numero(0.0, 100.0)
        self.posicion_b = self._numero(0.0, 100.0)
        self.panel_perpendicular = QWidget()
        formulario_perpendicular = QFormLayout(self.panel_perpendicular)
        formulario_perpendicular.setContentsMargins(0, 0, 0, 0)
        formulario_perpendicular.addRow("Separación entre apoyos / pórticos:", self.par_porticos)
        formulario_perpendicular.addRow("Tramo receptor en el primer apoyo:", self.tramo_a)
        formulario_perpendicular.addRow("Posición local en el primer tramo (m):", self.posicion_a)
        formulario_perpendicular.addRow("Tramo receptor en el segundo apoyo:", self.tramo_b)
        formulario_perpendicular.addRow("Posición local en el segundo tramo (m):", self.posicion_b)
        self.previsualizacion_muro = QLabel()
        self.previsualizacion_muro.setWordWrap(True)
        formulario_perpendicular.addRow("Reacciones estimadas:", self.previsualizacion_muro)

        for nombre in estructura:
            self.portico.addItem(nombre, nombre)
        if portico_inicial in estructura:
            self.portico.setCurrentIndex(self.portico.findData(portico_inicial))
        for primero, segundo, separacion in self.pares_porticos:
            self.par_porticos.addItem(
                f"{primero} — {segundo} ({separacion:.3f} m)", (primero, segundo, separacion)
            )

        nombre_superficie = "cubierta" if elemento.get("tipo") == "cubierta" else "losa"
        self.ancho_modo = QComboBox()
        self.ancho_modo.addItem(f"Media luz de la {nombre_superficie}", "media_luz")
        self.ancho_modo.addItem(f"Luz completa de la {nombre_superficie}", "luz_completa")
        self.ancho_modo.addItem(
            "Ancho tributario del pórtico (medias separaciones)", "entre_porticos"
        )
        self.ancho_modo.addItem("Ancho tributario personalizado", "manual")
        self.ancho_manual = self._numero(0.01, 100.0)
        self.previsualizacion_ancho = QLabel()
        self.modo_legacy = QComboBox()
        self.modo_legacy.addItem("Sumar; conservar cargas D/L de P00", "sumar")
        self.modo_legacy.addItem("Reemplazar D/L distribuidas antiguas de P00", "reemplazar")
        modo_existente = next((a.get("modo_cargas_previas") for a in self.aplicaciones
                               if a.get("modo_cargas_previas")), "sumar")
        self.modo_legacy.setCurrentIndex(max(0, self.modo_legacy.findData(modo_existente)))
        self.panel_modo_legacy = QWidget()
        fila_modo_legacy = QHBoxLayout(self.panel_modo_legacy)
        fila_modo_legacy.setContentsMargins(0, 0, 0, 0)
        fila_modo_legacy.addWidget(QLabel("Cargas D/L previas de P00:"))
        fila_modo_legacy.addWidget(self.modo_legacy)
        self.boton_agregar = QPushButton("Añadir a la lista")
        self.boton_quitar = QPushButton("Quitar aplicación seleccionada")
        self.boton_agregar.clicked.connect(self._agregar)
        self.boton_quitar.clicked.connect(self._quitar)
        self.portico.currentIndexChanged.connect(self._cargar_tramos)
        self.tramo.currentIndexChanged.connect(self._actualizar_rango)
        self.tramo_a.currentIndexChanged.connect(
            lambda *_args: self._actualizar_rango_receptor(self.tramo_a, self.posicion_a)
        )
        self.tramo_b.currentIndexChanged.connect(
            lambda *_args: self._actualizar_rango_receptor(self.tramo_b, self.posicion_b)
        )
        self.par_porticos.currentIndexChanged.connect(self._cargar_tramos_reaccion)
        self.modo_aplicacion.currentIndexChanged.connect(self._actualizar_modo)
        self.ancho_modo.currentIndexChanged.connect(self._actualizar_ancho)
        self.portico.currentIndexChanged.connect(self._actualizar_ancho)
        self.portico.currentIndexChanged.connect(self._actualizar_pares_porticos)
        self._cargar_tramos()
        self._cargar_tramos_reaccion()
        self._actualizar_ancho()

        caja = QVBoxLayout(self)
        nota = QLabel(
            "Cada aplicación se agrega sin borrar las anteriores. Las posiciones lineales y puntuales "
            "se miden desde el inicio local del tramo. Para que una carga aparezca en Inicio, "
            "agregala a un tramo, guardá las aplicaciones y seleccioná allí el pórtico receptor; "
            "definirla en el catálogo, sin aplicarla, no la dibuja."
        )
        nota.setWordWrap(True)
        caja.addWidget(nota)
        caja.addWidget(self.tabla, 1)
        formulario = QFormLayout()
        if self.es_muro:
            formulario.addRow("Tipo de aplicación:", self.modo_aplicacion)
        formulario.addRow(self.panel_lineal)
        formulario.addRow(self.panel_perpendicular)
        if self.es_superficial:
            formulario.addRow("Criterio de ancho / luz:", self.ancho_modo)
            formulario.addRow("Ancho personalizado (m):", self.ancho_manual)
            formulario.addRow(self.previsualizacion_ancho)
        formulario.addRow(self.panel_modo_legacy)
        caja.addLayout(formulario)
        fila = QHBoxLayout()
        fila.addWidget(self.boton_agregar)
        fila.addWidget(self.boton_quitar)
        caja.addLayout(fila)
        botones = QDialogButtonBox(QDialogButtonBox.StandardButton.Save | QDialogButtonBox.StandardButton.Cancel)
        botones.button(QDialogButtonBox.StandardButton.Save).setText("Guardar aplicaciones")
        botones.accepted.connect(self.accept)
        botones.rejected.connect(self.reject)
        caja.addWidget(botones)
        self._refrescar_tabla()
        self._actualizar_modo()
        self.resize(900, 650 if self.es_muro else 520)

    @staticmethod
    def _numero(minimo: float, maximo: float) -> QDoubleSpinBox:
        control = QDoubleSpinBox()
        control.setRange(minimo, maximo)
        control.setDecimals(2)
        control.setSingleStep(0.1)
        return control

    def _pares_porticos(self, portico_referencia: str = "") -> list[tuple[str, str, float]]:
        referencia = self.estructura.get(portico_referencia, {})
        direccion = cargas.direccion_portico(referencia)
        posiciones = []
        for nombre, datos in self.estructura.items():
            if cargas.direccion_portico(datos) != direccion:
                continue
            try:
                posiciones.append((float(datos["posicion_planta_m"]), nombre))
            except (KeyError, TypeError, ValueError):
                return []
        cantidad = sum(
            cargas.direccion_portico(datos) == direccion
            for datos in self.estructura.values()
        )
        if (
            len(posiciones) != cantidad
            or len({posicion for posicion, _ in posiciones}) != len(posiciones)
        ):
            return []
        posiciones.sort()
        return [
            (primero[1], segundo[1], segundo[0] - primero[0])
            for primero, segundo in zip(posiciones, posiciones[1:])
            if segundo[0] > primero[0]
        ]

    def _actualizar_pares_porticos(self, indice: int) -> None:
        if indice < 0:
            return
        nombre_portico = str(self.portico.currentData() or "")
        pares = self._pares_porticos(nombre_portico)
        self.par_porticos.blockSignals(True)
        self.par_porticos.clear()
        for primero, segundo, separacion in pares:
            self.par_porticos.addItem(
                f"{primero} — {segundo} ({separacion:.3f} m)",
                (primero, segundo, separacion),
            )
        self.par_porticos.blockSignals(False)
        self._cargar_tramos_reaccion()

    def _cargar_tramos_de(self, selector: QComboBox, nombre: str) -> None:
        selector.blockSignals(True)
        selector.clear()
        portico = self.estructura.get(nombre, {})
        for viga in portico.get("vigas", {}).values():
            for tramo in viga.get("tramos", []):
                tid = tramo.get("id")
                longitud = float(tramo.get("longitud_m", 0))
                if tid and longitud > 0:
                    selector.addItem(f"{tid} ({longitud:.2f} m)", (tid, longitud))
        selector.blockSignals(False)

    def _cargar_tramos(self, *_args) -> None:
        self._cargar_tramos_de(self.tramo, self.portico.currentData())
        self._actualizar_rango()

    def _cargar_tramos_reaccion(self, *_args) -> None:
        par = self.par_porticos.currentData()
        if not par:
            self.tramo_a.clear()
            self.tramo_b.clear()
            self.previsualizacion_muro.setText(
                "Definí primero las coordenadas acumuladas de los pórticos."
            )
            return
        self._cargar_tramos_de(self.tramo_a, par[0])
        self._cargar_tramos_de(self.tramo_b, par[1])
        self._actualizar_rango_receptor(self.tramo_a, self.posicion_a)
        self._actualizar_rango_receptor(self.tramo_b, self.posicion_b)
        try:
            elementos = rutas.leer_json(rutas.MATERIALES, {}) or {}
            items = cargas.items_de_elemento(self.carga_id, self.elemento, elementos)
            q = sum(float(item["valor"]) for item in items if item["tipo"] == "D")
            reaccion = q * float(par[2]) / 2
            self.previsualizacion_muro.setText(
                f"q = {q:.2f} kN/m × {par[2]:.3f} m / 2 = {reaccion:.2f} kN en cada pórtico."
            )
        except (KeyError, TypeError, ValueError):
            self.previsualizacion_muro.setText("No se pudo calcular el peso lineal del muro.")

    def _actualizar_rango_receptor(self, selector: QComboBox, control: QDoubleSpinBox) -> None:
        tramo = selector.currentData()
        longitud = float(tramo[1]) if tramo else 0.0
        control.setRange(0.0, longitud)
        control.setValue(longitud / 2)

    def _actualizar_rango(self, *_args) -> None:
        tramo = self.tramo.currentData()
        longitud = float(tramo[1]) if tramo else 0.0
        self.desde.setRange(0.0, longitud)
        self.hasta.setRange(0.0, longitud)
        self.hasta.setValue(longitud)

    def _actualizar_ancho(self, *_args) -> None:
        modo = self.ancho_modo.currentData()
        self.ancho_manual.setEnabled(self.es_superficial and modo == "manual")
        nombre_superficie = "cubierta" if self.elemento.get("tipo") == "cubierta" else "losa"
        if modo == "entre_porticos":
            ancho = cargas.ancho_tributario_entre_porticos(
                self.estructura, str(self.portico.currentData() or "")
            )
            texto = (
                f"Ancho tributario calculado: {ancho:.3f} m, sumando medias separaciones "
                "hacia los pórticos vecinos. Supone carga continua sobre el pórtico; "
                "no identifica un módulo/paño individual."
                if ancho is not None else
                "Asigná coordenadas acumuladas a todos los pórticos para calcular el ancho."
            )
        elif modo == "media_luz":
            texto = f"Se aplicará la mitad de la luz transversal guardada para esta {nombre_superficie}."
        elif modo == "luz_completa":
            texto = f"Se aplicará la luz transversal completa guardada para esta {nombre_superficie}."
        elif modo == "manual":
            texto = "Se usará el ancho tributario ingresado manualmente."
        else:
            texto = ""
        self.previsualizacion_ancho.setText(texto)

    def _actualizar_modo(self, *_args) -> None:
        perpendicular = self.es_muro and self.modo_aplicacion.currentData() == "perpendicular"
        self.panel_lineal.setVisible(not perpendicular)
        self.panel_perpendicular.setVisible(perpendicular)
        self.panel_modo_legacy.setVisible(not perpendicular)

    def _id_aplicacion(self) -> str:
        numero = len(self.todos) + len(self.aplicaciones) + 1
        while any(a.get("id") == f"Aplicacion 0-{numero}" for a in self.todos + self.aplicaciones):
            numero += 1
        return f"Aplicacion 0-{numero}"

    def _agregar(self) -> None:
        perpendicular = self.es_muro and self.modo_aplicacion.currentData() == "perpendicular"
        if perpendicular:
            par = self.par_porticos.currentData()
            apoyo_a, apoyo_b = self.tramo_a.currentData(), self.tramo_b.currentData()
            if not par or not apoyo_a or not apoyo_b:
                QMessageBox.warning(
                    self, "Faltan posiciones o apoyos",
                    "Definí las posiciones de los pórticos y elegí un tramo receptor en cada uno.",
                )
                return
            reacciones = []
            for portico, tramo, posicion in (
                (par[0], apoyo_a, self.posicion_a.value()),
                (par[1], apoyo_b, self.posicion_b.value()),
            ):
                reacciones.append({
                    "portico": portico, "tramo_id": tramo[0],
                    "x_m": posicion, "influencia_m": float(par[2]) / 2,
                })
            self.aplicaciones.append({
                "id": self._id_aplicacion(), "carga_id": self.carga_id,
                "tipo_aplicacion": "muro_perpendicular",
                "reacciones": reacciones, "separacion_m": float(par[2]),
                "modo_cargas_previas": "sumar", "activa": True,
            })
            self._refrescar_tabla()
            return

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
            if modo in ("media_luz", "luz_completa") and luz <= 0:
                superficie = "cubierta" if self.elemento.get("tipo") == "cubierta" else "losa"
                QMessageBox.warning(
                    self, f"Falta la luz de la {superficie}",
                    f"Editá la composición e ingresá la luz transversal de la {superficie}.",
                )
                return
            if modo == "media_luz":
                ancho = luz / 2
            elif modo == "luz_completa":
                ancho = luz
            elif modo == "entre_porticos":
                ancho = cargas.ancho_tributario_entre_porticos(self.estructura, portico)
                if ancho is None:
                    QMessageBox.warning(
                        self, "Faltan posiciones",
                        "Configurá las coordenadas acumuladas de los pórticos antes de usar este ancho.",
                    )
                    return
            else:
                ancho = self.ancho_manual.value()
        else:
            ancho = None
        modos_mismo_tramo = {
            a.get("modo_cargas_previas", "reemplazar")
            for a in self.todos + self.aplicaciones
            if a.get("portico") == portico and a.get("tramo_id") == tramo[0]
        }
        modo_legacy = self.modo_legacy.currentData()
        if len(modos_mismo_tramo) > 1 or (modos_mismo_tramo and modo_legacy not in modos_mismo_tramo):
            QMessageBox.warning(
                self, "Cargas anteriores",
                "Todas las aplicaciones del mismo tramo deben usar el mismo criterio para las cargas de P00.",
            )
            return
        self.aplicaciones.append({
            "id": self._id_aplicacion(), "carga_id": self.carga_id,
            "tipo_aplicacion": "lineal",
            "portico": portico, "tramo_id": tramo[0],
            "x_inicio_m": x0, "x_fin_m": x1,
            "ancho_modo": modo, "ancho_tributario_m": ancho,
            "modo_cargas_previas": modo_legacy, "activa": True,
        })
        self._refrescar_tabla()

    def _quitar(self) -> None:
        fila = self.tabla.currentRow()
        if fila >= 0:
            self.aplicaciones.pop(fila)
            self._refrescar_tabla()

    def _refrescar_tabla(self) -> None:
        self.tabla.setRowCount(len(self.aplicaciones))
        for i, aplicacion in enumerate(self.aplicaciones):
            reacciones = aplicacion.get("reacciones", [])
            if reacciones:
                porticos = " / ".join(r.get("portico", "") for r in reacciones)
                tramos = " / ".join(r.get("tramo_id", "") for r in reacciones)
                posicion = " / ".join(f"{float(r.get('x_m', 0)):.2f}" for r in reacciones)
                valores = (porticos, "Muro perpendicular", tramos, posicion, f"L={aplicacion.get('separacion_m', 0)}")
            else:
                ancho = aplicacion.get("ancho_tributario_m", "—")
                if aplicacion.get("ancho_modo") == "entre_porticos":
                    ancho = cargas.ancho_tributario_entre_porticos(
                        self.estructura, str(aplicacion.get("portico", ""))
                    )
                valores = (
                    aplicacion.get("portico", ""),
                    "Lineal",
                    aplicacion.get("tramo_id", ""),
                    f"{float(aplicacion.get('x_inicio_m', 0)):.2f}–{float(aplicacion.get('x_fin_m', 0)):.2f}",
                    "—" if ancho is None else f"{float(ancho):.3f}",
                )
            for j, valor in enumerate(valores):
                self.tabla.setItem(i, j, QTableWidgetItem(str(valor)))
        self.tabla.resizeColumnsToContents()

    def resultado(self) -> list[dict]:
        return self.todos + self.aplicaciones


class PaginaCargas(QWidget):
    """Catálogo del proyecto y asignación de cargas a intervalos de tramos."""

    def __init__(self, al_guardar=None, portico_actual=None):
        super().__init__()
        self.al_guardar = al_guardar
        self.portico_actual = portico_actual or (lambda: "")
        self.etiqueta = QLabel()
        self.etiqueta.setWordWrap(True)
        self.viento_activo = QCheckBox("Aplicar viento horizontal general al pórtico")
        self.ancho_viento = QDoubleSpinBox()
        self.ancho_viento.setRange(0.0, 100.0)
        self.ancho_viento.setDecimals(2)
        self.ancho_viento.setSingleStep(0.1)
        catalogo_viento = rutas.leer_json(rutas.DATOS_GLOBAL / "viento_cirsoc_102_25.json", {}) or {}
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
        self.tabla = QTableWidget(0, 9)
        self.tabla.setHorizontalHeaderLabels(
            ("ID", "Descripción", "Aplicación", "Transferencia a apoyos", "Categoría",
             "D (kN/m²)", "L (kN/m²)", "Incluir", "Aplicaciones manuales")
        )
        self.tabla.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.tabla.setSelectionMode(QTableWidget.SelectionMode.SingleSelection)
        self.tabla.setAlternatingRowColors(True)
        self.tabla.verticalHeader().setVisible(False)
        self.tabla.horizontalHeader().setStretchLastSection(True)

        self.boton_nueva = QPushButton("Nueva carga")
        self.boton_nueva_losa_planta = QPushButton("Nueva losa desde ejes…")
        self.boton_guardar = QPushButton("Guardar cargas")
        self.boton_editar = QPushButton("Editar composición")
        self.boton_eliminar = QPushButton("Eliminar elemento…")
        self.boton_aplicaciones = QPushButton("Aplicar a tramos…")
        self.boton_editar_composicion = QPushButton("Editar composición del elemento")
        self.boton_quitar_aplicacion = QPushButton("Quitar asignación seleccionada")
        self.boton_editar_composicion.setEnabled(False)
        self.boton_quitar_aplicacion.setEnabled(False)
        self.selector_aplicaciones = QComboBox()
        self.contexto_aplicaciones = QLabel()
        self.contexto_aplicaciones.setWordWrap(True)
        self.tabla_aplicaciones = QTableWidget(0, 7)
        self.tabla_aplicaciones.setHorizontalHeaderLabels(
            ("Pórtico(s)", "Carga", "Tipo", "Barra / tramo", "Desde / x (m)", "Hasta / P (kN)", "Ancho / L (m)")
        )
        self.tabla_aplicaciones.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.tabla_aplicaciones.setSelectionMode(QTableWidget.SelectionMode.SingleSelection)
        self.tabla_aplicaciones.setAlternatingRowColors(True)
        self.tabla_aplicaciones.verticalHeader().setVisible(False)
        self.tabla_aplicaciones.horizontalHeader().setStretchLastSection(True)
        self.etiqueta_aplicaciones = QLabel()
        self.etiqueta_aplicaciones.setWordWrap(True)
        self.boton_disenar_losa = QPushButton("Calcular losa seleccionada")
        self.boton_disenar_losa.setEnabled(False)
        self.boton_editar_apoyos_losa = QPushButton("Editar paño y apoyos…")
        self.boton_editar_apoyos_losa.setEnabled(False)
        self.boton_informe = QPushButton("Generar informe TXT")
        self.boton_abrir = QPushButton("Abrir último informe")
        self.boton_abrir.setEnabled(False)
        self.ultimo_informe: Path | None = None

        self.pestanias_cargas = QTabWidget()
        pagina_elementos = QWidget()
        caja_elementos = QVBoxLayout(pagina_elementos)
        caja_elementos.addWidget(self.etiqueta)
        self.etiqueta_portico = QLabel()
        self.etiqueta_portico.setStyleSheet("font-weight: 600;")
        fila_viento = QHBoxLayout()
        fila_viento.addWidget(self.viento_activo)
        fila_viento.addWidget(QLabel("Ancho tributario del viento (m):"))
        fila_viento.addWidget(self.ancho_viento)
        fila_viento.addStretch(1)
        caja_elementos.addLayout(fila_viento)
        fila_referencia_viento = QHBoxLayout()
        fila_referencia_viento.addWidget(QLabel("Ubicación CIRSOC 102-25:"))
        fila_referencia_viento.addWidget(self.ciudad_viento)
        fila_referencia_viento.addWidget(QLabel("Categoría de riesgo:"))
        fila_referencia_viento.addWidget(self.categoria_riesgo_viento)
        fila_referencia_viento.addWidget(self.velocidad_viento_ref)
        fila_referencia_viento.addStretch(1)
        caja_elementos.addLayout(fila_referencia_viento)
        caja_elementos.addWidget(self.nota_viento_ref)
        caja_elementos.addWidget(self.tabla, 1)
        fila = QHBoxLayout()
        fila.addWidget(self.boton_nueva_losa_planta)
        fila.addWidget(self.boton_nueva)
        fila.addWidget(self.boton_editar)
        fila.addWidget(self.boton_eliminar)
        fila.addWidget(self.boton_guardar)
        caja_elementos.addLayout(fila)
        caja_elementos.addWidget(self.boton_editar_apoyos_losa)
        caja_elementos.addWidget(self.boton_disenar_losa)
        acciones = QHBoxLayout()
        acciones.addWidget(self.boton_informe)
        acciones.addWidget(self.boton_abrir)
        acciones.addStretch(1)
        caja_elementos.addLayout(acciones)
        self.pestanias_cargas.addTab(pagina_elementos, "Definir elementos")

        pagina_aplicaciones = QWidget()
        caja_aplicaciones = QVBoxLayout(pagina_aplicaciones)
        caja_aplicaciones.addWidget(QLabel(
            "Elegí un elemento del catálogo para asignarlo a una o varias barras. "
            "Los muros pueden aplicarse sobre una viga o perpendicularmente entre pórticos; "
            "las losas creadas desde ejes ya envían sus reacciones a las dos vigas de apoyo "
            "con ancho tributario igual a media luz; no las asignes manualmente para evitar duplicarlas."
        ))
        caja_aplicaciones.addWidget(self.etiqueta_portico)
        fila_seleccion = QHBoxLayout()
        fila_seleccion.addWidget(QLabel("Elemento de carga:"))
        fila_seleccion.addWidget(self.selector_aplicaciones, 1)
        fila_seleccion.addWidget(self.boton_editar_composicion)
        fila_seleccion.addWidget(self.boton_aplicaciones)
        caja_aplicaciones.addLayout(fila_seleccion)
        caja_aplicaciones.addWidget(self.contexto_aplicaciones)
        caja_aplicaciones.addWidget(self.etiqueta_aplicaciones)
        caja_aplicaciones.addWidget(self.tabla_aplicaciones, 1)
        acciones_aplicaciones = QHBoxLayout()
        acciones_aplicaciones.addWidget(self.boton_quitar_aplicacion)
        acciones_aplicaciones.addStretch(1)
        caja_aplicaciones.addLayout(acciones_aplicaciones)
        self.pestanias_cargas.addTab(pagina_aplicaciones, "Aplicar y revisar en barras")

        caja = QVBoxLayout(self)
        caja.addWidget(self.pestanias_cargas)

        self.tabla.itemChanged.connect(self._actualizar_estado_activo)
        self.tabla.currentCellChanged.connect(self._actualizar_boton_disenar_losa)
        self.tabla.currentCellChanged.connect(self._actualizar_boton_editar_paño)
        self.tabla_aplicaciones.cellClicked.connect(self._seleccionar_carga_aplicada)
        self.selector_aplicaciones.currentIndexChanged.connect(self._actualizar_boton_aplicaciones)
        self.tabla_aplicaciones.itemSelectionChanged.connect(self._actualizar_boton_quitar)
        self.boton_guardar.clicked.connect(self.guardar)
        self.boton_nueva.clicked.connect(self.nueva)
        self.boton_nueva_losa_planta.clicked.connect(self.nueva_losa_desde_ejes)
        self.boton_editar.clicked.connect(self.editar)
        self.boton_eliminar.clicked.connect(self.eliminar_elemento)
        self.boton_aplicaciones.clicked.connect(self.editar_aplicaciones)
        self.boton_editar_composicion.clicked.connect(self.editar_composicion_aplicaciones)
        self.boton_quitar_aplicacion.clicked.connect(self.quitar_aplicacion_seleccionada)
        self.boton_informe.clicked.connect(self.generar_informe)
        self.boton_abrir.clicked.connect(self.abrir_informe)
        self.boton_disenar_losa.clicked.connect(self.disenar_losa_seleccionada)
        self.boton_editar_apoyos_losa.clicked.connect(self.editar_paño_losa)
        self._actualizar_boton_aplicaciones()
        self.recargar()

    def ir_a_elementos(self) -> None:
        self.pestanias_cargas.setCurrentIndex(0)

    def ir_a_aplicaciones(self) -> None:
        self.pestanias_cargas.setCurrentIndex(1)

    def _actualizar_boton_aplicaciones(self, *_args) -> None:
        nombre = self.selector_aplicaciones.currentData()
        elementos = (cargas.datos_cargas().get("elementos", {}) or {})
        elemento = elementos.get(nombre, {}) if nombre else {}
        aplica_por_apoyos = bool(elemento.get("panel_ejes"))
        disponible = self.selector_aplicaciones.currentIndex() >= 0
        self.boton_aplicaciones.setEnabled(disponible and not aplica_por_apoyos)
        self.boton_aplicaciones.setToolTip(
            "La losa aplica automáticamente sus reacciones a las vigas de apoyo."
            if aplica_por_apoyos
            else "Aplicar el elemento seleccionado a tramos."
        )
        self.boton_editar_composicion.setEnabled(disponible)
        if hasattr(self, "contexto_aplicaciones"):
            portico = str(self.portico_actual() or "").strip() or "todos los pórticos"
            nombre = self.selector_aplicaciones.currentData()
            self.contexto_aplicaciones.setText(
                f"Pórtico activo: {portico}  ·  Elemento seleccionado: {nombre or 'ninguno'}"
            )

    def _actualizar_boton_quitar(self) -> None:
        self.boton_quitar_aplicacion.setEnabled(
            self.tabla_aplicaciones.currentRow() >= 0
        )

    def _seleccionar_carga_aplicada(self, fila: int, _columna: int) -> None:
        item = self.tabla_aplicaciones.item(fila, 1)
        if item is not None:
            indice = self.selector_aplicaciones.findData(
                item.data(Qt.ItemDataRole.UserRole)
            )
            if indice >= 0:
                self.selector_aplicaciones.setCurrentIndex(indice)

    def _cargar_tabla_aplicaciones(self) -> None:
        datos = cargas.datos_cargas()
        estructura = rutas.cargar_estructura()
        elementos = datos.get("elementos", {}) or {}
        nombres_por_id = {
            elemento.get("id"): nombre for nombre, elemento in elementos.items()
        }
        portico = str(self.portico_actual() or "").strip()
        aplicaciones = [
            (indice, aplicacion)
            for indice, aplicacion in enumerate(datos.get("aplicaciones", []))
            if (
                not portico
                or aplicacion.get("portico") == portico
                or any(r.get("portico") == portico for r in aplicacion.get("reacciones", []))
            )
        ]
        self.tabla_aplicaciones.setRowCount(len(aplicaciones))
        por_carga = set()
        for fila, (indice_origen, aplicacion) in enumerate(aplicaciones):
            carga_id = aplicacion.get("carga_id")
            nombre = nombres_por_id.get(carga_id, f"Carga desconocida ({carga_id or '?'})")
            elemento = elementos.get(nombre, {})
            por_carga.add(carga_id)
            reacciones = aplicacion.get("reacciones", [])
            if reacciones:
                relevantes = [
                    r for r in reacciones
                    if not portico or r.get("portico") == portico
                ]
                try:
                    biblioteca = rutas.leer_json(rutas.MATERIALES, {}) or {}
                    peso_lineal = sum(
                        float(item["valor"])
                        for item in cargas.items_de_elemento(nombre, elemento, biblioteca)
                        if item["tipo"] == "D"
                    )
                    reacciones_kN = [
                        peso_lineal * float(r.get("influencia_m", 0.0))
                        for r in relevantes
                    ]
                except (KeyError, TypeError, ValueError):
                    reacciones_kN = []
                valores = (
                    " / ".join(r.get("portico", "") for r in relevantes),
                    nombre,
                    "Muro perpendicular",
                    " / ".join(r.get("tramo_id", "") for r in relevantes),
                    " / ".join(f"{float(r.get('x_m', 0.0)):.2f}" for r in relevantes) or "—",
                    " / ".join(f"{valor:.2f}" for valor in reacciones_kN) or "—",
                    f"{float(aplicacion.get('separacion_m', 0.0)):.2f}" if relevantes else "—",
                )
            else:
                ancho = aplicacion.get("ancho_tributario_m")
                if aplicacion.get("ancho_modo") == "entre_porticos":
                    ancho = cargas.ancho_tributario_entre_porticos(
                        estructura, str(aplicacion.get("portico", ""))
                    )
                valores = (
                    aplicacion.get("portico", ""),
                    nombre,
                    elemento.get("tipo", ""),
                    aplicacion.get("tramo_id", ""),
                    f"{float(aplicacion.get('x_inicio_m', 0.0)):.2f}",
                    f"{float(aplicacion.get('x_fin_m', 0.0)):.2f}",
                    "—" if ancho is None else f"{float(ancho):.2f}",
                )
            for columna, valor in enumerate(valores):
                item = QTableWidgetItem(str(valor))
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                item.setData(Qt.ItemDataRole.UserRole, indice_origen)
                if not aplicacion.get("activa", True):
                    item.setForeground(Qt.GlobalColor.darkGray)
                    item.setToolTip("Aplicación desactivada; no entra en el cálculo.")
                self.tabla_aplicaciones.setItem(fila, columna, item)
        self.tabla_aplicaciones.resizeColumnsToContents()
        self.etiqueta_aplicaciones.setText(
            f"{len(aplicaciones)} aplicación(es) "
            f"{'en ' + portico if portico else 'en toda la obra'} · "
            f"{len(por_carga)} elemento(s) aplicado(s). "
            "Seleccioná una fila para editar sus aplicaciones; las filas grises están desactivadas."
        )
        self._actualizar_boton_quitar()
        self._actualizar_boton_aplicaciones()

    def _actualizar_boton_disenar_losa(self, *_args) -> None:
        nombre = self._fila_nombre()
        elemento = (cargas.datos_cargas().get("elementos", {}) or {}).get(nombre, {}) if nombre else {}
        tipologia = elemento.get("tipologia")
        try:
            disponible = (
                elemento.get("tipo") == "losa"
                and tipologia in ("alivianada", "maciza")
                and float(elemento.get("luz_transversal_m", 0.0) or 0.0) > 0
                and float(elemento.get("ancho_losa_m", 0.0) or 0.0) > 0
            )
        except (TypeError, ValueError):
            disponible = False
        self.boton_disenar_losa.setEnabled(disponible)
        if tipologia == "maciza":
            self.boton_disenar_losa.setText("Calcular solicitaciones de losa maciza")
            ayuda = (
                "Analiza un paño unidireccional simplemente apoyado. "
                "Requiere luz, ancho y hormigón armado con espesor. Si el paño se definió desde ejes, "
                "sus reacciones se transfieren automáticamente a las vigas de apoyo."
            )
        elif tipologia == "alivianada":
            self.boton_disenar_losa.setText("Calcular viguetas de losa alivianada")
            ayuda = "Dimensiona viguetas con esta misma composición."
        elif tipologia == "casetonada":
            self.boton_disenar_losa.setText("Análisis de losa casetonada pendiente")
            ayuda = (
                "La tipología casetonada todavía no tiene un motor de cálculo. Definí en la composición "
                "sus cargas permanentes; si el paño se definió desde ejes, las cargas se transfieren "
                "automáticamente a las vigas de apoyo."
            )
        else:
            self.boton_disenar_losa.setText("Calcular losa seleccionada")
            ayuda = "Elegí una losa y completá la luz y el ancho del paño."
        if not disponible and tipologia in ("alivianada", "maciza"):
            ayuda = "Completá la luz y el ancho del paño para habilitar el cálculo."
        self.boton_disenar_losa.setToolTip(ayuda)

    def _actualizar_boton_editar_paño(self) -> None:
        nombre = self._fila_nombre()
        elemento = (
            (cargas.datos_cargas().get("elementos", {}) or {}).get(nombre, {})
            if nombre else {}
        )
        disponible = elemento.get("tipo") == "losa"
        self.boton_editar_apoyos_losa.setEnabled(disponible)

    def nueva_losa_desde_ejes(self) -> None:
        self._guardar_datos()
        datos = cargas.datos_cargas()
        dialogo = EditorLosaDesdeEjes(
            cargar_ejes(), cargar_niveles(), rutas.cargar_estructura(), self
        )
        if dialogo.exec() != QDialog.DialogCode.Accepted:
            return
        try:
            nombre, geometria = dialogo.resultado()
        except ValueError as exc:
            QMessageBox.warning(self, "Geometría de losa incompleta", str(exc))
            return
        elementos = datos.setdefault("elementos", {})
        if nombre in elementos:
            QMessageBox.warning(
                self, "Nombre existente", "Ya existe un elemento con ese nombre."
            )
            return
        elemento = {
            "id": self._id_nuevo(elementos),
            "tipo": "losa",
            "activo": 1,
            "ancho_tributario_m": 1.0,
            "componentes": [],
            "sobrecarga": "vivienda",
            "viento_activo": 0,
            **geometria,
        }
        editor_cargas = EditorElemento(
            nombre, elemento, rutas.leer_json(rutas.MATERIALES, {}) or {}, self
        )
        editor_cargas.setWindowTitle(f"Composición de cargas — losa {nombre}")
        if editor_cargas.exec() != QDialog.DialogCode.Accepted:
            return
        nuevo_nombre, elemento = editor_cargas.resultado()
        if not nuevo_nombre:
            QMessageBox.warning(self, "Nombre requerido", "La losa necesita un nombre.")
            return
        if nuevo_nombre != nombre and nuevo_nombre in elementos:
            QMessageBox.warning(
                self, "Nombre existente", "Ya existe un elemento con ese nombre."
            )
            return
        elementos[nuevo_nombre] = elemento
        rutas.guardar_json(rutas.CARGAS, datos)
        self.recargar()
        fila = next(
            (
                indice for indice in range(self.tabla.rowCount())
                if self.tabla.item(indice, 1)
                and self.tabla.item(indice, 1).text() == nuevo_nombre
            ),
            -1,
        )
        if fila >= 0:
            self.tabla.selectRow(fila)
            self.tabla.scrollToItem(self.tabla.item(fila, 3))
        self.selector_aplicaciones.setCurrentIndex(
            self.selector_aplicaciones.findData(nuevo_nombre)
        )
        if self.al_guardar:
            self.al_guardar()

    def editar_paño_losa(self) -> None:
        nombre = self._fila_nombre()
        if not nombre:
            return
        self._guardar_datos()
        datos = cargas.datos_cargas()
        elemento = (datos.get("elementos", {}) or {}).get(nombre)
        if not elemento or elemento.get("tipo") != "losa":
            QMessageBox.information(
                self, "Editar paño", "Seleccioná una losa del catálogo."
            )
            return
        dialogo = EditorLosaDesdeEjes(
            cargar_ejes(), cargar_niveles(), rutas.cargar_estructura(),
            self, elemento, nombre,
        )
        if dialogo.exec() != QDialog.DialogCode.Accepted:
            return
        try:
            nuevo_nombre, geometria = dialogo.resultado()
        except ValueError as exc:
            QMessageBox.warning(self, "Geometría de losa incompleta", str(exc))
            return
        elementos = datos.get("elementos", {})
        if nuevo_nombre != nombre and nuevo_nombre in elementos:
            QMessageBox.warning(
                self, "Nombre existente", "Ya existe un elemento con ese nombre."
            )
            return
        elemento.update(geometria)
        if nuevo_nombre != nombre:
            elementos.pop(nombre)
        elementos[nuevo_nombre] = elemento
        rutas.guardar_json(rutas.CARGAS, datos)
        self.recargar()
        if self.al_guardar:
            self.al_guardar()

    def disenar_losa_seleccionada(self) -> None:
        nombre = self._fila_nombre()
        if not nombre:
            return
        datos = cargas.datos_cargas()
        elemento = (datos.get("elementos", {}) or {}).get(nombre)
        if not elemento:
            return
        try:
            if elemento.get("tipologia") == "maciza":
                resultado = losas_macizas.calcular_desde_elemento_carga(nombre, elemento)
                memoria, archivo = losas_macizas.guardar_resultado(resultado)
                detalle = (
                    f"Se calcularon las solicitaciones de {nombre}.\n\n"
                    f"Memoria: {memoria}\nResultado: {archivo}\n\n"
                    "El resultado no dimensiona armaduras ni contempla continuidad."
                )
            else:
                resultado = losas.calcular_desde_elemento_carga(nombre, elemento)
                memoria, computo, archivo = losas.guardar_resultado(resultado)
                detalle = (
                    f"Se calcularon las viguetas de {nombre}.\n\n"
                    f"Memoria: {memoria}\nCómputo: {computo}\nResultado: {archivo}"
                )
        except (KeyError, TypeError, ValueError, FileNotFoundError) as exc:
            QMessageBox.warning(self, "No se pudo calcular la losa", str(exc))
            return
        QMessageBox.information(
            self, "Cálculo de losa", detalle,
        )
        if self.al_guardar:
            self.al_guardar()

    def _actualizar_estado_activo(self, item: QTableWidgetItem) -> None:
        if item.column() == 6:
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
        self.actualizar_portico()
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
        biblioteca = rutas.leer_json(rutas.MATERIALES, {}) or {}
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
            q_d = q_l = "—"
            if tipo in ("losa", "cubierta"):
                try:
                    valores = cargas.valores_superficiales(elemento, biblioteca)
                    q_d = f"{valores['D_kNm2']:.2f}"
                    q_l = f"{valores['L_kNm2']:.2f}"
                except (KeyError, TypeError, ValueError):
                    pass
            aplicaciones_manuales = sum(
                1 for a in aplicaciones if a.get("carga_id") == elemento.get("id")
            )
            panel = elemento.get("panel_ejes")
            apoyos = elemento.get("apoya_en", {}) or {}
            destinos = [
                (lado, apoyos.get(lado) or {})
                for lado in ("izq", "der")
            ]
            destinos = [
                (lado, destino) for lado, destino in destinos
                if destino.get("portico") and destino.get("viga")
            ]
            if elemento.get("tipo") == "losa" and panel:
                resumen_apoyos = " + ".join(
                    f"{destino['portico']} / {destino['viga']}"
                    for _, destino in destinos
                )
                if not elemento.get("activo", True):
                    estado_apoyos = f"Inactiva · {resumen_apoyos or 'sin apoyos'}"
                    color_estado = Qt.GlobalColor.darkGray
                    ayuda_estado = (
                        "La losa está desactivada: sus cargas no se transfieren al análisis."
                    )
                elif len(destinos) == 2:
                    estado_apoyos = f"Aplicada automáticamente · {resumen_apoyos}"
                    color_estado = Qt.GlobalColor.darkGreen
                    ayuda_estado = (
                        "Paño asociado a sus dos vigas. Sus cargas D/L se transfieren "
                        "automáticamente con ancho tributario igual a media luz."
                    )
                else:
                    estado_apoyos = f"Incompleta · {resumen_apoyos or 'sin apoyos'}"
                    color_estado = Qt.GlobalColor.darkYellow
                    ayuda_estado = (
                        "El paño necesita dos apoyos válidos para transferir sus cargas "
                        "automáticamente a los pórticos."
                    )
            else:
                estado_apoyos = "—"
                color_estado = None
                ayuda_estado = ""
            valores = (
                elemento["id"], nombre, forma, estado_apoyos, categorias, q_d, q_l,
                "Sí" if elemento.get("activo") else "No",
                str(aplicaciones_manuales),
            )
            for columna, valor in enumerate(valores):
                item = QTableWidgetItem(str(valor))
                item.setData(Qt.ItemDataRole.UserRole, nombre)
                if columna == 3:
                    item.setToolTip(ayuda_estado)
                    if color_estado is not None:
                        item.setForeground(color_estado)
                    item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                elif columna == 7:
                    item.setToolTip(
                        "Marcada: esta carga entra en los informes y en las aplicaciones al pórtico. "
                        "Desmarcada: se conserva en la obra, pero se omite del cálculo."
                    )
                    item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
                    item.setCheckState(
                        Qt.CheckState.Checked if elemento.get("activo") else Qt.CheckState.Unchecked
                    )
                else:
                    item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self.tabla.setItem(fila, columna, item)
        self.etiqueta.setText(
            f"{len(elementos)} elementos de carga. "
            "Losas y cubiertas convierten kN/m² a kN/m con el ancho tributario de cada aplicación; "
            "sus D y L superficiales se calculan aunque no tengan aplicaciones. "
            "La columna «Transferencia a apoyos» marca las losas de ejes vinculadas: las activas "
            "transfieren D/L automáticamente a sus dos vigas, con ancho tributario igual a media luz. "
            "Muros y encadenados generan carga lineal. La casilla Incluir activa u omite "
            "el elemento en el análisis y en el pórtico. Los cambios se guardan con los botones de abajo; los TXT "
            "anteriores quedan como estaban. Una carga puede tener varias aplicaciones "
            "en distintos tramos. Las reacciones de muros perpendiculares se aplican como cargas puntuales "
            "en los tramos receptores."
        )
        self.tabla.resizeColumnsToContents()
        seleccion_previa = self.selector_aplicaciones.currentData()
        self.selector_aplicaciones.blockSignals(True)
        self.selector_aplicaciones.clear()
        for nombre, elemento in elementos.items():
            self.selector_aplicaciones.addItem(
                f"{nombre} · {elemento.get('tipo', 'carga')}", nombre
            )
        indice = self.selector_aplicaciones.findData(seleccion_previa)
        if indice >= 0:
            self.selector_aplicaciones.setCurrentIndex(indice)
        self.selector_aplicaciones.blockSignals(False)
        self._actualizar_boton_aplicaciones()
        self._cargar_tabla_aplicaciones()
        self._actualizar_boton_disenar_losa()
        self._actualizar_boton_editar_paño()

    def actualizar_portico(self) -> None:
        portico_actual = str(self.portico_actual() or "").strip()
        self.etiqueta_portico.setText(
            f"Pórtico seleccionado: {portico_actual or 'ninguno'}. "
            "Las aplicaciones nuevas se proponen para este pórtico; podés cambiarlo en el editor."
        )
        if hasattr(self, "tabla_aplicaciones"):
            self._cargar_tabla_aplicaciones()

    def _fila_nombre(self) -> str | None:
        fila = self.tabla.currentRow()
        item = self.tabla.item(fila, 0) if fila >= 0 else None
        return item.data(Qt.ItemDataRole.UserRole) if item else None

    def editar(self) -> None:
        nombre = self._fila_nombre()
        if not nombre:
            QMessageBox.information(self, "Cargas", "Elegí primero una carga de la lista.")
            return
        self._editar_elemento(nombre)

    def editar_composicion_aplicaciones(self) -> None:
        nombre = self.selector_aplicaciones.currentData()
        if not nombre:
            QMessageBox.information(self, "Cargas", "Elegí un elemento de carga.")
            return
        self._editar_elemento(nombre)

    def _editar_elemento(self, nombre: str) -> None:
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
        indice = self.selector_aplicaciones.findData(nuevo_nombre)
        if indice >= 0:
            self.selector_aplicaciones.setCurrentIndex(indice)
        if self.al_guardar:
            self.al_guardar()

    def eliminar_elemento(self) -> None:
        nombre = self._fila_nombre()
        if not nombre:
            QMessageBox.information(self, "Cargas", "Elegí primero un elemento del catálogo.")
            return
        datos = cargas.datos_cargas()
        elemento = (datos.get("elementos", {}) or {}).get(nombre)
        if not elemento:
            QMessageBox.warning(self, "Elemento inexistente", "Actualizá la lista e intentá de nuevo.")
            return
        carga_id = elemento.get("id")
        aplicaciones = datos.get("aplicaciones", [])
        cantidad = sum(1 for aplicacion in aplicaciones if aplicacion.get("carga_id") == carga_id)
        texto = f"¿Eliminar «{nombre}» del catálogo?"
        if cantidad:
            texto += (
                f"\n\nTambién se quitarán sus {cantidad} aplicación(es) "
                "de los pórticos."
            )
        respuesta = QMessageBox.question(
            self,
            "Eliminar elemento de carga",
            texto,
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        if respuesta != QMessageBox.StandardButton.Yes:
            return
        self._guardar_datos()
        datos = cargas.datos_cargas()
        datos.get("elementos", {}).pop(nombre, None)
        datos["aplicaciones"] = [
            aplicacion for aplicacion in datos.get("aplicaciones", [])
            if aplicacion.get("carga_id") != carga_id
        ]
        rutas.guardar_json(rutas.CARGAS, datos)
        self.recargar()
        if self.al_guardar:
            self.al_guardar()

    def quitar_aplicacion_seleccionada(self) -> None:
        fila = self.tabla_aplicaciones.currentRow()
        item = self.tabla_aplicaciones.item(fila, 0) if fila >= 0 else None
        if item is None:
            QMessageBox.information(
                self, "Quitar asignación", "Seleccioná primero una aplicación de la tabla."
            )
            return
        indice = item.data(Qt.ItemDataRole.UserRole)
        datos = cargas.datos_cargas()
        aplicaciones = datos.get("aplicaciones", [])
        if not isinstance(indice, int) or indice < 0 or indice >= len(aplicaciones):
            QMessageBox.warning(
                self, "Asignación no encontrada", "Actualizá la lista y volvé a seleccionarla."
            )
            return
        aplicacion = aplicaciones[indice]
        nombre = self.tabla_aplicaciones.item(fila, 1).text()
        portico = str(aplicacion.get("portico", ""))
        tramo = str(aplicacion.get("tramo_id", ""))
        respuesta = QMessageBox.question(
            self,
            "Quitar asignación",
            f"¿Quitar «{nombre}» de {portico} · {tramo}?\n\n"
            "El elemento seguirá en el catálogo y sus otras aplicaciones no cambiarán.",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        if respuesta != QMessageBox.StandardButton.Yes:
            return
        self._guardar_datos()
        datos = cargas.datos_cargas()
        aplicaciones = datos.get("aplicaciones", [])
        if indice >= len(aplicaciones):
            QMessageBox.warning(
                self, "Asignación modificada", "La lista cambió antes de guardar. Actualizá e intentá de nuevo."
            )
            self.recargar()
            return
        aplicaciones.pop(indice)
        rutas.guardar_json(rutas.CARGAS, datos)
        self.recargar()
        if self.al_guardar:
            self.al_guardar()

    def editar_aplicaciones(self, *_args) -> None:
        nombre = self.selector_aplicaciones.currentData()
        if not nombre:
            QMessageBox.information(self, "Aplicar cargas", "Elegí un elemento de carga.")
            return
        self._guardar_datos()
        datos = cargas.datos_cargas()
        elemento = datos.get("elementos", {}).get(nombre, {})
        carga_id = elemento.get("id")
        if not carga_id:
            QMessageBox.warning(self, "Falta ID", "No se pudo identificar la carga del proyecto.")
            return
        editor = EditorAplicaciones(
            carga_id, elemento, rutas.cargar_estructura(), datos.get("aplicaciones", []), self,
            portico_inicial=self.portico_actual(),
        )
        editor.setWindowTitle(
            f"Aplicaciones · {nombre} · pórtico activo: "
            f"{self.portico_actual() or 'ninguno'}"
        )
        if editor.exec() != QDialog.DialogCode.Accepted:
            return
        datos["aplicaciones"] = editor.resultado()
        rutas.guardar_json(rutas.CARGAS, datos)
        self.recargar()
        if self.al_guardar:
            self.al_guardar()

    def nueva(self) -> None:
        tipos = ("cubierta", "muro", "encadenado")
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
            activo = self.tabla.item(fila, 7)
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

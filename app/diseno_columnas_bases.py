"""Formularios de entrada para dimensionar columnas y bases desde la app."""

from __future__ import annotations

from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QLabel,
    QLineEdit,
    QVBoxLayout,
)

from calc import rutas


def respuestas_columnas(parent, portico: str, columnas: dict) -> str | None:
    """Devuelve las opciones por columna; el pórtico se pasa como argumento del proceso."""
    dialogo = QDialog(parent)
    dialogo.setWindowTitle(f"Datos de diseño de columnas · {portico}")
    dialogo.resize(520, 420)
    caja = QVBoxLayout(dialogo)
    caja.addWidget(QLabel(
        "Elegí el tipo de armadura y, para columnas con estribos, "
        "la disposición de las barras longitudinales. En 2 caras se colocan "
        "barras en dos caras opuestas; en 4 caras se distribuyen alrededor "
        "del perímetro, con barras en las esquinas."
    ))
    explicacion_rango = (
        "La disposición en 4 caras puede cambiar la armadura necesaria, pero "
        "no corrige una sección fuera del rango de los diagramas."
    )
    if "C1-8" in columnas:
        explicacion_rango += (
            " En este pórtico, C1-8 quedó fuera porque γ = 0,96 supera el máximo "
            "disponible de 0,90; hay que revisar sus dimensiones o ampliar/verificar "
            "los diagramas."
        )
    caja.addWidget(QLabel(explicacion_rango))
    formulario = QFormLayout()
    controles = {}
    catalogo = rutas.cargar_estructura().get("columnas", {}) or {}
    for columna_id in columnas:
        tipo = QComboBox()
        tipo.addItem("Rectangular con estribos", "1")
        tipo.addItem("Circular con sunchos", "2")
        forma_guardada = str(catalogo.get(columna_id, {}).get("forma_seccion", "rectangular")).lower()
        tipo.setCurrentIndex(1 if forma_guardada in ("circular", "circular con sunchos") else 0)
        seccion = QComboBox()
        seccion.addItem("Barras en 2 caras", "1")
        seccion.addItem("Barras en 4 caras", "2")
        tipo.currentIndexChanged.connect(
            lambda _indice, selector=tipo, caras=seccion:
                caras.setEnabled(selector.currentData() == "1")
        )
        fila = QFormLayout()
        fila.addRow("Armadura:", tipo)
        fila.addRow("Sección:", seccion)
        formulario.addRow(f"{columna_id}:", fila)
        controles[columna_id] = (tipo, seccion)
    caja.addLayout(formulario)
    botones = QDialogButtonBox(
        QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
    )
    botones.accepted.connect(dialogo.accept)
    botones.rejected.connect(dialogo.reject)
    caja.addWidget(botones)
    if dialogo.exec() != QDialog.DialogCode.Accepted:
        return None

    respuestas = []
    for tipo, seccion in controles.values():
        respuestas.append(str(tipo.currentData()))
        if tipo.currentData() == "1":
            respuestas.append(str(seccion.currentData()))
    return "\n".join(respuestas) + "\n"


def datos_bases(parent) -> dict | None:
    """Pide y persiste los parámetros geotécnicos usados por P05."""
    guardados = rutas.leer_json(rutas.TERRENO, {}) or {}
    dialogo = QDialog(parent)
    dialogo.setWindowTitle("Datos de fundación")
    caja = QVBoxLayout(dialogo)
    caja.addWidget(QLabel(
        "Ingresá los datos del estudio geotécnico. Los valores que aparecen "
        "inicialmente son solo referencias: verificá y reemplazalos por los "
        "del informe antes de calcular."
    ))
    formulario = QFormLayout()
    profundidad = QDoubleSpinBox()
    profundidad.setRange(0.1, 20.0)
    profundidad.setDecimals(2)
    profundidad.setSuffix(" m")
    profundidad.setValue(float(guardados.get("profundidad_fundacion_m", 0.8)))
    profundidad.setToolTip(
        "Distancia vertical desde el nivel del terreno hasta el fondo de la zapata."
    )
    q_adm = QLineEdit("" if guardados.get("q_adm_kPa") is None else str(guardados.get("q_adm_kPa")))
    q_adm.setPlaceholderText("Valor del informe geotÃ©cnico, en kPa")
    q_adm.setToolTip(
        "Presión admisible del suelo indicada por el estudio geotécnico, en kPa."
    )
    formulario.addRow("Profundidad hasta el fondo de zapata:", profundidad)
    formulario.addRow("Presión admisible del suelo (q_adm):", q_adm)
    caja.addLayout(formulario)
    caja.addWidget(QLabel(
        "Alcance actual: predimensiona la zapata por presión del suelo y propone "
        "armadura mínima. No verifica resistencia del hormigón, flexión ni "
        "punzonado; por eso todavía no solicita f'c."
    ))
    botones = QDialogButtonBox(
        QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
    )
    botones.accepted.connect(dialogo.accept)
    botones.rejected.connect(dialogo.reject)
    caja.addWidget(botones)
    if dialogo.exec() != QDialog.DialogCode.Accepted:
        return None

    datos = dict(guardados)
    try:
        q_adm_kpa = float(q_adm.text().strip().replace(",", "."))
    except ValueError:
        QMessageBox.warning(parent, "Falta q_adm", "IngresÃ¡ la presiÃ³n admisible del informe geotÃ©cnico.")
        return None
    if q_adm_kpa <= 0:
        QMessageBox.warning(parent, "q_adm invÃ¡lido", "q_adm debe ser mayor que cero.")
        return None
    datos.update({
        "profundidad_fundacion_m": profundidad.value(),
        "q_adm_kPa": q_adm_kpa,
    })
    rutas.guardar_json(rutas.TERRENO, datos)
    return datos

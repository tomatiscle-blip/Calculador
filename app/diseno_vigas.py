"""Entrada gráfica de secciones y ejecución de P02 sobre el motor actual."""

from __future__ import annotations

from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFormLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QVBoxLayout,
)

from calc import diseno_vigas, rutas


def dimensionar_desde_app(parent, portico: str):
    estructura = rutas.cargar_estructura()
    datos_portico = estructura.get(portico, {})
    vigas = datos_portico.get("vigas", {})
    if not vigas:
        QMessageBox.information(parent, "Vigas", f"{portico} no tiene vigas cargadas.")
        return None

    dialogo = QDialog(parent)
    dialogo.setWindowTitle(f"Datos de diseño de vigas · {portico}")
    dialogo.resize(440, 180)
    dialogo.setStyleSheet(
        "QLineEdit, QComboBox {"
        " color: #202124; background-color: #ffffff;"
        " border: 1px solid #8a8f98; border-radius: 4px; padding: 5px;"
        " selection-color: #202124; selection-background-color: #cfe3ff;"
        "}"
        "QLineEdit:focus, QComboBox:focus { border: 1px solid #2563eb; }"
        "QComboBox QAbstractItemView {"
        " color: #202124; background-color: #ffffff;"
        " selection-color: #202124; selection-background-color: #cfe3ff;"
        "}"
    )
    caja = QVBoxLayout(dialogo)
    caja.addWidget(QLabel("Ingresá el ancho y el hormigón para cada viga. La altura se predimensiona con el criterio actual de P02."))
    formulario = QFormLayout()
    controles = {}
    for viga_id, viga in vigas.items():
        ancho = QLineEdit("" if viga.get("b_cm") is None else str(viga["b_cm"]))
        ancho.setPlaceholderText("b en cm")
        fc = QComboBox()
        fc.addItem("Elegir…", None)
        for valor in (20, 25, 30):
            fc.addItem(f"H{valor}", valor)
        seleccionado = fc.findData(viga.get("fc_MPa"))
        if seleccionado >= 0:
            fc.setCurrentIndex(seleccionado)
        fila = QVBoxLayout()
        fila.addWidget(QLabel("Ancho b (cm)"))
        fila.addWidget(ancho)
        fila.addWidget(QLabel("Hormigón"))
        fila.addWidget(fc)
        formulario.addRow(f"{viga_id}:", fila)
        controles[viga_id] = (ancho, fc)
    caja.addLayout(formulario)
    botones = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
    botones.accepted.connect(dialogo.accept)
    botones.rejected.connect(dialogo.reject)
    caja.addWidget(botones)
    if dialogo.exec() != QDialog.DialogCode.Accepted:
        return None

    for viga_id, (ancho, fc) in controles.items():
        try:
            b_cm = float(ancho.text().replace(",", "."))
        except ValueError:
            QMessageBox.warning(parent, "Ancho inválido", f"Ingresá un ancho b válido para {viga_id}.")
            return None
        fc_mpa = fc.currentData()
        if b_cm <= 0 or fc_mpa not in (20, 25, 30):
            QMessageBox.warning(parent, "Faltan datos", f"Completá b y seleccioná H20, H25 o H30 para {viga_id}.")
            return None
        vigas[viga_id]["b_cm"] = b_cm
        vigas[viga_id]["fc_MPa"] = fc_mpa
    rutas.guardar_json(rutas.ESTRUCTURA, estructura)

    try:
        salida, avisos = diseno_vigas.dimensionar(portico)
    except Exception as exc:
        QMessageBox.critical(parent, "No se pudo dimensionar", str(exc))
        return None
    QMessageBox.information(
        parent,
        "Dimensionado de vigas",
        f"P02 terminó y guardó sus resultados en:\n{salida}"
        + ("\n\nA revisar:\n• " + "\n• ".join(avisos) if avisos else ""),
    )
    return salida, avisos

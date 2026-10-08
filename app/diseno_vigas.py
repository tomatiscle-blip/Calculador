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
    dialogo.resize(540, 560)
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
    caja.addWidget(QLabel(
        "Elegí el criterio de flecha inmediata y ajustá las secciones. "
        "Si no cumple, probá aumentar h y volver a dimensionar. La flecha de larga duración todavía no se verifica."
    ))
    formulario = QFormLayout()
    tipo_flecha = QComboBox()
    tipo_flecha.addItem("Viga de entrepiso · inmediata por sobrecarga de uso · δ ≤ ℓ/360", "piso")
    tipo_flecha.addItem("Viga de cubierta plana · inmediata por carga variable · δ ≤ ℓ/180", "cubierta")
    formulario.addRow("Flecha inmediata:", tipo_flecha)
    caja.addWidget(QLabel(
        "Para elementos no estructurales apoyados o unidos a la viga, CIRSOC distingue luego la flecha posterior a su colocación: "
        "ℓ/480 si son sensibles al daño y ℓ/240 si no lo son. Ese cálculo de larga duración todavía no está disponible."
    ))
    caja.addWidget(QLabel("Este cálculo corresponde a vigas; no usarlo para dimensionar losas."))
    caja.addLayout(formulario)
    formulario = QFormLayout()
    controles = {}
    for viga_id, viga in vigas.items():
        ancho = QLineEdit("" if viga.get("b_cm") is None else str(viga["b_cm"]))
        ancho.setPlaceholderText("b en cm")
        altura = QLineEdit("" if viga.get("h_cm") is None else str(viga["h_cm"]))
        altura.setPlaceholderText("vacío = predimensionado automático")
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
        fila.addWidget(QLabel("Peralte h (cm; vacío para automático)"))
        fila.addWidget(altura)
        fila.addWidget(QLabel("Hormigón"))
        fila.addWidget(fc)
        formulario.addRow(f"{viga_id}:", fila)
        controles[viga_id] = (ancho, altura, fc)
    caja.addLayout(formulario)
    botones = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
    botones.accepted.connect(dialogo.accept)
    botones.rejected.connect(dialogo.reject)
    caja.addWidget(botones)
    if dialogo.exec() != QDialog.DialogCode.Accepted:
        return None

    secciones = {}
    for viga_id, (ancho, altura, fc) in controles.items():
        try:
            b_cm = float(ancho.text().replace(",", "."))
        except ValueError:
            QMessageBox.warning(parent, "Ancho inválido", f"Ingresá un ancho b válido para {viga_id}.")
            return None
        fc_mpa = fc.currentData()
        if b_cm <= 0 or fc_mpa not in (20, 25, 30):
            QMessageBox.warning(parent, "Faltan datos", f"Completá b y seleccioná H20, H25 o H30 para {viga_id}.")
            return None
        try:
            h_cm = float(altura.text().replace(",", ".")) if altura.text().strip() else None
        except ValueError:
            QMessageBox.warning(parent, "Peralte inválido", f"Ingresá un h válido para {viga_id} o dejalo vacío.")
            return None
        if h_cm is not None and h_cm <= 0:
            QMessageBox.warning(parent, "Peralte inválido", f"El peralte h debe ser mayor que cero para {viga_id}.")
            return None
        secciones[viga_id] = {"b_cm": b_cm, "fc_MPa": fc_mpa}
        if h_cm is not None:
            secciones[viga_id]["h_cm"] = h_cm

    secciones_modificadas = False
    for viga_id, seccion in secciones.items():
        for clave in ("b_cm", "fc_MPa"):
            if vigas[viga_id].get(clave) != seccion[clave]:
                vigas[viga_id][clave] = seccion[clave]
                secciones_modificadas = True
        if "h_cm" in seccion:
            if vigas[viga_id].get("h_cm") != seccion["h_cm"]:
                vigas[viga_id]["h_cm"] = seccion["h_cm"]
                secciones_modificadas = True

    try:
        if secciones_modificadas:
            rutas.guardar_json(rutas.ESTRUCTURA, estructura)
        salida, avisos = diseno_vigas.dimensionar(portico, secciones, tipo_flecha.currentData())
    except Exception as exc:
        QMessageBox.critical(parent, "No se pudo dimensionar", str(exc))
        return None
    QMessageBox.information(
        parent,
        "Dimensionado de vigas",
        f"P02 terminó y guardó sus resultados en:\n{salida}\n\n"
        "Si la flecha inmediata no cumple, probá aumentar h y volver a dimensionar. Se actualizan flexión y armado. "
        "La flecha de larga duración y el efecto de cargas puntuales siguen pendientes."
        + ("\n\nA revisar:\n• " + "\n• ".join(avisos) if avisos else ""),
    )
    return salida, avisos

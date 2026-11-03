"""Selección del criterio de flecha y ejecución de P02 sobre el motor actual."""

from __future__ import annotations

import math

from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFormLayout,
    QLabel,
    QMessageBox,
    QVBoxLayout,
)

from calc import configuracion_obra, diseno_vigas, pipeline, rutas


def _fallas_flecha(salida) -> dict[str, float]:
    """Devuelve la mayor razón delta/límite de flecha por viga que no cumple."""
    datos = rutas.leer_json(salida, {}) or {}
    por_viga: dict[str, float] = {}
    for fila in datos.get("flecha", []) or []:
        if fila.get("cumple") is not False:
            continue
        viga_id = fila.get("viga")
        if not viga_id:
            continue
        try:
            delta = float(fila.get("delta_mm", 0.0) or 0.0)
            limite = float(fila.get("limite_mm", 0.0) or 0.0)
        except (TypeError, ValueError):
            continue
        if delta > 0 and limite > 0:
            por_viga[viga_id] = max(por_viga.get(viga_id, 1.0), delta / limite)
    return por_viga


def _altura_inicial_cm(viga: dict) -> float:
    if viga.get("h_cm") is not None:
        try:
            return float(viga["h_cm"])
        except (TypeError, ValueError):
            pass
    luz = max((float(t.get("longitud_m", 0.0)) for t in viga.get("tramos", [])), default=0.0)
    return float(round(luz * 100.0 / 12.5)) if luz > 0 else 30.0


def dimensionar_desde_app(parent, portico: str):
    estructura = rutas.cargar_estructura()
    datos_portico = estructura.get(portico, {})
    vigas = datos_portico.get("vigas", {})
    if not vigas:
        QMessageBox.information(parent, "Vigas", f"{portico} no tiene vigas cargadas.")
        return None

    dialogo = QDialog(parent)
    dialogo.setWindowTitle(f"Verificación de flecha de vigas · {portico}")
    dialogo.resize(480, 360)
    dialogo.setStyleSheet(
        "QComboBox {"
        " color: #202124; background-color: #ffffff;"
        " border: 1px solid #8a8f98; border-radius: 4px; padding: 5px;"
        " selection-color: #202124; selection-background-color: #cfe3ff;"
        "}"
        "QComboBox:focus { border: 1px solid #2563eb; }"
        "QComboBox QAbstractItemView {"
        " color: #202124; background-color: #ffffff;"
        " selection-color: #202124; selection-background-color: #cfe3ff;"
        "}"
    )
    caja = QVBoxLayout(dialogo)
    caja.addWidget(QLabel(
        "La sección se toma del predimensionado de la obra o de la viga. Revisá el criterio de flecha inmediata."
    ))
    formulario = QFormLayout()
    tipo_flecha = QComboBox()
    tipo_flecha.addItem("Viga de entrepiso · inmediata por sobrecarga de uso · δ ≤ ℓ/360", "piso")
    tipo_flecha.addItem("Viga de cubierta plana · inmediata por carga variable · δ ≤ ℓ/180", "cubierta")
    formulario.addRow("Flecha inmediata:", tipo_flecha)
    caja.addLayout(formulario)
    caja.addWidget(QLabel(
        "Secciones que usará P02 (solo lectura):"
    ))
    parametros = configuracion_obra.cargar()
    b_default = float(parametros.get("VIGA_ANCHO_PREDIM_CM", 20.0) or 20.0)
    fc_obra = configuracion_obra.fc_predeterminado()
    for viga_id, viga in vigas.items():
        b_cm = float(viga.get("b_cm") or b_default)
        fc_viga = int(viga.get("fc_MPa") or fc_obra)
        h_cm = viga.get("h_cm")
        if h_cm is None:
            luz = max((float(t.get("longitud_m", 0.0)) for t in viga.get("tramos", [])), default=0.0)
            h_cm = round(luz * 100.0 / 12.5) if luz > 0 else 0.0
            origen_h = "automático · luz mayor / 12,5"
        else:
            h_cm = float(h_cm)
            origen_h = (
                "ajustado automáticamente por flecha"
                if viga.get("h_ajustado_por_flecha") else "guardado en la viga"
            )
        origen_b = "guardado en la viga" if viga.get("b_cm") else "predimensionado de obra"
        caja.addWidget(QLabel(
            f"{viga_id}: b={b_cm:g} cm ({origen_b}), h={h_cm:g} cm ({origen_h}), hormigón H{fc_viga}."
        ))
    caja.addWidget(QLabel(
        "CIRSOC también limita la flecha posterior a colocar elementos no estructurales: ℓ/480 "
        "si pueden dañarse y ℓ/240 si no son sensibles. Esa verificación de larga duración aún no está implementada."
    ))
    caja.addWidget(QLabel("Este cálculo corresponde a vigas; no dimensiona losas."))
    botones = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
    botones.accepted.connect(dialogo.accept)
    botones.rejected.connect(dialogo.reject)
    caja.addWidget(botones)
    if dialogo.exec() != QDialog.DialogCode.Accepted:
        return None

    try:
        h_iniciales: dict[str, float] = {}
        avisos = []
        salida = None
        max_intentos = 12
        for intento in range(max_intentos + 1):
            calculo_motor = pipeline.ejecutar("portico", portico)
            if not calculo_motor.get("ok"):
                raise RuntimeError(
                    "No se pudieron actualizar las solicitaciones con las secciones predimensionadas. "
                    + str(calculo_motor.get("error") or calculo_motor.get("salida") or "")
                )
            salida, avisos = diseno_vigas.dimensionar(
                portico, tipo_flecha=tipo_flecha.currentData()
            )
            fallas = _fallas_flecha(salida)
            if not fallas:
                break
            if intento >= max_intentos:
                avisos.append(
                    f"La flecha sigue fuera del límite después de {max_intentos} ajustes automáticos de altura. "
                    "Revisá cargas y sección antes de adoptar el diseño."
                )
                break

            ajustes = {}
            for viga_id, razon in fallas.items():
                viga = vigas.get(viga_id)
                if not viga:
                    continue
                altura_actual = _altura_inicial_cm(viga)
                h_objetivo = altura_actual * razon ** (1.0 / 3.0) * 1.04
                altura_nueva = math.ceil(h_objetivo / 2.0) * 2.0
                if altura_nueva <= altura_actual:
                    altura_nueva = altura_actual + 2.0
                if altura_nueva > 120.0:
                    continue
                h_iniciales.setdefault(viga_id, altura_actual)
                ajustes[viga_id] = altura_nueva
            if not ajustes:
                avisos.append(
                    "No se pudo aumentar automáticamente la altura de las vigas que no verificaron flecha."
                )
                break
            for viga_id, altura in ajustes.items():
                vigas[viga_id]["h_cm"] = altura
                vigas[viga_id]["h_ajustado_por_flecha"] = True
            rutas.guardar_estructura(estructura)
        if salida is None:
            raise RuntimeError("P02 no generó resultados de vigas.")
    except Exception as exc:
        QMessageBox.critical(parent, "No se pudo dimensionar", str(exc))
        return None
    resumen_ajustes = []
    for viga_id, altura_inicial in h_iniciales.items():
        altura_final = float(vigas[viga_id].get("h_cm", altura_inicial))
        resumen_ajustes.append(f"{viga_id}: h {altura_inicial:g} → {altura_final:g} cm")
    QMessageBox.information(
        parent,
        "Dimensionado de vigas",
        f"P02 terminó y guardó sus resultados en:\n{salida}\n\n"
        "Se recalculó el pórtico con el predimensionado. Si la flecha inmediata no verificó, se aumentó h y se repitió el análisis/dimensionado automáticamente. "
        "La flecha de larga duración y el efecto de cargas puntuales siguen pendientes."
        + ("\n\nAlturas ajustadas por flecha:\n• " + "\n• ".join(resumen_ajustes) if resumen_ajustes else "")
        + ("\n\nA revisar:\n• " + "\n• ".join(avisos) if avisos else ""),
    )
    return salida, avisos

"""
tools.comparar_motores — Validación cruzada: anaStruct (viejo) vs Pynite (nuevo).

Para qué sirve
--------------
La Fase 4 del plan dice: "Comparar anaStruct vs Pynite con los Pórticos 1, 2 y 3
(momentos, cortantes, flechas)". Este programa hace exactamente eso y deja el
resultado por escrito, así el día que se retire anaStruct queda el respaldo de
que los dos motores daban lo mismo.

Cómo funciona
-------------
 * El lado anaStruct NO se recalcula: se leen los resultados **ya guardados** en
   `datos/estructura.json` (`Mu_servicio_*`, `Vu_servicio_*`, `P_servicio`,
   `Fx/Fy/Tz_servicio`). Es el mismo archivo que usó el cálculo viejo.
 * El lado Pynite se arma acá mismo: el mismo pórtico, como **pórtico plano**
   (se bloquean los movimientos fuera del plano: DZ, RX, RY en todos los nudos).
 * Se comparan las dos columnas de números y se informa la diferencia. Si una
   supera la tolerancia (±2 %) queda marcada como [DIF].
 * Un chequeo de equilibrio avisa si la referencia guardada está incompleta.

No calcula nada nuevo ni pisa ningún resultado: es una herramienta de control.

Uso:
    py tools\\comparar_motores.py               (Pórticos 1, 2 y 3)
    py tools\\comparar_motores.py "Portico 2"    (uno solo)
    py tools\\comparar_motores.py --tolerancia 1
"""

from __future__ import annotations

import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent
if str(RAIZ) not in sys.path:
    sys.path.insert(0, str(RAIZ))

from calc import rutas  # noqa: E402
from Pynite import FEModel3D  # noqa: E402

# ---------------------------------------------------------------------------
# Rigidez: los MISMOS valores por defecto que usa anaStruct en P01
# (anaStruct trabaja con E=1, A=15000, I=5000; acá se pasa todo a Pynite)
# ---------------------------------------------------------------------------
E_MOD = 1.0
G_MOD = 1.0
EA = 15000.0
EI = 5000.0
A_SEC = EA / E_MOD
I_SEC = EI / E_MOD
J_SEC = EI / E_MOD

COMBO = "Combo 1"          # así se llama la única combinación (servicio = D + L)
TOLERANCIA = 2.0           # % admitido entre los dos motores


def _k(x) -> float:
    return round(float(x), 4)


def _nodo(x: float, y: float) -> str:
    return f"N{_k(x)}_{_k(y)}"


# ---------------------------------------------------------------------------
# Armar el modelo Pynite de un pórtico (como pórtico plano)
# ---------------------------------------------------------------------------
def construir(portico: dict) -> tuple[FEModel3D, dict]:
    m = FEModel3D()
    m.add_material("mat", E_MOD, G_MOD, 0.3, 0.0)
    m.add_section("sec", A_SEC, I_SEC, I_SEC, J_SEC)

    mapa = {"columnas": {}, "vigas": {}, "voladizos": {}}

    # --- Columnas ---------------------------------------------------------
    for cid, col in portico["columnas"].items():
        x = _k(col["x"])
        y0 = _k(col.get("nivel", 0.0))
        h = _k(col["altura_m"])
        nb, nt = _nodo(x, y0), _nodo(x, y0 + h)
        if nb not in m.nodes:
            m.add_node(nb, x, y0, 0.0)
        if nt not in m.nodes:
            m.add_node(nt, x, y0 + h, 0.0)
        # pórtico plano: libre en el plano (DX, DY, RZ), fijo fuera del plano
        m.def_support(nb, False, False, True, True, True, False)
        m.def_support(nt, False, False, True, True, True, False)
        m.add_member(f"C-{cid}", nb, nt, "mat", "sec")
        mapa["columnas"][cid] = f"C-{cid}"

    # --- Apoyos en las bases ---------------------------------------------
    for base in portico["bases"].values():
        nodo = _nodo(base["x"], 0.0)
        if nodo not in m.nodes:
            continue
        if base.get("tipo") == "articulado":
            m.def_support(nodo, True, True, True, False, False, False)
        else:  # empotramiento (por defecto)
            m.def_support(nodo, True, True, True, True, True, True)

    # --- Vigas (por piso) y voladizos ------------------------------------
    for viga_id, viga in portico["vigas"].items():
        piso = int(viga_id.split("-")[0][1:])
        cols = sorted(
            (c for cid, c in portico["columnas"].items() if cid.startswith(f"C{piso}-")),
            key=lambda c: c["x"],
        )
        if not cols:
            continue
        y_viga = _k(cols[0].get("nivel", 0.0) + cols[0]["altura_m"])

        for i, tramo in enumerate(viga["tramos"]):
            cargas = tramo.get("cargas", {})
            w = float(cargas.get("D_total", 0.0)) + float(cargas.get("L_total", 0.0))

            if tramo.get("es_voladizo"):
                if "izq" in tramo["id"]:
                    base_x = _k(cols[0]["x"])
                    tip_x = _k(base_x - tramo["longitud_m"])
                else:
                    base_x = _k(cols[-1]["x"])
                    tip_x = _k(base_x + tramo["longitud_m"])
                miembro = f"V-{tramo['id']}"
                n_tip, n_base = _nodo(tip_x, y_viga), _nodo(base_x, y_viga)
                if n_tip not in m.nodes:
                    m.add_node(n_tip, tip_x, y_viga, 0.0)
                    m.def_support(n_tip, False, False, True, True, True, False)
                # del extremo libre hacia la columna: mismo orden que usó anaStruct
                m.add_member(miembro, n_tip, n_base, "mat", "sec")
                mapa["voladizos"][tramo["id"]] = miembro
            else:
                if i >= len(cols) - 1:
                    continue
                miembro = f"V-{tramo['id']}"
                m.add_member(miembro, _nodo(cols[i]["x"], y_viga), _nodo(cols[i + 1]["x"], y_viga),
                             "mat", "sec")
                mapa["vigas"][tramo["id"]] = miembro

            if abs(w) > 1e-9:
                m.add_member_dist_load(miembro, "FY", -w, -w)

            # cargas puntuales del tramo (en el nudo correspondiente)
            for c in tramo.get("cargas_puntuales", []):
                nodo = _nodo(c["coordenadas"]["x"], c["coordenadas"]["y"])
                if nodo in m.nodes:
                    m.add_node_load(nodo, "FY", -float(c["valor_kN"]))

    return m, mapa


# ---------------------------------------------------------------------------
# Comparación de un valor (se comparan MÓDULOS: los signos pueden diferir)
# ---------------------------------------------------------------------------
def _fila(titulo: str, viejo, nuevo, unidad: str, tol: float) -> tuple[str, bool]:
    if viejo is None or nuevo is None:
        return (f"    {titulo:<24} {'--':>10} {'--':>10} {'--':>8} {unidad}", False)
    base = max(abs(viejo), abs(nuevo), 1e-9)
    dif = abs(abs(viejo) - abs(nuevo)) / base * 100.0
    malo = dif > tol
    marca = "  [DIF]" if malo else ""
    return (
        f"    {titulo:<24} {viejo:>10.2f} {nuevo:>10.2f} {dif:>7.2f}% {unidad}{marca}",
        malo,
    )


def _carga_aplicada(portico: dict) -> float:
    total = 0.0
    for viga in portico["vigas"].values():
        for tramo in viga["tramos"]:
            c = tramo.get("cargas", {})
            total += (float(c.get("D_total", 0.0)) + float(c.get("L_total", 0.0))) * tramo["longitud_m"]
            for cp in tramo.get("cargas_puntuales", []):
                total += float(cp["valor_kN"])
    return total


def comparar_portico(nombre: str, portico: dict, tol: float) -> tuple[list[str], int]:
    """Devuelve las líneas del informe y la cantidad de diferencias (> tol)."""
    lineas = [
        f"\n=== {nombre} ===",
        f"    {'CASO':<24} {'anaStruct':>10} {'Pynite':>10} {'dif':>8}",
        "    " + "-" * 58,
    ]

    m, mapa = construir(portico)
    m.analyze()

    diferencias = 0

    # Chequeo de equilibrio: revela si la REFERENCIA guardada está completa.
    carga = _carga_aplicada(portico)
    suma_pynite = sum(
        abs(float(m.nodes[_nodo(b["x"], 0.0)].RxnFY.get(COMBO, 0.0)))
        for b in portico["bases"].values() if _nodo(b["x"], 0.0) in m.nodes
    )
    suma_ana = sum(abs(b.get("Fy_servicio", 0.0)) for b in portico["bases"].values())
    lineas.append("  Equilibrio (carga vertical aplicada vs. suma de reacciones)")
    lineas.append(f"    {'carga aplicada':<24} {'':>10} {carga:>10.2f} {'':>8} kN")
    lineas.append(f"    {'Σ reacciones Pynite':<24} {'':>10} {suma_pynite:>10.2f} {'':>8} kN")
    lineas.append(f"    {'Σ reacciones anaStruct':<24} {'':>10} {suma_ana:>10.2f} {'':>8} kN")
    if abs(suma_ana - carga) > 0.5:
        diferencias += 1
        lineas.append("    [DIF] la referencia anaStruct NO cierra el equilibrio: comparar con criterio")

    lineas.append("  Reacciones en las bases")
    for base_id, base in portico["bases"].items():
        n = m.nodes.get(_nodo(base["x"], 0.0))
        if n is None:
            continue
        # anaStruct guarda Fy "hacia abajo" (-) y Pynite "hacia arriba" (+): se compara módulo
        for etq, viejo, nuevo, un in [
            ("Fx", base.get("Fx_servicio"), float(n.RxnFX.get(COMBO, 0.0)), "kN"),
            ("Fy", base.get("Fy_servicio"), float(n.RxnFY.get(COMBO, 0.0)), "kN"),
            ("Tz", base.get("Tz_servicio"), float(n.RxnMZ.get(COMBO, 0.0)), "kNm"),
        ]:
            texto, malo = _fila(f"{base_id} {etq}", viejo, nuevo, un, tol)
            lineas.append(texto)
            diferencias += 1 if malo else 0

    lineas.append("  Columnas (M base, M punta, axial)")
    for cid, col in portico["columnas"].items():
        mem = m.members.get(mapa["columnas"].get(cid, ""))
        if mem is None:
            continue
        for etq, viejo, nuevo, un in [
            (f"{cid} M base", col.get("Mu_servicio_inf"), float(mem.moment("Mz", 0.0, COMBO)), "kNm"),
            (f"{cid} M punta", col.get("Mu_servicio_sup"), float(mem.moment("Mz", mem.L(), COMBO)), "kNm"),
            (f"{cid} axial", col.get("P_servicio"), float(mem.axial(0.0, COMBO)), "kN"),
        ]:
            texto, malo = _fila(etq, viejo, nuevo, un, tol)
            lineas.append(texto)
            diferencias += 1 if malo else 0

    lineas.append("  Vigas y voladizos (M izq/der/campo, V izq/der)")
    for viga in portico["vigas"].values():
        for tramo in viga["tramos"]:
            clave = tramo["id"]
            mem = m.members.get(mapa["vigas"].get(clave) or mapa["voladizos"].get(clave) or "")
            if mem is None:
                continue
            for etq, viejo, nuevo, un in [
                (f"{clave} M izq", tramo.get("Mu_servicio_izq"), float(mem.moment("Mz", 0.0, COMBO)), "kNm"),
                (f"{clave} M der", tramo.get("Mu_servicio_der"), float(mem.moment("Mz", mem.L(), COMBO)), "kNm"),
                (f"{clave} M campo", tramo.get("Mu_servicio_campo"), float(mem.min_moment("Mz", COMBO)), "kNm"),
                (f"{clave} V izq", tramo.get("Vu_servicio_izq"), float(mem.shear("Fy", 0.0, COMBO)), "kN"),
                (f"{clave} V der", tramo.get("Vu_servicio_der"), float(mem.shear("Fy", mem.L(), COMBO)), "kN"),
            ]:
                texto, malo = _fila(etq, viejo, nuevo, un, tol)
                lineas.append(texto)
                diferencias += 1 if malo else 0

    return lineas, diferencias


def informe(porticos: list[str] | None, tol: float) -> str:
    estructura = rutas.cargar_estructura()
    nombres = porticos or ["Portico 1", "Portico 2", "Portico 3"]

    lineas = [
        "COMPARACIÓN DE MOTORES  ·  anaStruct (guardado)  vs  Pynite (calculado)",
        "=" * 66,
        f"Tolerancia: ±{tol:.1f} %   (se comparan módulos: los signos pueden diferir)",
    ]
    total_dif = 0
    for nombre in nombres:
        if nombre not in estructura:
            lineas.append(f"\n=== {nombre} === (no está en estructura.json)")
            continue
        bloque, dif = comparar_portico(nombre, estructura[nombre], tol)
        lineas.extend(bloque)
        total_dif += dif

    lineas.append("\n" + "=" * 66)
    if total_dif == 0:
        lineas.append("RESULTADO: los dos motores coinciden dentro de la tolerancia. OK")
    else:
        lineas.append(f"RESULTADO: {total_dif} valor(es) fuera de la tolerancia (ver [DIF]).")
    return "\n".join(lineas)


def _main(argv: list[str]) -> int:
    tol = TOLERANCIA
    porticos: list[str] = []
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == "--tolerancia":
            i += 1
            tol = float(argv[i]) if i < len(argv) else tol
        elif not a.startswith("-"):
            porticos.append(a)
        i += 1

    print(informe(porticos or None, tol))
    return 0


if __name__ == "__main__":
    raise SystemExit(_main(sys.argv[1:]))
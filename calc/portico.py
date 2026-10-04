"""
calc.portico — El MOTOR de cálculo: resuelve el pórtico 2D y devuelve SOLO mecánica.

Qué es (y qué no es)
--------------------
Es la capa 2 de la arquitectura (ver ARQUITECTURA.md): el MOTOR. No sabe de
hormigón ni de acero, no dimensiona nada y no pide datos por teclado. Recibe el
pórtico (geometría + cargas) y devuelve las solicitaciones (M, V, N), las
reacciones en las bases y los desplazamientos, POR BARRA Y POR COMBINACIÓN. El
dimensionado (vigas, columnas, bases) se hace después leyendo esto.

Reemplaza a anaStruct (el viejo P01). Se validó que Pynite reproduce a anaStruct
en <0,1 % con los pórticos ya calculados (tools/comparar_motores.py), y esta
versión POR COMBINACIÓN dio idéntica a lo guardado en datos/estructura.json.

Para que los números coincidan se usan las MISMAS rigideces relativas que usaba
anaStruct (EA=15000, EI=5000): las fuerzas internas de un pórtico elástico no
dependen de E, solo del cociente entre la rigidez axial y la de flexión.

Convención de signos (la de Pynite): M en kN·m, V y N en kN. En las columnas
M_inf es el momento en la base y M_sup en la punta; en las vigas M_izq/M_der son
los extremos y M_campo el momento mínimo de la luz.

Uso:
    py -m calc.portico                        qué pórticos hay y qué está resuelto
    py -m calc.portico "Portico 1"            resuelve ese pórtico e imprime el informe
    py -m calc.portico "Portico 1" --guardar  además guarda salidas/solicitaciones
"""

from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

from Pynite import FEModel3D

from . import cargas, rutas

# ---------------------------------------------------------------------------
# Rigidez (los MISMOS valores con que se validó Pynite contra anaStruct)
# ---------------------------------------------------------------------------
E_MOD = 1.0
G_MOD = 1.0
EA = 15000.0
EI = 5000.0
A_SEC = EA / E_MOD
I_SEC = EI / E_MOD
J_SEC = EI / E_MOD
POISSON = 0.3


def combinaciones() -> dict[str, tuple[float, float, float, float]]:
    """
    Las combinaciones del pórtico: (fD, fL, fW, fP).

    Las 4 de CIRSOC salen de `cargas.FACTORES_COMBINACIONES` (una sola cuenta).
    'servicio' (D + L, sin viento) es el estado límite de servicio. Las cargas
    puntuales ('P') se aplican enteras en todas las combinaciones, igual que en
    el cálculo viejo (P01).
    """
    tabla: dict[str, tuple[float, float, float, float]] = {"servicio": (1.0, 1.0, 0.0, 1.0)}
    for nombre, (fD, fL, fW) in cargas.FACTORES_COMBINACIONES.items():
        tabla[nombre] = (fD, fL, fW, 1.0)
    return tabla


# ---------------------------------------------------------------------------
# Nudos y nombres
# ---------------------------------------------------------------------------
def _k(x) -> float:
    return round(float(x), 4)


def _nodo(x: float, y: float) -> str:
    return f"N{_k(x)}_{_k(y)}"


def _piso_de_viga(viga_id: str) -> int | None:
    """De 'V0-1' saca 0; de 'V12-2' saca 12. None si el nombre no sigue el patrón."""
    cabeza = viga_id.split("-")[0]
    if cabeza[:1].upper() == "V" and cabeza[1:].isdigit():
        return int(cabeza[1:])
    return None


# ---------------------------------------------------------------------------
# Armar el modelo (pórtico plano)
# ---------------------------------------------------------------------------
def construir(portico: dict) -> tuple[FEModel3D, dict, list[str]]:
    """
    Arma el modelo Pynite del pórtico como PÓRTICO PLANO (libre en el plano,
    fijo fuera de él: DZ, RX y RY bloqueados en todos los nudos). Aplica las
    cargas como casos base D, L, W y P y deja cargadas las combinaciones.

    Devuelve (modelo, mapa_de_barras, avisos). El mapa dice el nombre Pynite de
    cada barra: mapa['columnas'][cid], mapa['vigas'][tramo], mapa['voladizos'][tramo].
    """
    avisos: list[str] = []
    m = FEModel3D()
    m.add_material("mat", E_MOD, G_MOD, POISSON, 0.0)
    m.add_section("sec", A_SEC, I_SEC, I_SEC, J_SEC)

    mapa: dict[str, dict] = {"columnas": {}, "vigas": {}, "voladizos": {}}
    columnas = portico.get("columnas", {})

    # --- Columnas (de abajo hacia arriba: i-end = base, j-end = punta) -----
    for cid, col in columnas.items():
        x = _k(col["x"])
        y0 = _k(col.get("nivel", 0.0))
        h = _k(col["altura_m"])
        nb, nt = _nodo(x, y0), _nodo(x, y0 + h)
        if nb not in m.nodes:
            m.add_node(nb, x, y0, 0.0)
        if nt not in m.nodes:
            m.add_node(nt, x, y0 + h, 0.0)
        m.def_support(nb, False, False, True, True, True, False)
        m.def_support(nt, False, False, True, True, True, False)
        m.add_member(f"C-{cid}", nb, nt, "mat", "sec")
        mapa["columnas"][cid] = f"C-{cid}"

    # --- Bases (apoyos) ----------------------------------------------------
    for base_id, base in portico.get("bases", {}).items():
        nodo = _nodo(base["x"], 0.0)
        if nodo not in m.nodes:
            avisos.append(f"{base_id}: no hay columna en x={base['x']}, se ignora")
            continue
        if base.get("tipo") == "articulado":
            m.def_support(nodo, True, True, True, False, False, False)
        else:  # empotramiento (por defecto)
            m.def_support(nodo, True, True, True, True, True, True)

    # --- Vigas (por piso), voladizos y sus cargas --------------------------
    for viga_id, viga in portico.get("vigas", {}).items():
        piso = _piso_de_viga(viga_id)
        cols = sorted(
            (c for cid, c in columnas.items() if piso is not None and cid.startswith(f"C{piso}-")),
            key=lambda c: c["x"],
        )
        if not cols:
            avisos.append(f"{viga_id}: no encuentro sus columnas (piso {piso}); se ignora")
            continue
        y_viga = _k(cols[0].get("nivel", 0.0) + cols[0]["altura_m"])

        for i, tramo in enumerate(viga.get("tramos", [])):
            cargas_tramo = tramo.get("cargas", {})
            miembro_id = f"V-{tramo['id']}"

            if tramo.get("es_voladizo"):
                if "izq" in tramo["id"]:
                    tip_x, base_x = _k(cols[0]["x"] - tramo["longitud_m"]), _k(cols[0]["x"])
                else:
                    tip_x, base_x = _k(cols[-1]["x"] + tramo["longitud_m"]), _k(cols[-1]["x"])
                n_tip, n_base = _nodo(tip_x, y_viga), _nodo(base_x, y_viga)
                if n_tip not in m.nodes:
                    m.add_node(n_tip, tip_x, y_viga, 0.0)
                    m.def_support(n_tip, False, False, True, True, True, False)
                # del extremo libre hacia la columna: mismo orden que usó anaStruct
                m.add_member(miembro_id, n_tip, n_base, "mat", "sec")
                mapa["voladizos"][tramo["id"]] = miembro_id
            else:
                if i >= len(cols) - 1:
                    avisos.append(f"{tramo['id']}: faltan columnas para el tramo; se ignora")
                    continue
                m.add_member(miembro_id, _nodo(cols[i]["x"], y_viga), _nodo(cols[i + 1]["x"], y_viga), "mat", "sec")
                mapa["vigas"][tramo["id"]] = miembro_id

            # cargas repartidas (casos base D y L)
            wD = float(cargas_tramo.get("D_total", 0.0))
            wL = float(cargas_tramo.get("L_total", 0.0))
            if abs(wD) > 1e-9:
                m.add_member_dist_load(miembro_id, "FY", -wD, -wD, case="D")
            if abs(wL) > 1e-9:
                m.add_member_dist_load(miembro_id, "FY", -wL, -wL, case="L")

            # cargas puntuales (caso P, entero en toda combinación)
            for cp in tramo.get("cargas_puntuales", []):
                n_cp = _nodo(cp["coordenadas"]["x"], cp["coordenadas"]["y"])
                if n_cp in m.nodes:
                    m.add_node_load(n_cp, "FY", -float(cp["valor_kN"]), case="P")
                else:
                    avisos.append(f"{tramo['id']}: carga puntual fuera de nudo ({n_cp}); se ignora")

    # --- Viento (caso W): fuerza horizontal en la punta de cada columna ----
    for viga_id, viga in portico.get("vigas", {}).items():
        piso = _piso_de_viga(viga_id)
        if piso is None or not viga.get("tramos"):
            continue
        cols = [c for cid, c in columnas.items() if cid.startswith(f"C{piso}-")]
        if not cols:
            continue
        w_total = float(viga["tramos"][0].get("cargas", {}).get("W_total", 0.0))
        total = w_total * float(cols[0]["altura_m"])
        if abs(total) < 1e-9:
            continue
        por_col = total / len(cols)  # reparto igual entre las columnas del piso
        for cid, col in columnas.items():
            if not cid.startswith(f"C{piso}-"):
                continue
            n_punta = _nodo(col["x"], col.get("nivel", 0.0) + col["altura_m"])
            if n_punta in m.nodes:
                m.add_node_load(n_punta, "FX", por_col, case="W")

    # --- Combinaciones -----------------------------------------------------
    for nombre, (fD, fL, fW, fP) in combinaciones().items():
        factores: dict[str, float] = {}
        if fD:
            factores["D"] = fD
        if fL:
            factores["L"] = fL
        if fW:
            factores["W"] = fW
        if fP:
            factores["P"] = fP
        m.add_load_combo(nombre, factores)

    return m, mapa, avisos


# ---------------------------------------------------------------------------
# Resolver y extraer las solicitaciones
# ---------------------------------------------------------------------------
def _envolvente_vacia(mapa: dict) -> dict:
    return {
        "columnas": {cid: {} for cid in mapa["columnas"]},
        "vigas": {t: {} for t in mapa["vigas"]},
        "voladizos": {t: {} for t in mapa["voladizos"]},
        "bases": {},
    }


def _acumular(bolsa: dict, clave: str, valores: dict) -> None:
    """Guarda, por barra, el máximo EN MÓDULO de cada magnitud (la envolvente)."""
    actual = bolsa.setdefault(clave, {})
    for campo, valor in valores.items():
        v = abs(float(valor))
        if campo not in actual or v > actual[campo]:
            actual[campo] = v


def calcular(portico: dict, nombre: str = "") -> dict:
    """
    Resuelve el pórtico y devuelve las solicitaciones, reacciones y
    desplazamientos por barra y por combinación, más una ENVOLVENTE con los
    máximos en módulo (la que usarán los dimensionadores). No guarda nada: eso
    lo hace `guardar`.
    """
    m, mapa, avisos = construir(portico)
    m.analyze()

    combos = list(combinaciones())
    solicitaciones: dict[str, dict] = {}
    envolvente = _envolvente_vacia(mapa)

    for combo in combos:
        bloque: dict[str, dict] = {
            "columnas": {}, "vigas": {}, "voladizos": {}, "bases": {}, "desplazamientos": {},
        }

        for cid, miembro_id in mapa["columnas"].items():
            mem = m.members[miembro_id]
            col = portico["columnas"][cid]
            h = float(col.get("altura_m", 1.0)) or 1.0
            n_punta = _nodo(col["x"], col.get("nivel", 0.0) + col["altura_m"])
            datos = {
                "M_inf": float(mem.moment("Mz", 0.0, combo)),
                "M_sup": float(mem.moment("Mz", mem.L(), combo)),
                "N": float(mem.axial(0.0, combo)),
                "deriva": float(m.nodes[n_punta].DX[combo]) / h,
            }
            bloque["columnas"][cid] = datos
            _acumular(envolvente["columnas"], cid, datos)

        for grupo in ("vigas", "voladizos"):
            for clave, miembro_id in mapa[grupo].items():
                mem = m.members[miembro_id]
                datos = {
                    "M_izq": float(mem.moment("Mz", 0.0, combo)),
                    "M_der": float(mem.moment("Mz", mem.L(), combo)),
                    "M_campo": float(mem.min_moment("Mz", combo)),
                    "V_izq": float(mem.shear("Fy", 0.0, combo)),
                    "V_der": float(mem.shear("Fy", mem.L(), combo)),
                }
                bloque[grupo][clave] = datos
                _acumular(envolvente[grupo], clave, datos)

        for base_id, base in portico.get("bases", {}).items():
            nodo = _nodo(base["x"], 0.0)
            if nodo not in m.nodes:
                continue
            n = m.nodes[nodo]
            datos = {
                "Fx": float(n.RxnFX.get(combo, 0.0)),
                "Fy": float(n.RxnFY.get(combo, 0.0)),
                "Mz": float(n.RxnMZ.get(combo, 0.0)),
            }
            bloque["bases"][base_id] = datos
            _acumular(envolvente["bases"], base_id, datos)

        for nodo_id, nodo in m.nodes.items():
            bloque["desplazamientos"][nodo_id] = {
                "DX": float(nodo.DX[combo]),
                "DY": float(nodo.DY[combo]),
                "RZ": float(nodo.RZ[combo]),
            }

        solicitaciones[combo] = bloque

    return {
        "portico": nombre,
        "generado": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "rigidez": {"EA": EA, "EI": EI},
        "combinaciones": combos,
        "avisos": avisos,
        "solicitaciones": solicitaciones,
        "envolvente": envolvente,
    }


# ---------------------------------------------------------------------------
# Informe y guardado
# ---------------------------------------------------------------------------
def informe(resultado: dict) -> str:
    lineas = [
        f"SOLICITACIONES DEL PÓRTICO: {resultado['portico']}",
        "=" * 74,
        f"Motor Pynite (pórtico plano) · rigidez EA={EA:g} EI={EI:g} · {resultado['generado']}",
        "Unidades: M [kN·m], N y V [kN], deriva [m/m]. Signos según Pynite.",
    ]
    for aviso in resultado["avisos"]:
        lineas.append(f"[AVISO] {aviso}")

    tabla = resultado["solicitaciones"]
    for combo in resultado["combinaciones"]:
        bloque = tabla[combo]
        lineas.append(f"\n=== {combo} ===")

        if bloque["columnas"]:
            lineas.append("  Columnas:   barra        M_base     M_punta          N      deriva")
            for cid, d in bloque["columnas"].items():
                lineas.append(
                    f"    {cid:<10} {d['M_inf']:>10.2f} {d['M_sup']:>11.2f} "
                    f"{d['N']:>10.2f} {d['deriva']:>11.5f}"
                )

        for titulo, grupo in (("Vigas", "vigas"), ("Voladizos", "voladizos")):
            if not bloque[grupo]:
                continue
            lineas.append(f"  {titulo}:     tramo        M_izq      M_der    M_campo      V_izq      V_der")
            for clave, d in bloque[grupo].items():
                lineas.append(
                    f"    {clave:<12} {d['M_izq']:>8.2f} {d['M_der']:>10.2f} {d['M_campo']:>10.2f} "
                    f"{d['V_izq']:>10.2f} {d['V_der']:>10.2f}"
                )

        if bloque["bases"]:
            lineas.append("  Bases:      base             Fx         Fy         Mz")
            for base_id, d in bloque["bases"].items():
                lineas.append(f"    {base_id:<10} {d['Fx']:>10.2f} {d['Fy']:>10.2f} {d['Mz']:>10.2f}")

    lineas.append("\nENVOLVENTE (máximos en módulo, para el dimensionado)")
    lineas.append("-" * 74)
    for cid, d in resultado["envolvente"]["columnas"].items():
        if d:
            lineas.append(f"    {cid:<10} |M| máx = {max(d.get('M_inf', 0.0), d.get('M_sup', 0.0)):8.2f} kNm"
                          f"   |N| máx = {d.get('N', 0.0):8.2f} kN")
    for grupo in ("vigas", "voladizos"):
        for clave, d in resultado["envolvente"][grupo].items():
            if d:
                lineas.append(f"    {clave:<12} |M| máx = {d.get('M_campo', 0.0):8.2f} kNm"
                              f"   |V| máx = {d.get('V_izq', 0.0):8.2f} kN")
    return "\n".join(lineas)


def guardar(resultado: dict) -> Path:
    """Guarda las solicitaciones en salidas/solicitaciones/<portico>.json."""
    rutas.SAL_SOLICITACIONES.mkdir(parents=True, exist_ok=True)
    nombre = rutas.nombre_seguro(resultado["portico"]) or "portico"
    archivo = rutas.SAL_SOLICITACIONES / f"{nombre}.json"
    return rutas.guardar_json(archivo, resultado)


# ---------------------------------------------------------------------------
# Listado y consola
# ---------------------------------------------------------------------------
def estado_texto() -> str:
    estructura = rutas.cargar_estructura()
    lineas = ["PÓRTICOS DISPONIBLES (datos/estructura.json)", "=" * 74]
    if not estructura:
        lineas.append("(no hay pórticos cargados)")
    for nombre, p in estructura.items():
        n_col = len(p.get("columnas", {}))
        n_tra = sum(len(v.get("tramos", [])) for v in p.get("vigas", {}).values())
        lineas.append(f" - {nombre:<26} {n_col} columnas · {n_tra} tramos")

    archivos = rutas.listar(rutas.SAL_SOLICITACIONES, "*.json")
    lineas.append("")
    if archivos:
        lineas.append(f"Solicitaciones ya guardadas ({len(archivos)}):")
        for a in archivos:
            lineas.append(f"   {a.name}")
    else:
        lineas.append("(todavía no hay solicitaciones guardadas)")
    lineas.append("")
    lineas.append('Resolver uno:   py -m calc.portico "<nombre>" [--guardar]')
    return "\n".join(lineas)


def _main(argv: list[str]) -> int:
    opciones = [a for a in argv if a.startswith("-")]
    nombres = [a for a in argv if not a.startswith("-")]

    if not nombres:
        print(estado_texto())
        return 0

    estructura = rutas.cargar_estructura()
    codigo = 0
    for nombre in nombres:
        if nombre not in estructura:
            print(f"No existe el pórtico '{nombre}'. Ver: py -m calc.portico")
            codigo = 1
            continue
        resultado = calcular(estructura[nombre], nombre)
        print(informe(resultado))
        if "--guardar" in opciones:
            print(f"\nGuardado en: {guardar(resultado)}")
    return codigo


if __name__ == "__main__":
    raise SystemExit(_main(sys.argv[1:]))

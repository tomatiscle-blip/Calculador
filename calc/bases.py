import json
import os
from pathlib import Path

from . import rutas


def dimensionar_base_geotecnia(Fy, Mx=0.0, My=0.0, q_adm=100.0, Bx=None, By=None):
    """Dimensionado geotécnico biaxial para zapata rectangular con iteración."""
    N = abs(float(Fy or 0.0))
    Mx = float(Mx or 0.0)
    My = float(My or 0.0)
    q_adm = float(q_adm or 0.0)
    if N <= 0 or q_adm <= 0:
        return {"ok": False, "error": "N o q_adm inválidos"}

    if Bx is None or By is None:
        lado = (N / q_adm) ** 0.5
        Bx = By = max(0.6, round(lado * 1.1, 2))

    Bx = float(Bx)
    By = float(By)
    paso = 0.05
    limite = 3.0
    while True:
        A = Bx * By
        Ix = Bx * By**3 / 12.0
        Iy = By * Bx**3 / 12.0
        vertices = []
        for sx in (-1, 1):
            for sy in (-1, 1):
                x = sx * Bx / 2.0
                y = sy * By / 2.0
                p = N / A + Mx * y / Ix + My * x / Iy
                vertices.append({"x": x, "y": y, "p": p})
        ps = [v["p"] for v in vertices]
        pmax = max(ps)
        pmin = min(ps)
        despegue = pmin < 0
        ok = (pmax <= q_adm) and (not despegue)
        if ok or (Bx >= limite and By >= limite):
            return {
                "ok": ok,
                "Bx": round(Bx, 2),
                "By": round(By, 2),
                "area": round(A, 4),
                "pmax_kPa": pmax,
                "pmin_kPa": pmin,
                "despegue": despegue,
                "vertices": vertices,
                "q_adm_kPa": q_adm,
                "requiere_diseno_especial": (not ok),
            }
        Bx = round(min(Bx + paso, limite), 2)
        By = round(min(By + paso, limite), 2)
    entrada = os.path.join("salidas", "solicitaciones", f"{portico}.json")
    salida = os.path.join("salidas", "bases", f"resultados_{portico}_bases.json")

    with open(entrada, "r") as f:
        solicitaciones = json.load(f)

    resultados = []
    for bid, valores in solicitaciones.get("envolvente", {}).get("bases", {}).items():
        resultados.append({"id": bid, **valores})

    os.makedirs(os.path.dirname(salida), exist_ok=True)
    with open(salida, "w") as f:
        json.dump(resultados, f, indent=4)

    return resultados


def calcular_global() -> list[dict]:
    consolidado = rutas.leer_json(rutas.SOLICITACIONES_BASES, {}) or {}
    estructura = rutas.cargar_estructura()
    resultados = []
    for bid in sorted(consolidado):
        entradas = consolidado[bid] or []
        combinado: dict = {}
        # Cargas verticales/horizontales: se toma el máximo en módulo.
        for clave in ("Fx", "Fy"):
            valores = [float(e.get(clave, 0.0) or 0.0) for e in entradas]
            if valores:
                combinado[clave] = max(valores, key=abs)
        # Momentos por eje: Mz del pórtico se asigna a Mx o My según dirección.
        mx = 0.0
        my = 0.0
        for entrada in entradas:
            mz = float(entrada.get("Mz", 0.0) or 0.0)
            try:
                portico = rutas.pórtico(estructura, entrada.get("portico")) or {}
                direccion = (portico.get("referencia_planta", {}) or {}).get("direccion")
            except Exception:
                direccion = None
            if direccion == "x":
                if abs(mz) > abs(my):
                    my = mz
            else:
                if abs(mz) > abs(mx):
                    mx = mz
        combinado["Mx"] = mx
        combinado["My"] = my
        combinado["M"] = max(abs(mx), abs(my))
        resultados.append({
            "id": bid,
            "porticos": [e.get("portico") for e in entradas],
            "n_porticos": len(entradas),
            "momentos_otro_sentido": any(bool(e.get("momentos_otro_sentido")) for e in entradas),
            **combinado,
            "estado": "OK" if combinado.get("M", 0) < 200 and abs(combinado.get("Fy", 0)) < 400 else "Revisar",
        })

    salida = rutas.SAL_BASES / "resultados_bases_globales.json"
    rutas.SAL_BASES.mkdir(parents=True, exist_ok=True)
    rutas.guardar_json(salida, resultados)
    return resultados


def calcular(datos, portico="Portico 4"):
    entrada = os.path.join("salidas", "solicitaciones", f"{portico}.json")
    salida = os.path.join("salidas", "bases", f"resultados_{portico}_bases.json")

    with open(entrada, "r") as f:
        solicitaciones = json.load(f)

    resultados = []
    for bid, valores in solicitaciones.get("envolvente", {}).get("bases", {}).items():
        resultados.append({"id": bid, **valores})

    os.makedirs(os.path.dirname(salida), exist_ok=True)
    with open(salida, "w") as f:
        json.dump(resultados, f, indent=4)

    return resultados
def dimensionar_global(q_adm: float | None = None, Bx: float | None = None, By: float | None = None) -> list[dict]:
    """Lee resultados_bases_globales.json y aplica dimensionado biaxial a cada base."""
    from . import configuracion_obra
    if q_adm is None:
        datos_terreno = rutas.leer_json(rutas.TERRENO, {}) or {}
        valor_q = datos_terreno.get("q_adm_kPa")
        if valor_q in (None, "", 0):
            raise ValueError("Falta la presión admisible q_adm del estudio geotécnico.")
        q_adm = float(valor_q)
    parametros = configuracion_obra.cargar()
    if Bx is None and By is None:
        bx_guardado = parametros.get("BASE_BX_INICIAL_M")
        by_guardado = parametros.get("BASE_BY_INICIAL_M")
        if bx_guardado and by_guardado:
            Bx, By = float(bx_guardado), float(by_guardado)
    datos = rutas.leer_json(rutas.SAL_BASES / "resultados_bases_globales.json", []) or []
    resultados = []
    for b in datos:
        res = dimensionar_base_geotecnia(
            b.get("Fy", 0.0),
            b.get("Mx", 0.0),
            b.get("My", 0.0),
            q_adm,
            Bx,
            By,
        )
        resultados.append({
            "id": b.get("id"),
            "Fy": b.get("Fy"),
            "Mx": b.get("Mx"),
            "My": b.get("My"),
            "M": b.get("M"),
            "porticos": b.get("porticos", []),
            **res,
        })
    salida = rutas.SAL_BASES / "resultados_bases_globales.json"
    rutas.SAL_BASES.mkdir(parents=True, exist_ok=True)
    rutas.guardar_json(salida, resultados)
    return resultados
def optimizar_fundaciones(q_adm: float | None = None, Bx: float | None = None, By: float | None = None) -> list[dict]:
    """Orquesta: dimensionado inicial -> vigas de fundacion -> redistribucion -> redimensionado."""
    # 1. consolidado ya existe (resultados_bases_globales.json) con Mx/My
    # 2. dimensionado inicial
    bases_dim = dimensionar_global(q_adm, Bx, By)
    # 3. vigas de fundacion
    try:
        import importlib
        vf_mod = importlib.import_module("calc.vigas_fundacion")
        vf_mod.calcular()
    except Exception:
        pass
    # 4. redistribucion iterativa hasta que no haya despegue o max iteraciones
    # leer vigas para saber conexiones
    vigas = rutas.leer_json(rutas.SAL_BASES / "resultados_vigas_fundacion.json", []) or []
    vecinos = {}
    for v in vigas:
        b0, b1 = v.get("bases", [None, None])
        vecinos.setdefault(b0, set()).add(b1)
        vecinos.setdefault(b1, set()).add(b0)

    por_id = {b.get("id"): b for b in bases_dim}
    max_iter = 10
    pct = 0.3
    for it in range(max_iter):
        especiales = [b for b in bases_dim if b.get("requiere_diseno_especial")]
        if not especiales:
            break
        pct = min(pct + 0.1, 0.8)
        for b in especiales:
            for vec in (vecinos.get(b.get("id")) or []):
                vb = por_id.get(vec)
                if vb is None:
                    continue
                exceso_mx = b.get("Mx", 0.0) * pct
                exceso_my = b.get("My", 0.0) * pct
                b["_momento_inicial_mx"] = b.get("_momento_inicial_mx", b.get("Mx", 0.0) + exceso_mx)
                b["_momento_inicial_my"] = b.get("_momento_inicial_my", b.get("My", 0.0) + exceso_my)
                b["_momento_transferido_mx"] = b.get("_momento_transferido_mx", 0.0) + exceso_mx
                b["_momento_transferido_my"] = b.get("_momento_transferido_my", 0.0) + exceso_my
                b["Mx"] = b.get("Mx", 0.0) - exceso_mx
                b["My"] = b.get("My", 0.0) - exceso_my
                vb["Mx"] = vb.get("Mx", 0.0) + exceso_mx
                vb["My"] = vb.get("My", 0.0) + exceso_my
        datos_ajustados = []
        for b in bases_dim:
            datos_ajustados.append({
                "id": b.get("id"),
                "Fy": b.get("Fy"),
                "Mx": b.get("Mx"),
                "My": b.get("My"),
                "M": max(abs(b.get("Mx", 0.0)), abs(b.get("My", 0.0))),
                "porticos": b.get("porticos", []),
                "_momento_inicial_mx": b.get("_momento_inicial_mx", b.get("Mx")),
                "_momento_inicial_my": b.get("_momento_inicial_my", b.get("My")),
                "_momento_transferido_mx": b.get("_momento_transferido_mx", 0.0),
                "_momento_transferido_my": b.get("_momento_transferido_my", 0.0),
            })
        rutas.guardar_json(rutas.SAL_BASES / "resultados_bases_globales.json", datos_ajustados)
        bases_dim = dimensionar_global(q_adm, Bx, By)
        por_id = {b.get("id"): b for b in bases_dim}
        for b in bases_dim:
            old = next((x for x in datos_ajustados if x.get("id") == b.get("id")), {})
            b["_momento_inicial_mx"] = old.get("_momento_inicial_mx", b.get("Mx"))
            b["_momento_inicial_my"] = old.get("_momento_inicial_my", b.get("My"))
            b["_momento_transferido_mx"] = old.get("_momento_transferido_mx", 0.0)
            b["_momento_transferido_my"] = old.get("_momento_transferido_my", 0.0)
            b["_momento_final_mx"] = b.get("Mx")
            b["_momento_final_my"] = b.get("My")
    rutas.guardar_json(rutas.SAL_BASES / "resultados_bases_globales.json", bases_dim)
    return bases_dim
if __name__ == "__main__":
    print(calcular({}))

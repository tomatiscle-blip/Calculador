import io, json, sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from app.ejes import cargar_ejes, cargar_niveles
from calc import rutas

out = io.open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "_data.txt"), "w", encoding="utf-8")
ejes = cargar_ejes()
niveles = cargar_niveles()
out.write("EJES: " + json.dumps(ejes, ensure_ascii=False, indent=1) + "\n")
out.write("NIVELES: " + json.dumps(niveles, ensure_ascii=False, indent=1) + "\n")
for nombre in ("estructura.json", "losas.json", "cargas.json"):
    p = rutas.DATOS / nombre
    d = json.load(io.open(p, encoding="utf-8"))
    out.write(f"--- {nombre} keys: {list(d.keys())}\n")
    s = json.dumps(d, ensure_ascii=False, indent=1)
    out.write(s[:6000] + "\n")
out.close()

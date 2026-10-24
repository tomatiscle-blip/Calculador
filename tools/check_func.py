path = r"d:\03 Ingenieria\Calculador\calc\rutas.py"
src = open(path, encoding="utf-8").read()
i = src.find("def proximo_nombre_portico()")
j = src.find("def partes(", i)
func = src[i:j]
print("double Ã³ count:", func.count("Ã³"))
print("proper ó count:", func.count("ó"))
print("--- function ---")
print(func)

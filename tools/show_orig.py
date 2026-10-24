data = open(r"d:\03 Ingenieria\Calculador\tools\rutas_original.txt", "rb").read()
i = data.find(b"Pr\u00f3ximo".encode("utf-8"))
print("Próximo at:", i, repr(data[i-20:i+60])) if i != -1 else print("not found")
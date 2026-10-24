data = open(r"d:\03 Ingenieria\Calculador\calc\rutas.py", "rb").read()
idx = data.find(b"cat")
print("cat at:", idx, repr(data[idx-20:idx+60]))
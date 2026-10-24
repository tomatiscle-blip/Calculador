data = open(r"d:\03 Ingenieria\Calculador\calc\rutas.py", "rb").read()
i = data.find(b"Pr")
print("byte at Pr:", repr(data[i-8:i+30]))
# also check what utf-8 decode gives
text = data.decode("utf-8")
print("utf8 decoded:", repr(text[i-8:i+30]))
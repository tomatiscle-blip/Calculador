data = open(r"d:\03 Ingenieria\Calculador\tools\rutas_original.txt", "rb").read()
target = "Pr\u00f3ximo".encode("utf-8")
i = data.find(target)
print("Pr\u00f3ximo at:", i)
if i != -1:
    print(repr(data[i-20:i+60]))
else:
    # search for any 'Pr' followed by c3 b3
    j = data.find(b"Pr")
    print("first 'Pr' at:", j, repr(data[j-10:j+30]))
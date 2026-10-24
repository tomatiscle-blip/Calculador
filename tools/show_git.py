data = open(r"d:\03 Ingenieria\Calculador\tools\git_orig.txt", "rb").read()
print("size:", len(data))
target = "Pr\u00f3ximo".encode("utf-8")
i = data.find(target)
print("Pr\u00f3ximo at:", i)
if i != -1:
    print(repr(data[i-20:i+60]))
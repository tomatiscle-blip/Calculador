src = open("calc/portico.py", encoding="utf-8").read()
lines = src.splitlines()
print("TOTAL LINES:", len(lines))
# Print lines 543-700 (construir body)
for i in range(542, 700):
    print(f"{i+1:4d}| {lines[i]}")

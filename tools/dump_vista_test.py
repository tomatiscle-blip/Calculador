lines = open("tests/test_vista_planta.py", encoding="utf-8").read().splitlines()
for i in range(120, 216):
    print(f"{i+1}| {lines[i]}")

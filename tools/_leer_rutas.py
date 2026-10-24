with open('calc/rutas.py', 'r', encoding='utf-8') as f:
    lines = f.readlines()
total = len(lines)
print(f"Total lines: {total}")
for i in range(0, total, 100):
    print(f"\n===== LINES {i+1}-{min(i+100, total)} =====")
    for j in range(i, min(i+100, total)):
        print(f'{j+1}: {lines[j]}', end='')

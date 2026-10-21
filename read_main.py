#!/usr/bin/env python3
import linecache
path = r"d:\03 Ingenieria\Calculador\app\principal.py"
for i in range(2468, 2476):
    line = linecache.getline(path, i)
    print(f"{i}: {line}", end='')
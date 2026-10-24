import glob
import sys

for f in sorted(glob.glob("app/*.py")):
    print("===", f, "===")
    with open(f, encoding="utf-8") as fh:
        for i, linea in enumerate(fh, 1):
            if "estructura" in linea:
                print(f"{i}: {linea.rstrip()}")

import glob

FILES = sorted(
    glob.glob("app/*.py")
    + glob.glob("calc/*.py")
    + glob.glob("tests/*.py")
)


def main():
    for f in FILES:
        with open(f, encoding="utf-8") as fh:
            for i, linea in enumerate(fh, 1):
                linea = linea.rstrip()
                if "EditorPorticoPorNiveles" in linea:
                    print(f"{f}:{i}: {linea}")
                if "estructura.get" in linea and ".get(" in linea:
                    print(f"{f}:{i}: {linea}")


if __name__ == "__main__":
    main()

data = open("calc/rutas.py", "rb").read()
print("total bytes:", len(data))
try:
    text = data.decode("utf-8")
    print("utf-8 OK")
except UnicodeDecodeError as e:
    print("UTF-8 ERROR:", e)
    print(repr(data[max(0, e.start-40):e.start+40]))

# check for BOM
print("BOM:", data[:3] == b"\xef\xbb\xbf")

# find the 'def partes' position
idx = data.find(b"def partes")
print("def partes at byte:", idx)
print(repr(data[idx-200:idx+120]))

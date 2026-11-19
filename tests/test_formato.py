import csv
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from calc.formato import numero_local, parse_numero
from calc.losas import COLUMNAS_COMPUTO, _migrar_computo_csv


class FormatoRegionalTests(unittest.TestCase):
    def test_formatea_miles_y_decimales(self):
        self.assertEqual(numero_local(12345.678, 3), "12.345,678")
        self.assertEqual(numero_local(8, 3), "8,000")

    def test_lee_csv_historico_y_regional(self):
        self.assertEqual(parse_numero("12.345"), 12.345)
        self.assertEqual(parse_numero("12.345,67"), 12345.67)

    def test_migra_computo_historico_a_delimitador_excel(self):
        with TemporaryDirectory() as directorio:
            ruta = Path(directorio) / "computo.csv"
            fila = ["LX1", "3.0", "3.1", "A", "B", "2", "3", "0.125"]
            fila.extend(["4.5", "5.5", "12.0", "2.40x6", "10.0", "1", "4", "2", "8", "1.25", "8", "10.0"])
            with ruta.open("w", encoding="utf-8", newline="") as archivo:
                writer = csv.writer(archivo)
                writer.writerow(COLUMNAS_COMPUTO)
                writer.writerow(fila)

            _migrar_computo_csv(ruta)

            with ruta.open("r", encoding="utf-8-sig", newline="") as archivo:
                filas = list(csv.reader(archivo, delimiter=";"))
            self.assertEqual(tuple(filas[0]), COLUMNAS_COMPUTO)
            self.assertEqual(filas[1][1], "3,00")
            self.assertEqual(filas[1][7], "0,125")
            self.assertEqual(filas[1][8], "4,50")


if __name__ == "__main__":
    unittest.main()

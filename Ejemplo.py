import argparse

# 1. Definición de los argumentos de la linea de comandos
parser = argparse.ArgumentParser(
    description="Construye un resumen del regulon "
    "a partir de un archivo de interacciones"
)

parser.add_argument("input_file", help="Nombre del archivo de interacciones")

parser.add_argument(
    "output_file", help="Nombre del archivo de salida para el resumen del regulon."
)

## Leyendo los argumentos
args = parser.parse_args()
print(args)

filenames = args.input_file
output_filename = args.output_file

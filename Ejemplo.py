import argparse

# Crear el objeto ArgumentParser para manejar los argumentos de la línea de comandos
parser = argparse.ArgumentParser(
    description="Construye un resumen del regulón a partir de un archivo de interacciones."
)

# Definición de los argumentos de la línea de comandos
parser.add_argument("input_file", help="Nombre del archivo de interacciones")
parser.add_argument(
    "output_file", help="Nombre del archivo de salida para el resumen del regulón"
)

parser.add_argument(
    "-g",
    "--min_genes",
    type=int,
    default=1,
    required=True,
    help="Número mínimo de genes para incluir en el resumen.",
)

# Leyendo los argumentos
args = parser.parse_args()

filename = args.input_file
output_filename = args.output_file

# Imprimiendo los argumentos para verificar que se han leído correctamente
print(f"Archivo de entrada: {filename}")
print(f"Archivo de salida: {output_filename}")
print(f"Número mínimo de genes: {args.min_genes}")

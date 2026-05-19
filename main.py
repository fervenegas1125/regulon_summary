# ==========================================================
# Programa: regulon_summary.py
#
# Objetivo:
# Generar un resumen de regulones a partir de un archivo
# de interacciones de RegulonDB.
# ==========================================================

import os
import argparse


# ==========================================================
# Responsabilidad:
# Leer interacciones desde archivo
#
# Entrada:
#     filename (str)
#
# Salida:
#     interactions (list)
# ==========================================================
def load_interactions(filename):

    interactions = []

    # Verificar existencia del archivo
    if not os.path.exists(filename):
        print("Error: archivo no encontrado")
        exit(1)

    else:

        with open(filename, "r", encoding="utf-8") as f:

            for line in f:

                line = line.strip()

                # Ignorar líneas vacías
                if not line:
                    continue

                # Ignorar comentarios
                if line.startswith("#"):
                    continue

                # Ignorar encabezado
                if line.startswith("1)regulatorId"):
                    continue

                # Separar columnas
                fields = line.split("\t")

                # Validar número mínimo de columnas
                if len(fields) <= 6:
                    continue

                # Columnas utilizadas
                tf = fields[1]
                gene = fields[4]
                effect = fields[5]

                # Validar effect
                if effect not in ["+", "-", "-+", "+-"]:
                    continue

                # Guardar interacción
                interactions.append((tf, gene, effect))

    return interactions


# ==========================================================
# Responsabilidad:
# Construir regulón
#
# Entrada:
#     interactions (list)
#
# Salida:
#     regulon (dict)
# ==========================================================
def build_regulon(interactions):

    regulon = {}

    for tf, gene, effect in interactions:

        if tf not in regulon:

            regulon[tf] = {"genes": [], "activados": 0, "reprimidos": 0}

        # Agregar gen
        regulon[tf]["genes"].append(gene)

        # Contar activaciones y represiones
        if effect == "+":
            regulon[tf]["activados"] += 1

        elif effect == "-":
            regulon[tf]["reprimidos"] += 1

        elif effect in ["-+", "+-"]:

            regulon[tf]["activados"] += 1
            regulon[tf]["reprimidos"] += 1

    return regulon


# ==========================================================
# Responsabilidad:
# Determinar si TF es activador, represor o dual
#
# Entrada:
#     activados (int)
#     reprimidos (int)
#
# Salida:
#     tipo_tf (str)
# ==========================================================
def get_tf_type(activados, reprimidos):

    if activados > 0 and reprimidos == 0:
        tipo_tf = "Activador"

    elif activados == 0 and reprimidos > 0:
        tipo_tf = "Represor"

    elif activados > 0 and reprimidos > 0:
        tipo_tf = "Dual"

    else:
        tipo_tf = "Desconocido"

    return tipo_tf


# ==========================================================
# Responsabilidad:
# Filtrar regulones por número mínimo de genes
#
# Entrada:
#     regulon (dict)
#     min_genes (int)
#
# Salida:
#     regulon_filtrado (dict)
# ==========================================================
def filter_regulon_by_min_genes(regulon, min_genes):

    regulon_filtrado = {}

    for tf in regulon:

        total_genes = len(regulon[tf]["genes"])

        # Conservar TFs válidos
        if total_genes >= min_genes:

            regulon_filtrado[tf] = regulon[tf]

    return regulon_filtrado


# ==========================================================
# Responsabilidad:
# Escribir resumen del regulón
#
# Entrada:
#     regulon (dict)
#     output_filename (str)
#
# Salida:
#     archivo TSV
# ==========================================================
def write_summary(regulon, output_filename):

    with open(output_filename, "w", encoding="utf-8") as f:

        # Encabezado
        f.write(
            "TF\t"
            "Total_genes_regulados\t"
            "Activados\t"
            "Reprimidos\t"
            "Tipo_TF\t"
            "Lista_genes\n"
        )

        # Recorrer TFs
        for tf in sorted(regulon):

            total_genes = len(regulon[tf]["genes"])

            lista_genes = ",".join(regulon[tf]["genes"])

            activados = regulon[tf]["activados"]

            reprimidos = regulon[tf]["reprimidos"]

            tipo_tf = get_tf_type(activados, reprimidos)

            # Escribir línea
            f.write(
                f"{tf}\t"
                f"{total_genes}\t"
                f"{activados}\t"
                f"{reprimidos}\t"
                f"{tipo_tf}\t"
                f"{lista_genes}\n"
            )


# ==========================================================
# MAIN
# ==========================================================
def main():

    # ======================================================
    # Definición de argumentos
    # ======================================================
    parser = argparse.ArgumentParser(description="Construye resumen del regulón")

    parser.add_argument("input_file", help="Archivo TSV de entrada")

    parser.add_argument("output_file", help="Archivo TSV de salida")

    # Nuevo argumento opcional
    parser.add_argument(
        "--min_genes",
        type=int,
        default=0,
        help="Número mínimo de genes regulados por TF",
    )

    # Leer argumentos
    args = parser.parse_args()

    input_file = args.input_file

    output_file = args.output_file

    min_genes = args.min_genes

    # ======================================================
    # Cargar interacciones
    # ======================================================
    interactions = load_interactions(input_file)

    print(f"Interacciones cargadas: {len(interactions)}")

    # ======================================================
    # Construir regulón
    # ======================================================
    regulon = build_regulon(interactions)

    # ======================================================
    # Filtrar regulón
    # ======================================================
    regulon = filter_regulon_by_min_genes(regulon, min_genes)

    # ======================================================
    # Generar salida
    # ======================================================
    write_summary(regulon, output_file)

    print(f"Archivo de salida generado: {output_file}")


# ==========================================================
# Ejecutar programa
# ==========================================================
if __name__ == "__main__":
    main()

# Leer el archivo de entrada
# Contar # genes activados y reprimidos para cada TF
# Generar archivo de salida
# Calcular si un TF es activador o represor o dual

# Refactorizar el código para que sea más legible y eficiente
import os


# =========================================
# Lectura del archivo y construcción de interactions
# =========================================
# =========================================
# Responsabilidad: leer interacciones desde archivo
# Entrada: archivo
# Salida: lista de interactions
# =========================================
def load_interactions(filename):
    interactions = []
    if not os.path.exists(filename):
        print("Error: archivo no encontrado")
        exit(1)
    else:
        with open(filename) as f:
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

                fields = line.split("\t")

                # Validar número mínimo de columnas
                if len(fields) <= 6:
                    continue

                # columnas a utilizar
                TF = fields[1]
                gene = fields[4]
                effect = fields[5]

                # Validar effect
                if effect not in ["+", "-", "-+"]:
                    continue

                interactions.append((TF, gene, effect))
    return interactions


filename = "data/NetworkRegulatorGene.tsv"
interactions = load_interactions(filename)


# =========================================
# Generación de la salida
# imprimir en un archivo el resumen de cada TF
# =========================================
# Responsabilidad: generar un diccionario con información de cada TF (genes regulados, activados, reprimidos)
# Entrada: lista de interactions
# Salida: diccionario con información de cada TF (genes regulados, activados, reprimidos)
# ======
def build_regulon(interactions):
    regulon = {}  # diccionario con lista de genes
    for tf, gene, effect in interactions:
        if tf not in regulon:
            regulon[tf] = {"genes": [], "activados": 0, "reprimidos": 0}
        regulon[tf]["genes"].append(gene)

        # Contar activados y reprimidos
        if effect == "+":
            regulon[tf]["activados"] += 1
        elif effect == "-":
            regulon[tf]["reprimidos"] += 1
        elif effect == "-+":
            regulon[tf]["activados"] += 1
            regulon[tf]["reprimidos"] += 1
    return regulon


regulon = build_regulon(interactions)

# "AraC" {
#    "genes": [araC, araA, araB, araD],
#    "activados": 4,
#    "reprimidos": 0
# }


# =========
# Responsabilidad: determinar si un TF es activador, represor o dual
# Entrada: diccionario con información de cada TF (genes regulados, activados, reprimidos)
# Salida: tipo de TF (activador, represor o dual)
# =========
def get_tf_type(activados, reprimidos):

    # Determinar del TF es activador, represor o dual
    if activados > 0 and reprimidos == 0:
        tipo_tf = "Activador"
    elif activados == 0 and reprimidos > 0:
        tipo_tf = "Represor"
    elif activados > 0 and reprimidos > 0:
        tipo_tf = "Dual"
    else:
        tipo_tf = "Desconocido"
    return tipo_tf


# =========================================
# Generación de la salida
# imprimir en un archivo el resumen de cada TF
# =========================================
# =========================================
# Responsabilidad: generar un archivo con resumen de cada TF (total genes regulados, activados, reprimidos, tipo de TF, lista de genes)
# Entrada: diccionario con información de cada TF
# Salida: archivo con resumen de cada TF (total genes regulados, activados, reprimidos, tipo de TF, lista de genes)
# =========================================
def write_summary(regulon, output_filename):
    with open(output_filename, "w") as f:
        f.write(
            "TF\tTotal de genes regulados\tActivados\tReprimidos\tTipo TF\tLista de genes"
        )
        for tf in sorted(regulon):
            total_genes = len(regulon[tf]["genes"])
            lista_genes = ",".join(regulon[tf]["genes"])

            activados = regulon[tf]["activados"]
            reprimidos = regulon[tf]["reprimidos"]

            tipo_tf = get_tf_type(activados, reprimidos)

            f.write(
                f"{tf}\t{total_genes}\t{activados}\t{reprimidos}\t{tipo_tf}\t{lista_genes}\n"
            )


output_filename = "regulon_summary.tsv"
write_summary(regulon, output_filename)
print(f"Archivo de salida generado: {output_filename}")

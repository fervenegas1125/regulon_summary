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


def read_arguments():
    import argparse

    # Definición los argumentos de línea de comandos
    perser = argparse.Argumetparser(
        description="construye un resumen del regulon a partir de un archivo de interacciones."
    )
    parser.add_argument(
        "input_file", help="Nombre del archivo de entrada con las interacciones"
    )
    parser.add_argument(
        "output_file", help="Nombre del archivo de salida para el resumen del regulon"
    )

    parser.add_argument(
        "--min_genes",
        type=int,
        default=1,
        help="Numero minimo de genes regulados para incluir un regulador en el resumen (default:1)",
    )

    # ===========================================
    # main
    # ===========================================
    def main():
        # Cargar interacciones desde el archivo de TSV

        import os


import argparse

# =====================================================================
# PARTE 1: Análisis de Archivos FASTA (Casos de prueba)
# =====================================================================


def leer_fasta(ruta_archivo):
    """
    Lee un archivo en formato FASTA y recupera múltiples secuencias de ADN.

    Argumentos:
        ruta_archivo (str): Ruta del archivo FASTA a procesar.
    Regresa:
        dict: Diccionario con identificadores como llaves y secuencias como valores.
    """
    secuencias = {}
    id_actual = None
    with open(ruta_archivo, "r", encoding="utf-8") as archivo:
        for linea in archivo:
            linea = linea.strip()
            if not linea:
                continue
            if linea.startswith(">"):
                id_actual = linea[1:]
                secuencias[id_actual] = ""
            else:
                if id_actual is not None:
                    secuencias[id_actual] += linea
    return secuencias


def calcular_contenido_gc(secuencia):
    """
    Calcula el porcentaje de contenido GC de una secuencia de nucleótidos.

    Argumentos:
        secuencia (str): Cadena de caracteres con la secuencia de ADN.
    Regresa:
        float: Porcentaje de bases guanina (G) y citosina (C).
    """
    secuencia = secuencia.upper()
    total_bases = len(secuencia)
    if total_bases == 0:
        return 0.0
    conteo_g = secuencia.count("G")
    conteo_c = secuencia.count("C")
    return ((conteo_g + conteo_c) / total_bases) * 100


# =====================================================================
# PARTE 2: Procesamiento de RegulonDB (Fiel al código de clase)
# =====================================================================


def load_interactions(filename):
    """
    Responsabilidad: leer interacciones desde un archivo plano de RegulonDB.
    Entrada: ruta del archivo (str)
    Salida: lista de tuplas con interacciones (list)
    """
    interactions = []
    if not os.path.exists(filename):
        print(f"Error: archivo '{filename}' no encontrado")
        exit(1)
    else:
        with open(filename, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()

                if not line or line.startswith("#"):
                    continue

                if line.startswith("1) regulatorId") or line.startswith(
                    "1)regulatorId"
                ):
                    continue

                fields = line.split("\t")

                if len(fields) < 5:
                    continue

                TF = fields[1]
                gene = fields[3]
                effect = fields[4]

                if effect not in ["+", "-", "-+", "+-"]:
                    continue

                interactions.append((TF, gene, effect))
    return interactions


def build_regulon(interactions):
    """
    Responsabilidad: generar un diccionario agrupado con información de cada TF.
    Entrada: lista de interacciones (list)
    Salida: diccionario estructurado con los conteos por TF (dict)
    """
    regulon = {}
    for tf, gene, effect in interactions:
        if tf not in regulon:
            regulon[tf] = {"genes": [], "activados": 0, "reprimidos": 0}
        regulon[tf]["genes"].append(gene)

        if effect == "+":
            regulon[tf]["activados"] += 1
        elif effect == "-":
            regulon[tf]["reprimidos"] += 1
        elif effect in ["-+", "+-"]:
            regulon[tf]["activados"] += 1
            regulon[tf]["reprimidos"] += 1
    return regulon


def get_tf_type(activados, reprimidos):
    """
    Responsabilidad: determinar el tipo de TF (Activador, Represor o Dual).
    """
    if activados > 0 and reprimidos == 0:
        return "Activador"
    elif activados == 0 and reprimidos > 0:
        return "Represor"
    elif activados > 0 and reprimidos > 0:
        return "Dual"
    else:
        return "Desconocido"


def write_summary(regulon, output_filename):
    """
    Responsabilidad: exportar un archivo TSV con el resumen de la red.
    """
    with open(output_filename, "w", encoding="utf-8") as f:
        f.write(
            "TF\tTotal de genes regulados\tActivados\tReprimidos\tTipo TF\tLista de genes\n"
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


# =====================================================================
# COORDINADOR DEL PROGRAMA (MAIN CON ARGPARSE)
# =====================================================================


def main():
    parser = argparse.ArgumentParser(
        description="Procesador de archivos FASTA y Redes de RegulonDB."
    )
    parser.add_argument("fasta_file", help="Ruta al archivo FASTA de prueba")
    parser.add_argument(
        "input_file", help="Ruta al archivo de interacciones de RegulonDB"
    )
    parser.add_argument(
        "output_file", help="Nombre del archivo de salida para el resumen del regulon"
    )

    args = parser.parse_args()

    print("=== PARTE 1: EJECUTANDO ANÁLISIS DE CONTENIDO GC ===")
    datos_fasta = leer_fasta(args.fasta_file)
    for identificador, secuencia in datos_fasta.items():
        porcentaje = calcular_contenido_gc(secuencia)
        print(f"Secuencia: {identificador} -> Contenido GC: {porcentaje:.2f}%")

    print("\n=== PARTE 2: PROCESANDO REDES DE REGULONDB ===")
    interactions = load_interactions(args.input_file)
    print(f"Se cargaron {len(interactions)} interacciones válidas.")

    if interactions:
        regulon = build_regulon(interactions)
        write_summary(regulon, args.output_file)
        print(f"Archivo de salida generado con éxito: {args.output_file}")


if __name__ == "__main__":
    main()

    import os
import argparse

# =====================================================================
# PARTE 1: Análisis de Archivos FASTA (Casos de prueba)
# =====================================================================


def leer_fasta(ruta_archivo):
    """
    Lee un archivo en formato FASTA y recupera múltiples secuencias de ADN.

    Argumentos:
        ruta_archivo (str): Ruta del archivo FASTA a procesar.
    Regresa:
        dict: Diccionario con identificadores como llaves y secuencias como valores.
    """
    secuencias = {}
    id_actual = None
    with open(ruta_archivo, "r", encoding="utf-8") as archivo:
        for linea in archivo:
            linea = linea.strip()
            if not linea:
                continue
            if linea.startswith(">"):
                id_actual = linea[1:]
                secuencias[id_actual] = ""
            else:
                if id_actual is not None:
                    secuencias[id_actual] += linea
    return secuencias


def calcular_contenido_gc(secuencia):
    """
    Calcula el porcentaje de contenido GC de una secuencia de nucleótidos.

    Argumentos:
        secuencia (str): Cadena de caracteres con la secuencia de ADN.
    Regresa:
        float: Porcentaje de bases guanina (G) y citosina (C).
    """
    secuencia = secuencia.upper()
    total_bases = len(secuencia)
    if total_bases == 0:
        return 0.0
    conteo_g = secuencia.count("G")
    conteo_c = secuencia.count("C")
    return ((conteo_g + conteo_c) / total_bases) * 100


# =====================================================================
# PARTE 2: Procesamiento de RegulonDB (Fiel al código de clase)
# =====================================================================


def load_interactions(filename):
    """
    Responsabilidad: leer interacciones desde un archivo plano de RegulonDB.
    Entrada: ruta del archivo (str)
    Salida: lista de tuplas con interacciones (list)
    """
    interactions = []
    if not os.path.exists(filename):
        print(f"Error: archivo '{filename}' no encontrado")
        exit(1)
    else:
        with open(filename, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()

                if not line or line.startswith("#"):
                    continue

                if line.startswith("1) regulatorId") or line.startswith(
                    "1)regulatorId"
                ):
                    continue

                fields = line.split("\t")

                if len(fields) < 5:
                    continue

                TF = fields[1]
                gene = fields[3]
                effect = fields[4]

                if effect not in ["+", "-", "-+", "+-"]:
                    continue

                interactions.append((TF, gene, effect))
    return interactions


def build_regulon(interactions):
    """
    Responsabilidad: generar un diccionario agrupado con información de cada TF.
    Entrada: lista de interacciones (list)
    Salida: diccionario estructurado con los conteos por TF (dict)
    """
    regulon = {}
    for tf, gene, effect in interactions:
        if tf not in regulon:
            regulon[tf] = {"genes": [], "activados": 0, "reprimidos": 0}
        regulon[tf]["genes"].append(gene)

        if effect == "+":
            regulon[tf]["activados"] += 1
        elif effect == "-":
            regulon[tf]["reprimidos"] += 1
        elif effect in ["-+", "+-"]:
            regulon[tf]["activados"] += 1
            regulon[tf]["reprimidos"] += 1
    return regulon


def get_tf_type(activados, reprimidos):
    """
    Responsabilidad: determinar el tipo de TF (Activador, Represor o Dual).
    """
    if activados > 0 and reprimidos == 0:
        return "Activador"
    elif activados == 0 and reprimidos > 0:
        return "Represor"
    elif activados > 0 and reprimidos > 0:
        return "Dual"
    else:
        return "Desconocido"


def write_summary(regulon, output_filename):
    """
    Responsabilidad: exportar un archivo TSV con el resumen de la red.
    """
    with open(output_filename, "w", encoding="utf-8") as f:
        f.write(
            "TF\tTotal de genes regulados\tActivados\tReprimidos\tTipo TF\tLista de genes\n"
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


# =====================================================================
# COORDINADOR DEL PROGRAMA (MAIN CON ARGPARSE)
# =====================================================================


def main():
    parser = argparse.ArgumentParser(
        description="Procesador de archivos FASTA y Redes de RegulonDB."
    )
    parser.add_argument("fasta_file", help="Ruta al archivo FASTA de prueba")
    parser.add_argument(
        "input_file", help="Ruta al archivo de interacciones de RegulonDB"
    )
    parser.add_argument(
        "output_file", help="Nombre del archivo de salida para el resumen del regulon"
    )

    args = parser.parse_args()

    print("=== PARTE 1: EJECUTANDO ANÁLISIS DE CONTENIDO GC ===")
    datos_fasta = leer_fasta(args.fasta_file)
    for identificador, secuencia in datos_fasta.items():
        porcentaje = calcular_contenido_gc(secuencia)
        print(f"Secuencia: {identificador} -> Contenido GC: {porcentaje:.2f}%")

    print("\n=== PARTE 2: PROCESANDO REDES DE REGULONDB ===")
    interactions = load_interactions(args.input_file)
    print(f"Se cargaron {len(interactions)} interacciones válidas.")

    if interactions:
        regulon = build_regulon(interactions)
        write_summary(regulon, args.output_file)
        print(f"Archivo de salida generado con éxito: {args.output_file}")


if __name__ == "__main__":
    main()

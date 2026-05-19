# Analizador de Archivos Biológicos - FASTA y RegulonDB

Este repositorio contiene las herramientas desarrolladas en Python para automatizar el procesamiento de datos genómicos utilizando un diseño basado en funciones y manejo de argumentos por linea de comandos.

## Módulos del Proyecto

1. **Análisis FASTA (`main.py`):** Carga un archivo con múltiples secuencias de ADN, procesa los encabezados y calcula con exactitud el porcentaje de contenido GC de cada una.
2. **Procesamiento de Regulones:** Clasifica los factores de transcripción (TF) extraídos de RegulonDB en las categorías de *Activador*, *Represor* o *Dual*.
3. **Módulo de Argumentos (`ejemplo.py`):** Implementación fiel de la sesión de clase utilizando la librería `argparse` para capturar parámetros de entrada y salida desde la terminal.

## Requisitos de Ejecución
Este proyecto utiliza `uv` como gestor de entornos y paquetes.

Para ejecutar el programa principal y mapear los archivos biológicos, corre:
```bash
uv run python main.py secuencia.fasta data/NetworkRegulatorGene.tsv regulon_summary.tsv
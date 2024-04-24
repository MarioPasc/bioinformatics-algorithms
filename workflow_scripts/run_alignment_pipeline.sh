#!/bin/bash

# Rutas relativas
EXECUTABLE="./build/testSW"
PYTHON_SCRIPT="./Alignment/alignment_visualization.py"
OUTPUT_DIR="./workflow_scripts"
TEMP_FILE="$OUTPUT_DIR/temp_output.txt"

# Asegúrate de que las carpetas necesarias existen
mkdir -p "$OUTPUT_DIR"
mkdir -p "./images"

# Ejecutar el programa C++ y guardar la salida
$EXECUTABLE > "$TEMP_FILE"

# Verifica si el archivo de salida tiene contenido antes de continuar
if [ -s "$TEMP_FILE" ]; then
    # Usar la salida del programa C++ como entrada para el script Python
    python3 "$PYTHON_SCRIPT" < "$TEMP_FILE"
else
    echo "Error: Empty output file."
fi

# Limpieza (opcional)
rm "$TEMP_FILE"

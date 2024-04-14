cmake_minimum_required(VERSION 3.10)
project(Bioinformatics-Algorithms)

# Configuración de compilación
set(CMAKE_CXX_STANDARD 14)

# Directorios de inclusión
include_directories(Assembly/De_Brujin_Graphs/include)

# Archivos de origen
set(ALIGNMENT_SOURCES Assembly/De_Brujin_Graphs/src/assembly.cpp Assembly/De_Brujin_Graphs/src/graph.cpp)

# Ejecutable principal
add_executable(main main.cpp ${ALIGNMENT_SOURCES} ${SEQUENCING_SOURCES})

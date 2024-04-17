cmake_minimum_required(VERSION 3.10)
project(Bioinformatics-Algorithms)

# Configuración de compilación
set(CMAKE_CXX_STANDARD 14)
find_package(ZLIB REQUIRED)
# Directorios de inclusión para tu proyecto
include_directories(Assembly/De_Brujin_Graphs/include)

# Incluir los directorios de las cabeceras de kseq++
include_directories(external/kseqpp/include)

# Archivos de origen
set(SOURCES "main.cpp" "Assembly/De_Brujin_Graphs/src/graph.cpp")

# Ejecutable principal
add_executable(main ${SOURCES})
target_link_libraries(main PRIVATE ZLIB::ZLIB)


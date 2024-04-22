cmake_minimum_required(VERSION 3.10)
project(Bioinformatics-Algorithms)

# Configuración de compilación
set(CMAKE_CXX_STANDARD 14)
find_package(ZLIB REQUIRED)

# Directorios de inclusión para tu proyecto
include_directories(
    Assembly/De_Brujin_Graphs/include
    Alignment/NeedlemanWunsch/include # Añadir el nuevo módulo
    external/kseqpp/include
)

# Archivos de origen
set(SOURCES
    main.cpp
    Assembly/De_Brujin_Graphs/src/graph.cpp
    Alignment/NeedlemanWunsch/src/needleman_wunsch.cpp # Añadir el archivo de implementación
)

# Ejecutable principal
add_executable(main ${SOURCES})
target_link_libraries(main PRIVATE ZLIB::ZLIB)

# Crear un ejecutable para el test de NeedlemanWunsch
add_executable(testNW Alignment/NeedlemanWunsch/testWN/test_needleman.cpp)
target_sources(testNW PRIVATE Alignment/NeedlemanWunsch/src/needleman_wunsch.cpp)
cmake_minimum_required(VERSION 3.10)
project(Bioinformatics-Algorithms)

# Configuración de compilación
set(CMAKE_CXX_STANDARD 14)
find_package(ZLIB REQUIRED)

# Directorios de inclusión para tu proyecto
include_directories(
    Assembly/De_Brujin_Graphs/include
    Alignment/NeedlemanWunsch/include 
    Alignment/SmithWaterman/include 
    Alignment/MultipleSequenceAlignment/include 
    external/kseqpp/include
)

# Archivos de origen
set(SOURCES
    main.cpp
    Assembly/De_Brujin_Graphs/src/graph.cpp
    Alignment/NeedlemanWunsch/src/needleman_wunsch.cpp
    Alignment/SmithWaterman/src/smith_waterman.cpp
    Alignment/MultipleSequenceAlignment/src/neighbour_joining.cpp
)

# Ejecutable principal
add_executable(main ${SOURCES})
target_link_libraries(main PRIVATE ZLIB::ZLIB)

# Crear un ejecutable para el test de NeedlemanWunsch
add_executable(testNW Alignment/NeedlemanWunsch/testWN/test_needleman.cpp)
target_sources(testNW PRIVATE Alignment/NeedlemanWunsch/src/needleman_wunsch.cpp)

# Crear un ejecutable para el test de SmithWaterman
add_executable(testSW Alignment/SmithWaterman/testSW/test_smith.cpp)
target_sources(testSW PRIVATE Alignment/SmithWaterman/src/smith_waterman.cpp)

# Crear un ejecutable para el test de Neighbour Joining
add_executable(testNJ Alignment/MultipleSequenceAlignment/testNJ/test_neighbour_joining.cpp)
target_sources(testNJ PRIVATE Alignment/MultipleSequenceAlignment/src/neighbour_joining.cpp)

add_library(needleman_wunsch STATIC
Alignment/NeedlemanWunsch/src/needleman_wunsch.cpp
)

target_link_libraries(testNJ needleman_wunsch)
cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(2x2MCP_DY LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(ROOT CONFIG REQUIRED COMPONENTS Core RIO Tree Physics)
find_package(Threads REQUIRED)

find_path(PYTHIA8_INCLUDE_DIR
    NAMES Pythia8/Pythia.h
    HINTS $ENV{PYTHIA8}/include $ENV{CONDA_PREFIX}/include
)

find_library(PYTHIA8_LIBRARY
    NAMES pythia8
    HINTS
        $ENV{PYTHIA8}/lib
        $ENV{PYTHIA8}/lib64
        $ENV{CONDA_PREFIX}/lib
)

if(NOT PYTHIA8_INCLUDE_DIR OR NOT PYTHIA8_LIBRARY)
    message(FATAL_ERROR "Could not find Pythia8 in PYTHIA8 or CONDA_PREFIX")
endif()

add_library(Pythia8::Pythia8 UNKNOWN IMPORTED)
set_target_properties(Pythia8::Pythia8 PROPERTIES
    IMPORTED_LOCATION "${PYTHIA8_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${PYTHIA8_INCLUDE_DIR}"
)

add_executable(dy_mcp_production dy_mcp_production.cc)
target_link_libraries(dy_mcp_production PRIVATE
    ROOT::Core
    ROOT::RIO
    ROOT::Tree
    ROOT::Physics
    Threads::Threads
    Pythia8::Pythia8
)

configure_file(beam.config beam.config COPYONLY)
configure_file(momentum.config momentum.config COPYONLY)

message(STATUS "Using PYTHIA8 include dir: ${PYTHIA8_INCLUDE_DIR}")
message(STATUS "Using PYTHIA8 library    : ${PYTHIA8_LIBRARY}")

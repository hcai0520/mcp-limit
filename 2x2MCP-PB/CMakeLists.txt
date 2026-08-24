cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(2x2MCP_PB LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(ROOT CONFIG REQUIRED COMPONENTS Core RIO Tree Physics)

add_executable(pb_fww_production pb_fww_production.cc)
target_link_libraries(pb_fww_production PRIVATE
    ROOT::Core
    ROOT::RIO
    ROOT::Tree
    ROOT::Physics
)


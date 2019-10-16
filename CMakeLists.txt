cmake_minimum_required(VERSION 3.13)
project(memristor)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS "-Wall -Wpedantic -Werror=return-type -m64 -march=native -fmax-errors=1")
set(CMAKE_CXX_FLAGS_DEBUG "-Wextra -O0 -ggdb3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

file(GLOB_RECURSE MEMRISTOR_SRC "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
add_executable(memristor "${MEMRISTOR_SRC}")

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external/eslib")
target_link_libraries(memristor PRIVATE eslib)

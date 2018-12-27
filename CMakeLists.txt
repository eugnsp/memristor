cmake_minimum_required(VERSION 3.5)
project(memristor)

set(CMAKE_CXX_STANDARD 				17)
set(CMAKE_CXX_STANDARD_REQUIRED 	ON)
set(CMAKE_CXX_EXTENSIONS 			OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS 	ON)
set(CMAKE_CXX_COMPILER             	"/usr/bin/g++")
#set(CMAKE_CXX_COMPILER             	"/usr/bin/clang++")
set(CMAKE_CXX_FLAGS                	"-Wall -Wextra -Wpedantic -Wno-psabi -m64 -march=native")
set(CMAKE_CXX_FLAGS_DEBUG          	"-O0 -ggdb3 -D_DEBUG")
set(CMAKE_CXX_FLAGS_RELEASE        	"-O3 -DNDEBUG -flto")

#add_compile_options(-Wall -Wextra -Wpedantic -g -m64 -D_DEBUG -fmax-errors=5)
#add_compile_options(-Wall -Wextra -Wpedantic -m64 -O3 -DNDEBUG -DMKL_ILP64 -flto)

set(ESLIB_ROOT "/home/es/dev/projects")

add_library(es_fe STATIC IMPORTED)
set_target_properties(es_fe PROPERTIES
	IMPORTED_LOCATION "${ESLIB_ROOT}/es_fe/build/libes_fe.a"
	INTERFACE_INCLUDE_DIRECTORIES "${ESLIB_ROOT}/es_fe/include")

add_library(es_util INTERFACE IMPORTED)
set_target_properties(es_util PROPERTIES
	INTERFACE_INCLUDE_DIRECTORIES "${ESLIB_ROOT}/es_util/include")

add_library(es_la STATIC IMPORTED)
set_target_properties(es_la PROPERTIES
	IMPORTED_LOCATION "${ESLIB_ROOT}/es_la/build/libes_la.a"
	INTERFACE_INCLUDE_DIRECTORIES "${ESLIB_ROOT}/es_la/include")

add_library(es_math INTERFACE IMPORTED)
set_target_properties(es_math PROPERTIES
	INTERFACE_INCLUDE_DIRECTORIES "${ESLIB_ROOT}/es_math/include")

add_library(es_phys INTERFACE IMPORTED)
set_target_properties(es_phys PROPERTIES
	INTERFACE_INCLUDE_DIRECTORIES "${ESLIB_ROOT}/es_phys/include")

add_library(es_geom STATIC IMPORTED)
set_target_properties(es_geom PROPERTIES
	IMPORTED_LOCATION "${ESLIB_ROOT}/es_geom/build/libes_geom.a"
	INTERFACE_INCLUDE_DIRECTORIES "${ESLIB_ROOT}/es_geom/include")

set(MKLROOT $ENV{MKLROOT})
include_directories("${MKLROOT}/include")
link_directories("${MKLROOT}/lib/intel64")
add_definitions(-DMKL_ILP64)
set(MKL_LIBRARIES "-Wl,--start-group" libmkl_intel_ilp64.a libmkl_sequential.a libmkl_core.a "-Wl,--end-group" pthread m dl)

file(GLOB_RECURSE SRC ${PROJECT_SOURCE_DIR}/src/*.cpp)
add_executable(memristor ${SRC})
target_link_libraries(memristor es_util es_phys es_la es_fe es_geom es_math ${MKL_LIBRARIES})

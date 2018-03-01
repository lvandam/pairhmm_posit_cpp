cmake_minimum_required(VERSION 3.0)
project(pairhmm)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-DPOSIT_VERBOSE_OUTPUT")# -DDEBUG_VALUES -DDEBUG_VERBOSE #-DPOSIT_TRACE_DECODE -DPOSIT_TRACE_CONVERSION -DPOSIT_TRACE_ROUNDING")

set(EXTERNAL_INSTALL_LOCATION ${CMAKE_BINARY_DIR}/external)

include(ExternalProject)
ExternalProject_Add(universal
        GIT_REPOSITORY https://github.com/lvandam/universal
        GIT_TAG master
        PREFIX ${EXTERNAL_INSTALL_LOCATION}
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION}
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
        )

ExternalProject_Add(gnuplot-iostream
        GIT_REPOSITORY https://github.com/dstahlke/gnuplot-iostream
        GIT_TAG master
        PREFIX ${EXTERNAL_INSTALL_LOCATION}
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION}
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
        )

include_directories(${EXTERNAL_INSTALL_LOCATION}/src/universal)
include_directories(${EXTERNAL_INSTALL_LOCATION}/src/gnuplot-iostream)

file(GLOB pairhmm_SRC
        "src/pairhmm/*.hpp"
        "src/pairhmm/*.cpp"
        )

file(GLOB benchmarks_SRC
        "src/benchmarks/*.hpp"
        "src/benchmarks/*.cpp"
        )

add_executable(pairhmm ${pairhmm_SRC})
add_executable(phred_to_posit ${benchmarks_SRC})

target_link_libraries(phred_to_posit LINK_PUBLIC boost_iostreams)
target_link_libraries(phred_to_posit LINK_PUBLIC boost_system)
target_link_libraries(phred_to_posit LINK_PUBLIC boost_filesystem)

add_dependencies(pairhmm universal)
add_dependencies(phred_to_posit universal)
add_dependencies(phred_to_posit gnuplot-iostream)

cmake_minimum_required(VERSION 3.0)
project(pairhmm_simple)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-DPOSIT_VERBOSE_OUTPUT")# -DDEBUG_VALUES -DDEBUG_VERBOSE #-DPOSIT_TRACE_DECODE -DPOSIT_TRACE_CONVERSION -DPOSIT_TRACE_ROUNDING")

set(EXTERNAL_INSTALL_LOCATION ${CMAKE_BINARY_DIR}/external)

include(ExternalProject)
ExternalProject_Add(universal
        GIT_REPOSITORY https://github.com/stillwater-sc/universal
        GIT_TAG master
        PREFIX ${EXTERNAL_INSTALL_LOCATION}
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION}
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
        )

include_directories(${EXTERNAL_INSTALL_LOCATION}/src/universal/posit)

file(GLOB project_SRC
        "src/*.hpp"
        "src/*.cpp"
        )

add_executable(pairhmm_simple ${project_SRC})
add_dependencies(pairhmm_simple universal)
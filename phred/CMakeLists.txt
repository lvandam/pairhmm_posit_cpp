cmake_minimum_required(VERSION 3.0)
project(phred)

set(CMAKE_CXX_STANDARD 11)

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

include_directories(${EXTERNAL_INSTALL_LOCATION}/src/universal)

file(GLOB pairhmm_SRC
        "*.hpp"
        "main.cpp"
        )

# Host (CAPI)
add_executable(phred ${pairhmm_SRC})
add_dependencies(phred universal)

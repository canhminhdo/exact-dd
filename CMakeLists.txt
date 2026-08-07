cmake_minimum_required(VERSION 3.19...3.30)

project(exact-dd
    VERSION 0.1.0
    DESCRIPTION "Exact/weighted decision diagrams (D[w]) for quantum computation"
    HOMEPAGE_URL "https://github.com/canhminhdo/exact-dd"
    LANGUAGES CXX
)

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
    message(FATAL_ERROR "In-source builds are not allowed")
endif()

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# set c++ standard for the project
include(cmake/CompilerOptions.cmake)
# check external dependencies
include(cmake/ExternalDependencies.cmake)
# package add test macro
include(cmake/PackageAddTests.cmake)

# check whether this is the master project or consumed via add_subdirectory
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    set(EXACT_DD_MASTER_PROJECT ON)
else()
    set(EXACT_DD_MASTER_PROJECT OFF)
endif()

add_subdirectory(src)

option(BUILD_EXACT_DD_TESTS "Build tests for exact-dd" ${EXACT_DD_MASTER_PROJECT})
if(BUILD_EXACT_DD_TESTS)
    enable_testing()
    include(GoogleTest)
    add_subdirectory(test)
endif()

# Off by default even for standalone builds: exporting ExactDD requires the
# FetchContent-provided Boost.Multiprecision target and nlohmann_json to be part of the same
# export set, which install(EXPORT ...) does not support out of the box. The
# primary way to consume this library (add_subdirectory as a git submodule,
# as QRAT does) doesn't need install/find_package support at all; opt in with
# -DEXACT_DD_INSTALL=ON only if you've resolved that export limitation --
# which also means adding find_dependency(nlohmann_json) to
# cmake/exact-dd-config.cmake.in and recording the DD_EXACT_STATISTICS
# definition in the export, since it must match in every consumer.
option(EXACT_DD_INSTALL "Generate installation instructions for exact-dd" OFF)
if(EXACT_DD_INSTALL)
    include(cmake/InstallAndPackage.cmake)
endif()

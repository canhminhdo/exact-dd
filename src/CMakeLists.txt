add_library(
    ExactDD STATIC
    dd/exact/Dw.cpp
    dd/exact/DwPackage.cpp
    dd/exact/DwGateMatrixDefinitions.cpp
    dd/exact/ExactDDSimulation.cpp
    utility/StringHelper.cpp
)
add_library(ExactDD::Core ALIAS ExactDD)

include(GNUInstallDirs)
target_include_directories(
    ExactDD PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_features(ExactDD PUBLIC cxx_std_17)
target_link_libraries(ExactDD PUBLIC Boost::multiprecision)
# Warning suppression policy should not be forced on consumers
target_compile_options(ExactDD PRIVATE -Wno-deprecated-literal-operator)

# Opt-in GMP-backed Integer (dd::exact::Dw's arbitrary-precision type)
option(EXACT_DD_WITH_GMP "Whether Dw's Integer type is backed by GMP instead of Boost's cpp_int" ON)
if(EXACT_DD_WITH_GMP)
    if(APPLE AND CMAKE_OSX_ARCHITECTURES)
        list(LENGTH CMAKE_OSX_ARCHITECTURES EXACT_DD_NUM_OSX_ARCHS)
        if(EXACT_DD_NUM_OSX_ARCHS GREATER 1)
            # Homebrew's GMP is single-architecture, incompatible with a
            # universal (multi-arch) macOS build -- fail with a clear,
            # actionable message at configure time rather than a cryptic
            # "symbol(s) not found for architecture ..." link error.
            message(FATAL_ERROR "EXACT_DD_WITH_GMP=ON is incompatible with a universal "
                                 "(multi-arch) macOS build, since Homebrew's GMP is single-"
                                 "architecture. Reconfigure with a single "
                                 "-DCMAKE_OSX_ARCHITECTURES=<your-arch> (e.g. arm64 or x86_64) "
                                 "alongside -DEXACT_DD_WITH_GMP=ON.")
        endif()
    endif()
    find_package(GMP REQUIRED)
    target_compile_definitions(ExactDD PUBLIC DD_EXACT_WITH_GMP)
    target_link_libraries(ExactDD PUBLIC GMP::gmp GMP::gmpxx)
endif()

option(BUILD_EXACT_DD_EXAMPLES "Build example executables for exact-dd" ${EXACT_DD_MASTER_PROJECT})
if(BUILD_EXACT_DD_EXAMPLES)
    add_subdirectory(examples)
endif()

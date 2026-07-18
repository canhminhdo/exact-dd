# Declare all external dependencies and make sure that they are available.

include(FetchContent)
set(FETCH_PACKAGES "")

find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    # Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE AND NOT EXISTS "${PROJECT_SOURCE_DIR}/extern/googletest/CMakeLists.txt")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

option(USE_SYSTEM_BOOST "Whether to try to use the system Boost installation" OFF)
set(BOOST_MIN_VERSION
    1.80.0
    CACHE STRING "Minimum required Boost version")
if(USE_SYSTEM_BOOST)
    find_package(Boost ${BOOST_MIN_VERSION} CONFIG REQUIRED)
else()
    set(BOOST_VERSION
        1_89_0
        CACHE INTERNAL "Boost version")
    set(BOOST_URL
        https://github.com/boostorg/multiprecision/archive/refs/tags/Boost_${BOOST_VERSION}.tar.gz)
    if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.24)
        FetchContent_Declare(boost_mp URL ${BOOST_URL} FIND_PACKAGE_ARGS ${BOOST_MIN_VERSION} CONFIG
                                          NAMES boost_multiprecision)
        list(APPEND FETCH_PACKAGES boost_mp)
    else()
        find_package(boost_mp ${BOOST_MIN_VERSION} QUIET CONFIG NAMES boost_multiprecision)
        if(NOT boost_mp_FOUND)
            FetchContent_Declare(boost_mp URL ${BOOST_URL})
            list(APPEND FETCH_PACKAGES boost_mp)
        endif()
    endif()
endif()

# Make all declared dependencies available.
FetchContent_MakeAvailable(${FETCH_PACKAGES})

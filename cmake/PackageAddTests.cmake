macro(PACKAGE_ADD_TESTS testname linklibs)
    if(NOT TARGET ${testname})
        # create an executable in which the tests will be stored
        add_executable(${testname} ${ARGN})
        # link the Google test infrastructure and a default main function to the test executable.
        target_link_libraries(${testname} PRIVATE ${linklibs} gmock gtest_main)
        # discover tests
        gtest_discover_tests(
            ${testname}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" DISCOVERY_TIMEOUT 60)
        set_target_properties(${testname} PROPERTIES FOLDER test)
    endif()
endmacro()

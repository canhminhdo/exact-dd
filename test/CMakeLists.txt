# building googletest as a part of project
if(NOT TARGET gtest OR NOT TARGET gmock)
    # Prevent overriding the parent project's compiler/linker settings on Windows
    set(gtest_force_shared_crt # cmake-lint: disable=C0103
        ON
        CACHE BOOL "" FORCE)
    add_subdirectory("${PROJECT_SOURCE_DIR}/extern/googletest" "extern/googletest" EXCLUDE_FROM_ALL)
    set_target_properties(gtest gtest_main gmock gmock_main PROPERTIES FOLDER extern)
endif()

message("Building tests for exact-dd project")

package_add_tests(
    exact_dd_test "ExactDD::Core"
    test_dw.cpp
    test_dwpackage.cpp
    test_gate_builders.cpp
    test_gate_matrices.cpp
    test_exact_dd_simulation.cpp
    test_grover.cpp
    test_qft.cpp
    test_shor.cpp
)

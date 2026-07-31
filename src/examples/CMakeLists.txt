add_executable(ExactDDMain main.cpp)
target_link_libraries(ExactDDMain PUBLIC ExactDD::Core)

# Scaling benchmark: prints wall time / node counts per qubit count, so a
# performance change can be judged against a recorded baseline.
add_executable(ExactDDBench bench.cpp)
target_link_libraries(ExactDDBench PUBLIC ExactDD::Core)

# exact-dd

`exact-dd` is a C++17 library implementing `DwPackage`, an exact/weighted decision diagram (DD)
package for representing and manipulating quantum states and operators without floating-point
rounding error. Amplitudes are represented in the ring Q[w] (dyadic rationals adjoined with
w = e^{i*pi/4}) via the `dd::exact::Dw` scalar type, backed by Boost.Multiprecision (or GMP,
optionally).

Type hierarchy:

```
Dw (exact scalar)
  -> DwEdge<Node> (weighted edge)
  -> DwVNode / DwMNode (branching node, 2 or 4 children)
  -> DwPackage (hash-consing + pooled allocation + refcounted GC on top)
```

`exact-dd` is structurally parallel to [MQT Core](https://github.com/cda-tum/mqt-core)'s
`dd::Package`, but is an independent, from-scratch implementation — it does not depend on or link
against MQT Core.

See `specs/` for design notes on the memory manager, reference-counting GC scheme, and the
gate-DD construction algorithms.

## Build & test

Requires the `extern/googletest` submodule, CMake >= 3.19, and a C++17 compiler.

```shell
git submodule update --init --recursive
cmake -S . -B build -D CMAKE_BUILD_TYPE=Release
cmake --build build --config Release -j 8
ctest --test-dir build
```

Boost.Multiprecision is fetched automatically via `FetchContent` (see
`cmake/ExternalDependencies.cmake`); pass `-DUSE_SYSTEM_BOOST=ON` to use a system installation
instead.

### Options

- `BUILD_EXACT_DD_TESTS` (default `ON` when built standalone) — build the GoogleTest suite.
- `BUILD_EXACT_DD_EXAMPLES` (default `ON` when built standalone) — build the `ExactDDMain` demo
  executable (a Grover's-search circuit driver) under `src/examples/`.
- `EXACT_DD_WITH_GMP` (default `ON`) — back `Dw`'s arbitrary-precision integer type with GMP
  instead of Boost's `cpp_int`.
- `EXACT_DD_INSTALL` (default `ON` when built standalone) — generate install rules and a CMake
  package config for `find_package(exact-dd)`.

## Using from another CMake project

As a git submodule (recommended, mirrors how this library itself consumes `googletest`):

```cmake
add_subdirectory("${PROJECT_SOURCE_DIR}/extern/exact-dd" "extern/exact-dd" EXCLUDE_FROM_ALL)
target_link_libraries(your_target PUBLIC ExactDD::Core)
```

# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`exact-dd` is a C++17 library implementing `DwPackage`, an exact/weighted decision-diagram (DD)
package for representing and manipulating quantum states and operators over the Clifford+T gate
set with **no floating-point rounding error**. Amplitudes live in the ring Q[w] (dyadic rationals
adjoined with w = e^{i*pi/4}), backed by Boost.Multiprecision or GMP.

It is structurally similar to [MQT Core](https://github.com/cda-tum/mqt-core)'s `dd::Package`
but is an independent, from-scratch implementation —
it does not depend on or link against MQT Core.

## Build & test

Requires the `extern/googletest` submodule, CMake >= 3.19, and a C++17 compiler.

```shell
git submodule update --init --recursive
cmake -S . -B build -D CMAKE_BUILD_TYPE=Release
cmake --build build --config Release -j 8
ctest --test-dir build
```

Boost.Multiprecision and nlohmann_json are fetched automatically via `FetchContent`
(`cmake/ExternalDependencies.cmake`); pass `-DUSE_SYSTEM_BOOST=ON` to use a system Boost instead.
The nlohmann_json block is skipped when the target already exists. In-source builds are
rejected by CMakeLists.txt.

Run a single test binary or filter to one test case with GoogleTest's built-in flag:

```shell
./build/test/exact_dd_test --gtest_filter='Dw.OmegaToFourthIsMinusOne'
./build/test/exact_dd_test --gtest_filter='DwPackage.*'
```

(ctest also discovers each `TEST(...)` individually via `gtest_discover_tests`, so
`ctest --test-dir build -R OmegaToFourthIsMinusOne` works too.)

### CMake options

- `BUILD_EXACT_DD_TESTS` (default `ON` standalone) — build the GoogleTest suite (`test/`).
- `BUILD_EXACT_DD_EXAMPLES` (default `ON` standalone) — build the `ExactDDMain` demo executable
  (Grover's-search / Hadamard-transform / GHZ-entanglement driver) at `src/examples/main.cpp`,
  and `ExactDDBench`, the scaling benchmark at `src/examples/bench.cpp` (prints wall time and
  node counts per qubit count for `innerProduct`, `outerProduct`, a GHZ measure sweep, the
  n-fold Hadamard transform, and Grover — use it to check a performance change against a
  recorded baseline).
- `EXACT_DD_WITH_GMP` (default `ON`) — back `Dw`'s arbitrary-precision integer type with GMP
  instead of Boost's `cpp_int`. Incompatible with a universal (multi-arch) macOS build — pin
  `-DCMAKE_OSX_ARCHITECTURES` to a single arch alongside it. GMP wins for large integers, while `cpp_int` wins for small integers.
- `EXACT_DD_STATISTICS` (default `ON`) — collect DD package statistics (unique tables, memory
  managers, compute tables) and report them via `DwPackage::statistics()`. Sets the PUBLIC compile
  definition `DD_EXACT_STATISTICS`; when OFF, every tracking call compiles to nothing and
  `PackageStatistics::json()` emits a `statistics_disabled` marker instead of a tree of zeros.
  Turn it off when recording performance baselines with `ExactDDBench`, since the counters add an
  increment per table probe on the hottest path.
- `EXACT_DD_INSTALL` (default `OFF`) — generate install rules / CMake package config. Left off by
  default even standalone: the FetchContent-provided Boost and nlohmann_json targets can't
  currently be part of the same `install(EXPORT ...)` set. The primary consumption path
  (add_subdirectory as a submodule) doesn't need this.

When consumed via `add_subdirectory` from another project (`EXACT_DD_MASTER_PROJECT` becomes
`OFF`), tests and examples are not built unless explicitly turned on by the parent.

## Architecture

### Type hierarchy

```
Dw (exact scalar, include/dd/exact/Dw.hpp)
  -> DwEdge<Node> (weighted edge, include/dd/exact/DwNode.hpp)
  -> DwVNode / DwMNode (branching node, 2 or 4 children)
  -> DwPackage (hash-consing + pooled allocation + refcounted GC on top)
```

**`Dw`** (`include/dd/exact/Dw.hpp`, `src/dd/exact/Dw.cpp`) is the exact scalar type. A value is
stored as `(a + b*w + c*w^2 + d*w^3) / (sqrt(2)^k * e)` with `a,b,c,d` arbitrary-precision
integers, `k` a non-negative exponent, and `e` a positive odd integer denominator — this covers
D[w] (e=1) and, beyond it, all of Q[w] (needed for `inverse()`). Every instance is kept in
canonical form (minimal `k`, odd positive `e`, `gcd(a,b,c,d,e) == 1`) so `operator==`/`hash()` are
exact with zero floating-point tolerance; see `Dw::canonicalize()` and its two phases
(`reduceSqrt2Power`, `reduceRationalDenominator`). Notable non-obvious operations: `inverse()`
(field inverse in Q[w], paper's Algorithm 2), `gcd()`/`reduceAssociate()` (D[w]-only, Euclidean-ring
GCD via `quarticNorm()`, paper's Algorithm 3) — these back `DwPackage`'s node-weight normalization.

**`DwVNode`/`DwMNode`/`DwEdge<Node>`** (`include/dd/exact/DwNode.hpp`) are the DD nodes: `DwVNode`
(vector/state DDs) has 2 children indexed by a qubit's basis value; `DwMNode` (matrix/operator DDs)
has 4, indexed `[M00, M01, M10, M11]` by `(row bit, col bit)`. `var` is the qubit index; the
terminal is `p == nullptr` with the edge weight carrying the scalar value. `ref`/`next` exist only
for `MemoryManager`'s pooling/GC bookkeeping and never participate in edge equality/hashing.

**`MemoryManager<T>`** (`include/dd/exact/MemoryManager.hpp`, header-only) is a pooled bump
allocator ported from MQT Core: nodes are carved from growing chunks, freed nodes go on a
singly-linked free list threaded through `T::next`. Chunks are stored as `vector<vector<T>>` so
growth *moves* (never copies) existing chunks — pointers previously returned by `get()` stay valid.

**`DwPackage`** (`include/dd/exact/DwPackage.hpp`, `src/dd/exact/DwPackage.cpp`) is the central
class: fixed-qubit-count vector/matrix DD package with hash-consed unique tables
(`vUnique_`/`mUnique_`, keyed on `var` + children) and reference-counted, manual/opt-in garbage
collection (`incRef`/`decRef`/`garbageCollect`), mirroring MQT Core's `Package`/`UniqueTable`
design. Key things to know when touching it:
- **Node-weight normalization** (`NormalizationStrategy`: `None`/`Inverse`/`Gcd`, constructor
  parameter, default `Inverse`) controls whether `makeVEdge`/`makeMEdge` factor a common scalar out
  of a node's outgoing weights so structurally-scaled submatrices collapse to one shared node. Only
  node pointers are hash-consed/refcounted — edge weights (`Dw`) are not, since they're heavyweight
  arbitrary-precision values, not simple doubles.
- Gate-DD builders come in two shapes: single-control builders are top-down; the general
  multi-controlled builders (`makeControlledSingleQubitGateDD`/`makeControlledTwoQubitGateDD` taking
  a `vector<size_t>` of controls) are built bottom-up like MQT Core's `makeGateDD`.
- `makeControlledTwoQubitGateDD` has both a scalar-`control` overload and a `vector<size_t>`-
  `controls` overload; calling either with a bare `{}` resolves to the *scalar* overload (control=0,
  not "uncontrolled") — pass `std::vector<size_t>{}` explicitly for an uncontrolled gate.
- `kronecker()` requires both operands to belong to the *same* `DwPackage`; there is no cross-
  package / extended-qubit-count support.
- Measurement (`measureOneQubit`) deliberately does **not** renormalize the post-measurement state
  (dividing by `sqrt(probability)` generally isn't exact in D[w]); callers track probability as an
  exact `Dw` ratio and only convert to `double` at the final step (see `fidelity()`).
- `DwPackage` is move-only (copying would leave unique tables pointing at the original's pooled
  node storage).

**Statistics** (`include/dd/exact/statistics/`, `src/dd/exact/statistics/`) mirrors MQT Core's
`dd/statistics/` module: `Statistics` (JSON/`toString`/`operator<<` base) → `TableStatistics`
(hash-table counters) → `UniqueTableStatistics` (adds active entries + gc runs), plus the templated
`MemoryManagerStatistics<T>` and the `PackageStatistics` aggregate that `DwPackage::statistics()`
returns. Three things differ from MQT Core and are documented at length in the headers: the report
is a `const` member returning a by-value snapshot (DwPackage's members are private, unlike MQT's
`Package`); `collisions` is a current-state bucket measure computed at report time rather than MQT's
cumulative chain-step count (an insert-time probe would re-hash a `VKey`/`MKey` on the hot path);
and every MiB figure is a lower bound, since each `Dw`'s arbitrary-precision limbs live on the heap
where `sizeof` cannot see them. All tracking compiles away under `EXACT_DD_STATISTICS=OFF`.

**`ExactDDSimulation`** (`include/dd/exact/ExactDDSimulation.hpp`,
`src/dd/exact/ExactDDSimulation.cpp`) is a thin standalone driver on top of `DwPackage`: owns one
`pkg_` and one persistent, incRef'd `state_` root, and applies gates/measurements by name through
`DwPackage::applyOperation` (which incRef's the new state, decRef's the old, and opportunistically
garbage-collects). This is *not* wired into QRAT's own pipeline — see `src/examples/main.cpp` for
example circuits (Grover's search, Hadamard transform, GHZ entanglement) built directly on this
class.

**`DwGateMatrixDefinitions`** (`include/dd/exact/DwGateMatrixDefinitions.hpp`,
`src/dd/exact/DwGateMatrixDefinitions.cpp`) supplies exact D[w] matrices for the full Clifford+T
single-qubit set (`i,x,y,z,h,s,sdg,t,tdg,v,vdg,sx,sxdg` — all 13 are exactly representable in D[w])
plus four two-qubit gates (`swap,iswap,iswapdg,dcx`), looked up by name via `byName()`/
`twoQubitByName()`. Gate names match QRAT's DSL naming. `cx`/`cz` are deliberately absent from
`twoQubitByName` — build controlled gates via `DwPackage::makeControlledSingleQubitGateDD` instead.
Throws `std::invalid_argument` for any name outside the supported set rather than silently
approximating.

### Test layout

`test/` is one GoogleTest binary (`exact_dd_test`, registered via the `package_add_tests` macro in
`cmake/PackageAddTests.cmake`) covering: `Dw` arithmetic/canonicalization (`test_dw.cpp`),
`DwPackage` core operations (`test_dwpackage.cpp`), gate-DD builders (`test_gate_builders.cpp`),
gate matrix correctness (`test_gate_matrices.cpp`), the `ExactDDSimulation` driver
(`test_exact_dd_simulation.cpp`), and full algorithm-level circuits (`test_grover.cpp`,
`test_qft.cpp` with shared helpers in `qft_helper.hpp`, `test_shor.cpp`).

## Using from another CMake project

As a git submodule (the intended/only well-supported path, mirroring how this library itself
consumes `googletest`):

```cmake
add_subdirectory("${PROJECT_SOURCE_DIR}/extern/exact-dd" "extern/exact-dd" EXCLUDE_FROM_ALL)
target_link_libraries(your_target PUBLIC ExactDD::Core)
```

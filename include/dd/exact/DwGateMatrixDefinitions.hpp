#ifndef DD_EXACT_DW_GATE_MATRIX_DEFINITIONS_HPP
#define DD_EXACT_DW_GATE_MATRIX_DEFINITIONS_HPP

#include "dd/exact/Dw.hpp"

#include <array>
#include <string>

namespace dd::exact::gates {

/// Exact D[w] matrices for the Clifford+T gate set: entries are {m00, m01,
/// m10, m11}, i.e. row-major 2x2, matching DwPackage's edge convention.
[[nodiscard]] std::array<Dw, 4> i();
[[nodiscard]] std::array<Dw, 4> x();
[[nodiscard]] std::array<Dw, 4> y();
[[nodiscard]] std::array<Dw, 4> z();
[[nodiscard]] std::array<Dw, 4> h();
[[nodiscard]] std::array<Dw, 4> s();
[[nodiscard]] std::array<Dw, 4> sdg();
[[nodiscard]] std::array<Dw, 4> t();
[[nodiscard]] std::array<Dw, 4> tdg();
[[nodiscard]] std::array<Dw, 4> v();
[[nodiscard]] std::array<Dw, 4> vdg();
[[nodiscard]] std::array<Dw, 4> sx();
[[nodiscard]] std::array<Dw, 4> sxdg();

/// Exact D[w] matrices for two-qubit Clifford gates, as row-major 4x4
/// arrays whose row/column index is 2*bit(target0) + bit(target1) --
/// DwPackage::makeTwoQubitGateDD's convention, mirroring MQT Core's
/// TwoQubitGateMatrix definitions (SWAP_MAT/ISWAP_MAT/ISWAPDG_MAT/DCX_MAT).
[[nodiscard]] std::array<Dw, 16> swap();
[[nodiscard]] std::array<Dw, 16> iswap();
[[nodiscard]] std::array<Dw, 16> iswapdg();
[[nodiscard]] std::array<Dw, 16> dcx();

/// Looks up a gate matrix by the same lowercase names QRAT's DSL uses
/// (include/dd/DDSimulation.hpp's DEFINE_SINGLE_TARGET_OPERATION list:
/// i,x,y,z,h,s,sdg,t,tdg,v,vdg,sx,sxdg). All 13 are exactly representable in
/// D[w] -- v/vdg/sx/sxdg via invSqrt2PlusI = (1+i)/2 (k=2, exact); v()/sx()
/// (and vdg()/sxdg()) are literally the same matrix, V being sqrt(X). Any
/// name outside this list throws std::invalid_argument rather than silently
/// falling back to an inexact representation.
[[nodiscard]] std::array<Dw, 4> byName(const std::string &name);

/// Two-qubit analogue of byName(): swap, iswap, iswapdg, dcx (all exactly
/// representable in D[w]); any other name throws std::invalid_argument.
/// cx/cz are deliberately absent -- build those via
/// DwPackage::makeControlledSingleQubitGateDD instead.
[[nodiscard]] std::array<Dw, 16> twoQubitByName(const std::string &name);

} // namespace dd::exact::gates

#endif // DD_EXACT_DW_GATE_MATRIX_DEFINITIONS_HPP

#ifndef QRAT_TEST_DD_EXACT_QFT_HELPER_HPP
#define QRAT_TEST_DD_EXACT_QFT_HELPER_HPP

#include "dd/exact/DwGateMatrixDefinitions.hpp"
#include "dd/exact/DwPackage.hpp"

namespace qft_helper {

using dd::exact::DwPackage;
namespace gates = dd::exact::gates;

/// Standard 3-qubit QFT circuit: H(q2), controlled-S(q1->q2),
/// controlled-T(q0->q2), H(q1), controlled-S(q0->q1), H(q0), SWAP(q0,q2).
/// q0 is the least-significant qubit (matches this codebase's bits[i]
/// convention: value = sum_i bits[i]*2^i). Only R1=Z, R2=S, R3=T are needed
/// for 3 qubits -- all Clifford+T, exactly representable in D[w].
inline DwPackage::mEdge buildQFT3(DwPackage &pkg, std::size_t q0, std::size_t q1, std::size_t q2) {
    DwPackage::mEdge result = pkg.makeIdentity();
    const auto apply = [&](const DwPackage::mEdge &gate) { result = pkg.multiply(gate, result); };
    apply(pkg.makeSingleQubitGateDD(q2, gates::h()));
    apply(pkg.makeControlledSingleQubitGateDD(q1, q2, gates::s()));
    apply(pkg.makeControlledSingleQubitGateDD(q0, q2, gates::t()));
    apply(pkg.makeSingleQubitGateDD(q1, gates::h()));
    apply(pkg.makeControlledSingleQubitGateDD(q0, q1, gates::s()));
    apply(pkg.makeSingleQubitGateDD(q0, gates::h()));
    apply(pkg.makeTwoQubitGateDD(q0, q2, gates::swap()));
    return result;
}

/// Inverse of buildQFT3(): same gates in reverse order, S/T replaced by
/// their adjoints (H and SWAP are self-adjoint).
inline DwPackage::mEdge buildInverseQFT3(DwPackage &pkg, std::size_t q0, std::size_t q1, std::size_t q2) {
    DwPackage::mEdge result = pkg.makeIdentity();
    const auto apply = [&](const DwPackage::mEdge &gate) { result = pkg.multiply(gate, result); };
    apply(pkg.makeTwoQubitGateDD(q0, q2, gates::swap()));
    apply(pkg.makeSingleQubitGateDD(q0, gates::h()));
    apply(pkg.makeControlledSingleQubitGateDD(q0, q1, gates::sdg()));
    apply(pkg.makeSingleQubitGateDD(q1, gates::h()));
    apply(pkg.makeControlledSingleQubitGateDD(q0, q2, gates::tdg()));
    apply(pkg.makeControlledSingleQubitGateDD(q1, q2, gates::sdg()));
    apply(pkg.makeSingleQubitGateDD(q2, gates::h()));
    return result;
}

} // namespace qft_helper

#endif // QRAT_TEST_DD_EXACT_QFT_HELPER_HPP

#include "dd/exact/DwGateMatrixDefinitions.hpp"
#include "dd/exact/DwPackage.hpp"

#include <gtest/gtest.h>

#include <cstddef>
#include <functional>
#include <vector>

using dd::exact::Dw;
using dd::exact::DwPackage;
using dd::exact::NormalizationStrategy;
namespace gates = dd::exact::gates;

namespace {

std::vector<bool> bitsFromIndex(std::size_t index, std::size_t nqubits) {
    std::vector<bool> bits(nqubits);
    for (std::size_t q = 0; q < nqubits; ++q)
        bits[q] = ((index >> q) & 1U) != 0;
    return bits;
}

/// Asserts `op` acts on computational basis states as the bit permutation
/// `perm` (input basis index -> output basis index), with amplitude exactly 1.
void expectPermutation(DwPackage &pkg, const DwPackage::mEdge &op,
                       const std::function<std::size_t(std::size_t)> &perm) {
    const std::size_t dim = std::size_t{1} << pkg.numQubits();
    for (std::size_t in = 0; in < dim; ++in) {
        const auto v = pkg.makeBasisState(bitsFromIndex(in, pkg.numQubits()));
        const auto out = pkg.multiply(op, v);
        for (std::size_t idx = 0; idx < dim; ++idx) {
            const Dw amp = pkg.amplitude(out, bitsFromIndex(idx, pkg.numQubits()));
            if (idx == perm(in))
                EXPECT_EQ(amp, Dw::one()) << "input " << in << ": expected output " << idx;
            else
                EXPECT_TRUE(amp.isZero()) << "input " << in << ": spurious amplitude at " << idx;
        }
    }
}

std::size_t bit(std::size_t index, std::size_t q) { return (index >> q) & 1U; }

std::size_t flipBit(std::size_t index, std::size_t q) { return index ^ (std::size_t{1} << q); }

std::size_t swapBits(std::size_t index, std::size_t q0, std::size_t q1) {
    if (bit(index, q0) != bit(index, q1))
        return flipBit(flipBit(index, q0), q1);
    return index;
}

/// Operator product of a gate sequence in application order (first gate of
/// `seq` is applied first, i.e. the product is seq.back() * ... * seq.front()).
DwPackage::mEdge composeSequence(DwPackage &pkg, const std::vector<DwPackage::mEdge> &seq) {
    DwPackage::mEdge result = pkg.makeIdentity();
    for (const auto &gate : seq)
        result = pkg.multiply(gate, result);
    return result;
}

/// The N&C Fig. 4.9 Toffoli decomposition (same sequence as test_grover.cpp),
/// composed into a single operator DD.
DwPackage::mEdge toffoliViaDecomposition(DwPackage &pkg, std::size_t c1, std::size_t c2, std::size_t target) {
    const auto h = gates::h();
    const auto t = gates::t();
    const auto tdg = gates::tdg();
    const auto x = gates::x();
    return composeSequence(pkg, {
        pkg.makeSingleQubitGateDD(target, h),
        pkg.makeControlledSingleQubitGateDD(c2, target, x),
        pkg.makeSingleQubitGateDD(target, tdg),
        pkg.makeControlledSingleQubitGateDD(c1, target, x),
        pkg.makeSingleQubitGateDD(target, t),
        pkg.makeControlledSingleQubitGateDD(c2, target, x),
        pkg.makeSingleQubitGateDD(target, tdg),
        pkg.makeControlledSingleQubitGateDD(c1, target, x),
        pkg.makeSingleQubitGateDD(c2, t),
        pkg.makeSingleQubitGateDD(target, t),
        pkg.makeSingleQubitGateDD(target, h),
        pkg.makeControlledSingleQubitGateDD(c1, c2, x),
        pkg.makeSingleQubitGateDD(c1, t),
        pkg.makeSingleQubitGateDD(c2, tdg),
        pkg.makeControlledSingleQubitGateDD(c1, c2, x),
    });
}

} // namespace

// The default Inverse normalization is canonical (a fixed reduction +
// normalization rule with hash-consing), so equal operators built through
// different code paths must come out as the exact same edge. The tests
// below lean on that via EXPECT_EQ on mEdges directly. Two kinds of check
// share this file: MultiControlMatchesScalarControlOverload and
// EmptyControlsEqualsPlainSingleQubitGate compare convenience overloads that
// both delegate to the same underlying vector-controls builder today --
// they're delegation-consistency regression checks, not cross-checks against
// an independent algorithm. ToffoliMatchesCliffordTDecomposition,
// SwapMatchesCnotDecomposition, and DcxTargetOrderConvention remain genuine
// differential checks: their "decomposed" side composes single-qubit and
// single-control gates via multiply(), a code path distinct from the one-shot
// multi-control/two-qubit builders under test.

TEST(DwGateBuilders, MultiControlMatchesScalarControlOverload) {
    DwPackage pkg(3);
    const std::vector<std::pair<std::size_t, std::size_t>> placements{{2, 0}, {0, 2}};
    for (const auto &m : {gates::x(), gates::z(), gates::t(), gates::h()}) {
        for (const auto &[control, target] : placements) {
            const auto viaScalar = pkg.makeControlledSingleQubitGateDD(control, target, m);
            const auto viaVector = pkg.makeControlledSingleQubitGateDD(std::vector<std::size_t>{control}, target, m);
            EXPECT_EQ(viaVector, viaScalar) << "control=" << control << " target=" << target;
        }
    }
}

// Delegation-consistency check: makeSingleQubitGateDD forwards to this same
// builder with an empty controls vector, so both sides run identical code.
TEST(DwGateBuilders, EmptyControlsEqualsPlainSingleQubitGate) {
    DwPackage pkg(3);
    for (const auto &m : {gates::x(), gates::h(), gates::t()}) {
        EXPECT_EQ(pkg.makeControlledSingleQubitGateDD(std::vector<std::size_t>{}, 1, m),
                  pkg.makeSingleQubitGateDD(1, m));
    }
}

TEST(DwGateBuilders, ToffoliMatchesDenseSemantics) {
    DwPackage pkg(3);
    const auto ccx = pkg.makeControlledSingleQubitGateDD({1, 2}, 0, gates::x());
    expectPermutation(pkg, ccx, [](std::size_t in) {
        return (bit(in, 1) == 1 && bit(in, 2) == 1) ? flipBit(in, 0) : in;
    });
}

TEST(DwGateBuilders, ToffoliMatchesCliffordTDecomposition) {
    // Both sides are exact, so the native builder and the composed N&C
    // decomposition must produce the identical canonical DD. Cover all
    // control/target placements: target below, between, and above controls.
    DwPackage pkg(3);
    const std::size_t configs[][3] = {{1, 2, 0}, {0, 2, 1}, {0, 1, 2}}; // {c1, c2, target}
    for (const auto &cfg : configs) {
        const auto native = pkg.makeControlledSingleQubitGateDD({cfg[0], cfg[1]}, cfg[2], gates::x());
        const auto decomposed = toffoliViaDecomposition(pkg, cfg[0], cfg[1], cfg[2]);
        EXPECT_EQ(native, decomposed) << "c1=" << cfg[0] << " c2=" << cfg[1] << " target=" << cfg[2];
    }
}

TEST(DwGateBuilders, MultiControlZStraddlingTarget) {
    // 4 qubits, controls {0, 3} straddle target 1: phase -1 exactly when
    // qubits 0, 1, and 3 are all |1>.
    DwPackage pkg(4);
    const auto ccz = pkg.makeControlledSingleQubitGateDD({0, 3}, 1, gates::z());
    for (std::size_t in = 0; in < 16; ++in) {
        const auto v = pkg.makeBasisState(bitsFromIndex(in, 4));
        const auto out = pkg.multiply(ccz, v);
        const Dw expected = (bit(in, 0) == 1 && bit(in, 3) == 1 && bit(in, 1) == 1) ? -Dw::one() : Dw::one();
        EXPECT_EQ(pkg.amplitude(out, bitsFromIndex(in, 4)), expected) << "basis state " << in;
    }
}

TEST(DwGateBuilders, ControlsBelowTarget) {
    DwPackage pkg(4);
    const auto ccx = pkg.makeControlledSingleQubitGateDD({0, 1}, 3, gates::x());
    expectPermutation(pkg, ccx, [](std::size_t in) {
        return (bit(in, 0) == 1 && bit(in, 1) == 1) ? flipBit(in, 3) : in;
    });
}

TEST(DwGateBuilders, SwapMatchesCnotDecomposition) {
    // SWAP(a, b) == CX(a->b) CX(b->a) CX(a->b), adjacent and non-adjacent.
    const std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> cases{{0, 1, 2}, {0, 2, 3}};
    for (const auto &[t0, t1, nqubits] : cases) {
        DwPackage pkg(nqubits);
        const auto native = pkg.makeTwoQubitGateDD(t0, t1, gates::swap());
        const auto x = gates::x();
        const auto decomposed = composeSequence(pkg, {
            pkg.makeControlledSingleQubitGateDD(t0, t1, x),
            pkg.makeControlledSingleQubitGateDD(t1, t0, x),
            pkg.makeControlledSingleQubitGateDD(t0, t1, x),
        });
        EXPECT_EQ(native, decomposed) << "targets " << t0 << ", " << t1;
    }
}

TEST(DwGateBuilders, SwapIsSymmetricInTargetOrder) {
    DwPackage pkg(3);
    EXPECT_EQ(pkg.makeTwoQubitGateDD(0, 2, gates::swap()), pkg.makeTwoQubitGateDD(2, 0, gates::swap()));
}

TEST(DwGateBuilders, DcxTargetOrderConvention) {
    // DCX(t0, t1) == CX(t1->t0) * CX(t0->t1) (matrix MSB = t0); the
    // reversed argument order must give the transposed circuit instead.
    DwPackage pkg(2);
    const auto x = gates::x();
    const std::vector<std::pair<std::size_t, std::size_t>> orders{{0, 1}, {1, 0}};
    for (const auto &[t0, t1] : orders) {
        const auto native = pkg.makeTwoQubitGateDD(t0, t1, gates::dcx());
        const auto decomposed = composeSequence(pkg, {
            pkg.makeControlledSingleQubitGateDD(t0, t1, x),
            pkg.makeControlledSingleQubitGateDD(t1, t0, x),
        });
        EXPECT_EQ(native, decomposed) << "t0=" << t0 << " t1=" << t1;
    }
}

TEST(DwGateBuilders, IswapTimesIswapDagIsIdentity) {
    DwPackage pkg(2);
    const auto a = pkg.makeTwoQubitGateDD(0, 1, gates::iswap());
    const auto b = pkg.makeTwoQubitGateDD(0, 1, gates::iswapdg());
    EXPECT_EQ(pkg.multiply(a, b), pkg.makeIdentity());
    EXPECT_EQ(pkg.multiply(b, a), pkg.makeIdentity());
}

TEST(DwGateBuilders, FredkinAllControlPositions) {
    // Controlled-SWAP with the control above, between, and below the targets.
    const std::size_t configs[][3] = {{2, 0, 1}, {1, 0, 2}, {0, 1, 2}}; // {control, t0, t1}
    for (const auto &cfg : configs) {
        DwPackage pkg(3);
        const auto cswap = pkg.makeControlledTwoQubitGateDD({cfg[0]}, cfg[1], cfg[2], gates::swap());
        expectPermutation(pkg, cswap, [&cfg](std::size_t in) {
            return bit(in, cfg[0]) == 1 ? swapBits(in, cfg[1], cfg[2]) : in;
        });
    }
}

TEST(DwGateBuilders, DoublyControlledSwapStraddlingTargets) {
    // 4 qubits: controls {0, 3} below and above targets {1, 2}.
    DwPackage pkg(4);
    const auto ccswap = pkg.makeControlledTwoQubitGateDD({0, 3}, 1, 2, gates::swap());
    expectPermutation(pkg, ccswap, [](std::size_t in) {
        return (bit(in, 0) == 1 && bit(in, 3) == 1) ? swapBits(in, 1, 2) : in;
    });
}

TEST(DwGateBuilders, NormalizationStrategyInvariance) {
    // The builders only go through makeMEdge, so every strategy must yield
    // the same represented operator (if not the same DD shape).
    DwPackage reference(3, NormalizationStrategy::Inverse);
    const auto refToffoli = reference.makeControlledSingleQubitGateDD({1, 2}, 0, gates::x());
    const auto refCswap = reference.makeControlledTwoQubitGateDD({1}, 0, 2, gates::swap());
    for (const auto strategy :
         {NormalizationStrategy::None, NormalizationStrategy::Inverse, NormalizationStrategy::Gcd}) {
        DwPackage pkg(3, strategy);
        const auto toffoli = pkg.makeControlledSingleQubitGateDD({1, 2}, 0, gates::x());
        const auto cswap = pkg.makeControlledTwoQubitGateDD({1}, 0, 2, gates::swap());
        const std::size_t dim = 8;
        for (std::size_t row = 0; row < dim; ++row) {
            for (std::size_t col = 0; col < dim; ++col) {
                const auto rowBits = bitsFromIndex(row, 3);
                const auto colBits = bitsFromIndex(col, 3);
                EXPECT_EQ(pkg.matrixEntry(toffoli, rowBits, colBits),
                          reference.matrixEntry(refToffoli, rowBits, colBits));
                EXPECT_EQ(pkg.matrixEntry(cswap, rowBits, colBits),
                          reference.matrixEntry(refCswap, rowBits, colBits));
            }
        }
    }
}

TEST(DwGateBuilders, InvalidArgumentsThrow) {
    DwPackage pkg(3);
    const auto m = gates::x();
    const auto m2 = gates::swap();
    EXPECT_THROW((void)pkg.makeControlledSingleQubitGateDD(std::vector<std::size_t>{1, 1}, 0, m),
                 std::invalid_argument);
    EXPECT_THROW((void)pkg.makeControlledSingleQubitGateDD(std::vector<std::size_t>{0}, 0, m),
                 std::invalid_argument);
    EXPECT_THROW((void)pkg.makeControlledSingleQubitGateDD(std::vector<std::size_t>{3}, 0, m),
                 std::invalid_argument);
    EXPECT_THROW((void)pkg.makeControlledSingleQubitGateDD(std::vector<std::size_t>{}, 3, m),
                 std::invalid_argument);
    EXPECT_THROW((void)pkg.makeTwoQubitGateDD(1, 1, m2), std::invalid_argument);
    EXPECT_THROW((void)pkg.makeTwoQubitGateDD(0, 3, m2), std::invalid_argument);
    EXPECT_THROW((void)pkg.makeControlledTwoQubitGateDD({0}, 0, 1, m2), std::invalid_argument);
    EXPECT_THROW((void)pkg.makeControlledTwoQubitGateDD({2, 2}, 0, 1, m2), std::invalid_argument);
}

TEST(DwGateBuilders, SurvivesGarbageCollection) {
    DwPackage pkg(3);
    const auto ccx = pkg.makeControlledSingleQubitGateDD({1, 2}, 0, gates::x());
    pkg.incRef(ccx);
    pkg.garbageCollect(true);
    expectPermutation(pkg, ccx, [](std::size_t in) {
        return (bit(in, 1) == 1 && bit(in, 2) == 1) ? flipBit(in, 0) : in;
    });
    pkg.decRef(ccx);
}

TEST(DwGateBuilders, DcxMatchesDenseSemantics) {
    // DCX(t0=0, t1=1): |01> -> |11>, |11> -> |10>, |10> -> |01>, |00> fixed.
    DwPackage pkg(2);
    const auto dcx = pkg.makeTwoQubitGateDD(0, 1, gates::dcx());
    expectPermutation(pkg, dcx, [](std::size_t in) {
        switch (in) {
            case 0: return std::size_t{0};
            case 1: return std::size_t{2};
            case 2: return std::size_t{3};
            default: return std::size_t{1};
        }
    });
}

TEST(DwGateBuilders, IswapCarriesImaginaryPhaseOnSwappedBasisStates) {
    // ISWAP maps |01> -> i|10> and |10> -> i|01>, with amplitude exactly i
    // (not the bare permutation of SWAP), so this can't reuse expectPermutation.
    DwPackage pkg(2);
    const auto iswap = pkg.makeTwoQubitGateDD(0, 1, gates::iswap());
    const Dw iUnit = Dw::omega() * Dw::omega();

    const auto zeroOne = pkg.makeBasisState({false, true});
    const auto viaZeroOne = pkg.multiply(iswap, zeroOne);
    EXPECT_EQ(pkg.amplitude(viaZeroOne, {true, false}), iUnit);
    EXPECT_TRUE(pkg.amplitude(viaZeroOne, {false, true}).isZero());

    const auto oneZero = pkg.makeBasisState({true, false});
    const auto viaOneZero = pkg.multiply(iswap, oneZero);
    EXPECT_EQ(pkg.amplitude(viaOneZero, {false, true}), iUnit);
    EXPECT_TRUE(pkg.amplitude(viaOneZero, {true, false}).isZero());
}

TEST(DwGateBuilders, ControlledIswapAppliesPhaseOnlyWhenControlSet) {
    // Controlled-ISWAP on targets {1, 2}, control 0: acts as ISWAP on |1>
    // control states and as identity otherwise.
    DwPackage pkg(3);
    const auto ciswap = pkg.makeControlledTwoQubitGateDD({0}, 1, 2, gates::iswap());
    const Dw iUnit = Dw::omega() * Dw::omega();

    // control=0: identity on targets.
    const auto uncontrolled = pkg.makeBasisState({false, true, false});
    const auto viaUncontrolled = pkg.multiply(ciswap, uncontrolled);
    EXPECT_EQ(pkg.amplitude(viaUncontrolled, {false, true, false}), Dw::one());

    // control=1, targets |01>: ISWAP applies, giving i|10> on the targets.
    const auto controlled = pkg.makeBasisState({true, true, false});
    const auto viaControlled = pkg.multiply(ciswap, controlled);
    EXPECT_EQ(pkg.amplitude(viaControlled, {true, false, true}), iUnit);
    EXPECT_TRUE(pkg.amplitude(viaControlled, {true, true, false}).isZero());
}

TEST(DwGateBuilders, TwoQubitByNameLookup) {
    EXPECT_EQ(gates::twoQubitByName("SWAP"), gates::swap());
    EXPECT_EQ(gates::twoQubitByName("iswap"), gates::iswap());
    EXPECT_EQ(gates::twoQubitByName("iswapdg"), gates::iswapdg());
    EXPECT_EQ(gates::twoQubitByName("dcx"), gates::dcx());
    EXPECT_THROW((void)gates::twoQubitByName("cx"), std::invalid_argument);
    EXPECT_THROW((void)gates::twoQubitByName("rxx"), std::invalid_argument);
}

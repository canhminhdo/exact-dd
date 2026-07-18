#include "dd/exact/DwPackage.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <random>

using dd::exact::BasisState;
using dd::exact::Dw;
using dd::exact::DwPackage;

namespace {
const Dw kInvSqrt2(0, 1, 0, -1, 2); // 1/sqrt(2)

std::array<Dw, 4> matI() { return {Dw::one(), Dw::zero(), Dw::zero(), Dw::one()}; }
std::array<Dw, 4> matX() { return {Dw::zero(), Dw::one(), Dw::one(), Dw::zero()}; }
std::array<Dw, 4> matZ() { return {Dw::one(), Dw::zero(), Dw::zero(), -Dw::one()}; }
std::array<Dw, 4> matH() { return {kInvSqrt2, kInvSqrt2, kInvSqrt2, -kInvSqrt2}; }

std::vector<bool> bits(std::initializer_list<bool> b) { return std::vector<bool>(b); }
} // namespace

TEST(DwPackage, ZeroStateAmplitudes) {
    DwPackage pkg(2);
    const auto v = pkg.makeZeroState();
    EXPECT_EQ(pkg.amplitude(v, bits({false, false})), Dw::one());
    EXPECT_TRUE(pkg.amplitude(v, bits({false, true})).isZero());
    EXPECT_TRUE(pkg.amplitude(v, bits({true, false})).isZero());
    EXPECT_TRUE(pkg.amplitude(v, bits({true, true})).isZero());
}

TEST(DwPackage, IdentityLeavesStateUnchanged) {
    DwPackage pkg(2);
    const auto v = pkg.makeBasisState(bits({true, false}));
    const auto id = pkg.makeIdentity();
    const auto result = pkg.multiply(id, v);
    EXPECT_EQ(pkg.amplitude(result, bits({true, false})), Dw::one());
    EXPECT_TRUE(pkg.amplitude(result, bits({false, false})).isZero());
}

TEST(DwPackage, XGateFlipsBit) {
    DwPackage pkg(1);
    const auto v0 = pkg.makeZeroState();
    const auto x = pkg.makeSingleQubitGateDD(0, matX());
    const auto v1 = pkg.multiply(x, v0);
    EXPECT_TRUE(pkg.amplitude(v1, bits({false})).isZero());
    EXPECT_EQ(pkg.amplitude(v1, bits({true})), Dw::one());
}

TEST(DwPackage, HGateProducesEqualSuperposition) {
    DwPackage pkg(1);
    const auto v0 = pkg.makeZeroState();
    const auto h = pkg.makeSingleQubitGateDD(0, matH());
    const auto plus = pkg.multiply(h, v0);
    EXPECT_EQ(pkg.amplitude(plus, bits({false})), kInvSqrt2);
    EXPECT_EQ(pkg.amplitude(plus, bits({true})), kInvSqrt2);
}

TEST(DwPackage, HHIsIdentityOnState) {
    DwPackage pkg(1);
    const auto v0 = pkg.makeZeroState();
    const auto h = pkg.makeSingleQubitGateDD(0, matH());
    const auto hh = pkg.multiply(h, h);
    const auto result = pkg.multiply(hh, v0);
    EXPECT_EQ(pkg.amplitude(result, bits({false})), Dw::one());
    EXPECT_TRUE(pkg.amplitude(result, bits({true})).isZero());
}

TEST(DwPackage, MatrixMultiplyAssociativity) {
    DwPackage pkg(1);
    const auto v0 = pkg.makeZeroState();
    const auto h = pkg.makeSingleQubitGateDD(0, matH());
    const auto z = pkg.makeSingleQubitGateDD(0, matZ());
    const auto viaOperatorProduct = pkg.multiply(pkg.multiply(z, h), v0);
    const auto viaSequentialApplication = pkg.multiply(z, pkg.multiply(h, v0));
    EXPECT_EQ(pkg.amplitude(viaOperatorProduct, bits({false})), pkg.amplitude(viaSequentialApplication, bits({false})));
    EXPECT_EQ(pkg.amplitude(viaOperatorProduct, bits({true})), pkg.amplitude(viaSequentialApplication, bits({true})));
}

namespace {
// bits vectors are indexed by qubit number; build one by explicit qubit ->
// value assignment so tests don't have to reason about control/target order.
std::vector<bool> bitsAt(std::size_t nqubits, std::size_t q0, bool v0, std::size_t q1, bool v1) {
    std::vector<bool> b(nqubits, false);
    b[q0] = v0;
    b[q1] = v1;
    return b;
}
} // namespace

TEST(DwPackage, CnotControlHighTarget) {
    // 2 qubits, control = qubit 1 (higher index), target = qubit 0.
    DwPackage pkg(2);
    const auto cnot = pkg.makeControlledSingleQubitGateDD(1, 0, matX());
    for (bool c : {false, true}) {
        for (bool t : {false, true}) {
            const auto in = pkg.makeBasisState(bitsAt(2, 1, c, 0, t));
            const auto out = pkg.multiply(cnot, in);
            const bool expectedTarget = c ? !t : t;
            EXPECT_EQ(pkg.amplitude(out, bitsAt(2, 1, c, 0, expectedTarget)), Dw::one())
                << "control=" << c << " target(in)=" << t;
        }
    }
}

TEST(DwPackage, CnotControlLowTarget) {
    // 2 qubits, control = qubit 0 (lower index), target = qubit 1 (higher index,
    // hit before control while walking the DD top-down).
    DwPackage pkg(2);
    const auto cnot = pkg.makeControlledSingleQubitGateDD(0, 1, matX());
    for (bool c : {false, true}) {
        for (bool t : {false, true}) {
            const auto in = pkg.makeBasisState(bitsAt(2, 0, c, 1, t));
            const auto out = pkg.multiply(cnot, in);
            const bool expectedTarget = c ? !t : t;
            EXPECT_EQ(pkg.amplitude(out, bitsAt(2, 0, c, 1, expectedTarget)), Dw::one())
                << "control=" << c << " target(in)=" << t;
        }
    }
}

TEST(DwPackage, GhzStateAmplitudes) {
    // H on qubit 2, then CNOT(2->1), CNOT(2->0): (|000> + |111>)/sqrt(2).
    DwPackage pkg(3);
    auto v = pkg.makeZeroState();
    v = pkg.multiply(pkg.makeSingleQubitGateDD(2, matH()), v);
    v = pkg.multiply(pkg.makeControlledSingleQubitGateDD(2, 1, matX()), v);
    v = pkg.multiply(pkg.makeControlledSingleQubitGateDD(2, 0, matX()), v);

    EXPECT_EQ(pkg.amplitude(v, bits({false, false, false})), kInvSqrt2);
    EXPECT_EQ(pkg.amplitude(v, bits({true, true, true})), kInvSqrt2);
    for (auto b : {bits({true, false, false}), bits({false, true, false}), bits({false, false, true}),
                   bits({true, true, false}), bits({true, false, true}), bits({false, true, true})}) {
        EXPECT_TRUE(pkg.amplitude(v, b).isZero());
    }
}

TEST(DwPackage, NodeSharingReusesIdenticalStates) {
    DwPackage pkg(2);
    const auto v1 = pkg.makeBasisState(bits({true, false}));
    const auto countAfterFirst = pkg.vNodeCount();
    const auto v2 = pkg.makeBasisState(bits({true, false}));
    EXPECT_EQ(v1.p, v2.p);
    EXPECT_EQ(pkg.vNodeCount(), countAfterFirst);
}

TEST(DwPackage, AddOfOrthogonalBasisStates) {
    DwPackage pkg(1);
    const auto v0 = pkg.makeBasisState(bits({false}));
    const auto v1 = pkg.makeBasisState(bits({true}));
    const auto sum = pkg.add(v0, v1);
    EXPECT_EQ(pkg.amplitude(sum, bits({false})), Dw::one());
    EXPECT_EQ(pkg.amplitude(sum, bits({true})), Dw::one());
}

TEST(DwPackage, VectorDiagramToStringShowsStoredNodes) {
    DwPackage pkg(2);
    const auto v = pkg.makeBasisState(bits({true, false}));
    const std::string one = Dw::one().toString();
    // Zero-weight children are skipped entirely by appendVectorTree (unlike
    // appendMatrixTree, which prints all 4 children unconditionally) -- no
    // "-> ZERO" line is emitted for either branch here.
    EXPECT_EQ(pkg.vectorDiagramToString(v),
              "root -> @0 * " + one + "\n"
              "  node @0 var=1\n"
              "    0 -> @1 * " + one + "\n"
              "      node @1 var=0\n"
              "        1 -> T(" + one + ")\n");
}

TEST(DwPackage, VectorDiagramToStringShowsReducedTerminalRoot) {
    DwPackage pkg(1);
    const auto plus = pkg.multiply(pkg.makeSingleQubitGateDD(0, matH()), pkg.makeZeroState());
    EXPECT_EQ(pkg.vectorDiagramToString(plus), "root -> T(" + kInvSqrt2.toString() + ")\n");
}

TEST(DwPackage, MatrixDiagramToStringShowsStoredNodes) {
    DwPackage pkg(1);
    const auto x = pkg.makeSingleQubitGateDD(0, matX());
    const std::string one = Dw::one().toString();
    EXPECT_EQ(pkg.matrixDiagramToString(x),
              "root -> @0 * " + one + "\n"
              "  node @0 var=0\n"
              "    00 -> ZERO\n"
              "    01 -> T(" + one + ")\n"
              "    10 -> T(" + one + ")\n"
              "    11 -> ZERO\n");
}

TEST(DwPackage, MatrixDiagramToStringShowsReducedIdentityRoot) {
    DwPackage pkg(2);
    EXPECT_EQ(pkg.matrixDiagramToString(pkg.makeIdentity()), "root -> T(" + Dw::one().toString() + ")\n");
}

TEST(DwPackage, WindowedMakeBasisStateMatchesFullWidthOverload) {
    DwPackage pkg(2);
    const auto viaWindow = pkg.makeBasisState(2, bits({true, false}), 0);
    const auto viaFullWidth = pkg.makeBasisState(bits({true, false}));
    EXPECT_EQ(viaWindow.p, viaFullWidth.p);
    EXPECT_EQ(viaWindow.w, viaFullWidth.w);
}

TEST(DwPackage, WindowedMakeBasisStateBuildsNarrowerThanNumQubits) {
    // n=1 at start=1 within a 3-qubit package: topVar == 1, qubits 0 and 2
    // are left untouched by the resulting edge.
    DwPackage pkg(3);
    const auto y = pkg.makeBasisState(1, bits({true}), 1);
    EXPECT_EQ(pkg.amplitude(y, bitsAt(3, 1, true, 0, false)), Dw::one());
    EXPECT_TRUE(pkg.amplitude(y, bitsAt(3, 1, false, 0, false)).isZero());
    // qubit 0 and qubit 2 don't affect the result (pass-through).
    EXPECT_EQ(pkg.amplitude(y, bitsAt(3, 1, true, 0, true)), Dw::one());
    EXPECT_EQ(pkg.amplitude(y, bitsAt(3, 1, true, 2, true)), Dw::one());
}

TEST(DwPackage, WindowedMakeBasisStateComposesWithKronecker) {
    // y = |1> on qubit 0 only (narrow), x = |1> on qubit 1 (narrow gate-free
    // vEdge isn't directly available, so use a single-qubit gate DD applied
    // to a windowed basis state at qubit 0, then kronecker-shift it).
    DwPackage pkg(2);
    const auto y = pkg.makeBasisState(1, bits({true}), 0);
    const auto xNarrow = pkg.makeBasisState(1, bits({true}), 0); // built at qubit 0, will be shifted
    const auto combined = pkg.kronecker(xNarrow, y, 1);          // shifts xNarrow's qubit 0 -> qubit 1
    EXPECT_EQ(pkg.amplitude(combined, bits({true, true})), Dw::one());
    EXPECT_TRUE(pkg.amplitude(combined, bits({false, true})).isZero());
    EXPECT_TRUE(pkg.amplitude(combined, bits({true, false})).isZero());
    EXPECT_TRUE(pkg.amplitude(combined, bits({false, false})).isZero());
}

TEST(DwPackage, WindowedMakeBasisStateRejectsWindowExceedingNumQubits) {
    DwPackage pkg(2);
    EXPECT_THROW((void)pkg.makeBasisState(2, bits({true, false}), 1), std::invalid_argument);
}

TEST(DwPackage, WindowedMakeBasisStateRejectsMismatchedStateSize) {
    DwPackage pkg(2);
    EXPECT_THROW((void)pkg.makeBasisState(1, bits({true, false}), 0), std::invalid_argument);
}

TEST(DwPackage, BasisStateZeroOneMatchesBoolVectorOverload) {
    DwPackage pkg(1);
    const auto viaBasisState =
        pkg.makeBasisState(1, std::vector<BasisState>{BasisState::One}, 0);
    const auto viaBool = pkg.makeBasisState(1, bits({true}), 0);
    EXPECT_EQ(viaBasisState.p, viaBool.p);
    EXPECT_EQ(viaBasisState.w, viaBool.w);
}

TEST(DwPackage, BasisStatePlusGivesEqualSuperposition) {
    DwPackage pkg(1);
    const auto plus = pkg.makeBasisState(1, std::vector<BasisState>{BasisState::Plus}, 0);
    EXPECT_EQ(pkg.amplitude(plus, bits({false})), kInvSqrt2);
    EXPECT_EQ(pkg.amplitude(plus, bits({true})), kInvSqrt2);
}

TEST(DwPackage, BasisStateMinusNegatesOneAmplitude) {
    DwPackage pkg(1);
    const auto minus = pkg.makeBasisState(1, std::vector<BasisState>{BasisState::Minus}, 0);
    EXPECT_EQ(pkg.amplitude(minus, bits({false})), kInvSqrt2);
    EXPECT_EQ(pkg.amplitude(minus, bits({true})), -kInvSqrt2);
}

TEST(DwPackage, BasisStateRightGivesImaginaryOneAmplitude) {
    DwPackage pkg(1);
    const auto right = pkg.makeBasisState(1, std::vector<BasisState>{BasisState::Right}, 0);
    const Dw iOverSqrt2 = Dw::omega() * Dw::omega() * kInvSqrt2;
    EXPECT_EQ(pkg.amplitude(right, bits({false})), kInvSqrt2);
    EXPECT_EQ(pkg.amplitude(right, bits({true})), iOverSqrt2);
}

TEST(DwPackage, BasisStateLeftGivesNegativeImaginaryOneAmplitude) {
    DwPackage pkg(1);
    const auto left = pkg.makeBasisState(1, std::vector<BasisState>{BasisState::Left}, 0);
    const Dw iOverSqrt2 = Dw::omega() * Dw::omega() * kInvSqrt2;
    EXPECT_EQ(pkg.amplitude(left, bits({false})), kInvSqrt2);
    EXPECT_EQ(pkg.amplitude(left, bits({true})), -iOverSqrt2);
}

TEST(DwPackage, BasisStateTwoChainedPlusQubitsGivesJointOneHalfAmplitude) {
    // Regression test for the risk identified while designing this: a
    // literal port of MQT Core's loop (scaling each child by a bare
    // constant rather than the running edge's own accumulated weight)
    // would give kInvSqrt2 here instead of the correct kInvSqrt2*kInvSqrt2.
    DwPackage pkg(2);
    const auto plusPlus = pkg.makeBasisState(
        2, std::vector<BasisState>{BasisState::Plus, BasisState::Plus}, 0);
    const Dw expected = kInvSqrt2 * kInvSqrt2;
    for (auto b : {bits({false, false}), bits({true, false}), bits({false, true}), bits({true, true})})
        EXPECT_EQ(pkg.amplitude(plusPlus, b), expected);
}

TEST(DwPackage, BasisStateWindowedComposesWithKronecker) {
    DwPackage pkg(2);
    const auto y = pkg.makeBasisState(1, std::vector<BasisState>{BasisState::One}, 0);
    const auto xNarrow = pkg.makeBasisState(1, std::vector<BasisState>{BasisState::One}, 0);
    const auto combined = pkg.kronecker(xNarrow, y, 1);
    EXPECT_EQ(pkg.amplitude(combined, bits({true, true})), Dw::one());
    EXPECT_TRUE(pkg.amplitude(combined, bits({false, true})).isZero());
}

TEST(DwPackage, BasisStateRejectsWindowExceedingNumQubits) {
    DwPackage pkg(2);
    EXPECT_THROW((void)pkg.makeBasisState(
                     2, std::vector<BasisState>{BasisState::Plus, BasisState::One},
                     1),
                 std::invalid_argument);
}

TEST(DwPackage, BasisStateRejectsMismatchedStateSize) {
    DwPackage pkg(2);
    EXPECT_THROW((void)pkg.makeBasisState(2, std::vector<BasisState>{BasisState::Plus}, 0),
                 std::invalid_argument);
}

TEST(DwPackage, RandomSingleQubitStateIsExactlyNormalized) {
    DwPackage pkg(1);
    std::mt19937 rng(42);
    const auto state = pkg.makeRandomSingleQubitState(rng);
    EXPECT_EQ(pkg.innerProduct(state, state), Dw::one());
}

TEST(DwPackage, RandomSingleQubitStateZeroDepthIsZeroState) {
    DwPackage pkg(1);
    std::mt19937 rng(42);
    const auto state = pkg.makeRandomSingleQubitState(rng, 0);
    EXPECT_EQ(pkg.amplitude(state, bits({false})), Dw::one());
    EXPECT_TRUE(pkg.amplitude(state, bits({true})).isZero());
}

TEST(DwPackage, RandomSingleQubitStateVariesAcrossSeeds) {
    DwPackage pkg(1);
    std::vector<Dw> amplitudesAtOne;
    for (unsigned seed = 0; seed < 10; ++seed) {
        std::mt19937 rng(seed);
        const auto state = pkg.makeRandomSingleQubitState(rng);
        EXPECT_EQ(pkg.innerProduct(state, state), Dw::one());
        amplitudesAtOne.push_back(pkg.amplitude(state, bits({true})));
    }
    const bool allSame =
        std::all_of(amplitudesAtOne.begin(), amplitudesAtOne.end(),
                    [&](const Dw &a) { return a == amplitudesAtOne.front(); });
    EXPECT_FALSE(allSame);
}

TEST(DwPackage, RandomSingleQubitStateWindowedStaysNarrow) {
    // qubit=1 within a 3-qubit package: topVar should never exceed 1, and
    // qubits 0/2 don't affect the result (pass-through). innerProduct()
    // isn't meaningful here (it sums over all numQubits() levels, and a
    // narrow, windowed vEdge isn't a normalized state of the FULL
    // register -- pass-through at the untouched qubits would double-count
    // per skipped qubit), so normalization is checked directly via the
    // two amplitudes at the target qubit instead.
    DwPackage pkg(3);
    std::mt19937 rng(7);
    const auto state = pkg.makeRandomSingleQubitState(rng, 10, 1);
    const Dw a = pkg.amplitude(state, bitsAt(3, 1, false, 0, false));
    const Dw b = pkg.amplitude(state, bitsAt(3, 1, true, 0, false));
    EXPECT_EQ(a.normSquared() + b.normSquared(), Dw::one());
    EXPECT_EQ(pkg.amplitude(state, bitsAt(3, 1, false, 0, false)),
              pkg.amplitude(state, bitsAt(3, 1, false, 0, true)));
    EXPECT_EQ(pkg.amplitude(state, bitsAt(3, 1, false, 2, false)),
              pkg.amplitude(state, bitsAt(3, 1, false, 2, true)));
}

TEST(DwPackage, MakeStateFromVectorMatchesZeroState) {
    DwPackage pkg(2);
    const auto viaVector = pkg.makeStateFromVector({Dw::one(), Dw::zero(), Dw::zero(), Dw::zero()});
    const auto viaBasis = pkg.makeZeroState();
    EXPECT_EQ(viaVector.p, viaBasis.p);
    EXPECT_EQ(viaVector.w, viaBasis.w);
}

TEST(DwPackage, MakeStateFromVectorMatchesArbitraryBasisState) {
    DwPackage pkg(2);
    // index 0b10 = qubit1=1, qubit0=0 -> bits({false, true}) (qubit0, qubit1)
    const auto viaVector = pkg.makeStateFromVector({Dw::zero(), Dw::zero(), Dw::one(), Dw::zero()});
    const auto viaBasis = pkg.makeBasisState(bits({false, true}));
    EXPECT_EQ(viaVector.p, viaBasis.p);
    EXPECT_EQ(viaVector.w, viaBasis.w);
}

TEST(DwPackage, MakeStateFromVectorSuperposition) {
    DwPackage pkg(1);
    const auto v = pkg.makeStateFromVector({kInvSqrt2, kInvSqrt2});
    EXPECT_EQ(pkg.amplitude(v, bits({false})), kInvSqrt2);
    EXPECT_EQ(pkg.amplitude(v, bits({true})), kInvSqrt2);
}

TEST(DwPackage, MakeStateFromVectorRejectsWrongSize) {
    DwPackage pkg(2);
    EXPECT_THROW((void)pkg.makeStateFromVector({Dw::one(), Dw::zero()}), std::invalid_argument);
}

TEST(DwPackage, MakeDDFromMatrixMatchesIdentity) {
    DwPackage pkg(2);
    const std::vector<std::vector<Dw>> id = {{Dw::one(), Dw::zero(), Dw::zero(), Dw::zero()},
                                              {Dw::zero(), Dw::one(), Dw::zero(), Dw::zero()},
                                              {Dw::zero(), Dw::zero(), Dw::one(), Dw::zero()},
                                              {Dw::zero(), Dw::zero(), Dw::zero(), Dw::one()}};
    const auto viaMatrix = pkg.makeDDFromMatrix(id);
    const auto viaIdentity = pkg.makeIdentity();
    EXPECT_EQ(viaMatrix.p, viaIdentity.p);
    EXPECT_EQ(viaMatrix.w, viaIdentity.w);
}

TEST(DwPackage, MakeDDFromMatrixMatchesSingleQubitGate) {
    DwPackage pkg(1);
    const std::vector<std::vector<Dw>> x = {{Dw::zero(), Dw::one()}, {Dw::one(), Dw::zero()}};
    const auto viaMatrix = pkg.makeDDFromMatrix(x);
    const auto viaGate = pkg.makeSingleQubitGateDD(0, matX());
    EXPECT_EQ(viaMatrix.p, viaGate.p);
    EXPECT_EQ(viaMatrix.w, viaGate.w);
    for (bool row : {false, true}) {
        for (bool col : {false, true}) {
            EXPECT_EQ(pkg.matrixEntry(viaMatrix, bits({row}), bits({col})),
                      pkg.matrixEntry(viaGate, bits({row}), bits({col})));
        }
    }
}

TEST(DwPackage, MakeDDFromMatrixRejectsWrongSize) {
    DwPackage pkg(1);
    EXPECT_THROW((void)pkg.makeDDFromMatrix({{Dw::one(), Dw::zero(), Dw::zero()}}), std::invalid_argument);
}

TEST(DwPackage, MakeDDFromMatrixRejectsNonSquare) {
    DwPackage pkg(1);
    EXPECT_THROW((void)pkg.makeDDFromMatrix({{Dw::one(), Dw::zero()}, {Dw::zero()}}), std::invalid_argument);
}

TEST(DwPackage, KroneckerMatrixCombinesDisjointSingleQubitGates) {
    // X and Z are each naturally narrow (topVar == 0) since
    // makeSingleQubitGateDD never pads above its target.
    DwPackage pkg(2);
    const auto x = pkg.makeSingleQubitGateDD(0, matX());
    const auto z = pkg.makeSingleQubitGateDD(0, matZ());
    const auto combined = pkg.kronecker(x, z, 1); // shifts x to qubit 1
    const auto reference = pkg.multiply(pkg.makeSingleQubitGateDD(1, matX()), pkg.makeSingleQubitGateDD(0, matZ()));
    for (bool row1 : {false, true}) {
        for (bool row0 : {false, true}) {
            for (bool col1 : {false, true}) {
                for (bool col0 : {false, true}) {
                    EXPECT_EQ(pkg.matrixEntry(combined, bitsAt(2, 1, row1, 0, row0), bitsAt(2, 1, col1, 0, col0)),
                              pkg.matrixEntry(reference, bitsAt(2, 1, row1, 0, row0), bitsAt(2, 1, col1, 0, col0)));
                }
            }
        }
    }
}

TEST(DwPackage, KroneckerVectorGraftsYWhenXIsTrivial) {
    DwPackage pkg(3);
    const auto y = pkg.makeBasisState(bits({true, false, true}));
    const auto x = DwPackage::vEdge::terminal(kInvSqrt2);
    const auto combined = pkg.kronecker(x, y, 3);
    EXPECT_EQ(combined.p, y.p);
    EXPECT_EQ(combined.w, kInvSqrt2 * y.w);
}

TEST(DwPackage, KroneckerVectorScalesXWhenYIsTrivial) {
    DwPackage pkg(2);
    const auto x = pkg.makeBasisState(bits({true, false}));
    const auto trivialY = DwPackage::vEdge::terminal(Dw::one());
    const auto combined = pkg.kronecker(x, trivialY, 0);
    EXPECT_EQ(combined.p, x.p);
    EXPECT_EQ(combined.w, x.w);
    for (auto b : {bits({false, false}), bits({true, false}), bits({false, true}), bits({true, true})})
        EXPECT_EQ(pkg.amplitude(combined, b), pkg.amplitude(x, b));
}

TEST(DwPackage, KroneckerRejectsShiftExceedingNumQubits) {
    DwPackage pkg(2);
    const auto x = pkg.makeSingleQubitGateDD(1, matX()); // topVar == 1
    const auto y = pkg.makeSingleQubitGateDD(0, matZ());
    EXPECT_THROW((void)pkg.kronecker(x, y, 1), std::invalid_argument);
}

namespace {
std::vector<bool> idxBits(std::size_t nqubits, std::size_t index) {
    std::vector<bool> b(nqubits);
    for (std::size_t q = 0; q < nqubits; ++q)
        b[q] = ((index >> q) & 1U) != 0;
    return b;
}
} // namespace

TEST(DwPackage, KroneckerMatrixHandlesSharedDagStructure) {
    // A 2-controlled gate shares several of wrapWithControl's hash-consed
    // nodes across sibling positions (matX()'s repeated 0/1 entries collapse
    // to the same nodes at multiple slots), giving x a genuine DAG -- the
    // same node reachable via more than one path -- rather than a tree.
    // Regression check that memoizing kroneckerRec on x.p still produces
    // the same result an unmemoized walk would.
    DwPackage pkg(4);
    const auto x = pkg.makeControlledSingleQubitGateDD({0, 1}, 2, matX()); // topVar == 2
    const auto y = pkg.makeSingleQubitGateDD(0, matZ());                   // topVar == 0
    const auto combined = pkg.kronecker(x, y, 1); // shifts x's qubits {0,1,2} -> {1,2,3}
    const auto reference =
        pkg.multiply(pkg.makeControlledSingleQubitGateDD({1, 2}, 3, matX()), pkg.makeSingleQubitGateDD(0, matZ()));
    const std::size_t dim = std::size_t{1} << 4;
    for (std::size_t row = 0; row < dim; ++row) {
        for (std::size_t col = 0; col < dim; ++col) {
            EXPECT_EQ(pkg.matrixEntry(combined, idxBits(4, row), idxBits(4, col)),
                      pkg.matrixEntry(reference, idxBits(4, row), idxBits(4, col)));
        }
    }
}

TEST(DwPackage, OuterProductOfZeroStates) {
    DwPackage pkg(1);
    const auto zero = pkg.makeZeroState();
    const auto op = pkg.outerProduct(zero, zero);
    EXPECT_EQ(pkg.matrixEntry(op, bits({false}), bits({false})), Dw::one());
    EXPECT_TRUE(pkg.matrixEntry(op, bits({false}), bits({true})).isZero());
    EXPECT_TRUE(pkg.matrixEntry(op, bits({true}), bits({false})).isZero());
    EXPECT_TRUE(pkg.matrixEntry(op, bits({true}), bits({true})).isZero());
}

TEST(DwPackage, OuterProductOfPlusStates) {
    DwPackage pkg(1);
    const auto plus = pkg.multiply(pkg.makeSingleQubitGateDD(0, matH()), pkg.makeZeroState());
    const auto op = pkg.outerProduct(plus, plus);
    for (bool row : {false, true}) {
        for (bool col : {false, true}) {
            EXPECT_EQ(pkg.matrixEntry(op, bits({row}), bits({col})), kInvSqrt2 * kInvSqrt2);
        }
    }
}

TEST(DwPackage, OuterProductConjugatesYSide) {
    // |0><1|: only M[row=0][col=1] is nonzero.
    DwPackage pkg(1);
    const auto zero = pkg.makeBasisState(bits({false}));
    const auto one = pkg.makeBasisState(bits({true}));
    const auto op = pkg.outerProduct(zero, one);
    EXPECT_TRUE(pkg.matrixEntry(op, bits({false}), bits({false})).isZero());
    EXPECT_EQ(pkg.matrixEntry(op, bits({false}), bits({true})), Dw::one());
    EXPECT_TRUE(pkg.matrixEntry(op, bits({true}), bits({false})).isZero());
    EXPECT_TRUE(pkg.matrixEntry(op, bits({true}), bits({true})).isZero());
}

TEST(DwPackage, OuterProductOfGhzStateWithItselfMatchesAmplitudeFormula) {
    // Same construction as GhzStateAmplitudes: (|000> + |111>)/sqrt2, a
    // genuinely multi-level, non-collapsing state -- a general correctness
    // check on the restructured outerProductRec (not hand-picked like the
    // 1-qubit tests above), and since x and y are the same edge, this
    // naturally exercises the memo hitting repeatedly on matching
    // (x.p, y.p, level) triples.
    DwPackage pkg(3);
    auto psi = pkg.makeZeroState();
    psi = pkg.multiply(pkg.makeSingleQubitGateDD(2, matH()), psi);
    psi = pkg.multiply(pkg.makeControlledSingleQubitGateDD(2, 1, matX()), psi);
    psi = pkg.multiply(pkg.makeControlledSingleQubitGateDD(2, 0, matX()), psi);

    const auto op = pkg.outerProduct(psi, psi);
    const std::size_t dim = std::size_t{1} << 3;
    for (std::size_t row = 0; row < dim; ++row) {
        for (std::size_t col = 0; col < dim; ++col) {
            EXPECT_EQ(pkg.matrixEntry(op, idxBits(3, row), idxBits(3, col)),
                      pkg.amplitude(psi, idxBits(3, row)) * pkg.amplitude(psi, idxBits(3, col)).conjugate());
        }
    }
}

TEST(DwPackage, MeasureOneQubitOnPlusStateGivesEqualProbabilities) {
    DwPackage pkg(1);
    const auto plus = pkg.multiply(pkg.makeSingleQubitGateDD(0, matH()), pkg.makeZeroState());
    const auto measureFalse = pkg.measureOneQubit(plus, 0, false);
    const auto measureTrue = pkg.measureOneQubit(plus, 0, true);
    EXPECT_EQ(measureFalse.probability, kInvSqrt2 * kInvSqrt2);
    EXPECT_EQ(measureTrue.probability, kInvSqrt2 * kInvSqrt2);
    EXPECT_EQ(pkg.amplitude(measureFalse.state, bits({false})), kInvSqrt2);
    EXPECT_TRUE(pkg.amplitude(measureFalse.state, bits({true})).isZero());
    EXPECT_TRUE(pkg.amplitude(measureTrue.state, bits({false})).isZero());
    EXPECT_EQ(pkg.amplitude(measureTrue.state, bits({true})), kInvSqrt2);
}

TEST(DwPackage, FidelityOfStateWithItselfIsOne) {
    DwPackage pkg(1);
    const auto plus = pkg.multiply(pkg.makeSingleQubitGateDD(0, matH()), pkg.makeZeroState());
    EXPECT_DOUBLE_EQ(pkg.fidelity(plus, plus), 1.0);
}

TEST(DwPackage, FidelityOfOrthogonalStatesIsZero) {
    DwPackage pkg(1);
    const auto zero = pkg.makeBasisState(bits({false}));
    const auto one = pkg.makeBasisState(bits({true}));
    EXPECT_DOUBLE_EQ(pkg.fidelity(zero, one), 0.0);
}

TEST(DwPackage, FidelityUnaffectedByUnnormalizedScale) {
    // Projecting |+> onto outcome 0 collapses it to the (unnormalized,
    // scaled-by-kInvSqrt2) |0> ray -- a DIFFERENT ray from |+>, but the
    // SAME ray as a plain |0> basis state, so fidelity against |0> should
    // still be 1.0 despite the scale mismatch (numerator/denominator
    // cancel it).
    DwPackage pkg(1);
    const auto plus = pkg.multiply(pkg.makeSingleQubitGateDD(0, matH()), pkg.makeZeroState());
    const auto unnormalized = pkg.measureOneQubit(plus, 0, false).state; // scaled by kInvSqrt2
    const auto zero = pkg.makeBasisState(bits({false}));
    EXPECT_DOUBLE_EQ(pkg.fidelity(unnormalized, zero), 1.0);
}


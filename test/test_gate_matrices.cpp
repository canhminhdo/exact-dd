#include "dd/exact/DwGateMatrixDefinitions.hpp"
#include "dd/exact/DwPackage.hpp"

#include <gtest/gtest.h>

using dd::exact::Dw;
using dd::exact::DwPackage;
namespace gates = dd::exact::gates;

namespace {
std::vector<bool> bit(bool b) { return {b}; }
} // namespace

TEST(GateMatrices, IdentityIsIdentity) {
    EXPECT_EQ(gates::i()[0], Dw::one());
    EXPECT_EQ(gates::i()[3], Dw::one());
    EXPECT_TRUE(gates::i()[1].isZero());
    EXPECT_TRUE(gates::i()[2].isZero());
}

TEST(GateMatrices, XFlipsZeroToOne) {
    DwPackage pkg(1);
    const auto v1 = pkg.multiply(pkg.makeSingleQubitGateDD(0, gates::x()), pkg.makeZeroState());
    EXPECT_EQ(pkg.amplitude(v1, bit(true)), Dw::one());
    EXPECT_TRUE(pkg.amplitude(v1, bit(false)).isZero());
}

TEST(GateMatrices, YOnZeroGivesIOnOne) {
    DwPackage pkg(1);
    const auto v = pkg.multiply(pkg.makeSingleQubitGateDD(0, gates::y()), pkg.makeZeroState());
    const Dw iUnit = Dw::omega() * Dw::omega();
    EXPECT_EQ(pkg.amplitude(v, bit(true)), iUnit);
    EXPECT_TRUE(pkg.amplitude(v, bit(false)).isZero());
}

TEST(GateMatrices, ZOnOneGivesMinusOne) {
    DwPackage pkg(1);
    const auto one = pkg.multiply(pkg.makeSingleQubitGateDD(0, gates::x()), pkg.makeZeroState());
    const auto v = pkg.multiply(pkg.makeSingleQubitGateDD(0, gates::z()), one);
    EXPECT_EQ(pkg.amplitude(v, bit(true)), -Dw::one());
    EXPECT_TRUE(pkg.amplitude(v, bit(false)).isZero());
}

TEST(GateMatrices, SOnOneGivesI) {
    DwPackage pkg(1);
    const auto one = pkg.multiply(pkg.makeSingleQubitGateDD(0, gates::x()), pkg.makeZeroState());
    const auto v = pkg.multiply(pkg.makeSingleQubitGateDD(0, gates::s()), one);
    const Dw iUnit = Dw::omega() * Dw::omega();
    EXPECT_EQ(pkg.amplitude(v, bit(true)), iUnit);
}

TEST(GateMatrices, TOnOneGivesOmega) {
    DwPackage pkg(1);
    const auto one = pkg.multiply(pkg.makeSingleQubitGateDD(0, gates::x()), pkg.makeZeroState());
    const auto v = pkg.multiply(pkg.makeSingleQubitGateDD(0, gates::t()), one);
    EXPECT_EQ(pkg.amplitude(v, bit(true)), Dw::omega());
}

TEST(GateMatrices, SThenSdgIsIdentity) {
    DwPackage pkg(1);
    const auto s = pkg.makeSingleQubitGateDD(0, gates::s());
    const auto sdg = pkg.makeSingleQubitGateDD(0, gates::sdg());
    const auto one = pkg.multiply(pkg.makeSingleQubitGateDD(0, gates::x()), pkg.makeZeroState());
    const auto v = pkg.multiply(sdg, pkg.multiply(s, one));
    EXPECT_EQ(pkg.amplitude(v, bit(true)), Dw::one());
    EXPECT_TRUE(pkg.amplitude(v, bit(false)).isZero());
}

TEST(GateMatrices, TThenTdgIsIdentity) {
    DwPackage pkg(1);
    const auto t = pkg.makeSingleQubitGateDD(0, gates::t());
    const auto tdg = pkg.makeSingleQubitGateDD(0, gates::tdg());
    const auto one = pkg.multiply(pkg.makeSingleQubitGateDD(0, gates::x()), pkg.makeZeroState());
    const auto v = pkg.multiply(tdg, pkg.multiply(t, one));
    EXPECT_EQ(pkg.amplitude(v, bit(true)), Dw::one());
}

TEST(GateMatrices, TTIsS) {
    DwPackage pkg(1);
    const auto tt = pkg.multiply(pkg.makeSingleQubitGateDD(0, gates::t()), pkg.makeSingleQubitGateDD(0, gates::t()));
    const auto sGate = pkg.makeSingleQubitGateDD(0, gates::s());
    const auto one = pkg.multiply(pkg.makeSingleQubitGateDD(0, gates::x()), pkg.makeZeroState());
    const auto viaTT = pkg.multiply(tt, one);
    const auto viaS = pkg.multiply(sGate, one);
    EXPECT_EQ(pkg.amplitude(viaTT, bit(true)), pkg.amplitude(viaS, bit(true)));
}

TEST(GateMatrices, ByNameResolvesCliffordT) {
    for (const std::string &name :
         {"i", "x", "y", "z", "h", "s", "sdg", "t", "tdg", "v", "vdg", "sx", "sxdg"}) {
        EXPECT_NO_THROW((void)gates::byName(name)) << name;
    }
}

TEST(GateMatrices, ByNameRejectsUnsupportedGates) {
    for (const std::string &name : {"rx", "not_a_gate"}) {
        EXPECT_THROW((void)gates::byName(name), std::invalid_argument) << name;
    }
}

TEST(GateMatrices, VThenVdgIsIdentity) {
    DwPackage pkg(1);
    const auto v = pkg.makeSingleQubitGateDD(0, gates::v());
    const auto vdg = pkg.makeSingleQubitGateDD(0, gates::vdg());
    const auto one = pkg.multiply(pkg.makeSingleQubitGateDD(0, gates::x()), pkg.makeZeroState());
    const auto result = pkg.multiply(vdg, pkg.multiply(v, one));
    EXPECT_EQ(pkg.amplitude(result, bit(true)), Dw::one());
    EXPECT_TRUE(pkg.amplitude(result, bit(false)).isZero());
}

TEST(GateMatrices, VSquaredIsMinusIX) {
    // MQT Core's V_MAT satisfies V^2 == -i*X (not X -- that's SX's identity).
    DwPackage pkg(1);
    const auto vv = pkg.multiply(pkg.makeSingleQubitGateDD(0, gates::v()), pkg.makeSingleQubitGateDD(0, gates::v()));
    const auto xGate = pkg.makeSingleQubitGateDD(0, gates::x());
    const auto zero = pkg.makeZeroState();
    const Dw iUnit = Dw::omega() * Dw::omega();
    const auto viaVV = pkg.multiply(vv, zero);
    const auto viaX = pkg.multiply(xGate, zero);
    EXPECT_EQ(pkg.amplitude(viaVV, bit(true)), -iUnit * pkg.amplitude(viaX, bit(true)));
    EXPECT_EQ(pkg.amplitude(viaVV, bit(false)), -iUnit * pkg.amplitude(viaX, bit(false)));
}

TEST(GateMatrices, SxSquaredIsX) {
    // SX_MAT satisfies SX^2 == X exactly (unlike V_MAT, which picks up a
    // global phase of w == e^{i*pi/4} relative to SX_MAT).
    DwPackage pkg(1);
    const auto ss = pkg.multiply(pkg.makeSingleQubitGateDD(0, gates::sx()), pkg.makeSingleQubitGateDD(0, gates::sx()));
    const auto xGate = pkg.makeSingleQubitGateDD(0, gates::x());
    const auto zero = pkg.makeZeroState();
    const auto viaSS = pkg.multiply(ss, zero);
    const auto viaX = pkg.multiply(xGate, zero);
    EXPECT_EQ(pkg.amplitude(viaSS, bit(true)), pkg.amplitude(viaX, bit(true)));
    EXPECT_EQ(pkg.amplitude(viaSS, bit(false)), pkg.amplitude(viaX, bit(false)));
}

TEST(GateMatrices, VDiffersFromSxByGlobalPhaseW) {
    // V_MAT == conjugate(w) * SX_MAT elementwise (equivalently SX_MAT == w *
    // V_MAT), matching MQT Core's dd::Operations.hpp dispatch of qc::V ->
    // V_MAT vs qc::SX -> SX_MAT as genuinely distinct matrices.
    const Dw omega = Dw::omega();
    const auto vMat = gates::v();
    const auto sxMat = gates::sx();
    for (std::size_t idx = 0; idx < 4; ++idx)
        EXPECT_EQ(vMat[idx], omega.conjugate() * sxMat[idx]) << "index " << idx;
}

TEST(GateMatrices, SwapMatrixIsPermutation) {
    // Row/col index is 2*bit(target0) + bit(target1); SWAP maps |01> <-> |10>
    // and leaves |00>, |11> fixed, with amplitude exactly 1 throughout.
    const auto m = gates::swap();
    const std::array<Dw, 16> expected{
        Dw::one(),  Dw::zero(), Dw::zero(), Dw::zero(),
        Dw::zero(), Dw::zero(), Dw::one(),  Dw::zero(),
        Dw::zero(), Dw::one(),  Dw::zero(), Dw::zero(),
        Dw::zero(), Dw::zero(), Dw::zero(), Dw::one(),
    };
    for (std::size_t idx = 0; idx < 16; ++idx)
        EXPECT_EQ(m[idx], expected[idx]) << "index " << idx;
}

TEST(GateMatrices, IswapMatrixHasImaginaryOffDiagonal) {
    // ISWAP maps |01> -> i|10> and |10> -> i|01>, fixes |00>/|11>.
    const Dw iUnit = Dw::omega() * Dw::omega();
    const auto m = gates::iswap();
    const std::array<Dw, 16> expected{
        Dw::one(),  Dw::zero(), Dw::zero(), Dw::zero(),
        Dw::zero(), Dw::zero(), iUnit,      Dw::zero(),
        Dw::zero(), iUnit,      Dw::zero(), Dw::zero(),
        Dw::zero(), Dw::zero(), Dw::zero(), Dw::one(),
    };
    for (std::size_t idx = 0; idx < 16; ++idx)
        EXPECT_EQ(m[idx], expected[idx]) << "index " << idx;
}

TEST(GateMatrices, IswapdgIsConjugateTransposeOfIswap) {
    // ISWAPDG's off-diagonal entries are -i, the conjugate of ISWAP's +i.
    const auto iswap = gates::iswap();
    const auto iswapdg = gates::iswapdg();
    for (std::size_t idx = 0; idx < 16; ++idx)
        EXPECT_EQ(iswapdg[idx], iswap[idx].conjugate()) << "index " << idx;
}

TEST(GateMatrices, DcxMatrixIsPermutation) {
    // DCX maps |01> -> |11>, |11> -> |10>, |10> -> |01>, fixes |00>.
    const auto m = gates::dcx();
    const std::array<Dw, 16> expected{
        Dw::one(),  Dw::zero(), Dw::zero(), Dw::zero(),
        Dw::zero(), Dw::zero(), Dw::one(),  Dw::zero(),
        Dw::zero(), Dw::zero(), Dw::zero(), Dw::one(),
        Dw::zero(), Dw::one(),  Dw::zero(), Dw::zero(),
    };
    for (std::size_t idx = 0; idx < 16; ++idx)
        EXPECT_EQ(m[idx], expected[idx]) << "index " << idx;
}

TEST(GateMatrices, TwoQubitByNameResolvesSupportedGates) {
    for (const std::string &name : {"swap", "iswap", "iswapdg", "dcx", "SWAP", "ISwap"})
        EXPECT_NO_THROW((void)gates::twoQubitByName(name)) << name;
}

TEST(GateMatrices, TwoQubitByNameRejectsUnsupportedGates) {
    for (const std::string &name : {"cx", "cz", "not_a_gate"})
        EXPECT_THROW((void)gates::twoQubitByName(name), std::invalid_argument) << name;
}

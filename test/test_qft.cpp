#include "qft_helper.hpp"

#include <gtest/gtest.h>

#include <cstddef>
#include <vector>

using dd::exact::Dw;
using dd::exact::DwPackage;
namespace gates = dd::exact::gates;

namespace {
std::vector<bool> bitsFromIndex(std::size_t index, std::size_t nqubits) {
    std::vector<bool> bits(nqubits);
    for (std::size_t q = 0; q < nqubits; ++q)
        bits[q] = ((index >> q) & 1U) != 0;
    return bits;
}
} // namespace

// The 3-qubit QFT's matrix kernel is exactly the 8-point DFT: QFT[k][j] =
// (1/sqrt(8)) * omega^(j*k mod 8), where omega = e^{2*pi*i/8} = e^{i*pi/4}
// is precisely Dw::omega() -- not an approximation, an exact algebraic
// identity available only because n=3 keeps every phase at pi/4 granularity.
TEST(QFT, MatchesDFTOnBasisStates) {
    DwPackage pkg(3);
    const auto qft3 = qft_helper::buildQFT3(pkg, 0, 1, 2);

    const Dw invSqrt2(0, 1, 0, -1, 2);
    const Dw invSqrt8 = invSqrt2 * invSqrt2 * invSqrt2;
    const Dw omega = Dw::omega();

    for (std::size_t j = 0; j < 8; ++j) {
        for (std::size_t k = 0; k < 8; ++k) {
            Dw omegaPow = Dw::one();
            for (std::size_t p = 0; p < (j * k) % 8; ++p)
                omegaPow = omegaPow * omega;
            const Dw expected = invSqrt8 * omegaPow;
            const auto rowBits = bitsFromIndex(k, 3); // row = k (output index)
            const auto colBits = bitsFromIndex(j, 3); // col = j (input index)
            EXPECT_EQ(pkg.matrixEntry(qft3, rowBits, colBits), expected) << "j=" << j << " k=" << k;
        }
    }
}

TEST(QFT, ThenInverseQFTIsIdentity) {
    DwPackage pkg(3);
    const auto qft3 = qft_helper::buildQFT3(pkg, 0, 1, 2);
    const auto inv3 = qft_helper::buildInverseQFT3(pkg, 0, 1, 2);
    EXPECT_EQ(pkg.multiply(inv3, qft3), pkg.makeIdentity());
    EXPECT_EQ(pkg.multiply(qft3, inv3), pkg.makeIdentity());
}

// A period-2 input (equal superposition over indices {0,2,4,6}, i.e. q0
// fixed at |0> with q1,q2 in superposition) must transform under QFT to a
// state supported only on frequencies that are multiples of 8/2=4 -- indices
// 0 and 4 -- each with amplitude exactly 1/sqrt(2), zero everywhere else.
// This is the exact period-finding fact Shor's algorithm relies on.
TEST(QFT, DetectsPeriodicity) {
    DwPackage pkg(3);
    auto state = pkg.makeZeroState();
    state = pkg.multiply(pkg.makeSingleQubitGateDD(1, gates::h()), state);
    state = pkg.multiply(pkg.makeSingleQubitGateDD(2, gates::h()), state);

    const auto qft3 = qft_helper::buildQFT3(pkg, 0, 1, 2);
    const auto result = pkg.multiply(qft3, state);

    const Dw invSqrt2(0, 1, 0, -1, 2);
    for (std::size_t k = 0; k < 8; ++k) {
        const Dw amp = pkg.amplitude(result, bitsFromIndex(k, 3));
        if (k == 0 || k == 4)
            EXPECT_EQ(amp, invSqrt2) << "k=" << k;
        else
            EXPECT_TRUE(amp.isZero()) << "k=" << k;
    }
}

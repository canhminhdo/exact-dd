#include "qft_helper.hpp"

#include <gtest/gtest.h>

#include <cstddef>
#include <vector>

using dd::exact::Dw;
using dd::exact::DwPackage;
namespace gates = dd::exact::gates;

// Shor's algorithm, factoring N=15 via order-finding with a=7. Chosen
// because the order of 7 mod 15 is r=4, which divides 2^3=8 exactly -- a
// 3-qubit counting register gives an exact, zero-error phase estimate using
// only R1=Z, R2=S, R3=T (all Clifford+T), and the modular exponentiation
// register is a pure permutation circuit (SWAP/X), also exactly
// representable. This is the standard "hello world" Shor demo.
namespace {
std::vector<bool> bitsFromIndex(std::size_t index, std::size_t nqubits) {
    std::vector<bool> bits(nqubits);
    for (std::size_t q = 0; q < nqubits; ++q)
        bits[q] = ((index >> q) & 1U) != 0;
    return bits;
}

// "Multiply by 7 mod 15" on a 4-qubit register (w0 = LSB .. w3 = MSB):
// cyclically rotate the 4 qubits' values down by one position (via 3
// adjacent SWAPs), then flip all 4 bits. Verified by hand on the reachable
// orbit {1,7,4,13}: 1=(1,0,0,0) -rotate-> (0,0,0,1) -flip-> (1,1,1,0)=7;
// 7=(1,1,1,0) -rotate-> (1,1,0,1) -flip-> (0,0,1,0)=4; 4=(0,0,1,0) -rotate->
// (0,1,0,0) -flip-> (1,0,1,1)=13; 13=(1,0,1,1) -rotate-> (0,1,1,1) -flip->
// (1,0,0,0)=1.
DwPackage::mEdge buildMultiplyBy7(DwPackage &pkg, std::size_t w0, std::size_t w1, std::size_t w2, std::size_t w3) {
    DwPackage::mEdge result = pkg.makeIdentity();
    const auto apply = [&](const DwPackage::mEdge &gate) { result = pkg.multiply(gate, result); };
    apply(pkg.makeTwoQubitGateDD(w0, w1, gates::swap()));
    apply(pkg.makeTwoQubitGateDD(w1, w2, gates::swap()));
    apply(pkg.makeTwoQubitGateDD(w2, w3, gates::swap()));
    apply(pkg.makeSingleQubitGateDD(w0, gates::x()));
    apply(pkg.makeSingleQubitGateDD(w1, gates::x()));
    apply(pkg.makeSingleQubitGateDD(w2, gates::x()));
    apply(pkg.makeSingleQubitGateDD(w3, gates::x()));
    return result;
}

// Same circuit with every elementary gate individually controlled by the
// same qubit: when control=0 every factor is identity (whole product is
// identity); when control=1 every factor becomes its active gate (whole
// product is buildMultiplyBy7's circuit). Standard technique for building a
// controlled multi-gate circuit from single-controlled primitives.
DwPackage::mEdge buildControlledMultiplyBy7(DwPackage &pkg, std::size_t control, std::size_t w0, std::size_t w1,
                                            std::size_t w2, std::size_t w3) {
    DwPackage::mEdge result = pkg.makeIdentity();
    const auto apply = [&](const DwPackage::mEdge &gate) { result = pkg.multiply(gate, result); };
    apply(pkg.makeControlledTwoQubitGateDD(control, w0, w1, gates::swap()));
    apply(pkg.makeControlledTwoQubitGateDD(control, w1, w2, gates::swap()));
    apply(pkg.makeControlledTwoQubitGateDD(control, w2, w3, gates::swap()));
    apply(pkg.makeControlledSingleQubitGateDD(control, w0, gates::x()));
    apply(pkg.makeControlledSingleQubitGateDD(control, w1, gates::x()));
    apply(pkg.makeControlledSingleQubitGateDD(control, w2, gates::x()));
    apply(pkg.makeControlledSingleQubitGateDD(control, w3, gates::x()));
    return result;
}

int modpow(int base, int exp, int mod) {
    int result = 1;
    for (int i = 0; i < exp; ++i)
        result = (result * base) % mod;
    return result;
}

int gcdOf(int a, int b) {
    while (b != 0) {
        const int t = b;
        b = a % b;
        a = t;
    }
    return a;
}
} // namespace

TEST(Shor, ControlledModExp7Mod15FollowsExpectedOrbit) {
    DwPackage pkg(4);
    const auto mult7 = buildMultiplyBy7(pkg, 0, 1, 2, 3);
    auto state = pkg.makeBasisState({true, false, false, false}); // value 1

    const std::vector<std::size_t> expectedOrbit{7, 4, 13, 1};
    for (const std::size_t expected : expectedOrbit) {
        state = pkg.multiply(mult7, state);
        for (std::size_t v = 0; v < 16; ++v) {
            const Dw amp = pkg.amplitude(state, bitsFromIndex(v, 4));
            if (v == expected)
                EXPECT_EQ(amp, Dw::one()) << "expected value " << expected;
            else
                EXPECT_TRUE(amp.isZero()) << "unexpected amplitude at " << v << " (expected " << expected << ")";
        }
    }
}

TEST(Shor, ShorFactors15UsingOrderFinding) {
    // 7 qubits: 0,1,2 = counting register; 3,4,5,6 = work register (values 0..15).
    DwPackage pkg(7);
    DwPackage::mEdge circuit = pkg.makeIdentity();
    const auto apply = [&](const DwPackage::mEdge &gate) { circuit = pkg.multiply(gate, circuit); };

    apply(pkg.makeSingleQubitGateDD(3, gates::x())); // init work register to |0001> = 1
    apply(pkg.makeSingleQubitGateDD(0, gates::h()));
    apply(pkg.makeSingleQubitGateDD(1, gates::h()));
    apply(pkg.makeSingleQubitGateDD(2, gates::h()));

    apply(buildControlledMultiplyBy7(pkg, 0, 3, 4, 5, 6)); // controlled-U_7^1
    apply(buildControlledMultiplyBy7(pkg, 1, 3, 4, 5, 6)); // controlled-U_7^2 (=U_4): apply twice
    apply(buildControlledMultiplyBy7(pkg, 1, 3, 4, 5, 6));
    // controlled-U_7^4 (=U_1=identity) from qubit 2: order is 4, so this is
    // a no-op and is skipped entirely.

    apply(qft_helper::buildInverseQFT3(pkg, 0, 1, 2));

    const auto finalState = pkg.multiply(circuit, pkg.makeZeroState());

    // Exact, zero-error distribution: marginal probability of the counting
    // register (summed over the, possibly entangled, work register) is
    // exactly 1/4 at each multiple of 2 (0,2,4,6) and exactly 0 at odd
    // outcomes, since r=4 divides 2^3=8 evenly.
    const Dw half(1, 0, 0, 0, 2); // 1/2
    const Dw quarter = half * half; // 1/4, computed via Dw arithmetic
    for (std::size_t counting = 0; counting < 8; ++counting) {
        Dw totalProb = Dw::zero();
        for (std::size_t work = 0; work < 16; ++work) {
            std::vector<bool> bits(7);
            bits[0] = (counting & 1U) != 0;
            bits[1] = (counting & 2U) != 0;
            bits[2] = (counting & 4U) != 0;
            bits[3] = (work & 1U) != 0;
            bits[4] = (work & 2U) != 0;
            bits[5] = (work & 4U) != 0;
            bits[6] = (work & 8U) != 0;
            totalProb = totalProb + pkg.amplitude(finalState, bits).normSquared();
        }
        if (counting % 2 == 0)
            EXPECT_EQ(totalProb, quarter) << "counting=" << counting;
        else
            EXPECT_TRUE(totalProb.isZero()) << "counting=" << counting;
    }

    // Classical post-processing (no DD involved): from the nonzero outcome
    // counting=2 (phase 2/8=1/4, continued-fraction denominator r=4), recover
    // the two nontrivial factors of 15.
    constexpr int a = 7;
    constexpr int n = 15;
    constexpr int r = 4;
    const int factor1 = gcdOf(modpow(a, r / 2, n) - 1, n);
    const int factor2 = gcdOf(modpow(a, r / 2, n) + 1, n);
    EXPECT_EQ(factor1, 3);
    EXPECT_EQ(factor2, 5);
    EXPECT_EQ(factor1 * factor2, n);
}

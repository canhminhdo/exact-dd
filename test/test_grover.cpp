#include "dd/exact/ExactDDSimulation.hpp"
#include "dd/exact/DwPackage.hpp"

#include <gtest/gtest.h>

#include <vector>

using dd::exact::Dw;
using dd::exact::ExactDDSimulation;

// End-to-end Grover's search over a 5-qubit register (qubits 0..4), marking
// the single basis state |01010> (qubit i holds the i-th character of
// "01010", i.e. q0=0, q1=1, q2=0, q3=1, q4=0).
//
// The first test deliberately avoids DwPackage's native multi-controlled
// gates: the oracle's and diffuser's 4-controlled-Z are built from Toffoli
// gates (themselves decomposed into H/T/T^dg/CNOT, Nielsen & Chuang Fig.
// 4.9) chained through 3 ancilla qubits (5,6,7) that compute the AND of the
// four control qubits and are always uncomputed back to |0> before the next
// gate. Because every primitive involved (H, X, T, T^dg, CNOT, Toffoli) is
// exactly representable in D[w], the whole circuit -- and therefore every
// amplitude checked below -- is exact, not an approximation of the usual
// sin^2 Grover formula. The second test runs the same search with a native
// 4-controlled Z (applyMultiControlledGate) on a plain 5-qubit register --
// no ancillas -- and must reach the exact same probabilities.
namespace {

void toffoli(ExactDDSimulation &sim, std::size_t c1, std::size_t c2, std::size_t target) {
    sim.applyGate("h", target);
    sim.applyControlledGate("x", c2, target);
    sim.applyGate("tdg", target);
    sim.applyControlledGate("x", c1, target);
    sim.applyGate("t", target);
    sim.applyControlledGate("x", c2, target);
    sim.applyGate("tdg", target);
    sim.applyControlledGate("x", c1, target);
    sim.applyGate("t", c2);
    sim.applyGate("t", target);
    sim.applyGate("h", target);
    sim.applyControlledGate("x", c1, c2);
    sim.applyGate("t", c1);
    sim.applyGate("tdg", c2);
    sim.applyControlledGate("x", c1, c2);
}

// Flips the phase of the all-ones state |11..1> on `controls` + `target`
// (a (controls.size()+1)-qubit multi-controlled Z), using
// controls.size()-1 ancillas to ladder the AND of the control bits, then
// uncomputing the ladder so the ancillas are restored to |0>.
void multiControlledZ(ExactDDSimulation &sim, const std::vector<std::size_t> &controls, std::size_t target,
                       const std::vector<std::size_t> &ancillas) {
    ASSERT_EQ(ancillas.size(), controls.size() - 1);
    toffoli(sim, controls[0], controls[1], ancillas[0]);
    for (std::size_t i = 2; i < controls.size(); ++i) {
        toffoli(sim, ancillas[i - 2], controls[i], ancillas[i - 1]);
    }
    sim.applyControlledGate("z", ancillas.back(), target);
    for (std::size_t i = controls.size() - 1; i > 1; --i) {
        toffoli(sim, ancillas[i - 2], controls[i], ancillas[i - 1]);
    }
    toffoli(sim, controls[0], controls[1], ancillas[0]);
}

// Phase-flips |mark> (a size-5 boolean pattern over qubits 0..4) by
// sandwiching multiControlledZ() with X gates on the qubits whose target
// bit is 0, mapping |mark> to |11111> and back.
void oracle(ExactDDSimulation &sim, const std::vector<bool> &mark, const std::vector<std::size_t> &ancillas) {
    for (std::size_t q = 0; q < mark.size(); ++q) {
        if (!mark[q])
            sim.applyGate("x", q);
    }
    multiControlledZ(sim, {0, 1, 2, 3}, 4, ancillas);
    for (std::size_t q = 0; q < mark.size(); ++q) {
        if (!mark[q])
            sim.applyGate("x", q);
    }
}

// D = H^5 (2|0><0| - I) H^5, up to the unobservable global phase -1 (built
// here as H^5 X^5 (C^4 Z) X^5 H^5, which realizes I - 2|0><0| instead).
void diffuser(ExactDDSimulation &sim, const std::vector<std::size_t> &ancillas) {
    for (std::size_t q = 0; q < 5; ++q)
        sim.applyGate("h", q);
    for (std::size_t q = 0; q < 5; ++q)
        sim.applyGate("x", q);
    multiControlledZ(sim, {0, 1, 2, 3}, 4, ancillas);
    for (std::size_t q = 0; q < 5; ++q)
        sim.applyGate("x", q);
    for (std::size_t q = 0; q < 5; ++q)
        sim.applyGate("h", q);
}

std::vector<bool> registerState(const std::vector<bool> &reg) {
    std::vector<bool> full = reg;
    full.resize(8, false); // ancillas (qubits 5,6,7) always uncomputed back to |0>
    return full;
}

} // namespace

TEST(Grover, FiveQubitSearchAmplifiesMarkedState) {
    const std::vector<bool> mark{false, true, false, true, false}; // |01010>
    const std::vector<std::size_t> ancillas{5, 6, 7};

    ExactDDSimulation sim(8);
    for (std::size_t q = 0; q < 5; ++q)
        sim.applyGate("h", q);

    // Optimal iteration count for N=32, M=1: round(pi/4*sqrt(32) - 1/2) == 4.
    const int iterations = 4;
    for (int i = 0; i < iterations; ++i) {
        oracle(sim, mark, ancillas);
        diffuser(sim, ancillas);
    }

    EXPECT_NEAR(sim.normSquared().toComplexDouble().real(), 1.0, 1e-12);

    const Dw markedAmp = sim.amplitude(registerState(mark));
    const double markedProb = markedAmp.normSquared().toComplexDouble().real();
    // Exact Grover formula: sin^2(9*theta), theta = asin(1/sqrt(32)).
    EXPECT_NEAR(markedProb, 0.9991823155432941, 1e-12);

    // Every ancilla must be uncomputed back to |0>, and the marked state
    // must dominate all other 31 basis states of the 5-qubit register.
    double otherProbSum = 0.0;
    for (std::size_t idx = 0; idx < 32; ++idx) {
        std::vector<bool> reg(5);
        for (std::size_t q = 0; q < 5; ++q)
            reg[q] = (idx >> q) & 1U;
        if (reg == mark)
            continue;
        const double p = sim.amplitude(registerState(reg)).normSquared().toComplexDouble().real();
        EXPECT_LT(p, 0.001) << "unmarked basis state " << idx << " has unexpectedly high probability";
        otherProbSum += p;
    }
    EXPECT_NEAR(markedProb + otherProbSum, 1.0, 1e-9);
}

TEST(Grover, NativeMultiControlledZMatchesDecomposedSearch) {
    const std::vector<bool> mark{false, true, false, true, false}; // |01010>

    ExactDDSimulation sim(5); // no ancillas needed
    for (std::size_t q = 0; q < 5; ++q)
        sim.applyGate("h", q);

    const int iterations = 4;
    for (int i = 0; i < iterations; ++i) {
        // Oracle: phase-flip |mark> via X-conjugated native 4-controlled Z.
        for (std::size_t q = 0; q < 5; ++q) {
            if (!mark[q])
                sim.applyGate("x", q);
        }
        sim.applyMultiControlledGate("z", {0, 1, 2, 3}, 4);
        for (std::size_t q = 0; q < 5; ++q) {
            if (!mark[q])
                sim.applyGate("x", q);
        }
        // Diffuser: H^5 X^5 (C^4 Z) X^5 H^5.
        for (std::size_t q = 0; q < 5; ++q)
            sim.applyGate("h", q);
        for (std::size_t q = 0; q < 5; ++q)
            sim.applyGate("x", q);
        sim.applyMultiControlledGate("z", {0, 1, 2, 3}, 4);
        for (std::size_t q = 0; q < 5; ++q)
            sim.applyGate("x", q);
        for (std::size_t q = 0; q < 5; ++q)
            sim.applyGate("h", q);
    }

    EXPECT_NEAR(sim.normSquared().toComplexDouble().real(), 1.0, 1e-12);
    const Dw markedAmp = sim.amplitude(mark);
    // Same exact probability as the ancilla-based decomposition above:
    // sin^2(9*theta), theta = asin(1/sqrt(32)).
    EXPECT_NEAR(markedAmp.normSquared().toComplexDouble().real(), 0.9991823155432941, 1e-12);
}

// The existing strategy coverage (DwGateBuilders.NormalizationStrategyInvariance)
// exercises the matrix gate builders on 3 qubits and nothing else -- it never
// reaches makeVEdge, multiply or add. Those are precisely the paths a real
// algorithm spends its time in, so a full Grover run under a non-default
// strategy would otherwise be running on code no test has driven.
//
// Every strategy must produce the same represented state (the DD shape may
// differ; the amplitudes may not). Amplitudes are compared exactly -- these
// are D[w] values, so there is no tolerance to allow for.
TEST(Grover, AllNormalizationStrategiesAgreeOnTheFullSearch) {
    constexpr std::size_t n = 5;
    const std::vector<std::size_t> controls{0, 1, 2, 3};

    const auto run = [&](dd::exact::NormalizationStrategy strategy) {
        ExactDDSimulation sim(n, strategy);
        for (std::size_t q = 0; q < n; ++q)
            sim.applyGate("h", q);
        // 4 iterations is the optimal count for n = 5.
        for (int it = 0; it < 4; ++it) {
            sim.applyMultiControlledGate("z", controls, n - 1); // oracle: mark |1...1>
            for (std::size_t q = 0; q < n; ++q)
                sim.applyGate("h", q);
            for (std::size_t q = 0; q < n; ++q)
                sim.applyGate("x", q);
            sim.applyMultiControlledGate("z", controls, n - 1);
            for (std::size_t q = 0; q < n; ++q)
                sim.applyGate("x", q);
            for (std::size_t q = 0; q < n; ++q)
                sim.applyGate("h", q);
        }
        std::vector<Dw> amplitudes;
        amplitudes.reserve(std::size_t{1} << n);
        for (std::size_t i = 0; i < (std::size_t{1} << n); ++i) {
            std::vector<bool> bits(n);
            for (std::size_t q = 0; q < n; ++q)
                bits[q] = ((i >> q) & 1U) != 0;
            amplitudes.push_back(sim.amplitude(bits));
        }
        return amplitudes;
    };

    const auto reference = run(dd::exact::NormalizationStrategy::Inverse);
    // Sanity: the search actually amplified something, so an all-equal vector
    // cannot pass this test vacuously.
    const Dw marked = reference.back(); // |1...1>
    EXPECT_GT(marked.normSquared().toComplexDouble().real(), 0.9);

    for (const auto strategy :
         {dd::exact::NormalizationStrategy::None, dd::exact::NormalizationStrategy::Gcd}) {
        const auto amplitudes = run(strategy);
        ASSERT_EQ(amplitudes.size(), reference.size());
        for (std::size_t i = 0; i < amplitudes.size(); ++i) {
            EXPECT_EQ(amplitudes[i], reference[i])
                << "strategy " << static_cast<int>(strategy) << " differs at basis state " << i;
        }
    }
}

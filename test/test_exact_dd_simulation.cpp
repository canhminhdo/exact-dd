#include "dd/exact/ExactDDSimulation.hpp"

#include <gtest/gtest.h>

using dd::exact::Dw;
using dd::exact::ExactDDSimulation;

namespace {
const Dw kHalf(1, 0, 0, 0, 2); // 1/2
std::vector<bool> bitsAt(std::size_t nqubits, std::size_t q0, bool v0, std::size_t q1, bool v1) {
    std::vector<bool> b(nqubits, false);
    b[q0] = v0;
    b[q1] = v1;
    return b;
}
std::vector<bool> bitsAt(std::size_t nqubits, std::size_t q0, bool v0, std::size_t q1, bool v1, std::size_t q2,
                          bool v2) {
    std::vector<bool> b = bitsAt(nqubits, q0, v0, q1, v1);
    b[q2] = v2;
    return b;
}
} // namespace

// --- module 4: basic driver behavior ---------------------------------

TEST(ExactDDSimulation, InitialStateIsZeroWithUnitNorm) {
    ExactDDSimulation sim(2);
    EXPECT_EQ(sim.amplitude({false, false}), Dw::one());
    EXPECT_TRUE(sim.amplitude({true, false}).isZero());
    EXPECT_EQ(sim.normSquared(), Dw::one());
}

TEST(ExactDDSimulation, HGateGivesEqualSuperposition) {
    ExactDDSimulation sim(1);
    sim.applyGate("h", 0);
    const Dw invSqrt2(0, 1, 0, -1, 2);
    EXPECT_EQ(sim.amplitude({false}), invSqrt2);
    EXPECT_EQ(sim.amplitude({true}), invSqrt2);
    EXPECT_EQ(sim.normSquared(), Dw::one());
}

TEST(ExactDDSimulation, MeasureProjectsAndTracksExactUnnormalizedProbability) {
    ExactDDSimulation sim(1);
    sim.applyGate("h", 0);
    const auto result = sim.measure(0, true);
    // Not renormalized: probability is exactly 1/2, and the surviving
    // amplitude is still 1/sqrt(2) (not rescaled to 1).
    EXPECT_EQ(result.probability, kHalf);
    EXPECT_EQ(sim.normSquared(), kHalf);
    const Dw invSqrt2(0, 1, 0, -1, 2);
    EXPECT_EQ(sim.amplitude({true}), invSqrt2);
    EXPECT_TRUE(sim.amplitude({false}).isZero());
}

TEST(ExactDDSimulation, RejectsUnknownGateName) {
    ExactDDSimulation sim(1);
    EXPECT_THROW(sim.applyGate("not_a_gate", 0), std::invalid_argument);
}

TEST(ExactDDSimulation, FidelityOfStateWithItselfIsOne) {
    ExactDDSimulation sim(2);
    sim.applyGate("h", 0);
    sim.applyControlledGate("x", 0, 1);
    EXPECT_NEAR(sim.fidelity(sim.state()), 1.0, 1e-12);
}

TEST(ExactDDSimulation, AppliesVAndSxGates) {
    ExactDDSimulation sim(1);
    sim.applyGate("v", 0);
    sim.applyGate("vdg", 0);
    EXPECT_EQ(sim.amplitude({false}), Dw::one());
    EXPECT_TRUE(sim.amplitude({true}).isZero());

    ExactDDSimulation sim2(1);
    sim2.applyGate("sx", 0);
    sim2.applyGate("sxdg", 0);
    EXPECT_EQ(sim2.amplitude({false}), Dw::one());
    EXPECT_TRUE(sim2.amplitude({true}).isZero());
}

// --- module 5: scenario tests (GHZ + Clifford+T teleportation) -------

TEST(ExactDDSimulationScenario, GhzState) {
    // H(q2); CNOT(2->1); CNOT(2->0) => (|000> + |111>)/sqrt(2), exactly.
    ExactDDSimulation sim(3);
    sim.applyGate("h", 2);
    sim.applyControlledGate("x", 2, 1);
    sim.applyControlledGate("x", 2, 0);

    const Dw invSqrt2(0, 1, 0, -1, 2);
    EXPECT_EQ(sim.amplitude({false, false, false}), invSqrt2);
    EXPECT_EQ(sim.amplitude({true, true, true}), invSqrt2);
    EXPECT_TRUE(sim.amplitude({true, false, false}).isZero());
    EXPECT_TRUE(sim.amplitude({false, true, false}).isZero());
    EXPECT_TRUE(sim.amplitude({false, false, true}).isZero());
    EXPECT_EQ(sim.normSquared(), Dw::one());

    // Measuring q2 collapses the entangled triple consistently, with exact
    // probability 1/2 for each outcome.
    ExactDDSimulation branch0(3);
    branch0.applyGate("h", 2);
    branch0.applyControlledGate("x", 2, 1);
    branch0.applyControlledGate("x", 2, 0);
    const auto r0 = branch0.measure(2, false);
    EXPECT_EQ(r0.probability, kHalf);
    // Unnormalized by design: the surviving amplitude stays 1/sqrt(2)
    // rather than being rescaled to 1.
    EXPECT_EQ(branch0.amplitude({false, false, false}), invSqrt2);
    EXPECT_TRUE(branch0.amplitude({true, true, true}).isZero());
}

// Clifford+T-restricted quantum teleportation: teleport the state H|0> (the
// |+> state, itself Clifford) from qubit 0 to qubit 2, using an EPR pair on
// qubits 1,2 and Pauli corrections (X^m1 Z^m0) conditioned on the two
// classical measurement bits from the Bell-basis measurement of qubits 0,1.
TEST(ExactDDSimulationScenario, CliffordTeleportationAllFourOutcomes) {
    const Dw eighth(1, 0, 0, 0, 6); // 1/8 == 1/sqrt(2)^6

    for (bool m0 : {false, true}) {
        for (bool m1 : {false, true}) {
            ExactDDSimulation sim(3);
            sim.applyGate("h", 0); // prepare |+> on qubit 0 (the state to teleport)
            sim.applyGate("h", 1); // EPR pair on qubits 1,2
            sim.applyControlledGate("x", 1, 2);
            sim.applyControlledGate("x", 0, 1); // Bell-basis measurement gates on qubits 0,1
            sim.applyGate("h", 0);
            sim.measure(0, m0);
            sim.measure(1, m1);
            if (m1)
                sim.applyGate("x", 2);
            if (m0)
                sim.applyGate("z", 2);

            // Joint outcome probability is exactly 1/4 (uniform, independent
            // of the teleported state), split evenly across qubit 2's two
            // basis values once correctly Pauli-corrected back to |+>: each
            // amplitude has |amp|^2 == 1/8, exactly.
            const Dw amp0 = sim.amplitude(bitsAt(3, 0, m0, 1, m1, 2, false));
            const Dw amp1 = sim.amplitude(bitsAt(3, 0, m0, 1, m1, 2, true));
            EXPECT_EQ(amp0.normSquared(), eighth) << "m0=" << m0 << " m1=" << m1;
            EXPECT_EQ(amp1.normSquared(), eighth) << "m0=" << m0 << " m1=" << m1;
            // |+> is a +1 eigenstate of X and Z|+> = |->, so after a
            // *correct* correction the two amplitudes must be equal
            // (not negated) regardless of m0 and m1.
            EXPECT_EQ(amp0, amp1) << "m0=" << m0 << " m1=" << m1;
        }
    }
}

// --- module 6: Hadamard transform, entanglement, QFT, Shor -----------

TEST(ExactDDSimulationScenario, HadamardTransformGivesEqualSuperposition) {
    // H^(x)3 |000> = equal superposition over all 8 basis states, each with
    // amplitude exactly (1/sqrt(2))^3 -- computed via Dw arithmetic (not
    // hand-transcribed) so the exact canonical form is whatever Dw actually
    // produces, matching the lesson from the diagram-string test fix.
    ExactDDSimulation sim(3);
    sim.applyGate("h", 0);
    sim.applyGate("h", 1);
    sim.applyGate("h", 2);

    const Dw invSqrt2(0, 1, 0, -1, 2);
    const Dw invSqrt2Cubed = invSqrt2 * invSqrt2 * invSqrt2;
    for (std::size_t idx = 0; idx < 8; ++idx) {
        const std::vector<bool> bits{(idx & 1U) != 0, (idx & 2U) != 0, (idx & 4U) != 0};
        EXPECT_EQ(sim.amplitude(bits), invSqrt2Cubed) << "basis index " << idx;
    }
    EXPECT_EQ(sim.normSquared(), Dw::one());
}

TEST(ExactDDSimulationScenario, BellStateAmplitudes) {
    // The canonical minimal entangled state: (|00> + |11>) / sqrt(2).
    ExactDDSimulation sim(2);
    sim.applyGate("h", 0);
    sim.applyControlledGate("x", 0, 1);

    const Dw invSqrt2(0, 1, 0, -1, 2);
    EXPECT_EQ(sim.amplitude({false, false}), invSqrt2);
    EXPECT_EQ(sim.amplitude({true, true}), invSqrt2);
    EXPECT_TRUE(sim.amplitude({true, false}).isZero());
    EXPECT_TRUE(sim.amplitude({false, true}).isZero());
    EXPECT_EQ(sim.normSquared(), Dw::one());
}

// Entanglement swapping (port of examples/entangleswap.qw): q0,q1 and q2,q3
// each start as independent Bell pairs; a Bell-basis measurement of q1,q2
// (the qubits that meet in the middle) leaves q0 and q3 entangled with each
// other, even though they never interacted directly. After the standard
// X^m2 Z^m1 correction on q3, the result is always the same Bell state
// (|00>+|11>)/sqrt(2) on (q0,q3), regardless of which of the 4 measurement
// outcomes occurred -- verified exactly for all 4 combinations.
TEST(ExactDDSimulationScenario, EntanglementSwappingProducesBellPair) {
    const Dw invSqrt2(0, 1, 0, -1, 2);
    const Dw invSqrt2Cubed = invSqrt2 * invSqrt2 * invSqrt2;

    for (bool m2 : {false, true}) {
        for (bool m1 : {false, true}) {
            ExactDDSimulation sim(4);
            sim.applyGate("h", 0);
            sim.applyControlledGate("x", 0, 1);
            sim.applyGate("h", 2);
            sim.applyControlledGate("x", 2, 3);
            sim.applyControlledGate("x", 1, 2);
            sim.applyGate("h", 1);
            sim.measure(2, m2);
            sim.measure(1, m1);
            if (m2)
                sim.applyGate("x", 3);
            if (m1)
                sim.applyGate("z", 3);

            const Dw amp00 = sim.amplitude({false, m1, m2, false});
            const Dw amp11 = sim.amplitude({true, m1, m2, true});
            const Dw amp01 = sim.amplitude({false, m1, m2, true});
            const Dw amp10 = sim.amplitude({true, m1, m2, false});
            EXPECT_EQ(amp00, invSqrt2Cubed) << "m1=" << m1 << " m2=" << m2;
            EXPECT_EQ(amp11, invSqrt2Cubed) << "m1=" << m1 << " m2=" << m2;
            EXPECT_TRUE(amp01.isZero()) << "m1=" << m1 << " m2=" << m2;
            EXPECT_TRUE(amp10.isZero()) << "m1=" << m1 << " m2=" << m2;
        }
    }
}

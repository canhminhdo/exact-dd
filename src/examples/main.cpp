//
// Created by CanhDo on 2026/07/06.
//

#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include "dd/exact/ExactDDSimulation.hpp"

using namespace dd::exact;

// Grover's search over an n-qubit register (qubits 0..n-1, n = mark.size()),
// marking the single basis state |mark> (qubit i holds mark[i]).
//
// The oracle's and diffuser's (n-1)-controlled-Z are applied directly via
// DwPackage/ExactDDSimulation's built-in multi-controlled gate support
// (ExactDDSimulation::applyMultiControlledGate), so no ancillas or manual
// Toffoli decomposition are needed. The whole circuit is driven purely by
// gate multiplication on the state vector (applyMultiControlledGate ->
// DwPackage::applyOperation -> multiply(mEdge, vEdge)) -- no
// amplitude/innerProduct-based reflection shortcuts.
namespace {

// Phase-flips |mark> by sandwiching a native (n-1)-controlled Z with X gates
// on the qubits whose target bit is 0, mapping |mark> to |11..1> and back.
void oracle(ExactDDSimulation &sim, const std::vector<bool> &mark, const std::vector<std::size_t> &controls,
            std::size_t target) {
    for (std::size_t q = 0; q < mark.size(); ++q) {
        if (!mark[q])
            sim.applyGate("x", q);
    }
    sim.applyMultiControlledGate("z", controls, target);
    for (std::size_t q = 0; q < mark.size(); ++q) {
        if (!mark[q])
            sim.applyGate("x", q);
    }
}

// D = H^n (2|0><0| - I) H^n, up to the unobservable global phase -1 (built
// here as H^n X^n (C^(n-1) Z) X^n H^n, which realizes I - 2|0><0| instead).
void diffuser(ExactDDSimulation &sim, std::size_t nqubits, const std::vector<std::size_t> &controls,
              std::size_t target) {
    for (std::size_t q = 0; q < nqubits; ++q)
        sim.applyGate("h", q);
    for (std::size_t q = 0; q < nqubits; ++q)
        sim.applyGate("x", q);
    sim.applyMultiControlledGate("z", controls, target);
    for (std::size_t q = 0; q < nqubits; ++q)
        sim.applyGate("x", q);
    for (std::size_t q = 0; q < nqubits; ++q)
        sim.applyGate("h", q);
}

// Renders a bit vector as a "0101..."-style string, qubit 0 first.
std::string toBinaryString(const std::vector<bool> &bits) {
    std::string s(bits.size(), '0');
    for (std::size_t q = 0; q < bits.size(); ++q)
        s[q] = bits[q] ? '1' : '0';
    return s;
}

// Randomly picks a marked basis state over `nqubits` qubits.
std::vector<bool> randomMark(std::size_t nqubits) {
    static std::mt19937 rng{std::random_device{}()};
    std::bernoulli_distribution coin(0.5);
    std::vector<bool> mark(nqubits);
    for (std::size_t q = 0; q < nqubits; ++q)
        mark[q] = coin(rng);
    return mark;
}

} // namespace



void runGroverSearch(std::size_t nqubits) {
    const std::vector<bool> mark = randomMark(nqubits);
    const std::string markStr = toBinaryString(mark);
    const std::size_t target = nqubits - 1;
    std::vector<std::size_t> controls(nqubits - 1);
    std::iota(controls.begin(), controls.end(), 0);

    ExactDDSimulation sim(nqubits); // no ancillas needed
    for (std::size_t q = 0; q < nqubits; ++q)
        sim.applyGate("h", q);

    // Optimal iteration count for N=2^n, M=1 marked state: round(pi/4*sqrt(N) - 1/2).
    const double n = static_cast<double>(1ULL << nqubits);
    const int iterations = static_cast<int>(std::lround(M_PI / 4.0 * std::sqrt(n) - 0.5));
    std::cout << "Running " << iterations << " Grover iteration(s) over " << nqubits << " qubits (N=" << n
               << ") for target |" << markStr << ">\n";

    for (int it = 0; it < iterations; ++it) {
        oracle(sim, mark, controls, target);
        diffuser(sim, nqubits, controls, target);
        if (it % 100 == 0 or it == iterations - 1) {
            auto markedProb = sim.amplitude(mark).normSquared().toComplexFloat().real();
            std::cout << "After iteration " << (it + 1) << ": P("<< markStr << ") = " << markedProb << "\n";
        }
    }
}

// H^(x)n |0...0> = equal superposition over all 2^n basis states, each with
// exact amplitude (1/sqrt(2))^n.
void runHadamardTransform(std::size_t n) {
    ExactDDSimulation sim(n);
    for (std::size_t q = 0; q < n; ++q)
        sim.applyGate("h", q);

    std::cout << "Hadamard transform over " << n << " qubits: equal superposition\n";
    const std::vector<bool> allZero(n, false);
    const std::vector<bool> allOne(n, true);
    std::cout << "  amplitude(0...0) = " << sim.amplitude(allZero).toString() << "\n";
    std::cout << "  amplitude(1...1) = " << sim.amplitude(allOne).toString() << "\n";
    // std::cout << "  normSquared() = " << sim.normSquared().toString() << "\n";
    sim.package().printVectorDiagram(sim.state());
}

// GHZ state: H on the top qubit, then CNOT from it to every other qubit,
// giving (|0...0> + |1...1>) / sqrt(2) -- all n qubits exactly entangled.
void runEntanglement(std::size_t n) {
    ExactDDSimulation sim(n);
    const std::size_t control = n - 1;
    sim.applyGate("h", control);
    for (std::size_t q = 0; q < control; ++q)
        sim.applyControlledGate("x", control, q);

    std::cout << "GHZ entanglement over " << n << " qubits\n";
    const std::vector<bool> allZero(n, false);
    const std::vector<bool> allOne(n, true);
    std::cout << "  amplitude(0...0) = " << sim.amplitude(allZero).toString() << "\n";
    std::cout << "  amplitude(1...1) = " << sim.amplitude(allOne).toString() << "\n";
    // std::cout << "  normSquared() = " << sim.normSquared().toString() << "\n";
    // sim.package().printVectorDiagram(sim.state());
}

int main() {
    runGroverSearch(15);
    runHadamardTransform(100);
    runEntanglement(100);
    return 0;
}

//
// Created by CanhDo on 2026/07/31.
//

//
// Scaling benchmark for the exact DD package. Prints wall time and node
// counts for the operations whose cost is expected to depend on the qubit
// count in a way the decision diagram itself does not justify, so that a
// change to DwPackage/Dw can be judged against a recorded baseline rather
// than against intuition.
//

#include "dd/exact/DwPackage.hpp"
#include "dd/exact/DwGateMatrixDefinitions.hpp"
#include "dd/exact/ExactDDSimulation.hpp"

#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using namespace dd::exact;

namespace {

using Clock = std::chrono::steady_clock;

double elapsedMs(const Clock::time_point start) {
    return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

void reportHeader(const std::string &title) {
    std::cout << "\n=== " << title << " ===\n";
    std::cout << std::left << std::setw(8) << "qubits" << std::right << std::setw(14) << "time (ms)"
              << std::setw(12) << "vNodes" << std::setw(12) << "mNodes" << "\n";
}

void reportRow(std::size_t nqubits, double ms, std::size_t vNodes, std::size_t mNodes) {
    std::cout << std::left << std::setw(8) << nqubits << std::right << std::setw(14) << std::fixed
              << std::setprecision(2) << ms << std::setw(12) << vNodes << std::setw(12) << mNodes << "\n";
    std::cout.unsetf(std::ios::floatfield);
}

// <+...+|+...+> == 1. The DD for |+>^(x)n is n nodes wide, so the only thing
// that can make this scale badly is innerProduct's own recursion.
void benchInnerProductUniform(const std::vector<std::size_t> &sizes) {
    reportHeader("innerProduct on |+>^(x)n  (expect exactly 1)");
    for (const std::size_t n : sizes) {
        DwPackage pkg(n);
        const std::vector<BasisState> plus(n, BasisState::Plus);
        const auto state = pkg.makeBasisState(n, plus, 0);

        const auto start = Clock::now();
        const Dw norm = pkg.innerProduct(state, state);
        const double ms = elapsedMs(start);

        reportRow(n, ms, pkg.vNodeCount(), pkg.mNodeCount());
        if (!norm.isOne())
            std::cout << "  !! expected 1, got " << norm.toString() << "\n";
    }
}

// GHZ = (|0...0> + |1...1>)/sqrt(2): builds an entangled state via
// n-1 CNOTs, then measures every qubit's outcome probability without
// collapsing. This is the shape of QRAT's own measurement path, where each
// measureOneQubit does a multiply plus an innerProduct.
void benchGhzMeasureSweep(const std::vector<std::size_t> &sizes) {
    reportHeader("GHZ build + per-qubit measureOneQubit sweep");
    for (const std::size_t n : sizes) {
        DwPackage pkg(n);
        const auto start = Clock::now();

        auto state = pkg.makeZeroState();
        pkg.incRef(state);
        const std::size_t control = n - 1;
        state = pkg.applyOperation(pkg.makeSingleQubitGateDD(control, gates::h()), state);
        for (std::size_t q = 0; q < control; ++q)
            state = pkg.applyOperation(pkg.makeControlledSingleQubitGateDD(control, q, gates::x()), state);

        for (std::size_t q = 0; q < n; ++q) {
            const auto r0 = pkg.measureOneQubit(state, q, false);
            const auto r1 = pkg.measureOneQubit(state, q, true);
            const Dw total = r0.probability + r1.probability;
            if (!total.isOne())
                std::cout << "  !! qubit " << q << ": probabilities sum to " << total.toString() << "\n";
        }

        const double ms = elapsedMs(start);
        reportRow(n, ms, pkg.vNodeCount(), pkg.mNodeCount());
    }
}

// |x><y| for two n-qubit states. outerProductRec already carries a per-call
// memo keyed on (x.p, y.p, level); this exercises it on both a maximally
// collapsed operand (|+>^(x)n reduces to a single terminal edge, so every
// level is a pass-through) and a genuinely entangled one (GHZ).
void benchOuterProduct(const std::vector<std::size_t> &sizes) {
    reportHeader("outerProduct |+>^(x)n  x  <+|^(x)n");
    for (const std::size_t n : sizes) {
        DwPackage pkg(n);
        const std::vector<BasisState> plus(n, BasisState::Plus);
        const auto state = pkg.makeBasisState(n, plus, 0);

        const auto start = Clock::now();
        const auto op = pkg.outerProduct(state, state);
        const double ms = elapsedMs(start);

        reportRow(n, ms, pkg.vNodeCount(), pkg.mNodeCount());
        if (op.w.isZero())
            std::cout << "  !! unexpected zero operator\n";
    }

    reportHeader("outerProduct GHZ x GHZ");
    for (const std::size_t n : sizes) {
        DwPackage pkg(n);
        auto state = pkg.makeZeroState();
        pkg.incRef(state);
        const std::size_t control = n - 1;
        state = pkg.applyOperation(pkg.makeSingleQubitGateDD(control, gates::h()), state);
        for (std::size_t q = 0; q < control; ++q)
            state = pkg.applyOperation(pkg.makeControlledSingleQubitGateDD(control, q, gates::x()), state);

        const auto start = Clock::now();
        const auto op = pkg.outerProduct(state, state);
        const double ms = elapsedMs(start);

        reportRow(n, ms, pkg.vNodeCount(), pkg.mNodeCount());
        if (op.w.isZero())
            std::cout << "  !! unexpected zero operator\n";
    }
}

// Phase-flips |1...1> with a native (n-1)-controlled Z.
void groverOracle(ExactDDSimulation &sim, std::size_t n, const std::vector<std::size_t> &controls) {
    sim.applyMultiControlledGate("z", controls, n - 1);
}

void groverDiffuser(ExactDDSimulation &sim, std::size_t n, const std::vector<std::size_t> &controls) {
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

// Grover over n qubits marking |1...1>, run for the optimal iteration count.
// Exercises multiply and add far harder than the state-construction cases:
// the diffuser's superposition-wide reflection is where add's recursion
// actually has shared subgraphs to re-walk.
void benchGrover(const std::vector<std::size_t> &sizes) {
    reportHeader("Grover (marking |1...1>, optimal iterations)");
    for (const std::size_t n : sizes) {
        std::vector<std::size_t> controls(n - 1);
        std::iota(controls.begin(), controls.end(), 0);

        const auto start = Clock::now();
        ExactDDSimulation sim(n);
        for (std::size_t q = 0; q < n; ++q)
            sim.applyGate("h", q);

        const double size = static_cast<double>(1ULL << n);
        const int iterations = static_cast<int>(std::lround((M_PI / 4.0 * std::sqrt(size)) - 0.5));
        for (int it = 0; it < iterations; ++it) {
            groverOracle(sim, n, controls);
            groverDiffuser(sim, n, controls);
        }
        const double ms = elapsedMs(start);

        const std::vector<bool> mark(n, true);
        const double prob = static_cast<double>(sim.amplitude(mark).normSquared().toComplexFloat().real());
        reportRow(n, ms, sim.package().vNodeCount(), sim.package().mNodeCount());
        std::cout << "  " << iterations << " iteration(s), P(|1...1>) = " << prob << "\n";
    }
}

} // namespace

// With no arguments, runs the whole suite at sizes that stay well under a
// second. Given qubit counts as arguments, runs only Grover at those sizes --
// Grover is the case whose coefficients grow fastest, so it is what the
// EXACT_DD_WITH_GMP choice hinges on, and the large sizes are far too slow to
// belong in the default run:
//     ExactDDBench 16 18 20
int main(int argc, char **argv) {
    std::cout << "exact-dd scaling benchmark\n";

    if (argc > 1) {
        std::vector<std::size_t> sizes;
        for (int i = 1; i < argc; ++i)
            sizes.push_back(static_cast<std::size_t>(std::stoul(argv[i])));
        benchGrover(sizes);
        std::cout << "\ndone\n";
        return 0;
    }

    benchInnerProductUniform({4, 8, 12, 16});
    benchGhzMeasureSweep({4, 8, 12, 16});
    benchOuterProduct({4, 8, 12, 16});
    benchGrover({6, 10, 14});

    std::cout << "\ndone\n";
    return 0;
}

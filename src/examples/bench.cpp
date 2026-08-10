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
// NOTE ON BASELINES: with EXACT_DD_STATISTICS=ON (the default) the timings
// below include one counter increment per unique-table and compute-table
// probe, i.e. on the hottest path in the package. Record and compare
// baselines at a fixed setting, and prefer -DEXACT_DD_STATISTICS=OFF for a
// clean measurement of a DwPackage/Dw change.
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

// H on every qubit of |0...0>: the n-fold Hadamard transform, giving the
// equal superposition in which all 2^n amplitudes are exactly 1/sqrt(2)^n.
// The state stays a product state throughout, so its DD is one node per
// level and no coefficient ever grows past a numerator of 1.
void benchHadamardTransform(const std::vector<std::size_t> &sizes) {
    reportHeader("Hadamard transform H^(x)n on |0...0>");
    for (const std::size_t n : sizes) {
        const auto start = Clock::now();
        ExactDDSimulation sim(n);
        for (std::size_t q = 0; q < n; ++q)
            sim.applyGate("h", q);
        const double ms = elapsedMs(start);

        reportRow(n, ms, sim.package().vNodeCount(), sim.package().mNodeCount());

        // Both extremes must equal 1/sqrt(2)^n on the nose. This is a
        // canonical-form check as much as a correctness one: a stray factor
        // of sqrt(2) left unreduced in k would show up here and nowhere else
        // in this bench, since every other quantity stays 1.
        const Dw expected(1, 0, 0, 0, n);
        const Dw first = sim.amplitude(std::vector<bool>(n, false));
        const Dw last = sim.amplitude(std::vector<bool>(n, true));
        if (first != expected)
            std::cout << "  !! amplitude(0...0) = " << first.toString() << ", expected "
                      << expected.toString() << "\n";
        if (last != expected)
            std::cout << "  !! amplitude(1...1) = " << last.toString() << ", expected "
                      << expected.toString() << "\n";
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
void benchGrover(const std::vector<std::size_t> &sizes, bool withStatistics = false,
                 NormalizationStrategy strategy = NormalizationStrategy::Inverse) {
    reportHeader("Grover (marking |1...1>, optimal iterations)");
    for (const std::size_t n : sizes) {
        std::vector<std::size_t> controls(n - 1);
        std::iota(controls.begin(), controls.end(), 0);

        const auto start = Clock::now();
        ExactDDSimulation sim(n, strategy);
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

        if (withStatistics) {
            sim.package().printStatistics();
            std::cout << "\n";
        }
    }
}

// One full statistics dump on a workload that touches all four tables: the
// GHZ build hits both unique tables and the matrix-vector cache, the
// outerProduct populates the matrix unique table further, and the forced
// collection gives the memory managers entries to hand back out. Printed
// once rather than per size -- the point is the shape of the report, not a
// scaling curve, so reportHeader/reportRow are deliberately left alone and
// the timing tables above stay comparable against recorded baselines.
void benchStatisticsReport(std::size_t n) {
    std::cout << "\n=== statistics report (" << n << " qubits) ===\n";

    DwPackage pkg(n);
    auto state = pkg.makeZeroState();
    pkg.incRef(state);
    const std::size_t control = n - 1;
    state = pkg.applyOperation(pkg.makeSingleQubitGateDD(control, gates::h()), state);
    for (std::size_t q = 0; q < control; ++q)
        state = pkg.applyOperation(pkg.makeControlledSingleQubitGateDD(control, q, gates::x()), state);

    const auto op = pkg.outerProduct(state, state);
    if (op.w.isZero())
        std::cout << "  !! unexpected zero operator\n";
    pkg.garbageCollect(true);

    pkg.printStatistics();
    std::cout << "\n";
}

} // namespace

// With no arguments, runs the whole suite at sizes that stay well under a
// second. Given qubit counts as arguments, runs only Grover at those sizes --
// Grover is the case whose coefficients grow fastest, so it is what the
// EXACT_DD_WITH_GMP choice hinges on, and the large sizes are far too slow to
// belong in the default run:
//     ExactDDBench 16 18 20
// The argv form additionally dumps a full statistics report per size, which is
// how "where does Grover's cost go" gets answered at sizes that matter. The
// no-argument form deliberately does not, so its output stays comparable
// against recorded baselines.
int main(int argc, char **argv) {
    std::cout << "exact-dd scaling benchmark\n";

    if (argc > 1) {
        std::vector<std::size_t> sizes;
        auto strategy = NormalizationStrategy::Inverse;
        std::string strategyName = "inverse";
        for (int i = 1; i < argc; ++i) {
            const std::string arg(argv[i]);
            if (arg.rfind("--strategy=", 0) == 0) {
                strategyName = arg.substr(11);
                if (strategyName == "none") {
                    strategy = NormalizationStrategy::None;
                } else if (strategyName == "gcd") {
                    strategy = NormalizationStrategy::Gcd;
                } else if (strategyName != "inverse") {
                    std::cerr << "unknown --strategy=" << strategyName
                              << " (expected none|inverse|gcd)\n";
                    return 1;
                }
                continue;
            }
            sizes.push_back(static_cast<std::size_t>(std::stoul(arg)));
        }
        std::cout << "normalization strategy: " << strategyName << "\n";
        benchGrover(sizes, true, strategy);
        std::cout << "\ndone\n";
        return 0;
    }

    benchInnerProductUniform({4, 8, 12, 16});
    benchGhzMeasureSweep({4, 8, 12, 16});
    benchOuterProduct({4, 8, 12, 16});
    benchHadamardTransform({64, 128, 256, 512});
    benchGrover({6, 10, 14});
    benchStatisticsReport(12);

    std::cout << "\ndone\n";
    return 0;
}

#ifndef DD_EXACT_EXACT_DD_SIMULATION_HPP
#define DD_EXACT_EXACT_DD_SIMULATION_HPP

#include "dd/exact/DwPackage.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace dd::exact {

/**
 * Standalone driver that exercises DwPackage end-to-end on Clifford+T
 * circuits: applies single- and two-qubit gates (with any number of
 * positive controls), performs computational-basis measurements, and
 * reports exact amplitudes and probabilities. Not wired into QRAT's
 * Interpreter/SearchGraph/DTMC pipeline -- that integration is a
 * separate, later plan.
 *
 * Normalization: after a measurement, the projected state is kept
 * unnormalized (dividing by sqrt(probability) is generally not exact in
 * D[w]). The state's exact squared norm is always recomputable via
 * DwPackage::innerProduct, so probabilities of further measurements and
 * fidelity are computed as ratios of exact quantities, converted to
 * double only at the very last step (see fidelity()).
 *
 * state_ is kept incRef'd with pkg_ as its persistent GC root: every gate
 * application/measurement goes through DwPackage::applyOperation(), which
 * incRef()s the new state, decRef()s the old one, and (opportunistically,
 * once pkg_'s node tables grow large enough) garbage-collects nodes that
 * became unreachable as a result.
 */
class ExactDDSimulation {
public:
    explicit ExactDDSimulation(std::size_t nqubits);

    [[nodiscard]] DwPackage &package() { return pkg_; }
    [[nodiscard]] const DwPackage::vEdge &state() const { return state_; }

    /// Resets to the |0...0> basis state.
    void reset();

    /// Applies a Clifford+T gate (see DwGateMatrixDefinitions::byName) to
    /// `target`. Throws std::invalid_argument for gates outside D[w].
    void applyGate(const std::string &gateName, std::size_t target);

    /// Applies a controlled Clifford+T gate; control may have either a
    /// higher or lower qubit index than target.
    void applyControlledGate(const std::string &gateName, std::size_t control, std::size_t target);

    /// Applies a Clifford+T gate to `target` under any number of positive
    /// controls (empty vector = uncontrolled).
    void applyMultiControlledGate(const std::string &gateName, const std::vector<std::size_t> &controls,
                                  std::size_t target);

    /// Applies a two-qubit gate (see DwGateMatrixDefinitions::twoQubitByName:
    /// swap, iswap, iswapdg, dcx) to targets (target0, target1), where
    /// target0 is the matrix's more significant bit.
    void applyTwoQubitGate(const std::string &gateName, std::size_t target0, std::size_t target1);

    /// Two-qubit gate under any number of positive controls.
    void applyControlledTwoQubitGate(const std::string &gateName, const std::vector<std::size_t> &controls,
                                     std::size_t target0, std::size_t target1);

    struct MeasurementResult {
        DwPackage::vEdge state; ///< projected, unnormalized post-measurement state
        Dw probability;         ///< exact probability of this outcome
    };

    /// Projects onto `qubit == outcome`, updates the (unnormalized)
    /// current state, and returns the exact outcome probability.
    MeasurementResult measure(std::size_t qubit, bool outcome);

    /// Exact squared norm of the current state (1 after |0...0> and any
    /// sequence of unitary gates alone; may drop below 1 after measure()).
    [[nodiscard]] Dw normSquared() const;

    [[nodiscard]] Dw amplitude(const std::vector<bool> &bits) const;

    /// |<state|other>|^2 / (<state|state> * <other|other>), computed
    /// exactly and converted to double only for the final division.
    [[nodiscard]] double fidelity(const DwPackage::vEdge &other) const;

private:
    DwPackage pkg_;
    DwPackage::vEdge state_;
};

} // namespace dd::exact

#endif // DD_EXACT_EXACT_DD_SIMULATION_HPP

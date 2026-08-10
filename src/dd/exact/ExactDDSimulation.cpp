#include "dd/exact/ExactDDSimulation.hpp"
#include "dd/exact/DwGateMatrixDefinitions.hpp"
#include <utility>

namespace dd::exact {

ExactDDSimulation::ExactDDSimulation(std::size_t nqubits, NormalizationStrategy strategy)
    : pkg_(nqubits, strategy) {
    reset();
}

void ExactDDSimulation::reset() {
    // decRef on the default-constructed (p == nullptr) initial state_ is a
    // safe no-op; this keeps state_ a properly incRef'd root of pkg_'s
    // reference-counted garbage collector across resets.
    pkg_.decRef(state_);
    state_ = pkg_.makeZeroState();
    pkg_.incRef(state_);
}

void ExactDDSimulation::applyGate(const std::string &gateName, std::size_t target) {
    const auto &matrix = gates::byName(gateName);
    state_ = pkg_.applyOperation(pkg_.makeSingleQubitGateDD(target, matrix), state_);
}

void ExactDDSimulation::applyControlledGate(const std::string &gateName, std::size_t control, std::size_t target) {
    const auto &matrix = gates::byName(gateName);
    state_ = pkg_.applyOperation(pkg_.makeControlledSingleQubitGateDD(control, target, matrix), state_);
}

void ExactDDSimulation::applyMultiControlledGate(const std::string &gateName,
                                                 const std::vector<std::size_t> &controls, std::size_t target) {
    const auto &matrix = gates::byName(gateName);
    state_ = pkg_.applyOperation(pkg_.makeControlledSingleQubitGateDD(controls, target, matrix), state_);
}

void ExactDDSimulation::applyTwoQubitGate(const std::string &gateName, std::size_t target0, std::size_t target1) {
    const auto &matrix = gates::twoQubitByName(gateName);
    state_ = pkg_.applyOperation(pkg_.makeTwoQubitGateDD(target0, target1, matrix), state_);
}

void ExactDDSimulation::applyControlledTwoQubitGate(const std::string &gateName,
                                                    const std::vector<std::size_t> &controls, std::size_t target0,
                                                    std::size_t target1) {
    const auto &matrix = gates::twoQubitByName(gateName);
    state_ = pkg_.applyOperation(pkg_.makeControlledTwoQubitGateDD(controls, target0, target1, matrix), state_);
}

ExactDDSimulation::MeasurementResult ExactDDSimulation::measure(std::size_t qubit, bool outcome) {
    auto result = pkg_.measureOneQubit(state_, qubit, outcome);
    pkg_.incRef(result.state);
    pkg_.decRef(state_);
    state_ = std::move(result.state);
    pkg_.garbageCollect();
    return {state_, std::move(result.probability)};
}

Dw ExactDDSimulation::normSquared() const { return pkg_.innerProduct(state_, state_); }

Dw ExactDDSimulation::amplitude(const std::vector<bool> &bits) const { return pkg_.amplitude(state_, bits); }

double ExactDDSimulation::fidelity(const DwPackage::vEdge &other) const { return pkg_.fidelity(state_, other); }

} // namespace dd::exact

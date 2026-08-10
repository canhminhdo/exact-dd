#include "dd/exact/DwPackage.hpp"
#include "dd/exact/DwGateMatrixDefinitions.hpp"
#include "dd/exact/statistics/StatisticsConfig.hpp"
#include "utility/HashUtil.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_set>
#include <utility>
#include <vector>

namespace dd::exact {

namespace {
/// GCD of a nonempty vector of nonzero Dw values (paper's Algorithm 3),
/// reduced pairwise via Dw::gcd().
Dw gcdOfMany(const std::vector<Dw> &weights) {
    Dw g = weights.front();
    for (std::size_t i = 1; i < weights.size(); ++i)
        g = g.gcd(weights[i]);
    return g;
}

/// Normalizes a node's nonzero outgoing edge weights by their (canonical-
/// associate) GCD, dividing each in place, and returns the GCD as the
/// factor to propagate to the incoming edge (paper's Algorithm 3).
Dw gcdNormalize(const std::vector<Dw *> &weights) {
    std::vector<Dw> values;
    values.reserve(weights.size());
    for (const Dw *w : weights)
        values.push_back(*w);
    const Dw eta = gcdOfMany(values).reduceAssociate();
    if (!eta.isOne()) {
        const Dw etaInv = eta.inverse();
        for (Dw *w : weights)
            *w = *w * etaInv;
    }
    return eta;
}
} // namespace

namespace detail {
std::size_t VKeyHash::operator()(const VKey &k) const noexcept {
    std::size_t h = std::hash<int>{}(k.var);
    for (const auto &edge : k.children) {
        h = HashUtil::combinedHash(h, std::hash<DwVNode *>{}(edge.p));
        h = HashUtil::combinedHash(h, edge.w.hash());
    }
    return h;
}

std::size_t MKeyHash::operator()(const MKey &k) const noexcept {
    std::size_t h = std::hash<int>{}(k.var);
    for (const auto &edge : k.children) {
        h = HashUtil::combinedHash(h, std::hash<DwMNode *>{}(edge.p));
        h = HashUtil::combinedHash(h, edge.w.hash());
    }
    return h;
}

std::size_t VVKeyHash::operator()(const VVKey &k) const noexcept {
    std::size_t h = std::hash<DwVNode *>{}(k.x);
    h = HashUtil::combinedHash(h, std::hash<DwVNode *>{}(k.y));
    h = HashUtil::combinedHash(h, std::hash<std::size_t>{}(k.level));
    return h;
}
} // namespace detail

DwPackage::DwPackage(std::size_t nqubits, NormalizationStrategy strategy)
    : nqubits_(nqubits), strategy_(strategy) {
    // Entry sizes are the only table statistic that is neither tracked on
    // the fly nor derivable from the map, so seed them once here. numBuckets
    // is deliberately NOT seeded: std::unordered_map rehashes as it grows,
    // so it is snapshotted at report time instead (TableStatistics).
    if constexpr (kStatisticsEnabled) {
        vUniqueStats_.entrySize = sizeof(decltype(vUnique_)::value_type);
        mUniqueStats_.entrySize = sizeof(decltype(mUnique_)::value_type);
        mvCacheStats_.entrySize = sizeof(decltype(mvCache_)::value_type);
        mmCacheStats_.entrySize = sizeof(decltype(mmCache_)::value_type);
    }
}

namespace {
int topVar(DwVNode *p) { return p != nullptr ? p->var : -1; }
int topVar(DwMNode *p) { return p != nullptr ? p->var : -1; }
} // namespace

DwPackage::vEdge DwPackage::makeVEdge(int var, std::array<vEdge, 2> children) {
    if (children[0].w.isZero())
        children[0] = vEdge::zero();
    if (children[1].w.isZero())
        children[1] = vEdge::zero();
    if (children[0].w.isZero() && children[1].w.isZero())
        return vEdge::zero();
    if (children[0] == children[1])
        return std::move(children[0]);

    // Node-weight normalization (see NormalizationStrategy): factor a value
    // eta out of the (about-to-be-created-or-shared) node's outgoing edge
    // weights, dividing each by it, and return eta as the incoming edge's
    // weight instead of always Dw::one(). Runs after the structural
    // collapse checks above, which remain valid regardless of the strategy.
    Dw eta = Dw::one();
    if (strategy_ == NormalizationStrategy::Inverse) {
        const bool leftmostIsSecond = children[0].w.isZero();
        eta = leftmostIsSecond ? children[1].w : children[0].w; // leftmost nonzero
        if (!eta.isOne()) {
            const Dw etaInv = eta.inverse();
            if (leftmostIsSecond) {
                children[1].w = Dw::one();
            } else {
                children[0].w = Dw::one();
                children[1].w = children[1].w * etaInv;
            }
        }
    } else if (strategy_ == NormalizationStrategy::Gcd) {
        std::vector<Dw *> nonzero;
        if (!children[0].w.isZero())
            nonzero.push_back(&children[0].w);
        if (!children[1].w.isZero())
            nonzero.push_back(&children[1].w);
        eta = gcdNormalize(nonzero);
    }

    detail::VKey key{var, children};
    vUniqueStats_.trackLookup();
    const auto [it, inserted] = vUnique_.try_emplace(std::move(key), nullptr);
    if (!inserted) {
        vUniqueStats_.trackHit();
        return {it->second, std::move(eta)};
    }

    DwVNode *raw = nullptr;
    raw = vMemory_.get();
    assert(raw != nullptr);
    raw->var = var;
    raw->e = std::move(children);
    it->second = raw;
    vUniqueStats_.trackInsert(vUnique_);
    return {raw, std::move(eta)};
}

DwPackage::mEdge DwPackage::makeMEdge(int var, std::array<mEdge, 4> children) {
    int nZero = 0;
    for (auto &c : children) {
        if (c.w.isZero()) {
            c = mEdge::zero();
            ++nZero;
        }
    }
    if (nZero == 4)
        return mEdge::zero();

    // Redundant-node elimination: identity-like pass-through (off-diagonal
    // zero, both diagonal children identical) doesn't depend on `var`.
    if (children[0] == children[3] && children[1].w.isZero() && children[2].w.isZero())
        return std::move(children[0]);

    // Node-weight normalization -- see the analogous block in makeVEdge().
    Dw eta = Dw::one();
    if (strategy_ == NormalizationStrategy::Inverse) {
        std::size_t idx = 0;
        while (children[idx].w.isZero())
            ++idx; // guaranteed to terminate: not all zero, checked above
        eta = children[idx].w;
        if (!eta.isOne()) {
            const Dw etaInv = eta.inverse();
            // children[idx] becomes exactly one by construction
            for (std::size_t i = idx; i < children.size(); ++i) {
                if (i == idx) {
                    children[i].w = Dw::one();
                } else {
                    children[i].w = children[i].w * etaInv;
                }
            }
        }
    } else if (strategy_ == NormalizationStrategy::Gcd) {
        std::vector<Dw *> nonzero;
        for (auto &c : children)
            if (!c.w.isZero())
                nonzero.push_back(&c.w);
        eta = gcdNormalize(nonzero);
    }

    detail::MKey key{var, children};
    mUniqueStats_.trackLookup();
    const auto [it, inserted] = mUnique_.try_emplace(std::move(key), nullptr);
    if (!inserted) {
        mUniqueStats_.trackHit();
        return {it->second, std::move(eta)};
    }

    DwMNode *raw = nullptr;
    raw = mMemory_.get();
    assert(raw != nullptr);
    raw->var = var;
    raw->e = std::move(children);
    it->second = raw;
    mUniqueStats_.trackInsert(mUnique_);
    return {raw, std::move(eta)};
}

// ---------------------------------------------------------------------
// State / operator construction
// ---------------------------------------------------------------------

DwPackage::vEdge DwPackage::makeBasisState(const std::vector<bool> &bits) {
    if (bits.size() != nqubits_)
        throw std::invalid_argument("makeBasisState: bits.size() must equal numQubits()");
    return makeBasisState(nqubits_, bits, 0);
}

DwPackage::vEdge DwPackage::makeBasisState(const std::size_t n, const std::vector<bool> &state,
                                                const std::size_t start) {
    if (start + n > nqubits_)
        throw std::invalid_argument("makeBasisState: start + n must not exceed numQubits()");
    if (state.size() != n)
        throw std::invalid_argument("makeBasisState: state.size() must equal n");
    vEdge e = vEdge::terminal(Dw::one());
    for (std::size_t i = 0; i < n; ++i) {
        const auto qubit = static_cast<int>(start + i);
        e = state[i] ? makeVEdge(qubit, {vEdge::zero(), std::move(e)})
                     : makeVEdge(qubit, {std::move(e), vEdge::zero()});
    }
    return e;
}

// Plus/Minus/Right/Left explicitly multiply the RUNNING edge's own
// accumulated weight (e.w) into each new child, rather than using a bare
// constant the way MQT Core's equivalent loop appears to (see the
// declaration's doc comment in DwPackage.hpp) -- hand-tracing a 2-qubit
// {Plus,Plus} chain with this approach gives the correct joint amplitude
// 1/2 at every index (invSqrt2 * invSqrt2, since e.w carries the prior
// qubit's invSqrt2 forward into the next multiplication), not just
// invSqrt2.
DwPackage::vEdge DwPackage::makeBasisState(const std::size_t n, const std::vector<BasisState> &state,
                                            const std::size_t start) {
    if (start + n > nqubits_)
        throw std::invalid_argument("makeBasisState: start + n must not exceed numQubits()");
    if (state.size() != n)
        throw std::invalid_argument("makeBasisState: state.size() must equal n");
    const Dw invSqrt2 = Dw::sqrt2().inverse();
    const Dw i = Dw::omega() * Dw::omega();
    vEdge e = vEdge::terminal(Dw::one());
    for (std::size_t idx = 0; idx < n; ++idx) {
        const std::size_t qubit = start + idx;
        switch (state[idx]) {
        case BasisState::Zero:
            e = makeVEdge(static_cast<int>(qubit), {std::move(e), vEdge::zero()});
            break;
        case BasisState::One:
            e = makeVEdge(static_cast<int>(qubit), {vEdge::zero(), std::move(e)});
            break;
        case BasisState::Plus:
            e = makeVEdge(static_cast<int>(qubit), {vEdge{e.p, e.w * invSqrt2}, vEdge{e.p, e.w * invSqrt2}});
            break;
        case BasisState::Minus:
            e = makeVEdge(static_cast<int>(qubit), {vEdge{e.p, e.w * invSqrt2}, vEdge{e.p, -(e.w * invSqrt2)}});
            break;
        case BasisState::Right:
            e = makeVEdge(static_cast<int>(qubit), {vEdge{e.p, e.w * invSqrt2}, vEdge{e.p, e.w * invSqrt2 * i}});
            break;
        case BasisState::Left:
            e = makeVEdge(static_cast<int>(qubit), {vEdge{e.p, e.w * invSqrt2}, vEdge{e.p, -(e.w * invSqrt2 * i)}});
            break;
        }
    }
    return e;
}

DwPackage::vEdge DwPackage::makeZeroState() {
    return makeBasisState(std::vector<bool>(nqubits_, false));
}

// Sidesteps the "exact state synthesis" problem (directly sampling a
// valid (a,b) with |a|^2+|b|^2=1 in D[w] in general needs Gaussian-
// integer-style factorization) by relying on unitarity instead: applying
// any sequence of exactly-representable single-qubit unitaries to |0>
// yields another exactly-normalized D[w] state, with zero new arithmetic
// needed. Starts from the windowed (narrow) |0> at `qubit`, so the result
// stays confined to `qubit` throughout (a qubit-local gate multiplied
// against a state already confined to that same qubit never touches
// other qubits), directly reusable as a kronecker() operand.
DwPackage::vEdge DwPackage::makeRandomSingleQubitState(std::mt19937 &rng, std::size_t depth, std::size_t qubit) {
    using GateMatrixFn = const std::array<Dw, 4> &(*)();
    static constexpr std::array<GateMatrixFn, 13> kGatePool{gates::i,   gates::x,  gates::y,  gates::z,
                                                             gates::h,   gates::s,  gates::sdg, gates::t,
                                                             gates::tdg, gates::v,  gates::vdg, gates::sx,
                                                             gates::sxdg};
    std::uniform_int_distribution<std::size_t> dist(0, kGatePool.size() - 1);

    vEdge e = makeBasisState(1, std::vector<BasisState>{BasisState::Zero}, qubit);
    for (std::size_t step = 0; step < depth; ++step) {
        const auto &matrix = kGatePool[dist(rng)]();
        e = multiply(makeSingleQubitGateDD(qubit, matrix), e);
    }
    return e;
}

DwPackage::mEdge DwPackage::makeIdentity() {
    mEdge e = mEdge::terminal(Dw::one());
    for (std::size_t qubit = 0; qubit < nqubits_; ++qubit)
        e = makeMEdge(static_cast<int>(qubit), {e, mEdge::zero(), mEdge::zero(), e});
    return e;
}

DwPackage::vEdge DwPackage::makeStateFromVectorRec(std::vector<Dw>::const_iterator begin,
                                                    std::vector<Dw>::const_iterator end) {
    if (end - begin == 1)
        return vEdge::terminal(*begin);
    const auto half = (end - begin) / 2;
    const int level = static_cast<int>(std::log2(end - begin)) - 1;
    vEdge e0 = makeStateFromVectorRec(begin, begin + half);
    vEdge e1 = makeStateFromVectorRec(begin + half, end);
    return makeVEdge(level, {std::move(e0), std::move(e1)});
}

DwPackage::vEdge DwPackage::makeStateFromVector(const std::vector<Dw> &amplitudes) {
    if (amplitudes.size() != (std::size_t{1} << nqubits_))
        throw std::invalid_argument("makeStateFromVector: amplitudes.size() must equal 2^numQubits()");
    return makeStateFromVectorRec(amplitudes.begin(), amplitudes.end());
}

DwPackage::mEdge DwPackage::makeDDFromMatrixRec(const std::vector<std::vector<Dw>> &matrix, std::size_t rowStart,
                                                 std::size_t rowEnd, std::size_t colStart, std::size_t colEnd) {
    if (rowEnd - rowStart == 1)
        return mEdge::terminal(matrix[rowStart][colStart]);
    const std::size_t rowMid = rowStart + (rowEnd - rowStart) / 2;
    const std::size_t colMid = colStart + (colEnd - colStart) / 2;
    const int level = static_cast<int>(std::log2(rowEnd - rowStart)) - 1;
    mEdge m00 = makeDDFromMatrixRec(matrix, rowStart, rowMid, colStart, colMid);
    mEdge m01 = makeDDFromMatrixRec(matrix, rowStart, rowMid, colMid, colEnd);
    mEdge m10 = makeDDFromMatrixRec(matrix, rowMid, rowEnd, colStart, colMid);
    mEdge m11 = makeDDFromMatrixRec(matrix, rowMid, rowEnd, colMid, colEnd);
    return makeMEdge(level, {std::move(m00), std::move(m01), std::move(m10), std::move(m11)});
}

DwPackage::mEdge DwPackage::makeDDFromMatrix(const std::vector<std::vector<Dw>> &matrix) {
    const std::size_t dim = std::size_t{1} << nqubits_;
    if (matrix.size() != dim)
        throw std::invalid_argument("makeDDFromMatrix: matrix.size() must equal 2^numQubits()");
    for (const auto &row : matrix) {
        if (row.size() != dim)
            throw std::invalid_argument("makeDDFromMatrix: matrix must be square with dimension 2^numQubits()");
    }
    return makeDDFromMatrixRec(matrix, 0, dim, 0, dim);
}

DwPackage::mEdge DwPackage::makeSingleQubitGateDD(std::size_t target, const std::array<Dw, 4> &matrix) {
    if (target >= nqubits_)
        throw std::invalid_argument("makeSingleQubitGateDD: target out of range");
    return makeControlledSingleQubitGateDD(std::vector<std::size_t>{}, target, matrix);
}

DwPackage::mEdge DwPackage::makeControlledSingleQubitGateDD(std::size_t control, std::size_t target,
                                                             const std::array<Dw, 4> &matrix) {
    if (control >= nqubits_ || target >= nqubits_ || control == target)
        throw std::invalid_argument("makeControlledSingleQubitGateDD: invalid control/target");
    return makeControlledSingleQubitGateDD(std::vector<std::size_t>{control}, target, matrix);
}

// ---------------------------------------------------------------------
// Bottom-up gate construction (multi-control, two-qubit) -- mirrors MQT
// Core's makeGateDD/makeTwoQubitGateDD. Unlike makeControlledRec above,
// no deferral state machine is needed: building starts from the gate's
// terminal matrix entries and walks the special qubits ascending, so by
// the time a control level is emitted its subtree is already known. See
// specs/multi-controlled-two-qubit-gate.md.
// ---------------------------------------------------------------------

std::vector<std::size_t> DwPackage::checkedControls(const std::vector<std::size_t> &controls, std::size_t target0,
                                                    std::size_t target1) const {
    std::vector<std::size_t> sorted(controls);
    std::sort(sorted.begin(), sorted.end());
    for (std::size_t i = 0; i < sorted.size(); ++i) {
        if (sorted[i] >= nqubits_)
            throw std::invalid_argument("controlled gate: control out of range");
        if (i > 0 && sorted[i] == sorted[i - 1])
            throw std::invalid_argument("controlled gate: duplicate control");
        if (sorted[i] == target0 || sorted[i] == target1)
            throw std::invalid_argument("controlled gate: control overlaps a target");
    }
    return sorted;
}

DwPackage::mEdge DwPackage::wrapWithControl(std::size_t ctrl, mEdge active, bool diagonal) {
    mEdge inactive = diagonal ? mEdge::terminal(Dw::one()) : mEdge::zero();
    return makeMEdge(static_cast<int>(ctrl), {std::move(inactive), mEdge::zero(), mEdge::zero(), std::move(active)});
}

DwPackage::mEdge DwPackage::makeControlledSingleQubitGateDD(const std::vector<std::size_t> &controls,
                                                             std::size_t target, const std::array<Dw, 4> &matrix) {
    if (target >= nqubits_)
        throw std::invalid_argument("makeControlledSingleQubitGateDD: target out of range");
    const auto sorted = checkedControls(controls, target, target);

    // The target's four block entries, starting as terminals (identity below
    // the target is implicit: a terminal edge is diagonal on all remaining
    // levels).
    std::array<mEdge, 4> em{mEdge::terminal(matrix[0]), mEdge::terminal(matrix[1]),
                            mEdge::terminal(matrix[2]), mEdge::terminal(matrix[3])};

    auto it = sorted.begin();
    for (; it != sorted.end() && *it < target; ++it) {
        for (std::size_t i = 0; i < 4; ++i)
            em[i] = wrapWithControl(*it, std::move(em[i]), i == 0 || i == 3);
    }
    mEdge e = makeMEdge(static_cast<int>(target), std::move(em));
    for (; it != sorted.end(); ++it)
        e = wrapWithControl(*it, std::move(e), true);
    return e;
}

DwPackage::mEdge DwPackage::makeControlledTwoQubitGateDD(const std::size_t control,
                                                          std::size_t target0, std::size_t target1,
                                                          const std::array<Dw, 16> &matrix) {
    if (control >= nqubits_ || target0 >= nqubits_ || target1 >= nqubits_ || control == target0 || control == target1 || target0 == target1)
        throw std::invalid_argument("makeControlledTwoQubitGateDD: invalid control/target0/target1");
    return makeControlledTwoQubitGateDD(std::vector<std::size_t>{control}, target0, target1, matrix);
}

DwPackage::mEdge DwPackage::makeControlledTwoQubitGateDD(const std::vector<std::size_t> &controls,
                                                          std::size_t target0, std::size_t target1,
                                                          const std::array<Dw, 16> &matrix) {
    if (target0 >= nqubits_ || target1 >= nqubits_)
        throw std::invalid_argument("makeControlledTwoQubitGateDD: target out of range");
    if (target0 == target1)
        throw std::invalid_argument("makeControlledTwoQubitGateDD: targets must be distinct");
    const auto sorted = checkedControls(controls, target0, target1);

    // 4x4 grid of terminal entries; row/col index = 2*bit(target0) +
    // bit(target1) (target0 is the MSB regardless of variable order).
    std::array<std::array<mEdge, 4>, 4> em{};
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t col = 0; col < 4; ++col)
            em[row][col] = mEdge::terminal(matrix[(row * 4) + col]);
    }

    const std::size_t lowTarget = std::min(target0, target1);
    const std::size_t highTarget = std::max(target0, target1);

    auto it = sorted.begin();
    for (; it != sorted.end() && *it < lowTarget; ++it) {
        for (std::size_t row = 0; row < 4; ++row) {
            for (std::size_t col = 0; col < 4; ++col)
                em[row][col] = wrapWithControl(*it, std::move(em[row][col]), row == col);
        }
    }

    // Collapse the lower target: one node per (row, col) bit pair of the
    // higher target. The higher target's bit is the matrix index's MSB when
    // it is target0 and its LSB otherwise.
    std::array<mEdge, 4> em0{};
    for (std::size_t row = 0; row < 2; ++row) {
        for (std::size_t col = 0; col < 2; ++col) {
            std::array<mEdge, 4> local{};
            for (std::size_t i = 0; i < 2; ++i) {
                for (std::size_t j = 0; j < 2; ++j) {
                    // (row, col, i, j) enumerates the 16 em entries bijectively,
                    // so each is read exactly once across the whole nest and can
                    // be moved from rather than copied.
                    local[(i * 2) + j] = target0 == highTarget ? std::move(em[(row * 2) + i][(col * 2) + j])
                                                                : std::move(em[(i * 2) + row][(j * 2) + col]);
                }
            }
            em0[(row * 2) + col] = makeMEdge(static_cast<int>(lowTarget), std::move(local));
        }
    }

    for (; it != sorted.end() && *it < highTarget; ++it) {
        for (std::size_t i = 0; i < 4; ++i)
            em0[i] = wrapWithControl(*it, std::move(em0[i]), i == 0 || i == 3);
    }

    mEdge e = makeMEdge(static_cast<int>(highTarget), std::move(em0));
    for (; it != sorted.end(); ++it)
        e = wrapWithControl(*it, std::move(e), true);
    return e;
}

DwPackage::mEdge DwPackage::makeTwoQubitGateDD(std::size_t target0, std::size_t target1,
                                                const std::array<Dw, 16> &matrix) {
    return makeControlledTwoQubitGateDD(std::vector<std::size_t>{}, target0, target1, matrix);
}

// ---------------------------------------------------------------------
// Multiply / add
// ---------------------------------------------------------------------

std::array<DwPackage::mEdge, 4> DwPackage::mChildrenAt(DwMNode *p, int var) {
    if (p != nullptr && p->var == var)
        return p->e;
    // Doesn't depend on `var`: identity-tensor expansion.
    return {mEdge{p, Dw::one()}, mEdge::zero(), mEdge::zero(), mEdge{p, Dw::one()}};
}

std::array<DwPackage::vEdge, 2> DwPackage::vChildrenAt(DwVNode *p, int var) {
    if (p != nullptr && p->var == var)
        return p->e;
    return {vEdge{p, Dw::one()}, vEdge{p, Dw::one()}};
}

DwPackage::vEdge DwPackage::multiplyRec(mEdge a, vEdge b) {
    if (a.w.isZero() || b.w.isZero())
        return vEdge::zero();
    const Dw scale = a.w * b.w;
    if (a.p == nullptr && b.p == nullptr)
        return {nullptr, scale};

    const std::pair<DwMNode *, DwVNode *> key{a.p, b.p};
    mvCacheStats_.trackLookup();
    if (auto it = mvCache_.find(key); it != mvCache_.end()) {
        mvCacheStats_.trackHit();
        return {it->second.p, scale * it->second.w};
    }

    const int var = std::max(topVar(a.p), topVar(b.p));
    const auto m = mChildrenAt(a.p, var);
    const auto v = vChildrenAt(b.p, var);
    const vEdge r0 = addRec(multiplyRec(m[0], v[0]), multiplyRec(m[1], v[1]));
    const vEdge r1 = addRec(multiplyRec(m[2], v[0]), multiplyRec(m[3], v[1]));
    const vEdge resultUnit = makeVEdge(var, {r0, r1});

    mvCache_.emplace(key, resultUnit);
    mvCacheStats_.trackInsert(mvCache_);
    return {resultUnit.p, scale * resultUnit.w};
}

DwPackage::vEdge DwPackage::multiply(const mEdge &m, const vEdge &v) {
    return multiplyRec(m, v);
}

DwPackage::mEdge DwPackage::multiplyRec(mEdge a, mEdge b) {
    if (a.w.isZero() || b.w.isZero())
        return mEdge::zero();
    const Dw scale = a.w * b.w;
    if (a.p == nullptr && b.p == nullptr)
        return {nullptr, scale};

    const std::pair<DwMNode *, DwMNode *> key{a.p, b.p};
    mmCacheStats_.trackLookup();
    if (auto it = mmCache_.find(key); it != mmCache_.end()) {
        mmCacheStats_.trackHit();
        return {it->second.p, scale * it->second.w};
    }

    const int var = std::max(topVar(a.p), topVar(b.p));
    const auto x = mChildrenAt(a.p, var); // 00,01,10,11
    const auto y = mChildrenAt(b.p, var);

    const mEdge r00 = addRec(multiplyRec(x[0], y[0]), multiplyRec(x[1], y[2]));
    const mEdge r01 = addRec(multiplyRec(x[0], y[1]), multiplyRec(x[1], y[3]));
    const mEdge r10 = addRec(multiplyRec(x[2], y[0]), multiplyRec(x[3], y[2]));
    const mEdge r11 = addRec(multiplyRec(x[2], y[1]), multiplyRec(x[3], y[3]));
    const mEdge resultUnit = makeMEdge(var, {r00, r01, r10, r11});

    mmCache_.emplace(key, resultUnit);
    mmCacheStats_.trackInsert(mmCache_);
    return {resultUnit.p, scale * resultUnit.w};
}

DwPackage::mEdge DwPackage::multiply(const mEdge &a, const mEdge &b) { return multiplyRec(a, b); }

DwPackage::vEdge DwPackage::addRec(vEdge a, vEdge b) {
    if (a.w.isZero())
        return b;
    if (b.w.isZero())
        return a;
    const int var = std::max(topVar(a.p), topVar(b.p));
    if (var < 0)
        return vEdge::terminal(a.w + b.w);

    const auto ac = vChildrenAt(a.p, var);
    const auto bc = vChildrenAt(b.p, var);
    vEdge c0 = addRec({ac[0].p, a.w * ac[0].w}, {bc[0].p, b.w * bc[0].w});
    vEdge c1 = addRec({ac[1].p, a.w * ac[1].w}, {bc[1].p, b.w * bc[1].w});
    return makeVEdge(var, {std::move(c0), std::move(c1)});
}

DwPackage::vEdge DwPackage::add(const vEdge &a, const vEdge &b) { return addRec(a, b); }

DwPackage::mEdge DwPackage::addRec(mEdge a, mEdge b) {
    if (a.w.isZero())
        return b;
    if (b.w.isZero())
        return a;
    const int var = std::max(topVar(a.p), topVar(b.p));
    if (var < 0)
        return mEdge::terminal(a.w + b.w);

    const auto ac = mChildrenAt(a.p, var);
    const auto bc = mChildrenAt(b.p, var);
    std::array<mEdge, 4> children{};
    for (std::size_t i = 0; i < 4; ++i)
        children[i] = addRec(mEdge{ac[i].p, a.w * ac[i].w}, mEdge{bc[i].p, b.w * bc[i].w});
    return makeMEdge(var, std::move(children));
}

// ---------------------------------------------------------------------
// Kronecker / outer product
// ---------------------------------------------------------------------

// x's own edge weight is applied exactly once, on return -- x's CHILDREN
// (x.p->e[i]) carry their own local weights and are recursed into
// directly, uncombined with x.w, exactly like mChildrenAt/vChildrenAt feed
// multiplyRec/addRec. y is threaded down unchanged and only actually
// combined (via its own weight y.w) at the point x bottoms out to a
// terminal, which is where its whole subtree gets grafted in. This is
// *not* a literal port of MQT Core's kronecker2: re-deriving that
// function's weight algebra by hand (multiplying xWeight*yWeight afresh
// at every recursion level, since MQT passes y unchanged at every depth
// too) shows it double-counts y's weight once per level of x's recursion
// depth. Verified against amplitudeRec's telescoping edge-weight
// semantics for a 1-qubit (X) 1-qubit example before settling on this
// version.
//
// `memo` is a map scoped to a single top-level kronecker() call (not a
// persistent DwPackage member like mvCache_/mmCache_): since x is
// hash-consed, it is generally a DAG, and without memoizing on x.p a
// shared x node would be re-walked and re-grafted once per root-to-leaf
// path that reaches it -- exponential in the worst case. The memo stores
// the "unit" result (ignoring the particular x.w that reached this node
// on any given call), rescaled by the caller's x.w on every use, exactly
// mirroring mvCache_'s rescale-on-hit convention.
DwPackage::vEdge DwPackage::kroneckerRec(vEdge x, vEdge y, std::size_t yNumQubits, bool incIdx,
                                          std::unordered_map<DwVNode *, vEdge> &memo) {
    if (x.w.isZero() || y.w.isZero())
        return vEdge::zero();
    if (x.p == nullptr)
        return {y.p, x.w * y.w};

    if (auto it = memo.find(x.p); it != memo.end())
        return {it->second.p, x.w * it->second.w};

    vEdge r0 = kroneckerRec(x.p->e[0], y, yNumQubits, incIdx, memo);
    vEdge r1 = kroneckerRec(x.p->e[1], y, yNumQubits, incIdx, memo);
    const int idx = x.p->var + (incIdx ? static_cast<int>(yNumQubits) : 0);
    const vEdge resultUnit = makeVEdge(idx, {std::move(r0), std::move(r1)});

    memo.emplace(x.p, resultUnit);
    return {resultUnit.p, x.w * resultUnit.w};
}

DwPackage::vEdge DwPackage::kronecker(const vEdge &x, const vEdge &y, std::size_t yNumQubits, bool incIdx) {
    const std::size_t shift = incIdx ? yNumQubits : 0U;
    if (static_cast<std::size_t>(topVar(x.p) + 1) + shift > nqubits_)
        throw std::invalid_argument("kronecker: combined qubit range exceeds numQubits()");
    std::unordered_map<DwVNode *, vEdge> memo;
    return kroneckerRec(x, y, yNumQubits, incIdx, memo);
}

DwPackage::mEdge DwPackage::kroneckerRec(mEdge x, mEdge y, std::size_t yNumQubits, bool incIdx,
                                          std::unordered_map<DwMNode *, mEdge> &memo) {
    if (x.w.isZero() || y.w.isZero())
        return mEdge::zero();
    if (x.p == nullptr)
        return {y.p, x.w * y.w};

    if (auto it = memo.find(x.p); it != memo.end())
        return {it->second.p, x.w * it->second.w};

    std::array<mEdge, 4> r{};
    for (std::size_t i = 0; i < 4; ++i)
        r[i] = kroneckerRec(x.p->e[i], y, yNumQubits, incIdx, memo);
    const int idx = x.p->var + (incIdx ? static_cast<int>(yNumQubits) : 0);
    const mEdge resultUnit = makeMEdge(idx, std::move(r));

    memo.emplace(x.p, resultUnit);
    return {resultUnit.p, x.w * resultUnit.w};
}

DwPackage::mEdge DwPackage::kronecker(const mEdge &x, const mEdge &y, std::size_t yNumQubits, bool incIdx) {
    const std::size_t shift = incIdx ? yNumQubits : 0U;
    if (static_cast<std::size_t>(topVar(x.p) + 1) + shift > nqubits_)
        throw std::invalid_argument("kronecker: combined qubit range exceeds numQubits()");
    std::unordered_map<DwMNode *, mEdge> memo;
    return kroneckerRec(x, y, yNumQubits, incIdx, memo);
}

// Unlike multiplyRec/addRec/kroneckerRec, this cannot stop recursing just
// because both x and y happen to have reduced to a terminal early: a
// terminal vEdge means "value independent of the remaining qubits", not
// "no remaining qubits" (same caveat documented on innerProductRec) --
// e.g. |+> reduces to a single terminal edge since its two branches are
// equal, but |+><+| is NOT diagonal, so collapsing outerProduct(|+>,|+>)
// to a scalar terminal the moment both operands look terminal would
// silently drop the off-diagonal cross terms. This walks all `level`
// levels explicitly, exactly like innerProductRec, rather than being
// driven by x.p/y.p's own (possibly-reduced) structure.
// scale = x.w * y.w.conjugate() is computed once per call, from the
// CURRENT call's own edge weights only, and reapplied exactly once on
// return -- xc[i]/yc[i] are passed to the recursive calls unweighted
// (their own local weight, not pre-multiplied by the parent), exactly
// mirroring multiplyRec's "cache the unit result, rescale on hit"
// convention. This is valid despite the conjugate because outer product
// is bilinear (linear in x, conjugate-linear in y): multiplying together
// each level's local conj(y.w) equals conj() of the accumulated product
// (conjugate is a multiplicative homomorphism), so factoring scale out
// per level instead of once at the leaf changes nothing mathematically --
// verified by hand against the |+>-outer-|+> example that originally
// motivated walking `level` instead of being driven by x.p/y.p directly.
//
// `memo` is local to one outerProduct() call (like kroneckerRec's fix),
// not a persistent DwPackage member. It's keyed on (x.p, y.p, level), not
// just (x.p, y.p): unlike multiplyRec's var (fully determined by its
// operands), level here is independent of x.p/y.p -- vChildrenAt's
// pass-through means the same node can legitimately recur at multiple
// different levels, so dropping level from the key would reuse a result
// computed for the wrong number of remaining qubits.
DwPackage::mEdge DwPackage::outerProductRec(vEdge x, vEdge y, std::size_t level,
                                             std::unordered_map<detail::VVKey, mEdge, detail::VVKeyHash> &memo) {
    if (x.w.isZero() || y.w.isZero())
        return mEdge::zero();
    const Dw scale = x.w * y.w.conjugate();
    if (level == 0)
        return mEdge::terminal(scale);

    const detail::VVKey key{x.p, y.p, level};
    if (auto it = memo.find(key); it != memo.end())
        return {it->second.p, scale * it->second.w};

    const int qubit = static_cast<int>(level - 1);
    const auto xc = vChildrenAt(x.p, qubit);
    const auto yc = vChildrenAt(y.p, qubit);
    mEdge m00 = outerProductRec(xc[0], yc[0], level - 1, memo);
    mEdge m01 = outerProductRec(xc[0], yc[1], level - 1, memo);
    mEdge m10 = outerProductRec(xc[1], yc[0], level - 1, memo);
    mEdge m11 = outerProductRec(xc[1], yc[1], level - 1, memo);
    const mEdge resultUnit = makeMEdge(qubit, {std::move(m00), std::move(m01), std::move(m10), std::move(m11)});

    memo.emplace(key, resultUnit);
    return {resultUnit.p, scale * resultUnit.w};
}

DwPackage::mEdge DwPackage::outerProduct(const vEdge &x, const vEdge &y) {
    std::unordered_map<detail::VVKey, mEdge, detail::VVKeyHash> memo;
    return outerProductRec(x, y, nqubits_, memo);
}

// ---------------------------------------------------------------------
// Measurement / fidelity
// ---------------------------------------------------------------------

DwPackage::MeasurementResult DwPackage::measureOneQubit(const vEdge &v, std::size_t qubit, bool outcome) {
    static const std::array<Dw, 4> projectorZero = {Dw::one(), Dw::zero(), Dw::zero(), Dw::zero()};
    static const std::array<Dw, 4> projectorOne = {Dw::zero(), Dw::zero(), Dw::zero(), Dw::one()};
    const auto &projector = outcome ? projectorOne : projectorZero;
    const auto projectorDD = makeSingleQubitGateDD(qubit, projector);
    vEdge projected = multiply(projectorDD, v);
    Dw probability = innerProduct(projected, projected);
    return {std::move(projected), std::move(probability)};
}

// fidelity(|a>,|b>) = |<a|b>|^2 / (<a|a><b|b>)
double DwPackage::fidelity(const vEdge &a, const vEdge &b) const {
    const Dw ip = innerProduct(a, b);
    const Dw num = ip.normSquared();
    const Dw na = innerProduct(a, a);
    const Dw nb = innerProduct(b, b);
    const double numD = num.toComplexDouble().real();
    const double denomD = (na * nb).toComplexDouble().real();
    return numD / denomD;
}

Dw DwPackage::amplitudeRec(vEdge e, const std::vector<bool> &bits, std::size_t level) const {
    if (level == 0)
        return e.w;
    const std::size_t qubit = level - 1;
    if (e.p != nullptr && e.p->var == static_cast<int>(qubit)) {
        const vEdge child = e.p->e[bits[qubit] ? 1 : 0];
        return e.w * amplitudeRec(child, bits, level - 1);
    }
    // Doesn't depend on this qubit: continue with the same node pointer
    // (or nullptr terminal) and unit weight for the remaining levels.
    return e.w * amplitudeRec(vEdge{e.p, Dw::one()}, bits, level - 1);
}

Dw DwPackage::amplitude(const vEdge &v, const std::vector<bool> &bits) const {
    if (bits.size() != nqubits_)
        throw std::invalid_argument("amplitude: bits.size() must equal numQubits()");
    return amplitudeRec(v, bits, nqubits_);
}

// Symmetric to amplitudeRec, but selects one of DwMNode::e's 4 children per
// level (row bit, col bit) rather than one of DwVNode::e's 2 (single bit).
// Unlike a vector node, a matrix node whose var doesn't match the current
// qubit represents an IMPLICIT IDENTITY BLOCK for that qubit (see
// mChildrenAt's identical "identity-tensor expansion"), not a value
// independent of it: the off-diagonal entries of that block are zero, only
// the diagonal (rowBits[qubit] == colBits[qubit]) carries e.p through.
Dw DwPackage::matrixEntryRec(mEdge e, const std::vector<bool> &rowBits, const std::vector<bool> &colBits,
                              std::size_t level) const {
    if (level == 0)
        return e.w;
    if (e.w.isZero())
        return Dw::zero();
    const std::size_t qubit = level - 1;
    if (e.p != nullptr && e.p->var == static_cast<int>(qubit)) {
        const std::size_t idx = (rowBits[qubit] ? 2U : 0U) + (colBits[qubit] ? 1U : 0U);
        const mEdge child = e.p->e[idx];
        return e.w * matrixEntryRec(child, rowBits, colBits, level - 1);
    }
    if (rowBits[qubit] != colBits[qubit])
        return Dw::zero();
    return e.w * matrixEntryRec(mEdge{e.p, Dw::one()}, rowBits, colBits, level - 1);
}

Dw DwPackage::matrixEntry(const mEdge &m, const std::vector<bool> &rowBits, const std::vector<bool> &colBits) const {
    if (rowBits.size() != nqubits_ || colBits.size() != nqubits_)
        throw std::invalid_argument("matrixEntry: rowBits/colBits.size() must equal numQubits()");
    return matrixEntryRec(m, rowBits, colBits, nqubits_);
}

// Unlike multiplyRec/addRec (which build a new, still-reduced DD and so
// don't need to know how many qubits remain -- amplitude() interprets that
// later), this aggregates a sum over all 2^level basis assignments. A
// terminal reached early (p == nullptr) only means "value independent of
// the remaining qubits", not "no remaining qubits" -- so, like
// amplitudeRec, this must walk all `level` levels explicitly rather than
// stopping as soon as both operands look terminal, or it would undercount
// by a factor of 2 per qubit that got reduced away this way (e.g. |+>,
// whose two branches collapse to the same terminal edge).
// `memo` is local to one innerProduct() call (like outerProductRec's and
// kroneckerRec's), not a persistent DwPackage member, which keeps this
// method const and makes GC invalidation a non-issue.
//
// The inner product is conjugate-bilinear:
//   <w_a*u | w_b*v> = conj(w_a) * w_b * <u|v>
// so the bracketed <u|v> depends only on (a.p, b.p, level) and every caller
// reaching that triple by a different path can share it.
Dw DwPackage::innerProductRec(vEdge a, vEdge b, std::size_t level,
                              std::unordered_map<detail::VVKey, Dw, detail::VVKeyHash> &memo) const {
    if (a.w.isZero() || b.w.isZero())
        return Dw::zero();
    const Dw scale = a.w.conjugate() * b.w;
    if (level == 0)
        return scale;

    const detail::VVKey key{a.p, b.p, level};
    if (auto it = memo.find(key); it != memo.end())
        return scale * it->second;

    const int qubit = static_cast<int>(level - 1);
    const auto ac = vChildrenAt(a.p, qubit);
    const auto bc = vChildrenAt(b.p, qubit);
    const Dw unit = innerProductRec(ac[0], bc[0], level - 1, memo) + innerProductRec(ac[1], bc[1], level - 1, memo);

    memo.emplace(key, unit);
    return scale * unit;
}

Dw DwPackage::innerProduct(const vEdge &a, const vEdge &b) const {
    std::unordered_map<detail::VVKey, Dw, detail::VVKeyHash> memo;
    return innerProductRec(a, b, nqubits_, memo);
}

// ---------------------------------------------------------------------
// Debug printing
// ---------------------------------------------------------------------

namespace {
/// bits[qubit] = (index >> qubit) & 1, i.e. bit `qubit` of `index`.
std::vector<bool> bitsFromIndex(std::size_t index, std::size_t nqubits) {
    std::vector<bool> bits(nqubits);
    for (std::size_t qubit = 0; qubit < nqubits; ++qubit)
        bits[qubit] = ((index >> qubit) & 1U) != 0;
    return bits;
}

/// `index`'s standard nqubits-bit binary representation, MSB (qubit
/// nqubits-1) first.
std::string binaryLabel(std::size_t index, std::size_t nqubits) {
    std::string label(nqubits, '0');
    for (std::size_t qubit = 0; qubit < nqubits; ++qubit) {
        if ((index >> qubit) & 1U)
            label[nqubits - 1 - qubit] = '1';
    }
    return label;
}

template <class Node> void assignNodeId(Node *p, std::unordered_map<Node *, std::size_t> &ids, std::vector<Node *> &order) {
    if (p == nullptr)
        return;
    if (ids.emplace(p, order.size()).second)
        order.push_back(p);
}

void appendIndent(std::ostringstream &os, std::size_t depth) { os << std::string(depth * 2, ' '); }

std::string formatEdge(const DwPackage::vEdge &e, const std::unordered_map<DwVNode *, std::size_t> &ids) {
    if (e.w.isZero())
        return "ZERO";
    if (e.p == nullptr)
        return "T(" + e.w.toString() + ")";
    std::ostringstream os;
    os << "@" << ids.at(e.p) << " * " << e.w.toString();
    return os.str();
}

std::string formatEdge(const DwPackage::mEdge &e, const std::unordered_map<DwMNode *, std::size_t> &ids) {
    if (e.w.isZero())
        return "ZERO";
    if (e.p == nullptr)
        return "T(" + e.w.toString() + ")";
    std::ostringstream os;
    os << "@" << ids.at(e.p) << " * " << e.w.toString();
    return os.str();
}

void appendVectorTree(std::ostringstream &os, const DwPackage::vEdge &e,
                      const std::unordered_map<DwVNode *, std::size_t> &ids,
                      std::unordered_set<DwVNode *> &expanded, std::size_t depth, const char *label) {
    appendIndent(os, depth);
    os << label << " -> " << formatEdge(e, ids);
    if (e.w.isZero() || e.p == nullptr) {
        os << "\n";
        return;
    }
    if (!expanded.emplace(e.p).second) {
        os << " [shared]\n";
        return;
    }
    os << "\n";
    appendIndent(os, depth + 1);
    os << "node @" << ids.at(e.p) << " var=" << e.p->var << "\n";
    if (!e.p->e[0].w.isZero())
        appendVectorTree(os, e.p->e[0], ids, expanded, depth + 2, "0");
    if (!e.p->e[1].w.isZero())
        appendVectorTree(os, e.p->e[1], ids, expanded, depth + 2, "1");
}

void appendMatrixTree(std::ostringstream &os, const DwPackage::mEdge &e,
                      const std::unordered_map<DwMNode *, std::size_t> &ids,
                      std::unordered_set<DwMNode *> &expanded, std::size_t depth, const char *label) {
    appendIndent(os, depth);
    os << label << " -> " << formatEdge(e, ids);
    if (e.w.isZero() || e.p == nullptr) {
        os << "\n";
        return;
    }
    if (!expanded.emplace(e.p).second) {
        os << " [shared]\n";
        return;
    }
    os << "\n";
    appendIndent(os, depth + 1);
    os << "node @" << ids.at(e.p) << " var=" << e.p->var << "\n";
    appendMatrixTree(os, e.p->e[0], ids, expanded, depth + 2, "00");
    appendMatrixTree(os, e.p->e[1], ids, expanded, depth + 2, "01");
    appendMatrixTree(os, e.p->e[2], ids, expanded, depth + 2, "10");
    appendMatrixTree(os, e.p->e[3], ids, expanded, depth + 2, "11");
}
} // namespace

std::string DwPackage::vectorToString(const vEdge &v) const {
    std::ostringstream os;
    const std::size_t dim = std::size_t{1} << nqubits_;
    for (std::size_t i = 0; i < dim; ++i) {
        const Dw amp = amplitude(v, bitsFromIndex(i, nqubits_));
        if (amp.isZero())
            continue;
        os << binaryLabel(i, nqubits_) << ": " << amp.toString() << "\n";
    }
    return os.str();
}

std::string DwPackage::matrixToString(const mEdge &m) const {
    std::ostringstream os;
    const std::size_t dim = std::size_t{1} << nqubits_;
    for (std::size_t row = 0; row < dim; ++row) {
        const auto rowBits = bitsFromIndex(row, nqubits_);
        for (std::size_t col = 0; col < dim; ++col) {
            if (col > 0)
                os << "  ";
            os << matrixEntry(m, rowBits, bitsFromIndex(col, nqubits_)).toString();
        }
        os << "\n";
    }
    return os.str();
}

std::string DwPackage::vectorDiagramToString(const vEdge &v) const {
    std::unordered_map<DwVNode *, std::size_t> ids;
    std::vector<DwVNode *> order;
    assignNodeId(v.p, ids, order);
    for (std::size_t i = 0; i < order.size(); ++i) {
        for (const auto &child : order[i]->e)
            assignNodeId(child.p, ids, order);
    }

    std::ostringstream os;
    std::unordered_set<DwVNode *> expanded;
    appendVectorTree(os, v, ids, expanded, 0, "root");
    return os.str();
}

std::string DwPackage::matrixDiagramToString(const mEdge &m) const {
    std::unordered_map<DwMNode *, std::size_t> ids;
    std::vector<DwMNode *> order;
    assignNodeId(m.p, ids, order);
    for (std::size_t i = 0; i < order.size(); ++i) {
        for (const auto &child : order[i]->e)
            assignNodeId(child.p, ids, order);
    }

    std::ostringstream os;
    std::unordered_set<DwMNode *> expanded;
    appendMatrixTree(os, m, ids, expanded, 0, "root");
    return os.str();
}

// ---------------------------------------------------------------------
// Reference counting / garbage collection
// ---------------------------------------------------------------------

namespace {
constexpr RefCount kMaxRefCount = std::numeric_limits<RefCount>::max();
} // namespace

void DwPackage::incRef(const vEdge &e) {
    if (e.p == nullptr || e.p->ref == kMaxRefCount)
        return;
    ++e.p->ref;
    if (e.p->ref == 1U) {
        // 0 -> 1: the node just became externally reachable, which tracks an active entry on
        // exactly this transition. Note that a node saturating at
        // kMaxRefCount above is never decremented again and so stays
        // counted as active -- garbageCollect()'s re-derivation bounds how
        // far that can drift (see UniqueTableStatistics::trackGcSweep).
        vUniqueStats_.trackActiveEntry();
        for (const auto &child : e.p->e)
            incRef(child);
    }
}

void DwPackage::incRef(const mEdge &e) {
    if (e.p == nullptr || e.p->ref == kMaxRefCount)
        return;
    ++e.p->ref;
    if (e.p->ref == 1U) {
        mUniqueStats_.trackActiveEntry();
        for (const auto &child : e.p->e)
            incRef(child);
    }
}

void DwPackage::decRef(const vEdge &e) {
    if (e.p == nullptr || e.p->ref == kMaxRefCount)
        return;
    assert(e.p->ref > 0 && "Dw decRef: unbalanced incRef/decRef (refcount already zero)");
    --e.p->ref;
    if (e.p->ref == 0U) {
        // 1 -> 0: the node is no longer an external root and is eligible
        // for the next sweep. The counterpart of the incRef tracking above.
        vUniqueStats_.untrackActiveEntry();
        for (const auto &child : e.p->e)
            decRef(child);
    }
}

void DwPackage::decRef(const mEdge &e) {
    if (e.p == nullptr || e.p->ref == kMaxRefCount)
        return;
    assert(e.p->ref > 0 && "Dw decRef: unbalanced incRef/decRef (refcount already zero)");
    --e.p->ref;
    if (e.p->ref == 0U) {
        mUniqueStats_.untrackActiveEntry();
        for (const auto &child : e.p->e)
            decRef(child);
    }
}

bool DwPackage::garbageCollect(bool force) {
    if (!force && vUnique_.size() < vGcLimit_ && mUnique_.size() < mGcLimit_)
        return false;

    std::size_t vFreed = 0;
    for (auto it = vUnique_.begin(); it != vUnique_.end();) {
        if (it->second->ref == 0U) {
            vMemory_.returnEntry(it->second);
            it = vUnique_.erase(it);
            ++vFreed;
        } else {
            ++it;
        }
    }
    // Placed after the early exit above, so gcRuns counts sweeps that
    // actually ran rather than every call.
    vUniqueStats_.trackGcSweep(vUnique_);

    std::size_t mFreed = 0;
    for (auto it = mUnique_.begin(); it != mUnique_.end();) {
        if (it->second->ref == 0U) {
            mMemory_.returnEntry(it->second);
            it = mUnique_.erase(it);
            ++mFreed;
        } else {
            ++it;
        }
    }
    mUniqueStats_.trackGcSweep(mUnique_);

    // A stale pointer to a freed node in either cache is unsafe regardless
    // of which table it came from, since mvCache_ is keyed on (mEdge*,
    // vEdge*) pairs -- mirrors MQT Core's rule that matrixVectorMultiplication
    // is cleared whenever either node type is collected.
    //
    // The corresponding compute-table statistics are deliberately NOT reset
    // here. MQT Core's ComputeTable::clear() calls stats.reset(), but that
    // only zeroes numEntries, which is snapshotted from the map anyway
    // (TableStatistics::snapshot), so the call would be a no-op; keeping
    // lookups/hits/inserts cumulative over the package's whole lifetime is
    // strictly more informative. The consequence, noted in
    // PackageStatistics, is that a compute table's numEntries reflects the
    // cache since the last collection while its other counters are lifetime
    // totals.
    if (vFreed > 0 || mFreed > 0)
        mvCache_.clear();
    if (mFreed > 0)
        mmCache_.clear();

    // Dynamic threshold adjustment, mirroring MQT Core: grow the limit if
    // still >90% full after the sweep, so a table that's genuinely large
    // doesn't trigger a full sweep on every subsequent allocation.
    if (vUnique_.size() > vGcLimit_ / 10 * 9)
        vGcLimit_ = vUnique_.size() + kInitialGcLimit;
    if (mUnique_.size() > mGcLimit_ / 10 * 9)
        mGcLimit_ = mUnique_.size() + kInitialGcLimit;

    return vFreed > 0 || mFreed > 0;
}

DwPackage::vEdge DwPackage::applyOperation(const mEdge &operation, const vEdge &e) {
    vEdge result = multiply(operation, e);
    incRef(result);
    decRef(e);
    garbageCollect();
    return result;
}

DwPackage::mEdge DwPackage::applyOperation(const mEdge &operation, const mEdge &e) {
    mEdge result = multiply(operation, e);
    incRef(result);
    decRef(e);
    garbageCollect();
    return result;
}

} // namespace dd::exact

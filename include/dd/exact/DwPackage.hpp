#ifndef DD_EXACT_DW_PACKAGE_HPP
#define DD_EXACT_DW_PACKAGE_HPP

#include "dd/exact/DwNode.hpp"
#include "dd/exact/MemoryManager.hpp"
#include "utility/HashUtil.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

namespace dd::exact {

namespace detail {
struct VKey {
    int var;
    std::array<DwVEdge, 2> children;
    [[nodiscard]] bool operator==(const VKey &o) const noexcept { return var == o.var && children == o.children; }
};
struct VKeyHash {
    [[nodiscard]] std::size_t operator()(const VKey &k) const noexcept;
};

struct MKey {
    int var;
    std::array<DwMEdge, 4> children;
    [[nodiscard]] bool operator==(const MKey &o) const noexcept { return var == o.var && children == o.children; }
};
struct MKeyHash {
    [[nodiscard]] std::size_t operator()(const MKey &k) const noexcept;
};

struct PtrPairHash {
    template <class A, class B> std::size_t operator()(const std::pair<A *, B *> &p) const noexcept {
        return HashUtil::combinedHash(std::hash<A *>{}(p.first), std::hash<B *>{}(p.second));
    }
};

/// Key for outerProductRec's per-call memo
struct VVKey {
    DwVNode *x;
    DwVNode *y;
    std::size_t level;
    [[nodiscard]] bool operator==(const VVKey &o) const noexcept {
        return x == o.x && y == o.y && level == o.level;
    }
};
struct VVKeyHash {
    [[nodiscard]] std::size_t operator()(const VVKey &k) const noexcept;
};
} // namespace detail

/**
 * Node-weight normalization strategy used by makeVEdge()/makeMEdge() to
 * factor a common value out of a node's outgoing edge weights (dividing
 * each by it) and propagate that factor to the incoming edge instead, so
 * that submatrices which differ only by a scalar factor collapse to one
 * shared node (see "Overcoming the Tradeoff Between Accuracy and
 * Compactness in Decision Diagrams for Quantum Computation", Niemann et
 * al., IEEE TCAD 2020, Section IV-B).
 *   - None:    no normalization (today's behavior; every DwPackage
 *              constructed without an explicit strategy keeps this).
 *   - Inverse: paper's Algorithm 2. The leftmost nonzero outgoing edge
 *              weight is used directly as the factor, using Dw::inverse()
 *              (exact division in the field Q[w]); the chosen weight
 *              becomes exactly Dw::one().
 *   - Gcd:     paper's Algorithm 3. The factor is the GCD (Dw::gcd(),
 *              canonicalized via Dw::reduceAssociate()) of all nonzero
 *              outgoing edge weights, staying within the ring D[w] (no
 *              field extension); unlike Inverse, no single weight is
 *              necessarily forced to exactly one. The paper's own
 *              evaluation found Inverse consistently achieves at least as
 *              much compaction as Gcd.
 */
enum class NormalizationStrategy { None, Inverse, Gcd };

/// Single-qubit basis states usable with the general makeBasisState()
/// overload: the computational basis (Zero/One) plus the Hadamard
/// (Plus/Minus) and circular (Right/Left) bases, mirroring MQT Core's
/// BasisStates enum (adapted to this repo's PascalCase enum convention).
enum class BasisState : std::uint8_t { Zero, One, Plus, Minus, Right, Left };

/**
 * A structurally parallel decision-diagram package to MQT Core's
 * dd::Package<Config>, weighted by the exact D[w] ring (dd::exact::Dw)
 * instead of double. Vector DDs represent quantum states, matrix DDs
 * represent (Clifford+T-only) operators, both over a fixed number of
 * qubits with a fixed variable order (var == qubit index, root = highest
 * qubit index, terminal == nullptr node).
 *
 * Nodes are hash-consed for exact structural equality (same var + same
 * children pointers/weights). With NormalizationStrategy::None (the
 * default), edge weights are not additionally factored out, so submatrices
 * that are scalar multiples of each other are not recognized as shareable.
 * Passing Inverse or Gcd to the constructor enables node-weight
 * normalization (see NormalizationStrategy), giving a more compact DD for
 * the same represented state/operator without losing exactness.
 *
 * Nodes are pooled via MemoryManager and reference-counted for garbage
 * collection, mirroring MQT Core's Package/UniqueTable/MemoryManager
 * design (node-level only -- unlike MQT Core, edge weights (Dw) are not
 * separately hash-consed/refcounted, since Dw is a heavyweight,
 * arbitrary-precision value rather than a simple double). Garbage
 * collection is manual/opt-in, exactly like MQT Core: nothing is ever
 * reclaimed unless a caller calls incRef()/decRef() to track root
 * reachability and (directly, or via applyOperation()) calls
 * garbageCollect(). A DwPackage that is never asked to collect behaves
 * exactly as before -- nodes simply accumulate for the package's lifetime.
 */
class DwPackage {
public:
    using vEdge = DwVEdge;
    using mEdge = DwMEdge;

    explicit DwPackage(std::size_t nqubits, NormalizationStrategy strategy = NormalizationStrategy::Inverse);

    // Node pointers stored in vUnique_/mUnique_ point into vMemory_/mMemory_'s
    // pooled storage; copying a DwPackage would leave the copy's unique
    // tables pointing at the original's nodes. Moving is safe (MemoryManager
    // only ever moves, never copies, its chunk storage).
    DwPackage(const DwPackage &) = delete;
    DwPackage &operator=(const DwPackage &) = delete;
    DwPackage(DwPackage &&) = default;
    DwPackage &operator=(DwPackage &&) = default;

    [[nodiscard]] std::size_t numQubits() const { return nqubits_; }

    /// Basis state |bits[0] bits[1] ... bits[nqubits_-1]>, built bottom-up:
    /// one node per qubit, glued ascending from qubit 0 to nqubits_-1.
    [[nodiscard]] vEdge makeBasisState(const std::vector<bool> &bits);

    /// Basis state over only qubits [start, start+n) (n <= numQubits() -
    /// start), built the same way as makeBasisState(bits) but leaving
    /// qubits outside that window untouched by the resulting edge -- reads
    /// at those qubits fall through vChildrenAt's existing "doesn't depend
    /// on this var" pass-through (weight-1 on both branches), the same
    /// mechanism that already makes narrow mEdge gate DDs work, applied
    /// here to vEdges. Lets a genuinely narrower-than-numQubits() state be
    /// built directly (e.g. for use as kronecker()'s y operand). Mirrors
    /// MQT Core's Package::makeBasisState(n, state, start).
    [[nodiscard]] vEdge makeBasisState(const std::size_t n, const std::vector<bool> &state,
                                        const std::size_t start = 0);

    /// General basis state over qubits [start, start+n) from a per-qubit
    /// BasisState list (computational, Hadamard, and circular bases).
    /// Mirrors MQT Core's Package::makeBasisState(n, vector<BasisStates>,
    /// start) -- but does NOT port its loop body literally: MQT's version
    /// scales each Plus/Minus/Right/Left child by a bare constant rather
    /// than the running edge's own accumulated weight, which (by hand
    /// tracing, see specs/ or DwPackage.cpp) looks like it would
    /// under-scale chains of more than one such qubit (see DwPackage.cpp).
    /// This overload
    /// instead threads the running edge's own weight into each new child
    /// explicitly, the same way makeBasisState(bits) already does.
    [[nodiscard]] vEdge makeBasisState(const std::size_t n, const std::vector<BasisState> &state,
                                        const std::size_t start = 0);

    [[nodiscard]] vEdge makeZeroState();
    [[nodiscard]] mEdge makeIdentity();

    /// Builds a random single-qubit state a|0> + b|1> at `qubit`, with a,
    /// b in Dw satisfying |a|^2 + |b|^2 = 1 exactly: applies `depth`
    /// randomly-chosen single-qubit Clifford+T gates to |0>, each exactly
    /// unitary in D[w], so the result is exactly normalized by
    /// construction -- this sidesteps the much harder problem of directly
    /// sampling a valid (a,b) pair (an instance of the "exact state
    /// synthesis" problem for Clifford+T, which in general needs
    /// Gaussian-integer-style factorization). DwPackage owns no RNG
    /// itself; the caller supplies and seeds `rng`.
    [[nodiscard]] vEdge makeRandomSingleQubitState(std::mt19937 &rng, std::size_t depth = 10, std::size_t qubit = 0);

    /// Builds the vEdge for an explicit 2^numQubits()-length amplitude
    /// vector, top-down: recursively splits the vector in half at each
    /// level (index bit `level` is the most significant bit of the range
    /// covered at that recursion depth, so "first half -> 0-child, second
    /// half -> 1-child" matches amplitude()'s bit convention directly).
    /// Mirrors MQT Core's Package::makeStateFromVector, adapted to Dw
    /// weights (no complex-to-exact conversion is performed -- amplitudes
    /// must already be exact Dw values).
    [[nodiscard]] vEdge makeStateFromVector(const std::vector<Dw> &amplitudes);

    /// Builds the mEdge for an explicit 2^numQubits() x 2^numQubits() dense
    /// matrix, top-down via quadrant splitting (row/col bit `level` is the
    /// most significant bit of the row/col range at that depth). Mirrors
    /// MQT Core's Package::makeDDFromMatrix, adapted to Dw weights.
    [[nodiscard]] mEdge makeDDFromMatrix(const std::vector<std::vector<Dw>> &matrix);

    /// Embeds a single-qubit gate matrix {m00,m01,m10,m11} acting on
    /// `target` into an nqubits-wide operator DD (identity elsewhere).
    [[nodiscard]] mEdge makeSingleQubitGateDD(std::size_t target, const std::array<Dw, 4> &matrix);

    /// Embeds a single-target gate controlled by `control` (control !=
    /// target, either may have the higher qubit index).
    [[nodiscard]] mEdge makeControlledSingleQubitGateDD(std::size_t control, std::size_t target,
                                                         const std::array<Dw, 4> &matrix);

    /// Embeds a single-target gate with any number of positive controls
    /// (fire on |1>; an empty vector gives the plain single-qubit gate).
    /// Controls may lie above and/or below `target` in the variable order.
    /// Built bottom-up like MQT Core's makeGateDD, in contrast to the
    /// top-down single-control builder above (see
    /// specs/multi-controlled-two-qubit-gate.md).
    [[nodiscard]] mEdge makeControlledSingleQubitGateDD(const std::vector<std::size_t> &controls,
                                                         std::size_t target, const std::array<Dw, 4> &matrix);

    /// Embeds a two-qubit gate given as a row-major 4x4 matrix whose row/
    /// column index is 2*bit(target0) + bit(target1) -- target0 is the more
    /// significant bit regardless of which target has the higher qubit
    /// index, mirroring MQT Core's TwoQubitGateMatrix convention.
    [[nodiscard]] mEdge makeTwoQubitGateDD(std::size_t target0, std::size_t target1,
                                            const std::array<Dw, 16> &matrix);

    /// Single-control convenience overload; forwards to the vector-controls
    /// overload below with {control}. CAUTION: with both overloads in scope,
    /// a bare `{}` argument resolves to *this* overload (an empty-list-to-
    /// scalar conversion ranks as identity, beating the vector overload's
    /// initializer_list constructor), value-initializing control to 0 --
    /// NOT "no control". Callers wanting an uncontrolled gate must call
    /// makeTwoQubitGateDD(), or spell std::vector<std::size_t>{} explicitly.
    [[nodiscard]] mEdge makeControlledTwoQubitGateDD(const std::size_t control,
                                                      std::size_t target0, std::size_t target1,
                                                      const std::array<Dw, 16> &matrix);

    /// Two-qubit gate with any number of positive controls; controls may
    /// lie below, between, or above the two targets. See the scalar overload
    /// above for a brace-init hazard when calling this with zero controls.
    [[nodiscard]] mEdge makeControlledTwoQubitGateDD(const std::vector<std::size_t> &controls,
                                                      std::size_t target0, std::size_t target1,
                                                      const std::array<Dw, 16> &matrix);

    [[nodiscard]] vEdge multiply(const mEdge &m, const vEdge &v);
    [[nodiscard]] mEdge multiply(const mEdge &a, const mEdge &b);
    [[nodiscard]] vEdge add(const vEdge &a, const vEdge &b);

    /// Exact amplitude of computational basis index `bits` (bits[i] for
    /// qubit i) under state `v`.
    [[nodiscard]] Dw amplitude(const vEdge &v, const std::vector<bool> &bits) const;

    /// Exact entry M[row][col] of operator `m`, where rowBits/colBits give
    /// each qubit's row/column bit (rowBits[i]/colBits[i] for qubit i).
    /// Symmetric to amplitude(), selecting one of DwMNode::e's 4 children
    /// per level instead of one of DwVNode::e's 2.
    [[nodiscard]] Dw matrixEntry(const mEdge &m, const std::vector<bool> &rowBits,
                                  const std::vector<bool> &colBits) const;

    /// Exact inner product <a|b> = sum over basis states of conj(amp_a) *
    /// amp_b, computed via DD recursion (no exponential basis-state
    /// enumeration). normSquared(v) == innerProduct(v, v).
    [[nodiscard]] Dw innerProduct(const vEdge &a, const vEdge &b) const;

    /// Tensor/Kronecker product x (X) y within this same package: y keeps
    /// its own variable numbering (occupying the low qubit range), x's
    /// node variables are shifted up by yNumQubits (occupying the high
    /// range) unless incIdx is false, in which case x's variables are left
    /// untouched (the caller is responsible for x/y already occupying
    /// disjoint ranges). x and y must belong to this same DwPackage --
    /// there is no cross-package/extended-qubit-count support. Throws
    /// std::invalid_argument if the shifted result would exceed
    /// numQubits(). Weight-threading is derived from scratch (see
    /// DwPackage.cpp) to avoid double-counting y's weight once per level
    /// of x's recursion depth.
    [[nodiscard]] vEdge kronecker(const vEdge &x, const vEdge &y, std::size_t yNumQubits, bool incIdx = true);
    [[nodiscard]] mEdge kronecker(const mEdge &x, const mEdge &y, std::size_t yNumQubits, bool incIdx = true);

    /// |x><y|, computed directly via DD recursion for arbitrary-width x/y.
    /// y's weight is conjugated, matching innerProduct's bra convention.
    [[nodiscard]] mEdge outerProduct(const vEdge &x, const vEdge &y);

    struct MeasurementResult {
        vEdge state;    ///< projected, unnormalized post-measurement state
        Dw probability; ///< exact probability of this outcome
    };

    /// Projects v onto qubit == outcome and returns the (unnormalized)
    /// projected state plus its exact probability -- dividing by
    /// sqrt(probability) is generally not exact in D[w], so (unlike MQT
    /// Core's measureOneQubit) this deliberately does not renormalize,
    /// matching ExactDDSimulation::measure's existing convention.
    [[nodiscard]] MeasurementResult measureOneQubit(const vEdge &v, std::size_t qubit, bool outcome);

    /// |<a|b>|^2 / (<a|a> * <b|b>), computed as a ratio of exact Dw
    /// quantities and converted to double only for the final division.
    [[nodiscard]] double fidelity(const vEdge &a, const vEdge &b) const;

    /// Debug helpers: render the full 2^nqubits_ state vector / 2^nqubits_ x
    /// 2^nqubits_ operator matrix that `v`/`m` represents, using each
    /// entry's exact Dw::toString() form. vectorToString() skips zero
    /// entries (one "<n-bit binary index>: <value>" line per nonzero
    /// amplitude); matrixToString() prints every entry as a dense,
    /// space-separated, row-major grid. Intended for small circuits only
    /// (output is exponential in nqubits_), mirroring MQT Core's
    /// Edge::printVector()/printMatrix().
    [[nodiscard]] std::string vectorToString(const vEdge &v) const;
    [[nodiscard]] std::string matrixToString(const mEdge &m) const;
    void printVector(const vEdge &v, std::ostream &os = std::cout) const { os << vectorToString(v); }
    void printMatrix(const mEdge &m, std::ostream &os = std::cout) const { os << matrixToString(m); }

    /// Render the stored reduced DD itself rather than its dense semantic
    /// vector/matrix as an indented tree: each edge shows either ZERO,
    /// T(weight), or @id * weight, and shared nonterminal nodes are only
    /// expanded once (later references are marked shared).
    [[nodiscard]] std::string vectorDiagramToString(const vEdge &v) const;
    [[nodiscard]] std::string matrixDiagramToString(const mEdge &m) const;
    void printVectorDiagram(const vEdge &v, std::ostream &os = std::cout) const { os << vectorDiagramToString(v); }
    void printMatrixDiagram(const mEdge &m, std::ostream &os = std::cout) const { os << matrixDiagramToString(m); }

    [[nodiscard]] std::size_t vNodeCount() const { return vUnique_.size(); }
    [[nodiscard]] std::size_t mNodeCount() const { return mUnique_.size(); }

    /// Increments the reference count of e.p (and, only the first time it
    /// becomes reachable, recursively of its children), marking it and
    /// everything below it as an external root that garbageCollect() must
    /// not reclaim. Mirrors MQT Core's Package::incRef -- but, per this
    /// module's scope, only refcounts nodes, not edge weights.
    void incRef(const vEdge &e);
    void incRef(const mEdge &e);

    /// Decrements the reference count of e.p (and, only when it drops to
    /// zero, recursively of its children). Asserts against decrementing an
    /// already-zero refcount (an unbalanced incRef/decRef call).
    void decRef(const vEdge &e);
    void decRef(const mEdge &e);

    /// Mark-and-sweep garbage collection: unless `force` is set, this is a
    /// no-op until either node table has grown past its current threshold.
    /// Otherwise, erases every hash-consed node with ref == 0 from its
    /// unique table, returns it to its MemoryManager, clears whichever
    /// compute-table caches could reference a freed node, and grows the
    /// relevant threshold(s) if still nearly full afterward. Mirrors MQT
    /// Core's Package::garbageCollect.
    bool garbageCollect(bool force = false);

    /// Convenience method mirroring MQT Core's Package::applyOperation:
    /// multiplies `operation` into `e`, incRef()s the result, decRef()s
    /// `e`, calls garbageCollect(), and returns the result. Intended for
    /// callers holding a persistent root (e.g. a simulator's current state)
    /// that they want replaced with the product.
    [[nodiscard]] vEdge applyOperation(const mEdge &operation, const vEdge &e);
    [[nodiscard]] mEdge applyOperation(const mEdge &operation, const mEdge &e);

private:
    std::size_t nqubits_;
    NormalizationStrategy strategy_;

    [[nodiscard]] vEdge makeVEdge(int var, std::array<vEdge, 2> children);
    [[nodiscard]] mEdge makeMEdge(int var, std::array<mEdge, 4> children);

    /// Validates a control set against nqubits_ and the given target(s)
    /// (pass the same qubit twice for a single-target gate) and returns it
    /// sorted ascending; throws std::invalid_argument on an out-of-range or
    /// duplicate control, or a control that overlaps a target.
    [[nodiscard]] std::vector<std::size_t> checkedControls(const std::vector<std::size_t> &controls,
                                                           std::size_t target0, std::size_t target1) const;

    /// Bottom-up building block shared by the multi-controlled builders:
    /// wraps `active` in a positive-control node at `ctrl`, whose control=0
    /// child is the corresponding entry of the identity operator --
    /// T(Dw::one()) if the entry being wrapped sits on the block diagonal
    /// (`diagonal`), zero otherwise.
    [[nodiscard]] mEdge wrapWithControl(std::size_t ctrl, const mEdge &active, bool diagonal);

    [[nodiscard]] static std::array<mEdge, 4> mChildrenAt(DwMNode *p, int var);
    [[nodiscard]] static std::array<vEdge, 2> vChildrenAt(DwVNode *p, int var);

    [[nodiscard]] vEdge makeStateFromVectorRec(std::vector<Dw>::const_iterator begin,
                                                std::vector<Dw>::const_iterator end);
    [[nodiscard]] mEdge makeDDFromMatrixRec(const std::vector<std::vector<Dw>> &matrix, std::size_t rowStart,
                                             std::size_t rowEnd, std::size_t colStart, std::size_t colEnd);

    [[nodiscard]] vEdge multiplyRec(mEdge a, vEdge b);
    [[nodiscard]] mEdge multiplyRec(mEdge a, mEdge b);
    [[nodiscard]] vEdge addRec(vEdge a, vEdge b);
    [[nodiscard]] mEdge addRec(mEdge a, mEdge b);

    [[nodiscard]] vEdge kroneckerRec(vEdge x, vEdge y, std::size_t yNumQubits, bool incIdx,
                                      std::unordered_map<DwVNode *, vEdge> &memo);
    [[nodiscard]] mEdge kroneckerRec(mEdge x, mEdge y, std::size_t yNumQubits, bool incIdx,
                                      std::unordered_map<DwMNode *, mEdge> &memo);
    [[nodiscard]] mEdge outerProductRec(vEdge x, vEdge y, std::size_t level,
                                         std::unordered_map<detail::VVKey, mEdge, detail::VVKeyHash> &memo);

    [[nodiscard]] Dw amplitudeRec(vEdge e, const std::vector<bool> &bits, std::size_t level) const;
    [[nodiscard]] Dw matrixEntryRec(mEdge e, const std::vector<bool> &rowBits, const std::vector<bool> &colBits,
                                     std::size_t level) const;
    [[nodiscard]] Dw innerProductRec(vEdge a, vEdge b, std::size_t level) const;

    // pooled node storage (MemoryManager) + unique tables (hash-consing)
    MemoryManager<DwVNode> vMemory_;
    MemoryManager<DwMNode> mMemory_;
    std::unordered_map<detail::VKey, DwVNode *, detail::VKeyHash> vUnique_;
    std::unordered_map<detail::MKey, DwMNode *, detail::MKeyHash> mUnique_;

    // GC trigger thresholds (grown after a sweep if still nearly full,
    // mirroring MQT Core's dynamic gcLimit adjustment)
    static constexpr std::size_t kInitialGcLimit = 4096U;
    std::size_t vGcLimit_{kInitialGcLimit};
    std::size_t mGcLimit_{kInitialGcLimit};

    // compute tables (memoization), keyed on unit-weight node pointer pairs
    std::unordered_map<std::pair<DwMNode *, DwVNode *>, vEdge, detail::PtrPairHash> mvCache_;
    std::unordered_map<std::pair<DwMNode *, DwMNode *>, mEdge, detail::PtrPairHash> mmCache_;
};

} // namespace dd::exact

#endif // DD_EXACT_DW_PACKAGE_HPP

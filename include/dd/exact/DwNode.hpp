#ifndef DD_EXACT_DW_NODE_HPP
#define DD_EXACT_DW_NODE_HPP

#include "dd/exact/Dw.hpp"

#include <array>
#include <cstdint>

namespace dd::exact {

/// Saturating reference count used by DwPackage's garbage collector: once a
/// node's ref reaches the maximum value it is treated as permanently
/// reachable (never incremented/decremented further), matching MQT Core's
/// dd::RefCount semantics.
using RefCount = std::uint32_t;

/// An edge to a DD node, weighted by an exact D[w] value. The terminal is
/// represented by p == nullptr; a terminal edge's weight is the scalar
/// value for all remaining (already-decided) qubits.
template <class Node> struct DwEdge {
    Node *p = nullptr;
    Dw w = Dw::zero();

    [[nodiscard]] bool operator==(const DwEdge &other) const noexcept {
        return p == other.p && w == other.w;
    }
    [[nodiscard]] bool operator!=(const DwEdge &other) const noexcept { return !(*this == other); }

    [[nodiscard]] static DwEdge zero() { return {}; }
    [[nodiscard]] static DwEdge terminal(Dw weight) { return {nullptr, std::move(weight)}; }
};

/// Vector (state) DD node: qubit `var`, two children indexed by that
/// qubit's basis value (e[0] = |0>-branch, e[1] = |1>-branch).
/// `ref`/`next` support DwPackage's MemoryManager-based pooling and
/// reference-counted garbage collection; they play no part in DwEdge
/// equality/hashing (which only ever compares node pointers and weights).
struct DwVNode {
    int var = -1;
    std::array<DwEdge<DwVNode>, 2> e{};
    RefCount ref{0};
    DwVNode *next{nullptr};
};

/// Matrix (operator) DD node: qubit `var`, four children -- e[0]=M00,
/// e[1]=M01, e[2]=M10, e[3]=M11, indexed by (row bit, col bit) of `var`.
/// See DwVNode's doc comment regarding `ref`/`next`.
struct DwMNode {
    int var = -1;
    std::array<DwEdge<DwMNode>, 4> e{};
    RefCount ref{0};
    DwMNode *next{nullptr};
};

using DwVEdge = DwEdge<DwVNode>;
using DwMEdge = DwEdge<DwMNode>;

} // namespace dd::exact

#endif // DD_EXACT_DW_NODE_HPP

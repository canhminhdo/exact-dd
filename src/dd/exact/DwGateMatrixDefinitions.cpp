#include "dd/exact/DwGateMatrixDefinitions.hpp"
#include "utility/StringHelper.hpp"

#include <stdexcept>
#include <unordered_map>

namespace dd::exact::gates {

// Each gate matrix is a function-local static built once and handed back by
// reference. Returning by value rebuilt it on every gate application -- four
// (or sixteen) Dw copy constructions, each up to five mpz_init_set -- and for
// the matrices whose entries are derived rather than named (y, z, h, tdg, v,
// vdg, sx, sxdg) it also redid the negation, conjugation or multiplication
// every time. The values are constants; only their construction was repeated.
//
// The namespace-scope constants below are initialised before main, so they are
// always ready by the time a function-local static first runs.
namespace {
const Dw invSqrt2(0, 1, 0, -1, 2); // 1/sqrt(2)
const Dw invSqrt2PlusI{1, 0, 1, 0, 2}; // (1 + i)/2
const Dw k0 = Dw::zero();
const Dw k1 = Dw::one();
const Dw omega = Dw::omega();
const Dw iUnit = omega * omega;
} // namespace

const std::array<Dw, 4> &i() {
    static const std::array<Dw, 4> m{k1, k0, k0, k1};
    return m;
}

const std::array<Dw, 4> &x() {
    static const std::array<Dw, 4> m{k0, k1, k1, k0};
    return m;
}

const std::array<Dw, 4> &y() {
    // i == w^2.
    static const std::array<Dw, 4> m{k0, -iUnit, iUnit, k0};
    return m;
}

const std::array<Dw, 4> &z() {
    static const std::array<Dw, 4> m{k1, k0, k0, -k1};
    return m;
}

const std::array<Dw, 4> &h() {
    static const std::array<Dw, 4> m{invSqrt2, invSqrt2, invSqrt2, -invSqrt2};
    return m;
}

const std::array<Dw, 4> &s() {
    static const std::array<Dw, 4> m{k1, k0, k0, iUnit};
    return m;
}

const std::array<Dw, 4> &sdg() {
    static const std::array<Dw, 4> m{k1, k0, k0, -iUnit};
    return m;
}

const std::array<Dw, 4> &t() {
    static const std::array<Dw, 4> m{k1, k0, k0, omega};
    return m;
}

const std::array<Dw, 4> &tdg() {
    // w^-1 == conjugate(w) since |w| == 1.
    static const std::array<Dw, 4> m{k1, k0, k0, omega.conjugate()};
    return m;
}

const std::array<Dw, 4> &v() {
    static const Dw offDiag = -iUnit * invSqrt2;
    static const std::array<Dw, 4> m{invSqrt2, offDiag, offDiag, invSqrt2};
    return m;
}

const std::array<Dw, 4> &vdg() {
    static const Dw offDiag = iUnit * invSqrt2;
    static const std::array<Dw, 4> m{invSqrt2, offDiag, offDiag, invSqrt2};
    return m;
}

const std::array<Dw, 4> &sx() {
    static const std::array<Dw, 4> m{invSqrt2PlusI, invSqrt2PlusI.conjugate(), invSqrt2PlusI.conjugate(), invSqrt2PlusI};
    return m;
}

const std::array<Dw, 4> &sxdg() {
    static const std::array<Dw, 4> m{invSqrt2PlusI.conjugate(), invSqrt2PlusI, invSqrt2PlusI, invSqrt2PlusI.conjugate()};
    return m;
}

const std::array<Dw, 16> &swap() {
    static const std::array<Dw, 16> m{k1, k0, k0, k0,
            k0, k0, k1, k0,
            k0, k1, k0, k0,
            k0, k0, k0, k1};
    return m;
}

const std::array<Dw, 16> &iswap() {
    static const std::array<Dw, 16> m{k1, k0, k0,    k0,
            k0, k0, iUnit, k0,
            k0, iUnit, k0, k0,
            k0, k0, k0,    k1};
    return m;
}

const std::array<Dw, 16> &iswapdg() {
    static const std::array<Dw, 16> m{k1, k0, k0,     k0,
            k0, k0, -iUnit, k0,
            k0, -iUnit, k0, k0,
            k0, k0, k0,     k1};
    return m;
}

const std::array<Dw, 16> &dcx() {
    static const std::array<Dw, 16> m{k1, k0, k0, k0,
            k0, k0, k1, k0,
            k0, k0, k0, k1,
            k0, k1, k0, k0};
    return m;
}

const std::array<Dw, 4> &byName(const std::string &name) {
    static const std::unordered_map<std::string, const std::array<Dw, 4> &(*)()> kCliffordT{
        {"i", i}, {"x", x}, {"y", y}, {"z", z}, {"h", h}, {"s", s}, {"sdg", sdg}, {"t", t}, {"tdg", tdg},
        {"v", v}, {"vdg", vdg}, {"sx", sx}, {"sxdg", sxdg}
    };
    if (auto it = kCliffordT.find(StringHelper::toLowerCase(name)); it != kCliffordT.end())
        return it->second();
    throw std::invalid_argument("gate '" + name +
                                 "' is not exactly representable in D[w] (only Clifford+T "
                                 "gates i,x,y,z,h,s,sdg,t,tdg,v,vdg,sx,sxdg are supported in exact mode)");
}

const std::array<Dw, 16> &twoQubitByName(const std::string &name) {
    static const std::unordered_map<std::string, const std::array<Dw, 16> &(*)()> kTwoQubit{
        {"swap", swap}, {"iswap", iswap}, {"iswapdg", iswapdg}, {"dcx", dcx}
    };
    if (auto it = kTwoQubit.find(StringHelper::toLowerCase(name)); it != kTwoQubit.end())
        return it->second();
    throw std::invalid_argument("two-qubit gate '" + name +
                                 "' is not exactly representable in D[w] (only swap,iswap,iswapdg,dcx "
                                 "are supported in exact mode; build cx/cz as controlled gates)");
}

} // namespace dd::exact::gates

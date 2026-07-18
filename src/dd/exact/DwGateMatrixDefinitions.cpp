#include "dd/exact/DwGateMatrixDefinitions.hpp"
#include "utility/StringHelper.hpp"

#include <stdexcept>
#include <unordered_map>

namespace dd::exact::gates {

namespace {
const Dw invSqrt2(0, 1, 0, -1, 2); // 1/sqrt(2)
const Dw invSqrt2PlusI{1, 0, 1, 0, 2}; // (1 + i)/2
const Dw k0 = Dw::zero();
const Dw k1 = Dw::one();
const Dw omega = Dw::omega();
const Dw iUnit = omega * omega;
} // namespace

std::array<Dw, 4> i() { return {k1, k0, k0, k1}; }

std::array<Dw, 4> x() { return {k0, k1, k1, k0}; }

std::array<Dw, 4> y() {
    // i == w^2.
    return {k0, -iUnit, iUnit, k0};
}

std::array<Dw, 4> z() { return {k1, k0, k0, -k1}; }

std::array<Dw, 4> h() { return {invSqrt2, invSqrt2, invSqrt2, -invSqrt2}; }

std::array<Dw, 4> s() {
    return {k1, k0, k0, iUnit};
}

std::array<Dw, 4> sdg() {
    return {k1, k0, k0, -iUnit};
}

std::array<Dw, 4> t() { return {k1, k0, k0, omega}; }

std::array<Dw, 4> tdg() {
    // w^-1 == conjugate(w) since |w| == 1.
    return {k1, k0, k0, omega.conjugate()};
}

std::array<Dw, 4> v() {
    const Dw offDiag = -iUnit * invSqrt2;
    return {invSqrt2, offDiag, offDiag, invSqrt2};
}

std::array<Dw, 4> vdg() {
    const Dw offDiag = iUnit * invSqrt2;
    return {invSqrt2, offDiag, offDiag, invSqrt2};
}

std::array<Dw, 4> sx() {
    return {invSqrt2PlusI, invSqrt2PlusI.conjugate(), invSqrt2PlusI.conjugate(), invSqrt2PlusI};
}

std::array<Dw, 4> sxdg() {
    return {invSqrt2PlusI.conjugate(), invSqrt2PlusI, invSqrt2PlusI, invSqrt2PlusI.conjugate()};
}

std::array<Dw, 16> swap() {
    return {k1, k0, k0, k0,
            k0, k0, k1, k0,
            k0, k1, k0, k0,
            k0, k0, k0, k1};
}

std::array<Dw, 16> iswap() {
    return {k1, k0, k0,    k0,
            k0, k0, iUnit, k0,
            k0, iUnit, k0, k0,
            k0, k0, k0,    k1};
}

std::array<Dw, 16> iswapdg() {
    return {k1, k0, k0,     k0,
            k0, k0, -iUnit, k0,
            k0, -iUnit, k0, k0,
            k0, k0, k0,     k1};
}

std::array<Dw, 16> dcx() {
    return {k1, k0, k0, k0,
            k0, k0, k1, k0,
            k0, k0, k0, k1,
            k0, k1, k0, k0};
}

std::array<Dw, 4> byName(const std::string &name) {
    static const std::unordered_map<std::string, std::array<Dw, 4> (*)()> kCliffordT{
        {"i", i}, {"x", x}, {"y", y}, {"z", z}, {"h", h}, {"s", s}, {"sdg", sdg}, {"t", t}, {"tdg", tdg},
        {"v", v}, {"vdg", vdg}, {"sx", sx}, {"sxdg", sxdg}
    };
    if (auto it = kCliffordT.find(StringHelper::toLowerCase(name)); it != kCliffordT.end())
        return it->second();
    throw std::invalid_argument("gate '" + name +
                                 "' is not exactly representable in D[w] (only Clifford+T "
                                 "gates i,x,y,z,h,s,sdg,t,tdg,v,vdg,sx,sxdg are supported in exact mode)");
}

std::array<Dw, 16> twoQubitByName(const std::string &name) {
    static const std::unordered_map<std::string, std::array<Dw, 16> (*)()> kTwoQubit{
        {"swap", swap}, {"iswap", iswap}, {"iswapdg", iswapdg}, {"dcx", dcx}
    };
    if (auto it = kTwoQubit.find(StringHelper::toLowerCase(name)); it != kTwoQubit.end())
        return it->second();
    throw std::invalid_argument("two-qubit gate '" + name +
                                 "' is not exactly representable in D[w] (only swap,iswap,iswapdg,dcx "
                                 "are supported in exact mode; build cx/cz as controlled gates)");
}

} // namespace dd::exact::gates

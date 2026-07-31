#ifndef DD_EXACT_DW_HPP
#define DD_EXACT_DW_HPP

#include <boost/multiprecision/number.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

#if defined(DD_EXACT_WITH_GMP)
// GMP-backed types (by default), which handle coefficients grow large effectively
#include <boost/multiprecision/gmp.hpp>
#else
// Boost's native C++ types, which keep small values inline rather than
// heap-allocating every one. Faster only while coefficients stay small.
#include <boost/multiprecision/cpp_int.hpp>
#endif
#include <complex>
#include <cstddef>
#include <functional>
#include <iostream>
#include <string>

namespace dd::exact {

#if defined(DD_EXACT_WITH_GMP)
using Integer = boost::multiprecision::mpz_int;
using Float   = boost::multiprecision::mpf_float_50;
#else
using Integer = boost::multiprecision::cpp_int;
using Float   = boost::multiprecision::cpp_dec_float_50;
#endif

/**
 * Exact representation of a value in the field Q[w], i.e. dyadic rationals
 * adjoined with w = e^{i*pi/4} (w^4 = -1, w^8 = 1), further extended with
 * arbitrary odd-integer denominators. This covers D[w] (the ring of dyadic
 * rationals adjoined with w) as the special case where the odd denominator
 * is 1, and all of Q[w] beyond it (needed to represent multiplicative
 * inverses of non-unit D[w] elements, which generally require such a
 * denominator).
 *
 * A value is stored as (a + b*w + c*w^2 + d*w^3) / (sqrt(2)^k * e), with
 * a,b,c,d arbitrary-precision integers over the Z-basis {1,w,w^2,w^3} of
 * Z[w], k a non-negative exponent, and e a positive odd integer. Every Dw
 * instance is kept in canonical form:
 *   - k is minimal, i.e. (a,b,c,d) is not divisible by sqrt(2) as a Z[w]
 *     element unless k is already 0 (same invariant as before);
 *   - e is always odd and positive (any factor of 2 discovered in e is
 *     folded into k instead, since (1/sqrt(2))^2 == 1/2);
 *   - gcd(a,b,c,d,e) == 1 (the rational part of the denominator is reduced
 *     as far as possible).
 * This makes equality and hashing exact, with no floating-point tolerance
 * involved -- see Dw::canonicalize. The default e = 1 on every constructor
 * means all pre-existing D[w]-only usage (gate matrices, plain arithmetic)
 * is completely unaffected by this extension.
 */
class Dw {
public:
    Dw() = default;
    Dw(Integer a, Integer b, Integer c, Integer d, std::size_t k = 0, Integer e = 1);

    /// Implicit conversion from a plain integer (k = 0, b = c = d = 0, e = 1).
    // NOLINTNEXTLINE(google-explicit-constructor)
    Dw(long value);

    [[nodiscard]] static Dw zero() { return {}; }
    [[nodiscard]] static Dw one() { return Dw{1}; }
    /// The ring element sqrt(2) = w - w^3.
    [[nodiscard]] static Dw sqrt2();
    /// The ring element w = e^{i*pi/4}.
    [[nodiscard]] static Dw omega();

    [[nodiscard]] const Integer &a() const { return a_; }
    [[nodiscard]] const Integer &b() const { return b_; }
    [[nodiscard]] const Integer &c() const { return c_; }
    [[nodiscard]] const Integer &d() const { return d_; }
    [[nodiscard]] std::size_t k() const { return k_; }
    /// The odd, positive integer denominator (1 for a plain D[w] value).
    [[nodiscard]] const Integer &e() const { return e_; }

    [[nodiscard]] bool isZero() const noexcept;
    [[nodiscard]] bool isOne() const noexcept;

    [[nodiscard]] bool operator==(const Dw &other) const noexcept;
    [[nodiscard]] bool operator!=(const Dw &other) const noexcept { return !(*this == other); }

    [[nodiscard]] Dw operator+(const Dw &other) const;
    [[nodiscard]] Dw operator-(const Dw &other) const;
    [[nodiscard]] Dw operator*(const Dw &other) const;
    [[nodiscard]] Dw operator-() const;

    Dw &operator+=(const Dw &other) { return *this = *this + other; }
    Dw &operator-=(const Dw &other) { return *this = *this - other; }
    Dw &operator*=(const Dw &other) { return *this = *this * other; }

    /// Complex conjugate: w -> -w^3, w^2 -> -w^2, i.e. (a,b,c,d) -> (a,-d,-c,-b).
    [[nodiscard]] Dw conjugate() const;

    /// |this|^2 == this * this->conjugate(). Always a "real" element of
    /// D[w], i.e. one with c() == 0 and d() == -b() (of the form
    /// x + y*sqrt(2) for integers x,y, scaled by sqrt(2)^-k / e).
    [[nodiscard]] Dw normSquared() const { return *this * conjugate(); }

    /// Multiplicative inverse in Q[w] (paper's Algorithm 2). Precondition:
    /// !isZero(). See Dw.cpp for the derivation (matches the paper's
    /// worked Example 8: (1+i*sqrt(2))^-1 == (1-i*sqrt(2))/3).
    [[nodiscard]] Dw inverse() const;

    /// Greatest common divisor of *this and other in the Euclidean ring
    /// D[w] (paper's Algorithm 3), via repeated rounded division. Result is
    /// not yet a canonical associate -- see reduceAssociate(). Precondition:
    /// e() == 1 && other.e() == 1 (genuine D[w] elements; GCDs are a D[w]-
    /// only notion, unlike inverse() which needs the full field Q[w]).
    [[nodiscard]] Dw gcd(const Dw &other) const;

    /// Picks the associate (via multiplication by the unit 1+sqrt(2) and
    /// its inverse) with minimal quarticNorm(), as a step towards a
    /// canonical GCD representative (paper's Algorithm 3, reduceNorm step).
    [[nodiscard]] Dw reduceNorm() const;

    /// Among the 4 rotations of *this by w, picks the one whose
    /// (|a|,|b|,|c|,|d|) tuple is lexicographically smallest, then forces
    /// d() positive by negation (paper's Algorithm 3, lexicographicalMinimum
    /// step).
    [[nodiscard]] Dw lexicographicalMinimum() const;

    /// reduceNorm() followed by lexicographicalMinimum(): a full canonical-
    /// associate selection for a GCD result (paper's Algorithm 3).
    [[nodiscard]] Dw reduceAssociate() const { return reduceNorm().lexicographicalMinimum(); }

    /// E(z) = |(a^2+b^2+c^2+d^2)^2 - 2*(ab+bc+cd+da)^2|, the quartic norm
    /// used by the paper to show D[w] is Euclidean (E(r) <= (9/16)*E(z2) for
    /// the remainder of gcd()'s rounded division). Defined directly on the
    /// stored (a,b,c,d) tuple; primarily test/reduceNorm() support.
    [[nodiscard]] Integer quarticNorm() const;

    /// Approximate value as a std::complex<double>, for debugging/display
    /// and for converting exact quantities to double at integration seams.
    [[nodiscard]] std::complex<double> toComplexDouble() const;

    [[nodiscard]] std::complex<Float> toComplexFloat() const;

    [[nodiscard]] std::string toString() const;

    [[nodiscard]] std::size_t hash() const noexcept;

private:
    Integer a_{0};
    Integer b_{0};
    Integer c_{0};
    Integer d_{0};
    std::size_t k_{0};
    Integer e_{1};

    /// Restores the full canonical-form invariant (minimal k, e odd and
    /// positive, gcd(a,b,c,d,e) == 1) after construction/arithmetic.
    void canonicalize();

    /// Phase 1 of canonicalize(): reduce k by dividing (a,b,c,d) by sqrt(2)
    /// (as a Z[w] element) for as long as that division is exact. Purely a
    /// numerator-side reduction, independent of e_.
    void reduceSqrt2Power();

    /// Phase 2 of canonicalize(): first fold all factors of 2 from e_ into k_
    /// (since (1/sqrt(2))^2 == 1/2), making e_ odd; then reduce
    /// gcd(a,b,c,d,e_) to 1. Because that gcd is odd at this point, dividing
    /// (a,b,c,d) by it does not change the parity conditions used by
    /// reduceSqrt2Power().
    void reduceRationalDenominator();
};

std::ostream &operator<<(std::ostream &os, const Dw &v);

} // namespace dd::exact

namespace std {
template <> struct hash<dd::exact::Dw> {
    std::size_t operator()(const dd::exact::Dw &v) const noexcept { return v.hash(); }
};
} // namespace std

#endif // DD_EXACT_DW_HPP

#include "dd/exact/Dw.hpp"
#include "utility/HashUtil.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace dd::exact {

namespace {
/// Negacyclic convolution of two degree-<4 polynomials in w, reduced modulo
/// w^4 = -1. This is Z[w] multiplication expressed on the {1,w,w^2,w^3}
/// coefficient tuple, independent of any sqrt(2)^k/e scaling.
std::array<Integer, 4> mulTuple(const std::array<Integer, 4> &x, const std::array<Integer, 4> &y) {
    const auto &[a0, a1, a2, a3] = x;
    const auto &[b0, b1, b2, b3] = y;
    return {
        a0 * b0 - a1 * b3 - a2 * b2 - a3 * b1,
        a0 * b1 + a1 * b0 - a2 * b3 - a3 * b2,
        a0 * b2 + a1 * b1 + a2 * b0 - a3 * b3,
        a0 * b3 + a1 * b2 + a2 * b1 + a3 * b0,
    };
}

const std::array<Integer, 4> kSqrt2Tuple{0, 1, 0, -1};

/// Scale an integer tuple up by multiplying by sqrt(2) `steps` times, so
/// that it can be expressed over a larger, common denominator exponent (or,
/// used standalone, to compute an exact multiplication of the underlying
/// value by sqrt(2)^steps).
std::array<Integer, 4> scaleUp(std::array<Integer, 4> tuple, std::size_t steps) {
    if (steps == 0) {
        return tuple;
    }
    const std::size_t half = steps / 2;
    if (half != 0) {
        const Integer twoPowHalf = Integer(1) << half; // 2^half exactly
        for (auto &x : tuple) {
            x *= twoPowHalf;
        }
        // tuple = mulTuple(tuple, {twoPowHalf, 0, 0, 0});
    }
    if ((steps & 1U) != 0U) {
        const auto &[a, b, c, d] = tuple;
        tuple = {b - d, a + c, b + d, c - a}; // multiply by sqrt(2) = w - w^3
        // tuple = mulTuple(tuple, kSqrt2Tuple);
    }
    // for (std::size_t i = 0; i < steps; ++i) {
    //     tuple = mulTuple(tuple, kSqrt2Tuple);
    // }
    return tuple;
}

/// Round num/den (den > 0) to the nearest integer, ties away from zero.
Integer roundRatio(const Integer &num, const Integer &den) {
    Integer q = num / den; // truncates toward zero
    Integer r = num % den; // same sign as num
    if (r == 0)
        return q;
    const Integer twiceAbsR = 2 * (r < 0 ? -r : r);
    if (twiceAbsR >= den) {
        q += (num < 0 ? -1 : 1);
    }
    return q;
}

/// Rounds z1/z2 (an exact Q[w] quotient) to the nearest Z[w] point by
/// rounding each of its 4 basis coefficients independently to the nearest
/// integer -- exactly the paper's Algorithm 3 division-with-rounding step.
Dw roundedDivide(const Dw &z1, const Dw &z2) {
    const Dw quotient = z1 * z2.inverse();
    // quotient = (a,b,c,d) / (sqrt2^k * e). Rationalize by scaling the
    // tuple up by sqrt2^k (turning the sqrt(2) part of the denominator into
    // a plain integer 2^k), giving a single plain-integer denominator
    // 2^k * e shared by all 4 components.
    const std::array<Integer, 4> rationalized =
        scaleUp({quotient.a(), quotient.b(), quotient.c(), quotient.d()}, quotient.k());
    const Integer den = (Integer(1) << quotient.k()) * quotient.e();
    return {roundRatio(rationalized[0], den), roundRatio(rationalized[1], den), roundRatio(rationalized[2], den),
            roundRatio(rationalized[3], den)};
}
} // namespace

Dw::Dw(Integer a, Integer b, Integer c, Integer d, std::size_t k, Integer e)
    : a_(std::move(a)), b_(std::move(b)), c_(std::move(c)), d_(std::move(d)), k_(k), e_(std::move(e)) {
    canonicalize();
}

Dw::Dw(long value) : a_(value) {}

Dw Dw::sqrt2() { return {0, 1, 0, -1, 0}; }

Dw Dw::omega() { return {0, 1, 0, 0, 0}; }

void Dw::canonicalize() {
    if (a_ == 0 && b_ == 0 && c_ == 0 && d_ == 0) {
        k_ = 0;
        e_ = 1;
        return;
    }
    reduceRationalDenominator();
    reduceSqrt2Power();
}

// Paper's Algorithm 1, "Compute Minimal D[w] Representation", translated
// into this class's basis order. The paper writes a value as
// (a*w^3 + b*w^2 + c*w + d)/sqrt(2)^k -- i.e. its a is our d_, its b is our
// c_, its c is our b_, its d is our a_ (the reverse of this class's
// a_=coeff(1), b_=coeff(w), c_=coeff(w^2), d_=coeff(w^3)). Substituting that
// mapping into the paper's derivation (sqrt(2) = w - w^3, so
// alpha = alpha*sqrt(2)/sqrt(2), expanded and regrouped by power of w)
// gives, entirely in this class's own basis:
//   a' = (b-d)/2, b' = (a+c)/2, c' = (b+d)/2, d' = (c-a)/2, k' = k-1
// valid (a',b',c',d' all integers) exactly when a=c (mod 2) and b=d
// (mod 2) -- the paper's stated minimality criterion, translated the same
// way. Applying this while the criterion holds finds the minimal k
// directly, via the same total arithmetic the paper's closed-form
// substitution performs, rather than (equivalently, but at the cost of a
// full Z[w] ring multiplication per step) multiplying by sqrt(2) and
// checking 4-way divisibility by 2.
void Dw::reduceSqrt2Power() {
    while (k_ > 0 && (a_ - c_) % 2 == 0 && (b_ - d_) % 2 == 0) {
        const Integer newA = (b_ - d_) / 2;
        const Integer newB = (a_ + c_) / 2;
        const Integer newC = (b_ + d_) / 2;
        const Integer newD = (c_ - a_) / 2;
        a_ = newA;
        b_ = newB;
        c_ = newC;
        d_ = newD;
        --k_;
    }
}

void Dw::reduceRationalDenominator() {
    // Phase 1: Extract all powers of 2 from e into k exponent.
    // This ensures e_ becomes odd before we compute gcd.
    while (e_ % 2 == 0) {
        e_ /= 2;
        k_ += 2;
    }

    // Phase 2: Reduce gcd(a,b,c,d,e_) to 1.
    // Since e_ is now odd, gcd contains only odd prime factors.
    // Dividing (a,b,c,d) by an odd factor preserves the parity of (a-c)
    // and (b-d), so this phase alone can never expose a new sqrt(2) common
    // factor -- only phase 1 above can.
    if (e_ != 1) {
        const Integer g = boost::multiprecision::gcd(
            boost::multiprecision::gcd(boost::multiprecision::gcd(a_, b_), boost::multiprecision::gcd(c_, d_)), e_);
        if (g > 1) {
            a_ /= g;
            b_ /= g;
            c_ /= g;
            d_ /= g;
            e_ /= g;
        }
    }
}

bool Dw::isZero() const noexcept { return a_ == 0 && b_ == 0 && c_ == 0 && d_ == 0; }

bool Dw::isOne() const noexcept { return k_ == 0 && e_ == 1 && a_ == 1 && b_ == 0 && c_ == 0 && d_ == 0; }

bool Dw::operator==(const Dw &other) const noexcept {
    return k_ == other.k_ && e_ == other.e_ && a_ == other.a_ && b_ == other.b_ && c_ == other.c_ &&
           d_ == other.d_;
}

Dw Dw::operator+(const Dw &other) const {
    const std::size_t k = std::max(k_, other.k_);
    auto lhs = scaleUp({a_, b_, c_, d_}, k - k_);
    auto rhs = scaleUp({other.a_, other.b_, other.c_, other.d_}, k - other.k_);
    // Common (not necessarily minimal) denominator e_ * other.e_;
    // canonicalize() reduces it back down.
    for (auto &x : lhs)
        x *= other.e_;
    for (auto &x : rhs)
        x *= e_;
    return {lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2], lhs[3] + rhs[3], k, e_ * other.e_};
}

Dw Dw::operator-() const { return {-a_, -b_, -c_, -d_, k_, e_}; }

Dw Dw::operator-(const Dw &other) const { return *this + (-other); }

Dw Dw::operator*(const Dw &other) const {
    if (isOne()) return other;
    if (other.isOne()) return *this;
    if (isZero() || other.isZero()) return Dw::zero();
    const auto product = mulTuple({a_, b_, c_, d_}, {other.a_, other.b_, other.c_, other.d_});
    return {product[0], product[1], product[2], product[3], k_ + other.k_, e_ * other.e_};
}

Dw Dw::conjugate() const { return {a_, -d_, -c_, -b_, k_, e_}; }

Dw Dw::inverse() const {
    if (isZero()) {
        throw std::domain_error("Dw::inverse: division by zero");
    }
    const std::array<Integer, 4> conjN{a_, -d_, -c_, -b_};
    const std::array<Integer, 4> normTuple = mulTuple({a_, b_, c_, d_}, conjN);
    // normTuple is always of the form (x, y, 0, -y) -- see normSquared()'s
    // doc comment (a mathematical fact of Z[w] conjugation).
    const Integer &x = normTuple[0];
    const Integer &y = normTuple[1];

    // 1/(x + y*sqrt(2)) = (x - y*sqrt(2)) / (x^2 - 2*y^2), exact integers.
    Integer d = x * x - 2 * y * y;
    // d == 0 would require sqrt(2) == +-x/y (rational), impossible for
    // integers not both zero; (x,y)==(0,0) would mean *this == 0, already
    // excluded above. Unreachable for a valid nonzero input.
    assert(d != 0);

    // (x - y*sqrt(2)) as a Z[w] tuple: sqrt2 = w - w^3, so x - y*sqrt(2) =
    // x - y*w + y*w^3 -> (x, -y, 0, y).
    const std::array<Integer, 4> invNormNumer{x, -y, 0, y};
    std::array<Integer, 4> numerator = mulTuple(conjN, invNormNumer);

    // Fold in this->k_: multiplying the value by sqrt(2)^k_ (as opposed to
    // dividing by it) is done by scaling the raw tuple by sqrt(2)^k_
    // directly, leaving the exponent field at 0 for this contribution.
    numerator = scaleUp(numerator, k_);

    // Fold in this->e_ (a plain integer multiplier: dividing by 1/e_ ==
    // multiplying by e_).
    for (auto &v : numerator)
        v *= e_;

    // Normalize d's sign into the numerator (Dw's e must stay positive).
    if (d < 0) {
        d = -d;
        for (auto &v : numerator)
            v = -v;
    }

    return {numerator[0], numerator[1], numerator[2], numerator[3], 0, d};
}

Dw Dw::gcd(const Dw &other) const {
    assert(e_ == 1 && other.e_ == 1); // gcd() is a D[w]-only notion (see header doc)
    Dw z1 = *this;
    Dw z2 = other;
    while (!z2.isZero()) {
        const Dw q = roundedDivide(z1, z2);
        const Dw r = z1 - q * z2;
        z1 = z2;
        z2 = r;
    }
    return z1;
}

Dw Dw::reduceNorm() const {
    if (isZero())
        return *this;
    const Dw unit = Dw::sqrt2() + Dw::one(); // 1 + sqrt(2)
    const Dw unitInv = unit.inverse();       // sqrt(2) - 1
    Dw candidate = *this;
    Integer bestNorm = candidate.quarticNorm();
    bool improved = true;
    while (improved) {
        improved = false;
        const Dw up = candidate * unit;
        const Dw down = candidate * unitInv;
        const Integer upNorm = up.quarticNorm();
        const Integer downNorm = down.quarticNorm();
        if (upNorm < bestNorm && upNorm <= downNorm) {
            candidate = up;
            bestNorm = upNorm;
            improved = true;
        } else if (downNorm < bestNorm) {
            candidate = down;
            bestNorm = downNorm;
            improved = true;
        }
    }
    return candidate;
}

Dw Dw::lexicographicalMinimum() const {
    auto absTuple = [](const Dw &v) {
        auto abs1 = [](const Integer &x) { return x < 0 ? -x : x; };
        return std::array<Integer, 4>{abs1(v.a()), abs1(v.b()), abs1(v.c()), abs1(v.d())};
    };
    Dw best = *this;
    Dw rotated = *this;
    for (int i = 1; i < 4; ++i) {
        rotated = rotated * Dw::omega();
        if (absTuple(rotated) < absTuple(best)) {
            best = rotated;
        }
    }
    if (best.d() < 0) {
        best = -best;
    }
    return best;
}

Integer Dw::quarticNorm() const {
    const Integer s = a_ * a_ + b_ * b_ + c_ * c_ + d_ * d_;
    const Integer t = a_ * b_ + b_ * c_ + c_ * d_ + d_ * a_;
    const Integer result = s * s - 2 * t * t;
    return result < 0 ? -result : result;
}

std::complex<double> Dw::toComplexDouble() const {
    const std::complex<Float> value = toComplexFloat();
    return {static_cast<double>(value.real()), static_cast<double>(value.imag())};
}

std::complex<Float> Dw::toComplexFloat() const {
    // 1. High-precision constants
    static const Float sqrt2 = boost::multiprecision::sqrt(Float(2.0));
    static const std::complex<Float> omegaF{sqrt2 / 2.0, sqrt2 / 2.0};

    // Precompute powers to optimize performance
    static const std::complex<Float> omegaF2 = omegaF * omegaF;
    static const std::complex<Float> omegaF3 = omegaF2 * omegaF;

    // 2. Use .template convert_to<Float>() to eliminate ambiguity errors
    //    for both cpp_int and mpz_int backends.
    const std::complex<Float> value = Float(a_)
                                    + Float(b_) * omegaF
                                    + Float(c_) * omegaF2
                                    + Float(d_) * omegaF3;

    // 3. Exponentiate with a native int to bypass template deduction issues
    Float scale = boost::multiprecision::pow(Float(2.0), static_cast<int>(k_ / 2));
    if (k_ % 2 != 0) {
        scale *= sqrt2;
    }

    // 4. Return final safely divided complex number
    return value / (scale * Float(e_));
}

std::string Dw::toString() const {
    std::ostringstream os;
    os << "(" << a_ << " + " << b_ << "*w + " << c_ << "*w^2 + " << d_ << "*w^3) / (sqrt(2)^" << k_ << " * " << e_
       << ")";
    return os.str();
}

std::size_t Dw::hash() const noexcept {
    std::size_t h = hash_value(a_);
    for (const Integer *v : {&b_, &c_, &d_, &e_}) {
        h = HashUtil::combinedHash(h, hash_value(*v));
    }
    h = HashUtil::combinedHash(h, std::hash<std::size_t>{}(k_));
    return h;
}

std::ostream &operator<<(std::ostream &os, const Dw &v) { return os << v.toString(); }

} // namespace dd::exact

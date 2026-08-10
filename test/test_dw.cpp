#include "dd/exact/Dw.hpp"

#include <boost/multiprecision/number.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <cstddef>
#include <random>
#include <stdexcept>
#include <utility>

using dd::exact::Dw;
using dd::exact::Integer;

TEST(Dw, ZeroAndOne) {
    EXPECT_TRUE(Dw::zero().isZero());
    EXPECT_FALSE(Dw::zero().isOne());
    EXPECT_TRUE(Dw::one().isOne());
    EXPECT_FALSE(Dw::one().isZero());
    EXPECT_EQ(Dw::zero(), Dw(0));
    EXPECT_EQ(Dw::one(), Dw(1));
}

TEST(Dw, Equality) {
    EXPECT_EQ(Dw(3), Dw(3));
    EXPECT_NE(Dw(3), Dw(4));
    EXPECT_EQ(Dw::omega(), Dw(0, 1, 0, 0));
}

TEST(Dw, Sqrt2SquaredIsTwo) {
    const Dw sqrt2 = Dw::sqrt2();
    EXPECT_EQ(sqrt2 * sqrt2, Dw(2));
}

TEST(Dw, OmegaToFourthIsMinusOne) {
    const Dw w = Dw::omega();
    const Dw w2 = w * w;
    const Dw w4 = w2 * w2;
    EXPECT_EQ(w4, Dw(-1));
}

TEST(Dw, OmegaToEighthIsOne) {
    const Dw w = Dw::omega();
    Dw w8 = Dw::one();
    for (int i = 0; i < 8; ++i) {
        w8 = w8 * w;
    }
    EXPECT_EQ(w8, Dw::one());
}

// 1/sqrt(2) is exactly representable and, once canonicalized, has a trivial
// numerator: (1,0,0,0)/sqrt(2)^1.
TEST(Dw, ReciprocalSqrt2CanonicalForm) {
    const Dw sqrt2 = Dw::sqrt2();
    const Dw invSqrt2(0, 1, 0, -1, 2); // (w - w^3) / sqrt(2)^2 == 1/sqrt(2)
    EXPECT_EQ(invSqrt2.a(), 1);
    EXPECT_EQ(invSqrt2.b(), 0);
    EXPECT_EQ(invSqrt2.c(), 0);
    EXPECT_EQ(invSqrt2.d(), 0);
    EXPECT_EQ(invSqrt2.k(), 1U);
    // sqrt(2) * (1/sqrt(2)) == 1
    EXPECT_EQ(sqrt2 * invSqrt2, Dw::one());
}

TEST(Dw, AdditionAlignsDenominators) {
    const Dw invSqrt2(0, 1, 0, -1, 2); // == 1/sqrt(2), k reduces to 1
    const Dw sum = invSqrt2 + invSqrt2; // == 2/sqrt(2) == sqrt(2)
    EXPECT_EQ(sum, Dw::sqrt2());
}

TEST(Dw, ArithmeticIdentities) {
    const Dw x(3, -2, 1, 5);
    const Dw y(-4, 7, 0, 2);
    EXPECT_EQ(x + y, y + x);
    EXPECT_EQ(x * y, y * x);
    EXPECT_EQ(x - x, Dw::zero());
    EXPECT_EQ(x * Dw::one(), x);
    EXPECT_EQ(x + Dw::zero(), x);
}

TEST(Dw, ConjugateOfRealIsItself) {
    const Dw real(7);
    EXPECT_EQ(real.conjugate(), real);
}

TEST(Dw, ConjugateOfOmega) {
    // conjugate(w) == -w^3
    EXPECT_EQ(Dw::omega().conjugate(), Dw(0, 0, 0, -1));
}

TEST(Dw, NormSquaredOfOmegaIsOne) {
    EXPECT_EQ(Dw::omega().normSquared(), Dw::one());
}

TEST(Dw, NormSquaredOfInvSqrt2IsOneHalf) {
    const Dw invSqrt2(0, 1, 0, -1, 2); // 1/sqrt(2)
    const Dw expected(1, 0, 0, 0, 2);  // 1/2
    EXPECT_EQ(invSqrt2.normSquared(), expected);
}

TEST(Dw, HashConsistentWithEquality) {
    const Dw x(0, 1, 0, -1, 2);
    const Dw y(1, 0, 0, 0, 1); // same value as x after canonicalization
    EXPECT_EQ(x, y);
    EXPECT_EQ(std::hash<Dw>{}(x), std::hash<Dw>{}(y));
}

TEST(Dw, ToComplexDoubleApproximatesOmega) {
    const auto v = Dw::omega().toComplexDouble();
    EXPECT_NEAR(v.real(), std::sqrt(2.0) / 2.0, 1e-9);
    EXPECT_NEAR(v.imag(), std::sqrt(2.0) / 2.0, 1e-9);
}

TEST(Dw, ToComplexDoubleApproximatesInvSqrt2) {
    const Dw invSqrt2(0, 1, 0, -1, 2);
    const auto v = invSqrt2.toComplexDouble();
    EXPECT_NEAR(v.real(), 1.0 / std::sqrt(2.0), 1e-9);
    EXPECT_NEAR(v.imag(), 0.0, 1e-9);
}

// ---------------------------------------------------------------------
// inverse()
//
// inverse() had no direct coverage before these. It also carries the two
// closed-form specialisations in Dw.cpp (normXY, mulByConjNormNumer) that
// replace generic 16-product convolutions, and a wrong identity there would
// produce plausible-but-wrong amplitudes rather than a crash -- so the
// defining property is asserted over random inputs, not just examples.
// ---------------------------------------------------------------------

TEST(Dw, InverseOfZeroThrows) { EXPECT_THROW((void)Dw::zero().inverse(), std::domain_error); }

TEST(Dw, InverseTimesValueIsOneOnKnownCases) {
    const Dw cases[] = {
        Dw::one(),  Dw::omega(),         Dw::sqrt2(),
        Dw(3),      Dw(0, 1, 0, -1, 2), // 1/sqrt(2)
        Dw(1, 1, 0, 0),                 // 1 + w
        Dw(1, 0, 0, 1),                 // 1 + w^3 (exercises the -ad term in y)
        Dw(-5, 2, -7, 3, 5),            // large k
        Dw(4, -3, 2, -1, 3, 9),         // e != 1
    };
    for (const auto &z : cases) {
        EXPECT_EQ(z * z.inverse(), Dw::one()) << "failed for " << z.toString();
    }
}

// The specialisation in normXY() relies on z * conj(z) always having the shape
// (x, y, 0, -y). canonicalize() preserves that shape (dividing all four
// components by a common factor obviously does; the sqrt(2) reduction maps
// (x,y,0,-y) to (y, x/2, 0, -x/2)), so it can be asserted on the canonical
// product built from the *generic* operator* -- i.e. independently of the
// specialised code path it justifies.
TEST(Dw, NormTupleHasExpectedShape) {
    std::mt19937_64 rng(20260808);
    std::uniform_int_distribution<int> coeff(-40, 40);
    std::uniform_int_distribution<std::size_t> kDist(0, 6);
    for (int i = 0; i < 400; ++i) {
        const Dw z(coeff(rng), coeff(rng), coeff(rng), coeff(rng), kDist(rng));
        if (z.isZero()) {
            continue;
        }
        const Dw norm = z * z.conjugate();
        EXPECT_EQ(norm.c(), 0) << "failed for " << z.toString();
        EXPECT_EQ(norm.d(), -norm.b()) << "failed for " << z.toString();
    }
}

TEST(Dw, InverseTimesValueIsOneRandomised) {
    std::mt19937_64 rng(20260809);
    std::uniform_int_distribution<int> coeff(-60, 60);
    std::uniform_int_distribution<std::size_t> kDist(0, 8);
    std::uniform_int_distribution<int> oddE(0, 12); // e = 2*n + 1, kept odd and positive
    int checked = 0;
    for (int i = 0; i < 600; ++i) {
        const Dw z(coeff(rng), coeff(rng), coeff(rng), coeff(rng), kDist(rng), 2 * oddE(rng) + 1);
        if (z.isZero()) {
            continue;
        }
        const Dw inv = z.inverse();
        EXPECT_EQ(z * inv, Dw::one()) << "failed for " << z.toString();
        EXPECT_EQ(inv * z, Dw::one()) << "failed for " << z.toString();
        ++checked;
    }
    EXPECT_GT(checked, 500); // guard against the generator degenerating to zero
}

// operator+ has three paths (see the test below for which). They must all
// compute the same value, and each shortcut skips work the general path does,
// so an error in one would not show up in the others. Note that a "typical"
// random sweep is not enough here: under Inverse normalization most operands
// carry a non-unit denominator (measured at 83% of left operands on Grover),
// so the interesting cases are the ones a naive generator reaches least often.

namespace {
/// Random Dw over a range wide enough to produce both k mismatches and
/// non-unit denominators. `denomWidth` of 0 forces e == 1 (a genuine D[w]
/// value, the fast path's domain).
Dw randomDw(std::mt19937_64 &rng, int denomWidth, std::size_t kMax) {
    std::uniform_int_distribution<int> coeff(-60, 60);
    std::uniform_int_distribution<std::size_t> kDist(0, kMax);
    std::uniform_int_distribution<int> oddE(0, denomWidth);
    return Dw(coeff(rng), coeff(rng), coeff(rng), coeff(rng), kDist(rng), 2 * oddE(rng) + 1);
}

/// A random Dw whose canonical form has exactly the requested `k` and odd
/// denominator `oddE`. Both are pinned on purpose, so that a test can select
/// which branch of operator+ it exercises instead of waiting for the right
/// combination to come up by chance:
///   - a == 1 forces gcd(a,b,c,d,e) == 1, so reduceRationalDenominator()
///     cannot divide the denominator away;
///   - an even third coefficient breaks the a == c (mod 2) criterion, so
///     reduceSqrt2Power() cannot lower k.
Dw randomDwWith(std::mt19937_64 &rng, int oddE, std::size_t k) {
    std::uniform_int_distribution<int> coeff(-30, 30);
    return Dw(1, coeff(rng), 2 * coeff(rng), coeff(rng), k, oddE);
}
} // namespace

TEST(Dw, AdditionObeysGroupAxiomsAcrossAllPaths) {
    // operator+ picks one of three paths, on two independent conditions --
    // whether the sqrt(2) exponents already agree, and whether the denominators
    // are equal:
    //   1. k aligned, e_ == other.e_  -> coefficient-wise add, no tuple built;
    //   2. k differs, e_ == other.e_  -> scale one side, add, keep the shared e;
    //   3. e_ != other.e_             -> cross-multiply onto e_ * other.e_.
    //
    // Paths 1 and 2 with a *non-unit* shared denominator are the ones worth
    // guarding: they return the shared e directly rather than inflating the
    // denominator to e^2 and leaving canonicalize()'s gcd to reduce it back
    // (which is the same value by a more expensive route). A purely random
    // sweep would reach that case only when two independent draws happened to
    // land on the same denominator, so the four shapes are constructed instead
    // of waited for, and each is counted.
    std::mt19937_64 rng(20260810);

    // The pinning helper is load-bearing for the coverage claim below, so
    // check it does what it says before relying on it.
    const Dw pinned = randomDwWith(rng, 7, 4);
    ASSERT_EQ(pinned.e(), 7) << "randomDwWith failed to pin the denominator: " << pinned.toString();
    ASSERT_EQ(pinned.k(), 4U) << "randomDwWith failed to pin the exponent: " << pinned.toString();

    int alignedSameE = 0;
    int scaledSameE = 0;
    int differingE = 0;
    int sameEnonUnit = 0;

    for (int i = 0; i < 400; ++i) {
        const auto pair = [&]() -> std::pair<Dw, Dw> {
            switch (i % 4) {
            case 0: // path 1, e == 1
                return {randomDwWith(rng, 1, 3), randomDwWith(rng, 1, 3)};
            case 1: // path 1, shared non-unit denominator
                return {randomDwWith(rng, 7, 3), randomDwWith(rng, 7, 3)};
            case 2: // path 2, shared non-unit denominator, exponents differ
                return {randomDwWith(rng, 9, 2), randomDwWith(rng, 9, 5)};
            default: // path 3, denominators differ
                return {randomDwWith(rng, 7, 3), randomDwWith(rng, 9, 4)};
            }
        }();
        const Dw &x = pair.first;
        const Dw &y = pair.second;
        const Dw z = randomDw(rng, 12, 6);

        // Mirrors the implementation's scaleL/scaleR/sameE exactly.
        const bool sameE = x.e() == y.e();
        if (sameE && x.k() == y.k()) {
            ++alignedSameE;
        } else if (sameE) {
            ++scaledSameE;
        } else {
            ++differingE;
        }
        if (sameE && x.e() != 1) {
            ++sameEnonUnit;
        }

        EXPECT_EQ(x + y, y + x) << "commutativity failed for " << x.toString() << " + " << y.toString();
        EXPECT_EQ((x + y) + z, x + (y + z))
            << "associativity failed for " << x.toString() << ", " << y.toString() << ", " << z.toString();
        EXPECT_EQ(x + Dw::zero(), x) << "additive identity failed for " << x.toString();
        EXPECT_TRUE((x + (-x)).isZero()) << "additive inverse failed for " << x.toString();
    }

    // A path that is never reached is a path this test says nothing about.
    EXPECT_GT(alignedSameE, 50);
    EXPECT_GT(scaledSameE, 50);
    EXPECT_GT(differingE, 50);
    EXPECT_GT(sameEnonUnit, 50) << "the shared-non-unit-denominator shortcut was never exercised";
}

TEST(Dw, MultiplicationDistributesOverAdditionRandomised) {
    // The strongest single check on the pair: distributivity forces operator*
    // and operator+ to agree with each other, so an error in either the
    // coefficient convolution or the scale/denominator handling shows up here
    // even though each operator alone might look self-consistent.
    std::mt19937_64 rng(20260811);
    for (int i = 0; i < 400; ++i) {
        const int width = (i % 3 == 0) ? 0 : 12;
        const Dw x = randomDw(rng, width, 6);
        const Dw y = randomDw(rng, width, 6);
        const Dw z = randomDw(rng, width, 6);
        EXPECT_EQ(x * (y + z), (x * y) + (x * z))
            << "distributivity failed for " << x.toString() << ", " << y.toString() << ", " << z.toString();
    }
}

TEST(Dw, ArithmeticAgreesWithComplexEvaluation) {
    // Cross-check against a path that shares none of the integer tuple
    // arithmetic: evaluate each operand as a complex number and compare. This
    // catches a coefficient-order or sign error that stays self-consistent
    // across the algebraic identities above.
    std::mt19937_64 rng(20260812);
    for (int i = 0; i < 200; ++i) {
        const int width = (i % 2 == 0) ? 0 : 12;
        const Dw x = randomDw(rng, width, 5);
        const Dw y = randomDw(rng, width, 5);

        const auto xv = x.toComplexDouble();
        const auto yv = y.toComplexDouble();
        const auto sum = (x + y).toComplexDouble();
        const auto product = (x * y).toComplexDouble();

        const double scale = 1.0 + std::abs(xv) + std::abs(yv) + std::abs(xv) * std::abs(yv);
        EXPECT_NEAR(sum.real(), xv.real() + yv.real(), 1e-9 * scale);
        EXPECT_NEAR(sum.imag(), xv.imag() + yv.imag(), 1e-9 * scale);
        EXPECT_NEAR(product.real(), (xv * yv).real(), 1e-9 * scale);
        EXPECT_NEAR(product.imag(), (xv * yv).imag(), 1e-9 * scale);
    }
}

// Canonical form requires e to be odd: any factor of 2 in the rational
// denominator must be folded into k, since (1/sqrt(2))^2 == 1/2. If only some
// of those factors are extracted, two representations of the same number
// survive with different (k, e), and operator==/hash() -- which compare the
// stored fields, not the value -- report them unequal. That silently breaks
// hash-consing, so the invariant is asserted directly here rather than being
// left to whichever higher-level test happens to construct such a value.
TEST(Dw, RationalDenominatorIsAlwaysOdd) {
    for (int twos = 0; twos <= 6; ++twos) {
        Integer e = 1;
        for (int j = 0; j < twos; ++j) {
            e *= 2;
        }
        const Dw viaE(1, 0, 0, 0, 0, e);
        EXPECT_NE(viaE.e() % 2, 0) << "e = " << viaE.e() << " is even for 2^" << twos;

        // The same number written with the power of two already in k.
        const Dw viaK(1, 0, 0, 0, static_cast<std::size_t>(2 * twos), 1);
        EXPECT_EQ(viaE, viaK) << "1/2^" << twos << " has two distinct canonical forms: k="
                              << viaE.k() << ",e=" << viaE.e() << " vs k=" << viaK.k()
                              << ",e=" << viaK.e();
        EXPECT_EQ(viaE.hash(), viaK.hash());
    }
}

// The canonical form requires gcd(a,b,c,d,e) == 1. How that factor is found
// (balanced tree vs. progressive fold) must not change the result: if any
// operand is skipped, an over-large e survives and the same number acquires
// two representations, which operator==/hash() -- comparing stored fields, not
// values -- report as different weights. Nothing crashes; hash-consing just
// silently stops consing. So the invariant is asserted directly.
TEST(Dw, ReducingFactorUsesEveryCoefficient) {
    std::mt19937_64 rng(20260813);
    std::uniform_int_distribution<int> coeff(-40, 40);
    std::uniform_int_distribution<std::size_t> kDist(0, 5);
    std::uniform_int_distribution<int> oddFactor(1, 8); // planted factor 2*n+1
    int checked = 0;
    for (int i = 0; i < 500; ++i) {
        // Plant a common odd factor in all four coefficients and in e, so a
        // fold that drops any one operand yields a strictly larger g.
        const Integer f = 2 * oddFactor(rng) + 1;
        const Dw z(f * coeff(rng), f * coeff(rng), f * coeff(rng), f * coeff(rng), kDist(rng), f * (2 * oddFactor(rng) + 1));
        if (z.isZero()) {
            continue;
        }
        const Integer g = boost::multiprecision::gcd(
            boost::multiprecision::gcd(boost::multiprecision::gcd(z.a(), z.b()), boost::multiprecision::gcd(z.c(), z.d())),
            z.e());
        EXPECT_EQ(g, 1) << "gcd(a,b,c,d,e) = " << g << " after canonicalization of " << z.toString();
        EXPECT_NE(z.e() % 2, 0) << "e stayed even: " << z.toString();
        ++checked;
    }
    EXPECT_GT(checked, 450);
}

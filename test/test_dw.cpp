#include "dd/exact/Dw.hpp"

#include <gtest/gtest.h>

using dd::exact::Dw;

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

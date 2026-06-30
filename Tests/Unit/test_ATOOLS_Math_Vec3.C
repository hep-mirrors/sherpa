#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <limits>
#include <sstream>

#include "ATOOLS/Math/Vector.H"

using ATOOLS::Vec3D;
using ATOOLS::Vec4D;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

// Tolerance for exact-in-principle floating-point results
static constexpr double TOL = 1.0e-15;

/// Component-wise approximate equality check for Vec3D.
static bool ApproxEqual(const Vec3D& a, const Vec3D& b,
                        double tol = TOL) {
  return std::abs(a[1] - b[1]) < tol &&
         std::abs(a[2] - b[2]) < tol &&
         std::abs(a[3] - b[3]) < tol;
}

// =========================================================================
// 1. Construction
// =========================================================================

TEST_CASE("Vec3 default construction initialises to zero", "[Vec3][construction]") {
  Vec3D v;
  REQUIRE(v[1] == 0.0);
  REQUIRE(v[2] == 0.0);
  REQUIRE(v[3] == 0.0);
}

TEST_CASE("Vec3 component construction", "[Vec3][construction]") {
  Vec3D v(1.0, -2.5, 3.7);
  REQUIRE(v[1] == 1.0);
  REQUIRE(v[2] == -2.5);
  REQUIRE(v[3] == 3.7);
}

TEST_CASE("Vec3 copy construction", "[Vec3][construction]") {
  Vec3D orig(1.0, 2.0, 3.0);
  Vec3D copy(orig);
  REQUIRE(copy[1] == orig[1]);
  REQUIRE(copy[2] == orig[2]);
  REQUIRE(copy[3] == orig[3]);

  // Verify deep copy: modifying orig must not affect copy.
  orig[1] = 99.0;
  REQUIRE(copy[1] == 1.0);
}

TEST_CASE("Vec3 construction from Vec4 extracts spatial part", "[Vec3][construction]") {
  // Vec4 stores (E, px, py, pz) with 0-based indexing.
  // Vec3 should extract the spatial components (indices 1,2,3 of Vec4).
  Vec4D p(100.0, 1.0, 2.0, 3.0);
  Vec3D v(p);
  REQUIRE(v[1] == 1.0);  // px
  REQUIRE(v[2] == 2.0);  // py
  REQUIRE(v[3] == 3.0);  // pz
}

// =========================================================================
// 2. Sqr and Abs (Euclidean norm)
// =========================================================================

TEST_CASE("Vec3 Sqr and Abs", "[Vec3][norm]") {
  SECTION("3-4-0 triangle: |v|^2 = 25, |v| = 5") {
    // Classic Pythagorean triple in 2D embedded in 3D.
    Vec3D v(3.0, 4.0, 0.0);
    REQUIRE_THAT(v.Sqr(), WithinAbs(25.0, TOL));
    REQUIRE_THAT(v.Abs(), WithinAbs(5.0, TOL));
  }

  SECTION("unit vector has |v| = 1") {
    Vec3D v(1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0));
    REQUIRE_THAT(v.Sqr(), WithinAbs(1.0, TOL));
    REQUIRE_THAT(v.Abs(), WithinAbs(1.0, TOL));
  }

  SECTION("zero vector") {
    Vec3D v;
    REQUIRE(v.Sqr() == 0.0);
    REQUIRE(v.Abs() == 0.0);
  }

  SECTION("single-component vectors") {
    REQUIRE_THAT(Vec3D(7.0, 0.0, 0.0).Abs(), WithinAbs(7.0, TOL));
    REQUIRE_THAT(Vec3D(0.0, -5.0, 0.0).Abs(), WithinAbs(5.0, TOL));
    REQUIRE_THAT(Vec3D(0.0, 0.0, 13.0).Abs(), WithinAbs(13.0, TOL));
  }
}

// =========================================================================
// 3. Dot product (operator*)
// =========================================================================

TEST_CASE("Vec3 Euclidean dot product", "[Vec3][dot]") {
  SECTION("dot product of orthogonal vectors is zero") {
    Vec3D a(1.0, 0.0, 0.0);
    Vec3D b(0.0, 1.0, 0.0);
    REQUIRE_THAT(a * b, WithinAbs(0.0, TOL));
  }

  SECTION("dot product of parallel vectors equals product of norms") {
    Vec3D a(2.0, 0.0, 0.0);
    Vec3D b(3.0, 0.0, 0.0);
    REQUIRE_THAT(a * b, WithinAbs(6.0, TOL));
  }

  SECTION("dot product of anti-parallel vectors is negative") {
    Vec3D a(1.0, 0.0, 0.0);
    Vec3D b(-4.0, 0.0, 0.0);
    REQUIRE_THAT(a * b, WithinAbs(-4.0, TOL));
  }

  SECTION("self dot product equals Sqr") {
    Vec3D v(1.5, -2.3, 4.7);
    REQUIRE_THAT(v * v, WithinAbs(v.Sqr(), TOL));
  }

  SECTION("general case: (1,2,3) . (4,-5,6) = 4 - 10 + 18 = 12") {
    Vec3D a(1.0, 2.0, 3.0);
    Vec3D b(4.0, -5.0, 6.0);
    REQUIRE_THAT(a * b, WithinAbs(12.0, TOL));
  }

  SECTION("commutativity: a.b == b.a") {
    Vec3D a(1.1, -2.2, 3.3);
    Vec3D b(-4.4, 5.5, -6.6);
    REQUIRE_THAT(a * b, WithinAbs(b * a, TOL));
  }
}

// =========================================================================
// 4. Arithmetic operators
// =========================================================================

TEST_CASE("Vec3 addition", "[Vec3][arithmetic]") {
  Vec3D a(1.0, 2.0, 3.0);
  Vec3D b(10.0, 20.0, 30.0);
  Vec3D c = a + b;
  REQUIRE(ApproxEqual(c, Vec3D(11.0, 22.0, 33.0)));
}

TEST_CASE("Vec3 subtraction", "[Vec3][arithmetic]") {
  Vec3D a(10.0, 20.0, 30.0);
  Vec3D b(1.0, 2.0, 3.0);
  Vec3D c = a - b;
  REQUIRE(ApproxEqual(c, Vec3D(9.0, 18.0, 27.0)));
}

TEST_CASE("Vec3 scalar multiplication", "[Vec3][arithmetic]") {
  Vec3D v(1.0, -2.0, 3.0);
  double s = 2.5;

  SECTION("vec * scalar") {
    Vec3D r = v * s;
    REQUIRE(ApproxEqual(r, Vec3D(2.5, -5.0, 7.5)));
  }

  SECTION("scalar * vec (commutativity via free function)") {
    Vec3D r = s * v;
    REQUIRE(ApproxEqual(r, Vec3D(2.5, -5.0, 7.5)));
  }

  SECTION("multiplication by zero yields zero vector") {
    Vec3D r = v * 0.0;
    REQUIRE(ApproxEqual(r, Vec3D(0.0, 0.0, 0.0)));
  }

  SECTION("multiplication by one is identity") {
    Vec3D r = v * 1.0;
    REQUIRE(ApproxEqual(r, v));
  }
}

TEST_CASE("Vec3 scalar division", "[Vec3][arithmetic]") {
  Vec3D v(4.0, -6.0, 8.0);
  Vec3D r = v / 2.0;
  REQUIRE(ApproxEqual(r, Vec3D(2.0, -3.0, 4.0)));
}

TEST_CASE("Vec3 unary negation", "[Vec3][arithmetic]") {
  Vec3D v(1.0, -2.0, 3.0);
  Vec3D neg = -v;
  REQUIRE(ApproxEqual(neg, Vec3D(-1.0, 2.0, -3.0)));

  // Double negation is identity.
  REQUIRE(ApproxEqual(-neg, v));
}

// =========================================================================
// 5. Compound assignment operators
// =========================================================================

TEST_CASE("Vec3 operator+=", "[Vec3][compound]") {
  Vec3D v(1.0, 2.0, 3.0);
  Vec3D w(10.0, 20.0, 30.0);
  v += w;
  REQUIRE(ApproxEqual(v, Vec3D(11.0, 22.0, 33.0)));
}

TEST_CASE("Vec3 operator-=", "[Vec3][compound]") {
  Vec3D v(10.0, 20.0, 30.0);
  Vec3D w(1.0, 2.0, 3.0);
  v -= w;
  REQUIRE(ApproxEqual(v, Vec3D(9.0, 18.0, 27.0)));
}

TEST_CASE("Vec3 operator*=", "[Vec3][compound]") {
  Vec3D v(1.0, -2.0, 3.0);
  v *= 3.0;
  REQUIRE(ApproxEqual(v, Vec3D(3.0, -6.0, 9.0)));
}

TEST_CASE("Vec3 compound operators return *this for chaining", "[Vec3][compound]") {
  Vec3D a(1.0, 1.0, 1.0);
  Vec3D b(2.0, 2.0, 2.0);
  // (a += b) should return a reference to a, allowing further use.
  Vec3D& ref = (a += b);
  REQUIRE(&ref == &a);
  REQUIRE(ApproxEqual(a, Vec3D(3.0, 3.0, 3.0)));
}

// =========================================================================
// 6. Cross product
// =========================================================================

TEST_CASE("Vec3 cross product", "[Vec3][cross]") {
  SECTION("basis vector identities: x cross y = z (cyclic)") {
    // e_x × e_y = e_z
    Vec3D cx = cross(Vec3D(1, 0, 0), Vec3D(0, 1, 0));
    REQUIRE(ApproxEqual(cx, Vec3D(0, 0, 1)));

    // e_y × e_z = e_x
    Vec3D cy = cross(Vec3D(0, 1, 0), Vec3D(0, 0, 1));
    REQUIRE(ApproxEqual(cy, Vec3D(1, 0, 0)));

    // e_z × e_x = e_y
    Vec3D cz = cross(Vec3D(0, 0, 1), Vec3D(1, 0, 0));
    REQUIRE(ApproxEqual(cz, Vec3D(0, 1, 0)));
  }

  SECTION("anti-cyclic: y cross x = -z") {
    Vec3D c = cross(Vec3D(0, 1, 0), Vec3D(1, 0, 0));
    REQUIRE(ApproxEqual(c, Vec3D(0, 0, -1)));
  }

  SECTION("cross product of parallel vectors is zero") {
    Vec3D a(1.0, 2.0, 3.0);
    Vec3D b(2.0, 4.0, 6.0);  // b = 2*a
    Vec3D c = cross(a, b);
    REQUIRE_THAT(c.Abs(), WithinAbs(0.0, TOL));
  }

  SECTION("self cross product is zero") {
    Vec3D v(3.14, 2.72, 1.41);
    Vec3D c = cross(v, v);
    REQUIRE_THAT(c.Abs(), WithinAbs(0.0, TOL));
  }

  SECTION("cross product is perpendicular to both operands") {
    Vec3D a(1.0, 2.0, 3.0);
    Vec3D b(4.0, -1.0, 2.0);
    Vec3D c = cross(a, b);
    REQUIRE_THAT(c * a, WithinAbs(0.0, TOL));
    REQUIRE_THAT(c * b, WithinAbs(0.0, TOL));
  }

  SECTION("|a x b|^2 = |a|^2 |b|^2 - (a.b)^2  (Lagrange identity)") {
    Vec3D a(1.0, 2.0, 3.0);
    Vec3D b(4.0, -1.0, 2.0);
    Vec3D c = cross(a, b);
    double lhs = c.Sqr();
    double rhs = a.Sqr() * b.Sqr() - (a * b) * (a * b);
    REQUIRE_THAT(lhs, WithinRel(rhs, TOL));
  }

  SECTION("anti-commutativity: a x b = -(b x a)") {
    Vec3D a(1.5, -3.2, 0.7);
    Vec3D b(-2.1, 0.4, 5.6);
    Vec3D c1 = cross(a, b);
    Vec3D c2 = cross(b, a);
    REQUIRE(ApproxEqual(c1, -c2));
  }

  SECTION("general numerical case: (1,2,3) x (4,5,6)") {
    // (2*6 - 3*5, 3*4 - 1*6, 1*5 - 2*4) = (-3, 6, -3)
    Vec3D c = cross(Vec3D(1, 2, 3), Vec3D(4, 5, 6));
    REQUIRE(ApproxEqual(c, Vec3D(-3, 6, -3)));
  }
}

// =========================================================================
// 7. Static basis vectors
// =========================================================================

TEST_CASE("Vec3 static basis vectors", "[Vec3][basis]") {
  REQUIRE(ApproxEqual(Vec3D::XVEC, Vec3D(1, 0, 0)));
  REQUIRE(ApproxEqual(Vec3D::YVEC, Vec3D(0, 1, 0)));
  REQUIRE(ApproxEqual(Vec3D::ZVEC, Vec3D(0, 0, 1)));

  // Orthonormality
  REQUIRE_THAT(Vec3D::XVEC * Vec3D::YVEC, WithinAbs(0.0, TOL));
  REQUIRE_THAT(Vec3D::XVEC * Vec3D::ZVEC, WithinAbs(0.0, TOL));
  REQUIRE_THAT(Vec3D::YVEC * Vec3D::ZVEC, WithinAbs(0.0, TOL));
  REQUIRE_THAT(Vec3D::XVEC.Abs(), WithinAbs(1.0, TOL));
  REQUIRE_THAT(Vec3D::YVEC.Abs(), WithinAbs(1.0, TOL));
  REQUIRE_THAT(Vec3D::ZVEC.Abs(), WithinAbs(1.0, TOL));
}

// =========================================================================
// 8. Nan() and IsZero() predicates
// =========================================================================

TEST_CASE("Vec3 Nan detection", "[Vec3][predicate]") {
  SECTION("well-formed vector is not NaN") {
    Vec3D v(1.0, 2.0, 3.0);
    REQUIRE_FALSE(v.Nan());
  }

  SECTION("zero vector is not NaN") {
    Vec3D v;
    REQUIRE_FALSE(v.Nan());
  }

  SECTION("vector with NaN in any component is detected") {
    double nan = std::numeric_limits<double>::quiet_NaN();
    REQUIRE(Vec3D(nan, 0.0, 0.0).Nan());
    REQUIRE(Vec3D(0.0, nan, 0.0).Nan());
    REQUIRE(Vec3D(0.0, 0.0, nan).Nan());
  }
}

TEST_CASE("Vec3 IsZero predicate", "[Vec3][predicate]") {
  // IsZero uses the framework's Accu() = 1e-12 per component.
  SECTION("zero vector") {
    Vec3D v;
    REQUIRE(v.IsZero());
  }

  SECTION("non-zero vector") {
    Vec3D v(1.0, 0.0, 0.0);
    REQUIRE_FALSE(v.IsZero());
  }

  SECTION("vector with components below Accu() threshold") {
    // All components smaller than 1e-12 → considered zero.
    Vec3D v(1.0e-13, -5.0e-14, 9.0e-13);
    REQUIRE(v.IsZero());
  }

  SECTION("vector with one component at threshold boundary") {
    // A component exactly at 1e-12 is not < 1e-12, so IsZero should be false
    // for at least that component.
    Vec3D v(ATOOLS::Accu(), 0.0, 0.0);
    REQUIRE_FALSE(v.IsZero());
  }
}

// =========================================================================
// 9. Stream output
// =========================================================================

TEST_CASE("Vec3 stream output", "[Vec3][io]") {
  Vec3D v(1.0, -2.0, 3.0);
  std::ostringstream os;
  os << v;
  std::string s = os.str();
  // Expect format: (x,y,z)
  REQUIRE(s.front() == '(');
  REQUIRE(s.back() == ')');
  // Must contain separator commas.
  REQUIRE(s.find(',') != std::string::npos);
}

// =========================================================================
// 10. Physics-consistency tests
// =========================================================================

TEST_CASE("Vec3 extraction from lightlike Vec4", "[Vec3][physics]") {
  // A massless particle along z: E = pz, px = py = 0.
  double E = 50.0;
  Vec4D p(E, 0.0, 0.0, E);
  Vec3D p3(p);

  REQUIRE_THAT(p3.Abs(), WithinAbs(E, TOL));
  // Direction should be along z.
  REQUIRE(ApproxEqual(p3 / p3.Abs(), Vec3D::ZVEC));
}

TEST_CASE("Vec3 momentum conservation in simple decay", "[Vec3][physics]") {
  // Isotropic back-to-back decay: p1 + p2 = 0 (rest frame).
  Vec3D p1(30.0, 40.0, 0.0);
  Vec3D p2 = -p1;
  Vec3D total = p1 + p2;
  REQUIRE_THAT(total.Abs(), WithinAbs(0.0, TOL));
}

TEST_CASE("Vec3 triple product for coplanar vectors vanishes", "[Vec3][physics]") {
  // If c lies in the plane of a and b, the scalar triple product
  // a . (b x c) must be zero.
  Vec3D a(1.0, 0.0, 0.0);
  Vec3D b(0.0, 1.0, 0.0);
  Vec3D c(1.0, 1.0, 0.0);  // in the x-y plane
  double triple = a * cross(b, c);
  REQUIRE_THAT(triple, WithinAbs(0.0, TOL));
}

TEST_CASE("Vec3 Jacobi identity for cross product", "[Vec3][physics]") {
  // a x (b x c) + b x (c x a) + c x (a x b) = 0
  Vec3D a(1.0, 2.0, 3.0);
  Vec3D b(4.0, -1.0, 2.0);
  Vec3D c(-1.0, 3.0, -2.0);
  Vec3D sum = cross(a, cross(b, c)) +
              cross(b, cross(c, a)) +
              cross(c, cross(a, b));
  REQUIRE_THAT(sum.Abs(), WithinAbs(0.0, TOL));
}

// =========================================================================
// 11. Numerical robustness
// =========================================================================

TEST_CASE("Vec3 norm of large vectors", "[Vec3][numerical]") {
  // Test with components near the sqrt of DBL_MAX to check against overflow
  // in Sqr().  (1e154)^2 = 1e308, within double range.
  double big = 1.0e154;
  Vec3D v(big, 0.0, 0.0);
  REQUIRE_THAT(v.Abs(), WithinRel(big, TOL));
}

TEST_CASE("Vec3 norm of small vectors", "[Vec3][numerical]") {
  // Components near sqrt(DBL_MIN) ≈ 1.5e-162.
  double small = 1.0e-160;
  Vec3D v(small, 0.0, 0.0);
  REQUIRE_THAT(v.Sqr(), WithinAbs(small * small, TOL));
  REQUIRE_THAT(v.Abs(), WithinAbs(small, TOL));
}

TEST_CASE("Vec3 dot product catastrophic cancellation", "[Vec3][numerical]") {
  // Two nearly-orthogonal vectors: the dot product should be close to zero
  // but floating-point cancellation can inflate relative error.  We check
  // that the absolute error remains bounded.
  double x = 1.0;
  double eps = 1.0e-8;
  Vec3D a(x, eps, 0.0);
  Vec3D b(-eps, x, 0.0);
  // Exact: a.b = -x*eps + eps*x = 0
  REQUIRE_THAT(a * b, WithinAbs(0.0, TOL));
}

TEST_CASE("Vec3 subtraction of nearly equal vectors", "[Vec3][numerical]") {
  // Catastrophic cancellation: a - b when a ≈ b.
  double base = 1.0e6;
  double delta = 1.0e-4;
  Vec3D a(base + delta, base + delta, base + delta);
  Vec3D b(base, base, base);
  Vec3D diff = a - b;
  REQUIRE_THAT(diff[1], WithinRel(delta, 1.0e-6));
  REQUIRE_THAT(diff[2], WithinRel(delta, 1.0e-6));
  REQUIRE_THAT(diff[3], WithinRel(delta, 1.0e-6));
}

// =========================================================================
// 12. Edge cases
// =========================================================================

TEST_CASE("Vec3 division by very small scalar", "[Vec3][edge]") {
  Vec3D v(1.0, 1.0, 1.0);
  double tiny = 1.0e-300;
  Vec3D r = v / tiny;
  // Result should be large but finite.
  REQUIRE_FALSE(r.Nan());
  REQUIRE_THAT(r[1], WithinRel(1.0 / tiny, TOL));
}

TEST_CASE("Vec3 operations with infinity", "[Vec3][edge]") {
  double inf = std::numeric_limits<double>::infinity();
  Vec3D v(inf, 0.0, 0.0);
  // Inf is not NaN.
  REQUIRE_FALSE(v.Nan());
  // But the vector is clearly not zero.
  REQUIRE_FALSE(v.IsZero());
}

TEST_CASE("Vec3 negative zero components", "[Vec3][edge]") {
  // IEEE 754: -0.0 == 0.0; the vector should still be considered zero.
  Vec3D v(-0.0, -0.0, -0.0);
  REQUIRE(v.IsZero());
  REQUIRE(v.Sqr() == 0.0);
}

TEST_CASE("Vec3 addition is commutative and associative", "[Vec3][algebraic]") {
  Vec3D a(1.23, -4.56, 7.89);
  Vec3D b(-0.12, 3.45, -6.78);
  Vec3D c(9.01, -2.34, 5.67);

  // Commutativity: a + b == b + a
  REQUIRE(ApproxEqual(a + b, b + a, TOL));

  // Associativity: (a + b) + c == a + (b + c)
  // May differ by ULP due to FP ordering, so use relaxed tolerance.
  REQUIRE(ApproxEqual((a + b) + c, a + (b + c), TOL*10.0));
}

TEST_CASE("Vec3 scalar multiplication distributes over addition", "[Vec3][algebraic]") {
  Vec3D a(1.0, 2.0, 3.0);
  Vec3D b(4.0, -1.0, 2.0);
  double s = 2.5;
  REQUIRE(ApproxEqual(s * (a + b), s * a + s * b, TOL));
}

TEST_CASE("Vec3 zero vector is additive identity", "[Vec3][algebraic]") {
  Vec3D a(3.14, 2.72, 1.41);
  Vec3D zero;
  REQUIRE(ApproxEqual(a + zero, a));
  REQUIRE(ApproxEqual(zero + a, a));
}

TEST_CASE("Vec3 additive inverse", "[Vec3][algebraic]") {
  Vec3D a(3.14, 2.72, 1.41);
  Vec3D sum = a + (-a);
  REQUIRE_THAT(sum.Abs(), WithinAbs(0.0, TOL));
}

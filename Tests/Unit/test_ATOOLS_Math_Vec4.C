#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Vec3.H"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <cmath>
#include <limits>
#include <sstream>

#include "ATOOLS/Math/Vec4.H"

// ---------------------------------------------------------------------------
// Convenience aliases and constants
// ---------------------------------------------------------------------------
using Vec4d = ATOOLS::Vec4<double>;
using Vec3d = ATOOLS::Vec3<double>;

// Relative-tolerance comparisons via Catch2::Approx
static constexpr double kTightRel = 1e-12; // for algebraically exact results
static constexpr double kLooseRel = 1e-9;  // where intermediate cancellation occurs

// Helper: make an on-shell massive four-momentum p = (E, px, py, pz)
// with E = sqrt(m^2 + px^2 + py^2 + pz^2).
Vec4d makeOnShell(double m, double px, double py, double pz) {
    double E = std::sqrt(m*m + px*px + py*py + pz*pz);
    return Vec4d(E, px, py, pz);
}

// Helper: massless (lightlike) four-momentum along a given 3-direction.
Vec4d makeMassless(double px, double py, double pz) {
    double E = std::sqrt(px*px + py*py + pz*pz);
    return Vec4d(E, px, py, pz);
}

// ---------------------------------------------------------------------------
// ==========================================================================
// TEST SECTION 1 — Construction and component access
// ==========================================================================
// ---------------------------------------------------------------------------

TEST_CASE("Vec4 construction", "[vec4][construction]") {

    SECTION("Default constructor initialises all components to zero") {
        Vec4d v;
        // All four components must be exactly zero
        CHECK(v[0] == 0.0);
        CHECK(v[1] == 0.0);
        CHECK(v[2] == 0.0);
        CHECK(v[3] == 0.0);
    }

    SECTION("Component constructor stores values in correct order") {
        // Order convention: (E, px, py, pz) i.e. index 0 = energy
        Vec4d v(10.0, 1.0, 2.0, 3.0);
        CHECK(v[0] == 10.0);
        CHECK(v[1] ==  1.0);
        CHECK(v[2] ==  2.0);
        CHECK(v[3] ==  3.0);
        CHECK(v.E()  == 10.0);
    }

    SECTION("Copy constructor (same type)") {
        Vec4d a(5.0, 1.0, 2.0, 3.0);
        Vec4d b(a);
        CHECK(b[0] == a[0]);
        CHECK(b[1] == a[1]);
        CHECK(b[2] == a[2]);
        CHECK(b[3] == a[3]);
    }

    SECTION("Non-const operator[] allows mutation") {
        Vec4d v(0.0, 0.0, 0.0, 0.0);
        v[0] = 7.0;
        v[3] = -3.0;
        CHECK(v[0] ==  7.0);
        CHECK(v[3] == -3.0);
    }

    SECTION("Vec4(E, Vec3) constructor maps Vec3 components correctly") {
        Vec3d v3(1.0, 2.0, 3.0);  // spatial components: x=1, y=2, z=3
        Vec4d v(10.0, v3);

        CHECK(v[0] == Catch::Approx(10.0));
        CHECK(v[1] == Catch::Approx(v3[1]));
        CHECK(v[2] == Catch::Approx(v3[2]));
        CHECK(v[3] == Catch::Approx(v3[3]));
    }
}

// ---------------------------------------------------------------------------
// ==========================================================================
// TEST SECTION 2 — Minkowski metric and arithmetic operators
// ==========================================================================
// ---------------------------------------------------------------------------

TEST_CASE("Vec4 Minkowski inner product", "[vec4][metric]") {

    SECTION("Metric signature is (+,-,-,-)") {
        // A purely timelike unit vector dotted with itself must give +1.
        Vec4d t(1.0, 0.0, 0.0, 0.0);
        CHECK((t*t) == Catch::Approx(1.0).epsilon(kTightRel));

        // A purely spacelike unit vector dotted with itself must give -1.
        Vec4d sx(0.0, 1.0, 0.0, 0.0);
        Vec4d sy(0.0, 0.0, 1.0, 0.0);
        Vec4d sz(0.0, 0.0, 0.0, 1.0);
        CHECK((sx*sx) == Catch::Approx(-1.0).epsilon(kTightRel));
        CHECK((sy*sy) == Catch::Approx(-1.0).epsilon(kTightRel));
        CHECK((sz*sz) == Catch::Approx(-1.0).epsilon(kTightRel));
    }

    SECTION("Inner product is symmetric") {
        Vec4d p(5.0, 1.0, 2.0, 3.0);
        Vec4d q(4.0, -1.0, 0.5, 2.0);
        CHECK((p*q) == Catch::Approx(q*p).epsilon(kTightRel));
    }

    SECTION("Inner product is bilinear (additivity)") {
        Vec4d p(3.0, 1.0, 0.0, 0.0);
        Vec4d q(2.0, 0.0, 1.0, 0.0);
        Vec4d r(1.0, 0.0, 0.0, 1.0);
        double alpha = 2.5;
        // (alpha*p + q) * r == alpha*(p*r) + (q*r)
        Vec4d lhs_v = (p*alpha) + q;
        double lhs = lhs_v * r;
        double rhs = alpha*(p*r) + (q*r);
        CHECK(lhs == Catch::Approx(rhs).epsilon(kTightRel));
    }

    SECTION("Massless (lightlike) vector has zero self-product") {
        // A photon with |p| along z: p = (E, 0, 0, E)
        Vec4d photon(5.0, 0.0, 0.0, 5.0);
        CHECK((photon*photon) == Catch::Approx(0.0).margin(1e-15));
    }

    SECTION("Abs2() is consistent with operator*") {
        // Abs2() uses a numerically improved factored form (E+pz)(E-pz).
        // It must agree with the direct metric contraction.
        Vec4d p(10.0, 3.0, 4.0, 0.0);
        double via_operator   = p * p;
        double via_abs2       = p.Abs2();
        CHECK(via_abs2 == Catch::Approx(via_operator).epsilon(kTightRel));
    }

    SECTION("Orthogonal basis vectors are mutually orthogonal") {
        Vec4d t(1.0, 0.0, 0.0, 0.0);
        Vec4d sx(0.0, 1.0, 0.0, 0.0);
        Vec4d sy(0.0, 0.0, 1.0, 0.0);
        Vec4d sz(0.0, 0.0, 0.0, 1.0);
        CHECK((t*sx)  == Catch::Approx(0.0).margin(1e-15));
        CHECK((t*sy)  == Catch::Approx(0.0).margin(1e-15));
        CHECK((t*sz)  == Catch::Approx(0.0).margin(1e-15));
        CHECK((sx*sy) == Catch::Approx(0.0).margin(1e-15));
        CHECK((sx*sz) == Catch::Approx(0.0).margin(1e-15));
        CHECK((sy*sz) == Catch::Approx(0.0).margin(1e-15));
    }
}

TEST_CASE("Vec4 arithmetic operators", "[vec4][operators]") {

    SECTION("Addition is commutative") {
        Vec4d a(1.0, 2.0, 3.0, 4.0);
        Vec4d b(5.0, 6.0, 7.0, 8.0);
        Vec4d apb = a + b;
        Vec4d bpa = b + a;
        for (int i = 0; i < 4; ++i)
            CHECK(apb[i] == Catch::Approx(bpa[i]).epsilon(kTightRel));
    }

    SECTION("Addition and subtraction are inverse operations") {
        Vec4d a(3.0, -1.0, 2.5, 0.0);
        Vec4d b(1.0,  2.0, -1.0, 4.0);
        Vec4d result = (a + b) - b;
        for (int i = 0; i < 4; ++i)
            CHECK(result[i] == Catch::Approx(a[i]).epsilon(kTightRel));
    }

    SECTION("Scalar multiplication distributes over addition") {
        Vec4d a(1.0, 2.0, 3.0, 4.0);
        Vec4d b(5.0, 6.0, 7.0, 8.0);
        double s = 3.0;
        Vec4d lhs = (a + b) * s;
        Vec4d rhs = a*s + b*s;
        for (int i = 0; i < 4; ++i)
            CHECK(lhs[i] == Catch::Approx(rhs[i]).epsilon(kTightRel));
    }

    SECTION("Division is consistent with multiplication by reciprocal") {
        Vec4d a(4.0, 8.0, 12.0, 16.0);
        double s = 4.0;
        Vec4d via_div  = a / s;
        Vec4d via_mult = a * (1.0/s);
        for (int i = 0; i < 4; ++i)
            CHECK(via_div[i] == Catch::Approx(via_mult[i]).epsilon(kTightRel));
    }

    SECTION("Unary negation negates all components") {
        Vec4d a(1.0, -2.0, 3.0, -4.0);
        Vec4d neg = -a;
        CHECK(neg[0] == Catch::Approx(-1.0));
        CHECK(neg[1] == Catch::Approx( 2.0));
        CHECK(neg[2] == Catch::Approx(-3.0));
        CHECK(neg[3] == Catch::Approx( 4.0));
    }

    SECTION("Compound assignment += and -= are consistent with + and -") {
        Vec4d a(1.0, 2.0, 3.0, 4.0);
        Vec4d b(5.0, 6.0, 7.0, 8.0);
        Vec4d c = a;
        c += b;
        Vec4d expected = a + b;
        for (int i = 0; i < 4; ++i)
            CHECK(c[i] == Catch::Approx(expected[i]).epsilon(kTightRel));
        c -= b;
        for (int i = 0; i < 4; ++i)
            CHECK(c[i] == Catch::Approx(a[i]).epsilon(kTightRel));
    }

    SECTION("Scalar commutativity: s*v == v*s") {
        Vec4d v(2.0, 3.0, 4.0, 5.0);
        double s = 7.0;
        Vec4d lhs = s * v;
        Vec4d rhs = v * s;
        for (int i = 0; i < 4; ++i)
            CHECK(lhs[i] == Catch::Approx(rhs[i]).epsilon(kTightRel));
    }
}

// ---------------------------------------------------------------------------
// ==========================================================================
// TEST SECTION 3 — Physics invariants and kinematic quantities
// ==========================================================================
// ---------------------------------------------------------------------------

TEST_CASE("Vec4 mass and Lorentz invariants", "[vec4][physics]") {

    SECTION("Massive particle: Abs2() equals m^2") {
        // Proton at rest: p = (m, 0, 0, 0), Abs2 = m^2
        double m = 0.938; // GeV, proton mass
        Vec4d proton_rest(m, 0.0, 0.0, 0.0);
        CHECK(proton_rest.Abs2() == Catch::Approx(m*m).epsilon(kTightRel));
        CHECK(proton_rest.Mass() == Catch::Approx(m).epsilon(kTightRel));
    }

    SECTION("Massive particle with spatial momentum: Abs2() still m^2") {
        double m = 0.938;
        Vec4d p = makeOnShell(m, 1.0, 2.0, 3.0);
        // Physics: E^2 - |p|^2 = m^2 regardless of the boost
        CHECK(p.Abs2() == Catch::Approx(m*m).epsilon(kTightRel));
        CHECK(p.Mass() == Catch::Approx(m).epsilon(kTightRel));
    }

    SECTION("Massless particle: Abs2() is zero") {
        Vec4d photon = makeMassless(3.0, 4.0, 0.0);  // |p| = 5
        CHECK(photon.Abs2() == Catch::Approx(0.0).margin(1e-14));
        CHECK(photon.E()    == Catch::Approx(5.0).epsilon(kTightRel));
    }

    SECTION("Mass() is always non-negative via Abs() of Abs2()") {
        // Spacelike: Abs2() < 0 (unphysical but can occur in intermediate steps)
        Vec4d spacelike(0.0, 1.0, 0.0, 0.0);
        // Mass() = sqrt(|Abs2()|) = sqrt(1) = 1 (not imaginary!)
        CHECK(spacelike.Mass() == Catch::Approx(1.0).epsilon(kTightRel));
    }

    SECTION("Four-momentum conservation: invariant mass of sum") {
        // Two-body decay: parent -> child1 + child2.
        // Total four-momentum must have correct invariant mass.
        double mParent = 91.2;  // Z boson (GeV)
        // Z at rest decays to two massless back-to-back particles
        double E = mParent / 2.0;
        Vec4d child1(E,  0.0, 0.0,  E);
        Vec4d child2(E,  0.0, 0.0, -E);
        Vec4d total = child1 + child2;
        CHECK(total.Abs2() == Catch::Approx(mParent*mParent).epsilon(kTightRel));
        CHECK(total.Mass() == Catch::Approx(mParent).epsilon(kTightRel));
    }

    SECTION("Invariant mass squared is Lorentz-invariant under analytic boost") {
        // Verify that Abs2 gives same result for a particle at rest
        // and after an analytic longitudinal boost (gamma=2, beta=sqrt(3)/2).
        double m = 5.0;
        Vec4d pRest(m, 0.0, 0.0, 0.0);

        double gamma = 2.0;
        double beta  = std::sqrt(3.0) / 2.0;  // gamma*beta = sqrt(3)
        double E_boosted  = gamma * m;
        double pz_boosted = gamma * beta * m;
        Vec4d pBoosted(E_boosted, 0.0, 0.0, pz_boosted);

        CHECK(pRest.Abs2()   == Catch::Approx(m*m).epsilon(kTightRel));
        CHECK(pBoosted.Abs2() == Catch::Approx(m*m).epsilon(kTightRel));
    }
}

TEST_CASE("Vec4 transverse and light-cone kinematics", "[vec4][physics]") {

    SECTION("PPerp2 and PPerp agree: PPerp == sqrt(PPerp2)") {
        Vec4d p(10.0, 3.0, 4.0, 5.0);
        double pt2 = p.PPerp2();
        double pt  = p.PPerp();
        CHECK(pt*pt == Catch::Approx(pt2).epsilon(kTightRel));
        CHECK(pt    == Catch::Approx(5.0).epsilon(kTightRel)); // sqrt(9+16)=5
    }

    SECTION("PPerp is zero for a purely longitudinal momentum") {
        Vec4d p(10.0, 0.0, 0.0, 7.0);
        CHECK(p.PPerp()  == Catch::Approx(0.0).margin(1e-15));
        CHECK(p.PPerp2() == Catch::Approx(0.0).margin(1e-15));
    }

    SECTION("PSpat2 and PSpat agree: PSpat == sqrt(PSpat2)") {
        Vec4d p(10.0, 2.0, 3.0, 6.0);
        double ps2 = p.PSpat2();
        double ps  = p.PSpat();
        CHECK(ps*ps == Catch::Approx(ps2).epsilon(kTightRel));
        CHECK(ps    == Catch::Approx(7.0).epsilon(kTightRel)); // sqrt(4+9+36)=7
    }

    SECTION("Light-cone components: PPlus and PMinus") {
        Vec4d p(5.0, 1.0, 2.0, 3.0);
        CHECK(p.PPlus()  == Catch::Approx(p[0] + p[3]).epsilon(kTightRel));
        CHECK(p.PMinus() == Catch::Approx(p[0] - p[3]).epsilon(kTightRel));
    }

    SECTION("MPerp2: MPerp2 = E^2 - pz^2 for zero transverse momentum") {
        // MPerp2 = (E+pz)(E-pz) when px=py=0
        Vec4d p(5.0, 0.0, 0.0, 3.0);
        double expected = (5.0+3.0)*(5.0-3.0); // = 16
        CHECK(p.MPerp2() == Catch::Approx(expected).epsilon(kTightRel));
        CHECK(p.MPerp()  == Catch::Approx(4.0).epsilon(kTightRel));
    }

    SECTION("Rapidity Y() is well-defined for timelike momentum with pz != E") {
        // p = (E, 0, 0, pz), Y = 0.5*ln((E+pz)/(E-pz))
        double E = 10.0, pz = 3.0;
        Vec4d p(E, 0.0, 0.0, pz);
        double expected_Y = 0.5 * std::log((E + pz)/(E - pz));
        CHECK(p.Y() == Catch::Approx(expected_Y).epsilon(kTightRel));
    }

    SECTION("Rapidity is zero for a particle at rest") {
        Vec4d p(5.0, 0.0, 0.0, 0.0);
        CHECK(p.Y() == Catch::Approx(0.0).margin(1e-15));
    }

    SECTION("RelAbs2 returns zero when Abs2 is zero (no division by zero)") {
        // Guard for the special case in RelAbs2(): if E=0 we could divide by 0.
        // The implementation guards with an explicit check on abs2, not on E.
        Vec4d zero;
        CHECK(zero.RelAbs2() == 0.0);

        // Massless photon: Abs2=0 -> RelAbs2 should also be 0
        Vec4d photon(5.0, 0.0, 0.0, 5.0);
        CHECK(photon.RelAbs2() == Catch::Approx(0.0).margin(1e-15));
    }

    SECTION("RelAbs2 is dimensionless and bounded for massive particles") {
        // For a particle with m << E (ultra-relativistic), RelAbs2 ~ (m/E)^2 << 1
        double m = 1.0, E = 1000.0;
        double pz = std::sqrt(E*E - m*m);
        Vec4d p(E, 0.0, 0.0, pz);
        double rel = p.RelAbs2();
        // RelAbs2 = m^2/E^2 = 1e-6
        CHECK(rel > 0.0);
        CHECK(rel < 1.0);
        CHECK(rel == Catch::Approx(m*m/(E*E)).epsilon(kLooseRel));
    }

    SECTION("P() delegates to PSpat()") {
        Vec4d p(10.0, 3.0, 4.0, 0.0);
        CHECK(p.P() == Catch::Approx(p.PSpat()).epsilon(kTightRel));
    }
}

// ---------------------------------------------------------------------------
// ==========================================================================
// TEST SECTION 4 — Light-cone decomposition completeness
// ==========================================================================
// ---------------------------------------------------------------------------

TEST_CASE("Vec4 light-cone decomposition", "[vec4][physics][decomposition]") {

    // A generic massive four-momentum
    Vec4d p = makeOnShell(1.0, 2.0, 3.0, 4.0);

    SECTION("Plus() + Minus() + Perp() reconstructs the full vector") {
        Vec4d pPlus  = p.Plus();
        Vec4d pMinus = p.Minus();
        Vec4d pPerp  = p.Perp();
        Vec4d reconstructed = pPlus + pMinus + pPerp;
        for (int i = 0; i < 4; ++i)
            CHECK(reconstructed[i] == Catch::Approx(p[i]).epsilon(kTightRel));
    }

    SECTION("Plus() is lightlike along +z direction") {
        // Plus() = (p+/2)(1,0,0,1), which should be lightlike.
        Vec4d pp = p.Plus();
        CHECK(pp.Abs2() == Catch::Approx(0.0).margin(1e-14));
    }

    SECTION("Minus() is lightlike along -z direction") {
        Vec4d pm = p.Minus();
        CHECK(pm.Abs2() == Catch::Approx(0.0).margin(1e-14));
    }

    SECTION("Perp() has zero energy and zero pz components") {
        Vec4d pe = p.Perp();
        CHECK(pe[0] == Catch::Approx(0.0).margin(1e-15));
        CHECK(pe[3] == Catch::Approx(0.0).margin(1e-15));
    }

    SECTION("Long() + Perp() reconstructs the full vector") {
        Vec4d pLong  = p.Long();
        Vec4d pPerp  = p.Perp();
        Vec4d reconstructed = pLong + pPerp;
        for (int i = 0; i < 4; ++i)
            CHECK(reconstructed[i] == Catch::Approx(p[i]).epsilon(kTightRel));
    }

    SECTION("Long() has zero transverse components") {
        Vec4d pl = p.Long();
        CHECK(pl[1] == Catch::Approx(0.0).margin(1e-15));
        CHECK(pl[2] == Catch::Approx(0.0).margin(1e-15));
    }
}

// ---------------------------------------------------------------------------
// ==========================================================================
// TEST SECTION 5 — Numerical robustness
// ==========================================================================
// ---------------------------------------------------------------------------

TEST_CASE("Vec4 numerical robustness of Abs2()", "[vec4][numerical]") {

    SECTION("Abs2() factored form vs naive form for nearly-lightlike vector") {
        // For p = (E, 0, 0, pz) with E \approx pz (nearly massless),
        // naive form E^2 - pz^2 suffers catastrophic cancellation.
        // The factored form (E+pz)(E-pz) used in Abs2() avoids this.
        //
        // We construct a vector with a known small mass and verify Abs2
        // is close to m^2 even when m << E.

        double m  = 1e-2;
        double pz = 1e2;
        double E  = std::sqrt(m*m + pz*pz);

        Vec4d p(E, 0.0, 0.0, pz);

        // The factored form should give m^2 with reasonable relative accuracy.
        double abs2 = p.Abs2();
        CHECK(abs2 == Catch::Approx(m*m).epsilon(kLooseRel));

        // Sanity: naive form E^2-pz^2 at these scales may not meet the tighter
        // tolerance — confirming why the factored form is needed.
        double naive = E*E - pz*pz;
        // We do not assert naive is inaccurate (compiler-dependent), but document:
        INFO("Naive E^2-pz^2 = " << naive << ", factored Abs2 = " << abs2
             << ", exact m^2 = " << m*m);
    }

    SECTION("Abs2() is stable for a massive vector at moderate boost") {
        double m  = 0.938;
        double pz = 1.0;
        Vec4d p   = makeOnShell(m, 0.0, 0.0, pz);
        CHECK(p.Abs2() == Catch::Approx(m*m).epsilon(kTightRel));
    }

    SECTION("Sum of back-to-back massless momenta gives correct invariant mass") {
        // Two photons: p1=(E,0,0,E), p2=(E,0,0,-E). Sum: (2E,0,0,0)
        // Invariant mass squared = 4E^2.
        double E = 45.6; // GeV (LEP beam energy)
        Vec4d p1(E, 0.0, 0.0,  E);
        Vec4d p2(E, 0.0, 0.0, -E);
        Vec4d sum = p1 + p2;
        CHECK(sum.Abs2() == Catch::Approx(4.0*E*E).epsilon(kTightRel));
        CHECK(sum.Mass() == Catch::Approx(2.0*E).epsilon(kTightRel));
    }

    SECTION("Abs() returns positive value for timelike vector") {
        Vec4d p(10.0, 1.0, 1.0, 1.0);
        double abs_val = p.Abs();
        CHECK(abs_val > 0.0);
        CHECK(abs_val*abs_val == Catch::Approx(p.Abs2()).epsilon(kTightRel));
    }

    SECTION("Mandelstam-s + Mandelstam-u + Mandelstam-t = sum of mass^2") {
        // For 2->2 scattering p1+p2->p3+p4 with all massless particles:
        // s + t + u = 0. This tests that Minkowski sums and inner products
        // are numerically consistent.
        Vec4d p1 = makeMassless(0.0,  0.0,  50.0);
        Vec4d p2 = makeMassless(0.0,  0.0, -50.0);
        // 2->2 with arbitrary angle theta
        double E  = 50.0;
        double th = M_PI / 4.0;
        Vec4d p3(E,  E*std::sin(th), 0.0,  E*std::cos(th));
        Vec4d p4(E, -E*std::sin(th), 0.0, -E*std::cos(th));

        double s = (p1+p2).Abs2();
        double t = (p1-p3).Abs2();
        double u = (p1-p4).Abs2();
        // All external masses are zero => s + t + u = 0
        CHECK(s + t + u == Catch::Approx(0.0).margin(1e-10));
    }
}

// ---------------------------------------------------------------------------
// ==========================================================================
// TEST SECTION 6 — Utility predicates: IsZero(), Nan()
// ==========================================================================
// ---------------------------------------------------------------------------

TEST_CASE("Vec4 utility predicates", "[vec4][utilities]") {

    SECTION("IsZero() returns true for the default-constructed zero vector") {
        Vec4d z;
        CHECK(z.IsZero());
    }

    SECTION("IsZero() returns false if any component is nonzero") {
        Vec4d v(0.0, 0.0, 0.0, ATOOLS::Accu());
        CHECK_FALSE(v.IsZero());
    }

    SECTION("Nan() returns false for a well-formed vector") {
        Vec4d p = makeOnShell(0.938, 1.0, 2.0, 3.0);
        CHECK_FALSE(p.Nan());
    }

    SECTION("Nan() returns true if any component is NaN") {
        double nan = std::numeric_limits<double>::quiet_NaN();
        {
            Vec4d v(nan, 0.0, 0.0, 0.0);
            CHECK(v.Nan());
        }
        {
            Vec4d v(1.0, 0.0, nan, 0.0);
            CHECK(v.Nan());
        }
    }
}

// ---------------------------------------------------------------------------
// ==========================================================================
// TEST SECTION 7 — Output stream operator
// ==========================================================================
// ---------------------------------------------------------------------------

TEST_CASE("Vec4 stream output", "[vec4][io]") {

    SECTION("operator<< produces a parseable string in (x0,x1,x2,x3) format") {
        Vec4d v(1.0, 2.0, 3.0, 4.0);
        std::ostringstream oss;
        oss << v;
        std::string s = oss.str();
        // Must start with '(' and end with ')'
        CHECK(s.front() == '(');
        CHECK(s.back()  == ')');
        // Must contain the energy component
        CHECK(s.find("1") != std::string::npos);
    }
}

// ---------------------------------------------------------------------------
// ==========================================================================
// TEST SECTION 8 — Free function: cross()
// ==========================================================================
// ---------------------------------------------------------------------------

TEST_CASE("Vec4 free function: cross()", "[vec4][cross]") {

    SECTION("cross() result is orthogonal to all three input vectors") {
        // The 4D generalised cross product cross(q,r,s) is orthogonal
        // to q, r, and s under the Minkowski metric.
        Vec4d q(1.0, 0.0, 0.0, 0.0);
        Vec4d r(0.0, 1.0, 0.0, 0.0);
        Vec4d s(0.0, 0.0, 1.0, 0.0);
        Vec4d c = ATOOLS::cross(q, r, s);
        CHECK((c*q) == Catch::Approx(0.0).margin(1e-14));
        CHECK((c*r) == Catch::Approx(0.0).margin(1e-14));
        CHECK((c*s) == Catch::Approx(0.0).margin(1e-14));
    }

    SECTION("cross() with two identical inputs gives zero vector") {
        // Antisymmetry: cross(q,r,r) must vanish componentwise.
        Vec4d q(1.0, 0.0, 0.0, 0.0);
        Vec4d r(0.0, 1.0, 2.0, 3.0);
        Vec4d c = ATOOLS::cross(q, r, r);
        for (int i = 0; i < 4; ++i)
            CHECK(c[i] == Catch::Approx(0.0).margin(1e-14));
    }

    SECTION("cross() is linear in each argument") {
        Vec4d q(1.0, 0.0, 0.0, 0.0);
        Vec4d r(0.0, 1.0, 0.0, 0.0);
        Vec4d s1(0.0, 0.0, 1.0, 0.0);
        Vec4d s2(0.0, 0.0, 0.0, 1.0);
        double alpha = 2.0;
        Vec4d c1     = ATOOLS::cross(q, r, s1);
        Vec4d c2     = ATOOLS::cross(q, r, s2);
        Vec4d cSum   = ATOOLS::cross(q, r, s1*alpha + s2);
        Vec4d cExpect = c1*alpha + c2;
        for (int i = 0; i < 4; ++i)
            CHECK(cSum[i] == Catch::Approx(cExpect[i]).epsilon(kTightRel));
    }

    SECTION("cross() basis result: cross(e0,e1,e2) should equal e3 (up to sign/metric)") {
        // cross(e_0, e_1, e_2): the remaining basis direction is e_3.
        // The exact sign and metric factor can be read off the formula.
        Vec4d e0(1.0, 0.0, 0.0, 0.0);
        Vec4d e1(0.0, 1.0, 0.0, 0.0);
        Vec4d e2(0.0, 0.0, 1.0, 0.0);
        Vec4d c = ATOOLS::cross(e0, e1, e2);
        // Orthogonal to e0, e1, e2 — verify orthogonality
        CHECK((c*e0) == Catch::Approx(0.0).margin(1e-14));
        CHECK((c*e1) == Catch::Approx(0.0).margin(1e-14));
        CHECK((c*e2) == Catch::Approx(0.0).margin(1e-14));
        // Only one component should be nonzero (e3 direction)
        CHECK(c[0] == Catch::Approx(0.0).margin(1e-14));
        CHECK(c[1] == Catch::Approx(0.0).margin(1e-14));
        CHECK(c[2] == Catch::Approx(0.0).margin(1e-14));
        CHECK(std::abs(c[3]) > 1e-10); // nonzero
    }
}

// ---------------------------------------------------------------------------
// ==========================================================================
// TEST SECTION 9 — Free function: CosPhi() (acoplanarity angle)
// ==========================================================================
// ---------------------------------------------------------------------------

TEST_CASE("Vec4 free function: CosPhi()", "[vec4][cosphi]") {

    SECTION("CosPhi() result is in [-1, 1]") {
        // This is a necessary condition; violations signal numerical issues.
        Vec4d pi = makeOnShell(0.0, 10.0, 0.0, 0.0);
        Vec4d pj = makeOnShell(0.0,  0.0, 10.0, 0.0);
        Vec4d pk = makeOnShell(0.0,  0.0, 0.0, 10.0);
        Vec4d pl = makeOnShell(0.0,  5.0, 5.0, 0.0);
        double cp = ATOOLS::CosPhi(pi, pj, pk, pl);
        CHECK(cp >= -1.0 - 1e-10);
        CHECK(cp <=  1.0 + 1e-10);
    }

    SECTION("CosPhi() is symmetric under swap of i<->j and simultaneous k<->l") {
        // The formula has a symmetry: swapping (i,j) and (k,l) simultaneously
        // should leave CosPhi invariant (both numerator and denominator are
        // symmetric under this combined exchange).
        Vec4d pi = makeOnShell(0.0, 4.0, 3.0, 0.0);
        Vec4d pj = makeOnShell(0.0, 0.0, 5.0, 0.0);
        Vec4d pk = makeOnShell(0.0, 2.0, 0.0, 4.0);
        Vec4d pl = makeOnShell(0.0, 1.0, 2.0, 3.0);
        double cp1 = ATOOLS::CosPhi(pi, pj, pk, pl);
        double cp2 = ATOOLS::CosPhi(pj, pi, pl, pk);
        CHECK(cp1 == Catch::Approx(cp2).epsilon(kLooseRel));
    }
}

// ---------------------------------------------------------------------------
// ==========================================================================
// TEST SECTION 10 — Edge cases and boundary conditions
// ==========================================================================
// ---------------------------------------------------------------------------

TEST_CASE("Vec4 edge cases", "[vec4][edge]") {

    SECTION("Zero vector: all derived quantities are zero or well-defined") {
        Vec4d z;
        CHECK(z.Abs2()   == 0.0);
        CHECK(z.PPerp2() == 0.0);
        CHECK(z.PSpat2() == 0.0);
        CHECK(z.PPlus()  == 0.0);
        CHECK(z.PMinus() == 0.0);
        CHECK(z.IsZero());
    }

    SECTION("Purely timelike vector: spatial momentum quantities are zero") {
        Vec4d p(10.0, 0.0, 0.0, 0.0);
        CHECK(p.PPerp()  == Catch::Approx(0.0).margin(1e-15));
        CHECK(p.PSpat()  == Catch::Approx(0.0).margin(1e-15));
        CHECK(p.Mass()   == Catch::Approx(10.0).epsilon(kTightRel));
    }

    SECTION("Purely transverse momentum: pz=0 simplifications") {
        Vec4d p(10.0, 3.0, 4.0, 0.0);
        CHECK(p.PPerp() == Catch::Approx(5.0).epsilon(kTightRel));
        CHECK(p.PPlus() == Catch::Approx(p.PMinus()).epsilon(kTightRel));
        // MPerp2 = (E+pz)(E-pz) = E^2 when pz=0
        CHECK(p.MPerp2() == Catch::Approx(100.0).epsilon(kTightRel));
    }

    SECTION("Exactly back-to-back massless momenta: s = 4E^2") {
        double E = 100.0;
        Vec4d p1(E, 0.0, 0.0,  E);
        Vec4d p2(E, 0.0, 0.0, -E);
        CHECK((p1+p2).Abs2() == Catch::Approx(4.0*E*E).epsilon(kTightRel));
    }

    SECTION("Very large energy values do not produce NaN or Inf") {
        double bigE = 1e15;
        Vec4d p = makeOnShell(0.938, 0.0, 0.0, bigE);
        CHECK_FALSE(p.Nan());
        CHECK(std::isfinite(p.Abs2()));
        CHECK(std::isfinite(p.Mass()));
    }

    SECTION("EPerp2() is dangerous when PSpat=0: documenting the division") {
        // EPerp2() = E^2 * PPerp2() / PSpat2().
        // If PSpat=0 this is 0/0 -> NaN.  This test documents the current
        // behaviour — a future fix might add a guard.
        Vec4d p(5.0, 0.0, 0.0, 0.0);  // purely timelike, PSpat=0
        double eperp2 = p.EPerp2();    // Will likely be NaN or 0*inf
        // Document (not assert correctness): caller must ensure PSpat != 0
        INFO("EPerp2() with PSpat=0 gives: " << eperp2);
        // For now we only verify the finite case works correctly:
        Vec4d q(10.0, 3.0, 4.0, 5.0);
        double expected_eperp2 = q[0]*q[0]*q.PPerp2()/q.PSpat2();
        CHECK(q.EPerp2() == Catch::Approx(expected_eperp2).epsilon(kTightRel));
    }
}

// Cancellation-free kinematics of the YFS Comix currents
// (METOOLS/Currents/YFS_Stable_Current.H). A photon nearly collinear to a
// massive leg makes p.k a small difference of large numbers in double; the
// exact light-cone dot product and the invariant table must not lose it.
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "METOOLS/Currents/YFS_Stable_Current.H"

#include <cmath>
#include <vector>

using ATOOLS::DDouble;
using ATOOLS::Vec4D;
using METOOLS::Invariant_Table;
using METOOLS::YFSExactLCDot;
typedef ATOOLS::Vec4<DDouble> Vec4DD;

namespace {

  Vec4DD Wide(const Vec4D &p)
  { return Vec4DD(DDouble(p[0]), DDouble(p[1]), DDouble(p[2]), DDouble(p[3])); }

  DDouble WideDot(const Vec4DD &p, const Vec4DD &q)
  { return p[0]*q[0] - p[1]*q[1] - p[2]*q[2] - p[3]*q[3]; }

  Vec4D OnShell(double m, double px, double py, double pz)
  { return Vec4D(std::sqrt(m*m + px*px + py*py + pz*pz), px, py, pz); }

}

TEST_CASE("Exact light-cone dot product for a collinear photon",
          "[METOOLS][YFS_Stable_Current]") {
  // electron along +z, photon at angle th to it
  const double m(0.000511), P(45.6), w(1.0), th(1e-7);
  const double E(std::sqrt(P*P + m*m));
  const Vec4D p(E, 0., 0., P), k(w, w*std::sin(th), 0., w*std::cos(th));
  // reference without cancellation: p.k = w (E - P cos th)
  //   = w (m^2/(E+P) + 2 P sin^2(th/2))
  const double s(std::sin(th/2.));
  const double ref(w*(m*m/(E + P) + 2.*P*s*s));
  const double exact((double)YFSExactLCDot(Wide(p), DDouble(m*m),
                                           Wide(k), DDouble(0.0)));
  CHECK(exact == Catch::Approx(ref).epsilon(1e-12));
  // and the naive double product really does fail here
  const double naive(p*k);
  CHECK(std::abs(naive/ref - 1.) > 1e-9);
}

TEST_CASE("Exact light-cone dot product agrees for generic momenta",
          "[METOOLS][YFS_Stable_Current]") {
  const Vec4D p(OnShell(0.105, 3., -1., 40.)), q(OnShell(0.105, -2., 1.5, -38.));
  const double exact((double)YFSExactLCDot(Wide(p), DDouble(0.105*0.105),
                                           Wide(q), DDouble(0.105*0.105)));
  CHECK(exact == Catch::Approx(p*q).epsilon(1e-12));
}

namespace {

  // four outgoing legs that sum to zero: two beams crossed to outgoing and a
  // massive final-state pair back to back in the CMS, all on shell
  Invariant_Table BalancedTable(std::vector<Vec4D> &legs, std::vector<double> &mass)
  {
    const double m(0.105), me(0.000511), E(45.6);
    const double pe(std::sqrt(E*E - me*me));
    const double pf(std::sqrt(E*E - m*m));
    const double n(std::sqrt(10.*10. + 5.*5. + 20.*20.));
    const Vec4D f3(E, pf*10./n, -pf*5./n, pf*20./n), f4(E, -f3[1], -f3[2], -f3[3]);
    legs = {Vec4D(-E, 0., 0., -pe), Vec4D(-E, 0., 0., pe), f3, f4};
    mass = {me, me, m, m};
    Invariant_Table t;
    t.m_n = legs.size();
    t.m_m2.assign(t.m_n, DDouble(0.0));
    t.m_s.assign(t.m_n*t.m_n, DDouble(0.0));
    t.m_p = legs;
    t.m_R = Vec4D();
    for (size_t a(0); a < t.m_n; ++a) {
      t.m_m2[a] = DDouble(mass[a]*mass[a]);
      t.m_R += legs[a];
      for (size_t b(a+1); b < t.m_n; ++b)
        t.m_s[a*t.m_n+b] = DDouble(2.0)*YFSExactLCDot(
          Wide(legs[a]), DDouble(mass[a]*mass[a]),
          Wide(legs[b]), DDouble(mass[b]*mass[b]));
    }
    t.m_tol = 1e-9*2.*E;
    t.m_valid = true;
    return t;
  }

}

TEST_CASE("Invariant table virtualities", "[METOOLS][YFS_Stable_Current]") {
  std::vector<Vec4D> legs;
  std::vector<double> mass;
  const Invariant_Table t(BalancedTable(legs, mass));
  REQUIRE(t.Balanced(Vec4D()));
  // every subset: the table's p^2 against the direct (sum p)^2
  for (size_t set(1); set < t.Full(); ++set) {
    Vec4D sum;
    for (size_t a(0); a < t.m_n; ++a) if (set & (size_t(1) << a)) sum += legs[a];
    const double direct(sum.Abs2()), table((double)t.SetVirtuality(set));
    CHECK(table == Catch::Approx(direct).epsilon(1e-9).margin(1e-9));
  }
  // a single leg is exactly its mass^2
  for (size_t a(0); a < t.m_n; ++a)
    CHECK(t.SetVirtuality(size_t(1) << a) == t.m_m2[a]);
}

TEST_CASE("Invariant table uses the complement only when the legs balance",
          "[METOOLS][YFS_Stable_Current]") {
  std::vector<Vec4D> legs;
  std::vector<double> mass;
  Invariant_Table t(BalancedTable(legs, mass));
  // three legs {0,1,2}: the complement {3} is smaller, so with balanced legs
  // the virtuality is leg 3's mass^2 exactly
  const size_t three(0b0111);
  CHECK(t.Virtuality(three) == t.m_m2[3]);
  // unbalanced legs: the direct sum over {0,1,2}
  t.m_R += Vec4D(1., 0., 0., 0.);
  REQUIRE_FALSE(t.Balanced(Vec4D()));
  CHECK(t.Virtuality(three) == t.SetVirtuality(three));
}

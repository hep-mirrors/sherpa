// Frame covariance of YFS::SoftRegion (YFS/Main/Soft_Region.H).
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "YFS/Main/Soft_Region.H"
#include "ATOOLS/Math/Poincare.H"

using ATOOLS::Vec4D;
using YFS::SoftRegion;

namespace {

  Vec4D OnShell(double m, double px, double py, double pz) {
    return Vec4D(std::sqrt(m*m + px*px + py*py + pz*pz), px, py, pz);
  }

  // an arbitrary boost: the rest frame of a moving vector, applied as BoostBack
  Vec4D Boosted(const ATOOLS::Poincare &L, Vec4D p) { L.BoostBack(p); return p; }

}

TEST_CASE("SoftRegion lab: energies are the given components", "[YFS][SoftRegion]") {
  const Vec4D p1(OnShell(0.105, 3., -1., 40.)), p2(OnShell(0.105, -2., 1.5, -38.));
  const SoftRegion r(SoftRegion::Lab(0.5));
  double p1p2, E1, E2;
  r.Legs(p1, p2, p1.Mass(), p2.Mass(), p1p2, E1, E2);
  CHECK(E1 == p1[0]);
  CHECK(E2 == p2[0]);
  CHECK(p1p2 == p1*p2);
  CHECK(r.PhotonEnergy(Vec4D(1., 0., 0., 1.)) == 1.);
}

TEST_CASE("SoftRegion rest_of: covariant under a common boost", "[YFS][SoftRegion]") {
  const Vec4D p1(OnShell(0.105, 3., -1., 40.)), p2(OnShell(0.105, -2., 1.5, -12.));
  const Vec4D V(p1 + p2);
  const ATOOLS::Poincare L(Vec4D(10., 2., -3., 6.));   // a generic boost
  const SoftRegion r(SoftRegion::RestOf(V, 1.));
  const SoftRegion rL(SoftRegion::RestOf(Boosted(L, V), 1.));
  double a, E1, E2, aL, E1L, E2L;
  r.Legs(p1, p2, p1.Mass(), p2.Mass(), a, E1, E2);
  rL.Legs(Boosted(L, p1), Boosted(L, p2), p1.Mass(), p2.Mass(), aL, E1L, E2L);
  CHECK(aL == Catch::Approx(a).epsilon(1e-12));
  CHECK(E1L == Catch::Approx(E1).epsilon(1e-12));
  CHECK(E2L == Catch::Approx(E2).epsilon(1e-12));
  // in V's rest frame the two legs of a pair are back to back with E = sqrt(s)/2
  CHECK(E1 + E2 == Catch::Approx(std::sqrt(V.Abs2())).epsilon(1e-12));
  // a photon's energy in the region's frame is k.V/sqrt(V^2), in any frame
  const Vec4D k(2., 0.6, -1.2, 1.5);
  CHECK(rL.PhotonEnergy(Boosted(L, k)) ==
        Catch::Approx(k*V/std::sqrt(V.Abs2())).epsilon(1e-12));
}

TEST_CASE("SoftRegion rest_of equals lab in the vector's rest frame", "[YFS][SoftRegion]") {
  const Vec4D p1(OnShell(0.000511, 1., 2., 45.)), p2(OnShell(0.000511, -1., -2., -30.));
  const Vec4D V(p1 + p2);
  ATOOLS::Poincare rest(V);
  Vec4D r1(p1), r2(p2);
  rest.Boost(r1); rest.Boost(r2);
  double a, E1, E2, b, F1, F2;
  SoftRegion::RestOf(V, 2.).Legs(p1, p2, p1.Mass(), p2.Mass(), a, E1, E2);
  SoftRegion::Lab(2.).Legs(r1, r2, p1.Mass(), p2.Mass(), b, F1, F2);
  CHECK(a == b);
  CHECK(E1 == F1);
  CHECK(E2 == F2);
}

TEST_CASE("SoftRegion pair_rest: two-body energies, threshold fallback", "[YFS][SoftRegion]") {
  const double s(91.1876*91.1876), m1(1.777), m2(0.105);
  double a, E1, E2;
  SoftRegion::PairRest(s, 1e-3, false).Legs(Vec4D(), Vec4D(), m1, m2, a, E1, E2);
  CHECK(E1 + E2 == Catch::Approx(std::sqrt(s)).epsilon(1e-14));
  CHECK(E1*E1 - m1*m1 == Catch::Approx(E2*E2 - m2*m2).epsilon(1e-12));
  // below threshold: sqrt(s)/2 each when the fallback is requested
  SoftRegion::PairRest(1.0, 1e-3, true).Legs(Vec4D(), Vec4D(), 1., 1., a, E1, E2);
  CHECK(E1 == 0.5);
  CHECK(E2 == 0.5);
}

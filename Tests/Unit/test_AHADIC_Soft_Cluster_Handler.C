#include <catch2/catch_all.hpp>

#include "AHADIC++/Tools/Soft_Cluster_Handler.H"
#include "ATOOLS/Math/Vector.H"

#include <cmath>
#include <vector>

using namespace AHADIC;
using namespace ATOOLS;

// SelectRescuePartner picks the hadron that takes the recoil when
// Soft_Cluster_Handler::Rescue turns a too-light 2-parton cluster into a single
// hadron. A diquark-antidiquark pair has no such hadron; Rescue used to build a
// kf_none particle for it whenever another hadron could take the recoil, and
// Gluon_Decayer threw "Couldn't deal with 2-parton singlet" when none could.
// The function only needs momenta and masses, so no particle table or hadpars
// singleton is set up here.
TEST_CASE("SelectRescuePartner", "[AHADIC::Soft_Cluster_Handler]") {
  // a cluster of mass 2.27 at rest, e.g. [su_0b, sd_0] just above its
  // constituent masses, and a hadron of mass 2 to turn it into
  const Vec4D  cluster(2.27, 0., 0., 0.);
  const double newmass(2.);
  const double mpi(0.14);
  const Vec4D  pion_at_rest(mpi, 0., 0., 0.);
  const Vec4D  pion_moving(sqrt(1. + mpi * mpi), 0., 0., 1.);

  SECTION("no hadron for the cluster to turn into") {
    CHECK(SelectRescuePartner(cluster, false, 0., {pion_moving}, {mpi}) == -1);
  }

  SECTION("no hadron to take the recoil") {
    CHECK(SelectRescuePartner(cluster, true, newmass, {}, {}) == -1);
  }

  SECTION("hadrons too light to take the recoil") {
    // (1 + 0.14)^2 / (2 + 0.14)^2 < 1
    const Vec4D light_cluster(1., 0., 0., 0.);
    CHECK(SelectRescuePartner(light_cluster, true, newmass,
                              {pion_at_rest}, {mpi}) == -1);
  }

  SECTION("picks the hadron with the largest invariant-mass ratio") {
    // ratios 1.27 (at rest) and 2.13 (moving)
    CHECK(SelectRescuePartner(cluster, true, newmass, {pion_at_rest}, {mpi})
          == 0);
    CHECK(SelectRescuePartner(cluster, true, newmass,
                              {pion_at_rest, pion_moving}, {mpi, mpi}) == 1);
  }
}

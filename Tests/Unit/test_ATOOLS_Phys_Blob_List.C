#include <catch2/catch_all.hpp>

#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Phys/KF_Table.H"

using namespace ATOOLS;

namespace {

  void Boot() {
    static bool booted = false;
    if (booted) return;
    if (!ATOOLS::msg) ATOOLS::msg = new ATOOLS::Message();
    Settings::InitializeMainSettings("");
    if (s_kftable.find(kf_d) == s_kftable.end())
      AddParticle(kf_d, 0.01, 0., 0., -1, 1, true, 1, "d", "d");
    booted = true;
  }

}

TEST_CASE("Shared particles are boosted exactly once",
          "[ATOOLS::Phys::Blob_List]") {
  // A particle can be outgoing in one blob and incoming in the next; the
  // treateds set must protect it against a double boost.
  Boot();
  Blob_List bloblist;
  Blob* blob1 = new Blob();
  Blob* blob2 = new Blob();
  const Vec4D pshared(5., 0., 0., 3.), pother(2., 1., 0., 0.5);
  Particle* shared = new Particle(-1, Flavour(kf_d), pshared);
  Particle* other  = new Particle(-1, Flavour(kf_d), pother);
  blob1->AddToOutParticles(shared);
  blob2->AddToInParticles(shared);
  blob2->AddToOutParticles(other);
  bloblist.push_back(blob1);
  bloblist.push_back(blob2);

  Poincare cms(Vec4D(10., 0., 0., 6.));
  std::set<Particle*> treateds;
  bloblist.Boost(cms, &treateds);

  CHECK(shared->Momentum() == cms * pshared);
  CHECK(other->Momentum() == cms * pother);
  bloblist.Clear();
  Blob::ResetCounter();
}

TEST_CASE("Boosting in and out restores momenta and positions",
          "[ATOOLS::Phys::Blob_List]") {
  Boot();
  Blob_List bloblist;
  Blob* blob1 = new Blob();
  Blob* blob2 = new Blob();
  const Vec4D pshared(5., 0.5, -0.3, 3.), pother(2., 1., 0., 0.5);
  const Vec4D xblob(0., 1., 2., 3.);
  Particle* shared = new Particle(-1, Flavour(kf_d), pshared);
  Particle* other  = new Particle(-1, Flavour(kf_d), pother);
  blob1->AddToOutParticles(shared);
  blob2->AddToInParticles(shared);
  blob2->AddToOutParticles(other);
  blob1->SetPosition(xblob);
  bloblist.push_back(blob1);
  bloblist.push_back(blob2);

  Poincare cms(Vec4D(12., 0., 0., -7.));
  std::set<Particle*> treateds;
  bloblist.Boost(cms, &treateds);
  Poincare lab(cms);
  lab.Invert();
  std::set<Particle*> treateds2;
  bloblist.Boost(lab, &treateds2);

  CHECK(shared->Momentum() == pshared);
  CHECK(other->Momentum() == pother);
  CHECK(blob1->Position() == xblob);
  bloblist.Clear();
  Blob::ResetCounter();
}

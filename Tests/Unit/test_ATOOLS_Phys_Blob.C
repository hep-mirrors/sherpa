#include <catch2/catch_all.hpp>
#include "ATOOLS/Phys/Blob.H"
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace ATOOLS;
using namespace std;

TEST_CASE("Blob getter", "[ATOOLS::Phys::Blob]") {
  Blob_Data<string> blob("test");
  CHECK(blob.Get() == "test");
}

TEST_CASE("Blob setter", "[ATOOLS::Phys::Blob]") {
  Blob_Data<string> blob("test");
  blob.Set("sherpa");
  CHECK(blob.Get() == "sherpa");
}

TEST_CASE("Blob streamer", "[ATOOLS::Phys::Blob]") {
  Blob_Data<string> blob("test");
  cout << blob << endl;
  const bool result = true;
  CHECK(result);
}

TEST_CASE("Blob clone", "[ATOOLS::Phys::Blob]") {
  Blob_Data<string> blob("test");
  Blob_Data_Base* blobPtr = blob.ClonePtr();
  CHECK(blobPtr);
  CHECK(dynamic_cast<Blob_Data<string>*>(blobPtr)->Get() == "test"s);
  delete blobPtr;
}

TEST_CASE("Blob base clone", "[ATOOLS::Phys::Blob]") {

  class TestBlob : public Blob_Data_Base {
    using Blob_Data_Base::Blob_Data_Base;
    ostream& operator>>(ostream& ostr) const { return ostr; };
  };

  TestBlob blob;
  CHECK(blob.ClonePtr() == NULL);
}

TEST_CASE("Blob base setter", "[ATOOLS::Phys::Blob]") {

  Blob_Data<string> blob("test");
  Blob_Data_Base* blobPtr = blob.ClonePtr();
  blobPtr->Set<string>("sherpa");
  CHECK(blobPtr->Get<string>() == "sherpa");
  delete blobPtr;
}

TEST_CASE("Blob status code logic", "[ATOOLS::Phys::Blob]") {

  vector<blob_status::code> codes{
    blob_status::inactive,
    blob_status::needs_signal,
    blob_status::needs_showers,
    blob_status::needs_harddecays,
    blob_status::needs_beams,
    blob_status::needs_softUE,
    blob_status::needs_reconnections,
    blob_status::needs_hadronization,
    blob_status::needs_hadrondecays,
    blob_status::needs_extraQED,
    blob_status::needs_minBias,
    blob_status::needs_beamRescatter,
    blob_status::needs_smearing,
    blob_status::internal_flag,
    blob_status::needs_yfs,
    blob_status::fully_active,
  };

  for (size_t i=0; i<codes.size()-1; ++i) {
    CHECK( (codes[i] | codes[i+1]) == (codes[i] + codes[i+1]) );
    CHECK( (codes[i] & codes[i]) == codes[i] );
    CHECK( (codes[i] & codes[i+1]) == 0 );
    if (i==0) cout << codes[i] << endl;
    cout << codes[i+1] << endl;
  }

}

TEST_CASE("Btp code logic", "[ATOOLS::Phys::Blob]") {

  vector<btp::code> codes{
    btp::Signal_Process,
    btp::Hard_Decay,
    btp::Hard_Collision,
    btp::Soft_Collision,
    btp::Shower,
    btp::QED_Radiation,
    btp::Beam,
    btp::Bunch,
    btp::Fragmentation,
    btp::Cluster_Formation,
    btp::Cluster_Decay,
    btp::Hadron_Decay,
    btp::Hadron_Mixing,
    btp::Hadron_To_Parton,
    btp::Elastic_Collision,
    btp::Soft_Diffractive_Collision,
    btp::Quasi_Elastic_Collision,
    btp::YFS_Initial,
    btp::Unspecified,
  };

  for (size_t i=0; i<codes.size()-1; ++i) {
    CHECK( (codes[i] | codes[i+1]) == (codes[i] + codes[i+1]) );
    CHECK( (codes[i] & codes[i]) == codes[i] );
    CHECK( (codes[i] & codes[i+1]) == 0 );
    if (i==0) cout << codes[i] << endl;
    cout << codes[i+1] << endl;
  }

}

TEST_CASE("Blob particles setter", "[ATOOLS::Phys::Blob]") {
  CHECK( Blob::Counter() == 0 );
  Blob blob;
  CHECK( Blob::Counter() == 1 );
  blob.ResetCounter();
  CHECK( Blob::Counter() == 0 );
  CHECK( blob.GetOutParticles().size() == 0 );
  CHECK( blob.GetInParticles().size() == 0 );
  CHECK( blob.OutParticles()->size() == 0 );
  CHECK( blob.InParticles()->size() == 0 );
  CHECK( blob.NOutP() == 0 );
  CHECK( blob.NInP() == 0 );
}

TEST_CASE("Blob particles status", "[ATOOLS::Phys::Blob]") {
  Blob blob;
  CHECK( blob.Status() == blob_status::inactive );

  blob.SetStatus(blob_status::needs_signal);
  CHECK( blob.Status() == blob_status::needs_signal );

  blob.AddStatus(blob_status::needs_showers);
  CHECK( blob.Status() == (blob_status::needs_signal | blob_status::needs_showers) );
  CHECK( blob.Has(blob_status::needs_signal) );
  CHECK( blob.Has(blob_status::needs_showers) );

  blob.UnsetStatus(blob_status::needs_signal);
  CHECK( !blob.Has(blob_status::needs_signal) );
  CHECK( blob.Has(blob_status::needs_showers) );

  blob.SetStatus(blob_status::needs_harddecays);
  CHECK( blob.Status() == blob_status::needs_harddecays );
}

TEST_CASE("Blob btp setter", "[ATOOLS::Phys::Blob]") {
  Blob blob;
  CHECK( blob.Type() == btp::Unspecified );

  blob.SetType(btp::Shower);
  CHECK( blob.Type() == btp::Shower );
}

TEST_CASE("Blob type spec setter", "[ATOOLS::Phys::Blob]") {
  Blob blob;
  CHECK( blob.TypeSpec() == "none" );

  blob.SetTypeSpec("test");
  CHECK( blob.TypeSpec() == "test" );
}

TEST_CASE("Blob beam setter", "[ATOOLS::Phys::Blob]") {
  Blob blob;
  blob.SetBeam(1234);
  CHECK(blob.Beam() == 1234);
}

TEST_CASE("Blob position setter", "[ATOOLS::Phys::Blob]") {
  Blob blob;
  CHECK( blob.Position().IsZero() );
  blob.SetPosition(Vec4D(0,1,2,3));
  Vec4D res = blob.Position();
  for (size_t i=0; i < 4; ++i) {
    CHECK( res[i] == i );
  }
}

TEST_CASE("Blob CMS setter", "[ATOOLS::Phys::Blob]") {
  Blob blob;
  CHECK( blob.CMS().IsZero() );
  blob.SetCMS(Vec4D(0,1,2,3));
  Vec4D res = blob.CMS();
  for (size_t i=0; i < 4; ++i) {
    CHECK( res[i] == i );
  }
}

TEST_CASE("Blob ID setter", "[ATOOLS::Phys::Blob]") {
  Blob blob;
  CHECK( blob.Id() == -1 );

  blob.SetId(1234); // could be any number
  CHECK( blob.Id() == 1 );

  blob.SetId(1234); // could be any number
  CHECK( blob.Id() == 2 );

  blob.Reset(1000);

  blob.SetId(1234); // could be any number
  CHECK( blob.Id() == 1001 );

  blob.SetId(-1000); // if negative, ID will be the absolute value
  CHECK( blob.Id() == 1000 );
}

TEST_CASE("Blob data accessor", "[ATOOLS::Phys::Blob]") {
  Blob blob;
  CHECK( blob["test"] == 0 );
}

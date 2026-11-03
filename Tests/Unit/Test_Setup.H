#ifndef TESTS_UNIT_TEST_SETUP_H
#define TESTS_UNIT_TEST_SETUP_H

// Global state for the unit tests. sherpa_tests does not run Sherpa's main, so
// the message and run-parameter singletons, the settings and the particle table
// start out empty, and Catch2 runs the test cases of all files in random order.
// The particle table is shared by all of them: RegisterParticles gives every
// test the same records, whichever registers first.

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/KF_Table.H"
#include "MODEL/Main/Model_Base.H"

namespace Test_Setup {

  // proton mass in GeV and radius in fm
  constexpr double proton_mass = 0.938272, proton_radius = 0.88;

  // the singletons that Sherpa's main creates
  inline void InitGlobals()
  {
    if (!ATOOLS::msg) ATOOLS::msg = new ATOOLS::Message();
    if (!ATOOLS::rpa) ATOOLS::rpa = new ATOOLS::Run_Parameter();
  }

  // the particles the unit tests use; AddParticle reads PARTICLE_DATA, so the
  // main settings must be initialised before
  inline void RegisterParticles()
  {
    struct Record {
      kf_code kf;
      double mass, radius;
      int icharge, spin; // in units of e/3 and 1/2
      const char *name;
    };
    static const Record records[] = {
      {kf_d, 0.01, 0., -1, 1, "d"},
      {kf_u, 0.005, 0., 2, 1, "u"},
      {kf_s, 0.2, 0., -1, 1, "s"},
      {kf_c, 1.42, 0., 2, 1, "c"},
      {kf_b, 4.92, 0., -1, 1, "b"},
      {kf_gluon, 0., 0., 0, 2, "G"},
      {kf_photon, 0., 0., 0, 2, "P"},
      {kf_e, 0.000511, 0., -3, 1, "e-"},
      {kf_mu, 0.1057, 0., -3, 1, "mu-"},
      {kf_pi, 0.135, 0., 0, 0, "pi"},
      {kf_pi_plus, 0.1396, 0., 3, 0, "pi+"},
      {kf_K_plus, 0.494, 0., 3, 0, "K+"},
      {kf_n, 0.9396, 0., 0, 1, "n"},
      {kf_p_plus, proton_mass, proton_radius, 3, 1, "P+"},
      {kf_Lambda, 1.116, 0., 0, 1, "Lambda"},
    };
    for (const auto &r : records)
      if (ATOOLS::s_kftable.find(r.kf) == ATOOLS::s_kftable.end())
        ATOOLS::AddParticle(r.kf, r.mass, r.radius, 0., r.icharge, r.spin,
                            true, 1, r.name, r.name);
  }

  // a model with the five quark flavours, enough for the coupling and colour
  // set-up of the subtraction kernels
  class Quark_Model : public MODEL::Model_Base {
  public:
    Quark_Model() : Model_Base(true)
    {
      for (kf_code kf : {kf_d, kf_u, kf_s, kf_c, kf_b}) {
        m_fls.push_back(ATOOLS::Flavour(kf));
        m_fls.push_back(ATOOLS::Flavour(kf).Bar());
      }
    }
    void InitVertices() override {}
    void ParticleInit() override {}
    bool ModelInit() override { return true; }
  };

}

#endif

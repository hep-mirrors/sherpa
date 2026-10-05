#include <catch2/catch_all.hpp>

#include "PHASIC++/Process/ME_Generator_Base.H"

#include <string>
#include <vector>

using PHASIC::Mass_Shift_Mapping;

namespace {

  // The selection reads only kf codes, so local particle records suffice and
  // the shared KF table stays untouched.
  ATOOLS::Flavour MakeFlavour(ATOOLS::Particle_Info& info, kf_code kf)
  {
    info.m_kfc = kf;
    return ATOOLS::Flavour(info);
  }

  struct Case {
    std::string name;
    ATOOLS::Flavour in0, in1;
    int pdfs, leptonbunches;
    Mass_Shift_Mapping::code mode;
    int fixed, pdf;
  };

}

TEST_CASE("Mass-shift mapping follows incoming flavours, PDFs and lepton beams",
          "[PHASIC::ME_Generator_Base]")
{
  ATOOLS::Particle_Info e_info, nu_info, u_info, g_info, a_info;
  const ATOOLS::Flavour e(MakeFlavour(e_info, kf_e)),
    nu(MakeFlavour(nu_info, kf_nue)), u(MakeFlavour(u_info, kf_u)),
    g(MakeFlavour(g_info, kf_gluon)), photon(MakeFlavour(a_info, kf_photon));
  const auto standard = Mass_Shift_Mapping::standard;
  const auto dis = Mass_Shift_Mapping::dis;
  const auto epa = Mass_Shift_Mapping::epa;

  // pdfs and leptonbunches carry bit i for beam i.
  const std::vector<Case> cases {
    {"ee without PDFs", e.Bar(), e, 0, 3, standard, -1, -1},
    {"ee with electron PDFs", e.Bar(), e, 3, 3, standard, -1, -1},
    {"pp", u, g, 3, 0, standard, -1, -1},
    {"ep, electron without PDF", e.Bar(), u, 2, 1, dis, 0, 1},
    {"ep, electron with PDF", e.Bar(), u, 3, 1, dis, 0, 1},
    {"pe, electron without PDF", u, e.Bar(), 1, 2, dis, 1, 0},
    {"pe, electron with PDF", u, e.Bar(), 3, 2, dis, 1, 0},
    {"neutrino DIS", nu.Bar(), u, 2, 1, dis, 0, 1},
    {"lepton from a proton PDF", e.Bar(), u, 3, 0, standard, -1, -1},
    {"direct photoproduction", photon, u, 2, 0, epa, 0, 1},
    {"direct photoproduction, photon on beam 1", u, photon, 1, 0, epa, 1, 0},
    {"resolved photoproduction", g, u, 3, 0, standard, -1, -1},
    {"direct photon-photon", photon, photon, 0, 0, standard, -1, -1},
    {"direct and resolved photon", photon, g, 2, 0, epa, 0, 1},
    {"lepton and direct photon", e.Bar(), photon, 0, 1, standard, -1, -1},
    {"lepton with PDF and direct photon", e.Bar(), photon, 1, 1, epa, 1, 0},
    {"two leptons, one with PDF", e.Bar(), e, 1, 3, standard, -1, -1},
  };
  for (const auto& c : cases) {
    CAPTURE(c.name);
    const auto mapping =
      PHASIC::SelectMassShiftMapping(c.in0, c.in1, c.pdfs, c.leptonbunches);
    CHECK(mapping.m_mode == c.mode);
    CHECK(mapping.m_fixed == c.fixed);
    CHECK(mapping.m_pdf == c.pdf);
  }
}

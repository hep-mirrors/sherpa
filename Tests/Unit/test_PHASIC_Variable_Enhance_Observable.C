#include <catch2/catch_all.hpp>

#include "Tests/Unit/PHASIC_Test_Utils.H"
#include "PHASIC++/Enhance/Enhance_Observable_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"

#include <cmath>
#include <memory>
#include <string>

TEST_CASE("Enhance-observable H_T tags sum over final-state particles only",
          "[PHASIC::Variable_Enhance_Observable]")
{
  PHASIC_Test::SeedGlobals();
  PHASIC_Test::ScopedParticle electron(kf_e), muon(kf_mu), gluon(kf_gluon);
  electron.Info().m_icharge = -3;
  muon.Info().m_icharge = -3;
  gluon.Info().m_strong = 8;
  const ATOOLS::Flavour_Vector flavours = {
      ATOOLS::Flavour(kf_e), ATOOLS::Flavour(kf_e, true),
      ATOOLS::Flavour(kf_mu), ATOOLS::Flavour(kf_mu, true),
      ATOOLS::Flavour(kf_gluon)};
  PHASIC_Test::Process process(flavours, 2);
  // e+e- -> mu+mu- g at sqrt(s) = 100 GeV: the colourless beams must not
  // enter the colourless sum of H_TMp2
  const double pz = std::sqrt(45. * 45. - 25.);
  const ATOOLS::Vec4D_Vector momenta = {
      ATOOLS::Vec4D(50., 0., 0., 50.), ATOOLS::Vec4D(50., 0., 0., -50.),
      ATOOLS::Vec4D(45., -3., -4., pz), ATOOLS::Vec4D(45., -3., -4., -pz),
      ATOOLS::Vec4D(10., 6., 8., 0.)};
  process.Integrator()->SetMomenta(momenta);
  const std::string tag = "VAR{H_TMp2}";
  std::unique_ptr<PHASIC::Enhance_Observable_Base> observable(
      PHASIC::Enhance_Observable_Base::Getter_Function::GetObject(
          tag, PHASIC::Enhance_Arguments(&process, tag)));
  REQUIRE(observable);
  const double htmp = (momenta[2] + momenta[3]).MPerp() + momenta[4].MPerp();
  CHECK_THAT((*observable)(&momenta.front(), &flavours.front(), momenta.size()),
             Catch::Matchers::WithinRel(htmp * htmp, 1.e-12));
}

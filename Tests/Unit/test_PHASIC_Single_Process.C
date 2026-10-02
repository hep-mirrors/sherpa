#include <catch2/catch_all.hpp>

#include "Tests/Unit/PHASIC_Selector_Test_Utils.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "PDF/Main/Intact.H"
#include "PDF/Main/ISR_Base.H"
#include "PDF/Main/PDF_Base.H"
#include "ATOOLS/Org/Scoped_Settings.H"

namespace {
class Process : public PHASIC::Single_Process {
public:
  // Select the cached mode in this stub, not via an immutable shared override.
  // This exercises the public conversion helper, not settings-string parsing.
  Process() { m_nfconvscheme = 2; }
  double Partonic(const ATOOLS::Vec4D_Vector&, ATOOLS::Variations_Mode,
                  int) override { return m_lastxs = m_lastbxs = 30.; }
};

// Six strong quarks, the gluon and the group of the first nlight quarks.
class Quarks {
  ATOOLS::Particle_Info m_info, m_particles[7];
  ATOOLS::Particle_Info* p_previous;
  ATOOLS::Particle_Info* p_particles[7];
public:
  explicit Quarks(int nlight = 6)
  {
    // Other unit fixtures use simplified quarks with zero strong charge.
    // Temporarily replace, rather than mutate or re-register, those entries.
    const double masses[] = {0., 0., .1, 1.5, 4.8, 173., 0.};
    for (int i = 0; i < 7; ++i) {
      const auto code = static_cast<kf_code>(i < 6 ? i + 1 : 21);
      auto old = ATOOLS::s_kftable.find(code);
      p_particles[i] = old == ATOOLS::s_kftable.end() ? nullptr : old->second;
      m_particles[i].m_kfc = code;
      m_particles[i].m_strong = i < 6 ? 3 : 8;
      m_particles[i].m_mass = m_particles[i].m_hmass = masses[i];
      m_particles[i].m_content.push_back(new ATOOLS::Flavour(m_particles[i]));
      ATOOLS::s_kftable[code] = &m_particles[i];
    }
    auto old = ATOOLS::s_kftable.find(kf_quark);
    p_previous = old == ATOOLS::s_kftable.end() ? nullptr : old->second;
    m_info.m_kfc = kf_quark;
    m_info.m_mass = -1.;
    m_info.SetIsGroup(true);
    for (int code = 1; code <= nlight; ++code) {
      ATOOLS::Flavour quark(static_cast<kf_code>(code));
      m_info.Add(quark);
      m_info.Add(quark.Bar());
    }
    ATOOLS::s_kftable[kf_quark] = &m_info;
  }
  ~Quarks()
  {
    if (p_previous) ATOOLS::s_kftable[kf_quark] = p_previous;
    else ATOOLS::s_kftable.erase(kf_quark);
    for (int i = 0; i < 7; ++i) {
      const auto code = static_cast<kf_code>(i < 6 ? i + 1 : 21);
      if (p_particles[i]) ATOOLS::s_kftable[code] = p_particles[i];
      else ATOOLS::s_kftable.erase(code);
    }
  }
};

// A PDF that only reports its flavour count, as an LHAPDF set's NumFlavors.
class FlavourCountPDF : public PDF::PDF_Base {
  void CalculateSpec(const double&, const double&) override {}
public:
  explicit FlavourCountPDF(int nf)
  {
    m_bunch = ATOOLS::Flavour(kf_gluon);
    m_partons.insert(m_bunch);
    m_asinfo.m_nf = nf;
  }
  double GetXPDF(const ATOOLS::Flavour&) override { return 1.; }
  double GetXPDF(const kf_code&, bool) override { return 1.; }
  PDF::PDF_Base* GetCopy() override { return new FlavourCountPDF(m_asinfo.m_nf); }
};

class ISR : public PDF::ISR_Base {
public:
  explicit ISR(int nf) : PDF::ISR_Base(new FlavourCountPDF(nf))
  { m_bunch = p_pdf->Bunch(); m_type = PDF::isrtype::hadron; }
  bool CalculateWeight(double, double, double, double, int) override { return true; }
  double Weight(ATOOLS::Flavour) override { return 1.; }
};

class ActiveAlphaS {
  MODEL::Running_AlphaS* p_previous;
public:
  explicit ActiveAlphaS(MODEL::Running_AlphaS& coupling)
    : p_previous(MODEL::as) { MODEL::as = &coupling; }
  ~ActiveAlphaS() { MODEL::as = p_previous; }
};
}

TEST_CASE("Single_Process takes the conversion flavour count from alphaS without PDFs",
          "[PHASIC::Single_Process]")
{
  PHASIC_Test::SeedGlobals();
  auto& settings = ATOOLS::Settings::GetMainSettings();
  settings["ASSOCIATED_CONTRIBUTIONS_VARIATIONS"].SetDefault(
    std::vector<std::vector<std::string>>{});

  Quarks quarks;
  MODEL::Running_AlphaS coupling(.118, 91.1876 * 91.1876, 0, 1);
  ActiveAlphaS active(coupling);
  REQUIRE(coupling.Nf(1.e20) == ATOOLS::Flavour(kf_quark).Size() / 2);

  const ATOOLS::Flavour none(kf_none);
  PDF::ISR_Handler isr({new PDF::Intact(none), new PDF::Intact(none)});
  Process process;
  process.Integrator()->SetISR(&isr);
  REQUIRE(isr.PDF(0) == nullptr);
  REQUIRE(isr.PDF(1) == nullptr);
  // Matching alphaS and matrix-element schemes need no conversion correction.
  double correction = 1.;
  REQUIRE_NOTHROW(correction = process.NfSchemeConversionTerms());
  CHECK_THAT(correction, Catch::Matchers::WithinAbs(0., 1.e-14));

  SECTION("Only mode cancels the nominal result when the schemes already match")
  {
    // Differential adds the returned correction to m_lastxs.
    const double nominal = process.Partonic(
      {}, ATOOLS::Variations_Mode::nominal_only, 0);
    CHECK_THAT(nominal + process.NfSchemeConversionTerms(),
               Catch::Matchers::WithinAbs(0., 1.e-14));
  }
}

TEST_CASE("Single_Process accepts a PDF flavour count below the alphaS thresholds",
          "[PHASIC::Single_Process]")
{
  PHASIC_Test::SeedGlobals();
  auto& settings = ATOOLS::Settings::GetMainSettings();
  settings["ASSOCIATED_CONTRIBUTIONS_VARIATIONS"].SetDefault(
    std::vector<std::vector<std::string>>{});

  // Five light quarks and a five-flavour PDF, as with PDF4LHC21 in the 5F
  // scheme. The top quark still sets an alphaS threshold, so alphaS counts
  // six flavours at high scales.
  Quarks quarks(5);
  MODEL::Running_AlphaS coupling(.118, 91.1876 * 91.1876, 0, 1);
  ActiveAlphaS active(coupling);
  REQUIRE(coupling.Nf(1.e20) == 6);
  REQUIRE(ATOOLS::Flavour(kf_quark).Size() / 2 == 5);

  PDF::ISR_Handler isr({new ISR(5), new ISR(5)});
  Process process;
  process.Integrator()->SetISR(&isr);
  REQUIRE(isr.PDF(0) != nullptr);
  // The schemes of PDF and matrix element match, so Only mode returns the
  // correction that cancels the nominal result.
  const double nominal = process.Partonic(
    {}, ATOOLS::Variations_Mode::nominal_only, 0);
  double correction = 0.;
  REQUIRE_NOTHROW(correction = process.NfSchemeConversionTerms());
  CHECK_THAT(nominal + correction, Catch::Matchers::WithinAbs(0., 1.e-14));
}

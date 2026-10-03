#include <catch2/catch_all.hpp>

#include "Tests/Unit/PHASIC_Selector_Test_Utils.H"
#include "PHASIC++/Scales/MINLO_Scale_Setter.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "PDF/Main/ISR_Base.H"

#include <memory>

namespace {
std::string ScaleFactor(const std::string& name)
{
  try { return ATOOLS::rpa->gen.Variable(name); }
  catch (const ATOOLS::fatal_error&) { return "1"; }
}

class Environment {
  ATOOLS::Particle_Info m_gluon;
  ATOOLS::Particle_Info* p_previous_gluon;
  MODEL::Running_AlphaS* p_previous_as;
  std::unique_ptr<MODEL::Running_AlphaS> m_as;
  double m_ecms;
  std::string m_rsf, m_fsf;
public:
  Environment()
  {
    PHASIC_Test::SeedGlobals();
    auto& settings = ATOOLS::Settings::GetMainSettings();
    settings["ASSOCIATED_CONTRIBUTIONS_VARIATIONS"].SetDefault(
      std::vector<std::vector<std::string>>{});
    settings["DIPOLES"]["NF_GSPLIT"].SetDefault(0);
    const auto it = ATOOLS::s_kftable.find(kf_gluon);
    p_previous_gluon = it == ATOOLS::s_kftable.end() ? nullptr : it->second;
    m_gluon.m_kfc = kf_gluon;
    m_gluon.m_strong = 8;
    ATOOLS::s_kftable[kf_gluon] = &m_gluon;
    p_previous_as = MODEL::as;
    m_as.reset(new MODEL::Running_AlphaS(.118, 91.1876 * 91.1876, 0, 1));
    MODEL::as = m_as.get();
    m_ecms = ATOOLS::rpa->gen.Ecms();
    // Seed missing nominal factors; preserve previously registered values.
    m_rsf = ScaleFactor("RENORMALIZATION_SCALE_FACTOR");
    m_fsf = ScaleFactor("FACTORIZATION_SCALE_FACTOR");
    ATOOLS::rpa->gen.SetEcms(1000.);
    ATOOLS::rpa->gen.SetVariable("RENORMALIZATION_SCALE_FACTOR", "1");
    ATOOLS::rpa->gen.SetVariable("FACTORIZATION_SCALE_FACTOR", "1");
  }
  ~Environment()
  {
    ATOOLS::rpa->gen.SetEcms(m_ecms);
    ATOOLS::rpa->gen.SetVariable("RENORMALIZATION_SCALE_FACTOR", m_rsf);
    ATOOLS::rpa->gen.SetVariable("FACTORIZATION_SCALE_FACTOR", m_fsf);
    MODEL::as = p_previous_as;
    m_as.reset();
    if (p_previous_gluon) ATOOLS::s_kftable[kf_gluon] = p_previous_gluon;
    else ATOOLS::s_kftable.erase(kf_gluon);
  }
};

// Only PDF bounds/membership are needed, not a physical convolution.
class BoundsPDF : public PDF::PDF_Base {
  void CalculateSpec(const double&, const double&) override {}
public:
  BoundsPDF()
  {
    m_bunch = ATOOLS::Flavour(kf_gluon);
    m_partons.insert(m_bunch);
    m_q2min = 1.;
    m_q2max = 1.e6;
    SetBounds();
  }
  double GetXPDF(const ATOOLS::Flavour&) override { return 1.; }
  double GetXPDF(const kf_code&, bool) override { return 1.; }
  PDF::PDF_Base* GetCopy() override { return new BoundsPDF(); }
};

class ISR : public PDF::ISR_Base {
public:
  ISR() : PDF::ISR_Base(new BoundsPDF())
  { m_bunch = p_pdf->Bunch(); m_type = PDF::isrtype::hadron; }
  bool CalculateWeight(double, double, double, double, int) override { return true; }
  double Weight(ATOOLS::Flavour) override { return 1.; }
};

class Process : public PHASIC::Single_Process {
public:
  Process()
  {
    m_flavs.assign(4, ATOOLS::Flavour(kf_gluon));
    m_nin = m_nout = 2;
    m_maxcpl = m_mincpl = {2., 0.};
  }
  double Partonic(const ATOOLS::Vec4D_Vector&, ATOOLS::Variations_Mode,
                  int) override { return 0.; }
};

std::unique_ptr<PHASIC::Scale_Setter_Base> MakeScale(Process& process)
{
  PHASIC::Scale_Setter_Arguments args(nullptr, "MINLO[VAR{100}]");
  args.p_proc = &process;
  return std::unique_ptr<PHASIC::Scale_Setter_Base>(
    PHASIC::Scale_Setter_Base::Scale_Getter_Function::GetObject("MINLO", args));
}

const ATOOLS::Vec4D_Vector momenta = {
  ATOOLS::Vec4D(10., 0., 0., 10.), ATOOLS::Vec4D(10., 0., 0., -10.),
  ATOOLS::Vec4D(10., 10., 0., 0.), ATOOLS::Vec4D(10., -10., 0., 0.)};
}

TEST_CASE("MINLO stops clustering at the configured minimum final multiplicity",
          "[PHASIC::MINLO_Scale_Setter]")
{
  Environment environment;
  PHASIC_Test::Jet_Group jets{ATOOLS::Flavour(kf_gluon)};
  PDF::ISR_Handler isr({new ISR(), new ISR()});
  Process process;
  process.Integrator()->SetISR(&isr);
  auto scale = MakeScale(process);
  REQUIRE(scale != nullptr);
  auto* minlo = dynamic_cast<PHASIC::MINLO_Scale_Setter*>(scale.get());
  REQUIRE(minlo != nullptr);
  CHECK_THAT(scale->CalculateScale(momenta), Catch::Matchers::WithinRel(100., 1.e-12));
  REQUIRE(minlo->Ampl() != nullptr);
  CHECK(minlo->Ampl()->Legs().size() == 4);
  CHECK(minlo->Ampl()->Next() == nullptr);
}

#ifndef TESTS_UNIT_PHASIC_SELECTOR_TEST_UTILS_H
#define TESTS_UNIT_PHASIC_SELECTOR_TEST_UTILS_H

#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/KF_Table.H"

namespace PHASIC_Test {

inline void SeedGlobals()
{
  if (!ATOOLS::msg) ATOOLS::msg = new ATOOLS::Message();
  if (!ATOOLS::rpa) ATOOLS::rpa = new ATOOLS::Run_Parameter();
  static bool seeded = false;
  if (!seeded) {
    ATOOLS::Settings::InitializeMainSettings(std::string{"{}"});
    seeded = true;
  }
  if (ATOOLS::s_kftable.find(kf_none) == ATOOLS::s_kftable.end())
    ATOOLS::AddParticle(kf_none, 0., 0., 0., 0, 0, true, 1, "none", "none");
}

// Supply external-leg metadata without constructing an ME generator. None of
// the event-generation stubs is used by the selector tests.
class Process : public PHASIC::Process_Base {
public:
  Process(const ATOOLS::Flavour_Vector& flavours, size_t nin)
  {
    m_flavs = flavours;
    m_nin = nin;
    m_nout = flavours.size() - nin;
  }
  void SetScale(const PHASIC::Scale_Setter_Arguments&) override {}
  void SetKFactor(const PHASIC::KFactor_Setter_Arguments&) override {}
  void InitializeTheReweighting(ATOOLS::Variations_Mode) override {}
  size_t Size() const override { return 1; }
  PHASIC::Process_Base* operator[](const size_t&) override { return this; }
  ATOOLS::Weight_Info* OneEvent(int, ATOOLS::Variations_Mode, int) override
  { return nullptr; }
  ATOOLS::Weights_Map Differential(const ATOOLS::Vec4D_Vector&,
                                   ATOOLS::Variations_Mode) override
  { return ATOOLS::Weights_Map{0.0}; }
  bool CalculateTotalXSec(const std::string&, bool) override { return false; }
  void SetLookUp(bool lookup) override { m_lookup = lookup; }
};

// Temporarily supply the jet group without altering any existing shared entry
// or leaving Flavours pointing at a destroyed Particle_Info after the test.
class Jet_Group {
  ATOOLS::Particle_Info m_jet;
  ATOOLS::Particle_Info* p_previous;
public:
  explicit Jet_Group(const ATOOLS::Flavour& constituent)
  {
    auto it = ATOOLS::s_kftable.find(kf_jet);
    p_previous = it == ATOOLS::s_kftable.end() ? nullptr : it->second;
    m_jet.m_kfc = kf_jet;
    m_jet.SetIsGroup(true);
    m_jet.Add(constituent);
    ATOOLS::s_kftable[kf_jet] = &m_jet;
  }
  ~Jet_Group()
  {
    if (p_previous) ATOOLS::s_kftable[kf_jet] = p_previous;
    else ATOOLS::s_kftable.erase(kf_jet);
  }
  Jet_Group(const Jet_Group&) = delete;
  Jet_Group& operator=(const Jet_Group&) = delete;
};

} // namespace PHASIC_Test

#endif

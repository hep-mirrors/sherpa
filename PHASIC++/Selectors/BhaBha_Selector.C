#ifndef ATOOLS_Phys_Standard_Selector_BHABHA_H
#define ATOOLS_Phys_Standard_Selector_BHABHA_H
#include "ATOOLS/Math/Vec4.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Flavour.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Selectors/Selector.H"

using namespace PHASIC;
using namespace ATOOLS;

namespace PHASIC {
class Angle_Selector : public Selector_Base {
  double m_angmin, m_angmax, m_emin;
  // The beam momenta the angle is measured against are lab-frame quantities,
  // while the Selector_List momenta are in the frame the phase space was
  // generated in. For symmetric beams these coincide; for asymmetric ones
  // (fixed target) they do not, so the two must be brought into the same frame
  // before an angular acceptance means anything.
  bool m_labframe;
  ATOOLS::Flavour m_flav;

public:
  Angle_Selector(Process_Base *const);
  ~Angle_Selector();
  void SetRange(ATOOLS::Flavour, double, double, double emin = 0.);
  inline void SetLabframe(const bool &lab) { m_labframe = lab; }
  bool Trigger(ATOOLS::Selector_List &);
  void BuildCuts(Cut_Data *);
};

class FourFermion_Selector : public Selector_Base {
  double m_angmin, m_angmax, m_emin;
  bool m_labframe;  // see Angle_Selector
  ATOOLS::Flavour m_flav1, m_flav2;

public:
  FourFermion_Selector(Process_Base *const);
  ~FourFermion_Selector();
  void SetRange(ATOOLS::Flavour, ATOOLS::Flavour, double, double,
                double emin = 0.);
  inline void SetLabframe(const bool &lab) { m_labframe = lab; }
  bool Trigger(ATOOLS::Selector_List &);
  void BuildCuts(Cut_Data *);
  bool OneBhabha(ATOOLS::Selector_List &);
  //! angular window + energy threshold for one leg, measured against \p beam
  bool PassesLeg(const ATOOLS::Vec4D &p, const ATOOLS::Vec4D &beam) const;
};

//   class IINEL_Selector : public Selector_Base {
//     double m_ymin, m_ymax;
//     ATOOLS::Flavour m_flav1,m_flav2;
//   public:
//     IINEL_Selector(Process_Base *const);
//     ~IINEL_Selector();
//     void     SetRange(ATOOLS::Flavour,ATOOLS::Flavour,double,double);
//     bool     Trigger(ATOOLS::Selector_List &);
//     void     BuildCuts(Cut_Data *);
//   };

} // namespace PHASIC

#endif

/*--------------------------------------------------------------------

  Bhabha angle Selector

  --------------------------------------------------------------------*/

Angle_Selector::Angle_Selector(Process_Base *const proc)
    : Selector_Base("BhaBha_Angle_Selector", proc), m_angmin(0.),
      m_angmax(10.0), m_emin(0.), m_labframe(false),
      m_flav(Flavour(kf_none)) {}

Angle_Selector::~Angle_Selector() {}

bool Angle_Selector::Trigger(Selector_List &sl) {
  if (!m_on)
    return true;
  const Flavour beamfl1 = p_proc->Integrator()->Beam()->GetBeam(0)->Beam();
  const Flavour beamfl2 = p_proc->Integrator()->Beam()->GetBeam(1)->Beam();
  if (beamfl1 != m_flav && beamfl2 != m_flav)
    return true;
  Vec4D beamMom;
  for (size_t i = m_nin; i < sl.size(); i++) {
    if (m_flav.Includes(sl[i].Flavour())) {
      Vec4D mom = sl[i].Momentum();
      if (m_flav == beamfl1)
        beamMom = p_proc->Integrator()->Beam()->GetBeam(0)->OutMomentum();
      else
        beamMom = p_proc->Integrator()->Beam()->GetBeam(1)->OutMomentum();
      // beamMom is lab frame; bring mom into it as well when asked.
      if (m_labframe) p_proc->Integrator()->Beam()->BoostBackLab(mom);
      const double ang = mom.Theta(beamMom);
      //   if(!m_use_radians) ang *= 180./M_PI;
      // Selector_Log counts the PASS condition (see Standard_Selector); this
      // used to be handed the failure condition, which swapped the passed and
      // rejected tallies and inverted the "phase space too tight" warning.
      if (!m_sel_log->CountingIdentity(((ang >= m_angmin) && (ang <= m_angmax) &&
                                        (mom[0] >= m_emin))))
        return false;
    }
  }
  return true;
}

void Angle_Selector::BuildCuts(Cut_Data *cuts) {
  // Without an energy cut nothing constrains Cut_Data: with YFS ISR the
  // selector-derived smin is what clamps vmax=1-smin/s, so a pure angular
  // window leaves the radiative return to s'->0 unregulated.
  if (!m_on || m_emin <= 0.) return;
  if (m_isnlo) {
    cuts->smin = Max(cuts->smin, m_smin);
    return;
  }
  for (int i = m_nin; i < m_n; i++)
    if (m_flav.Includes(p_fl[i]))
      cuts->energymin[i] =
          Max(Max(m_emin, p_fl[i].SelMass()), cuts->energymin[i]);
}

void Angle_Selector::SetRange(Flavour flav, double min, double max,
                              double emin) {
  m_flav = flav;
  m_angmin = min;
  m_angmax = max;
  m_emin = emin;
  m_smin = sqr(emin);
  m_on = true;
}

DECLARE_GETTER(Angle_Selector, "BhaBhaAngle", Selector_Base, Selector_Key);

Selector_Base *
ATOOLS::Getter<Selector_Base, Selector_Key, Angle_Selector>::operator()(
    const Selector_Key &key) const {
  Scoped_Settings s{key.m_settings};
  const auto parameters =
      s.SetDefault<std::string>({}).GetVector<std::string>();
  if (parameters.size() < 4 || parameters.size() > 6)
    THROW(critical_error, "Invalid syntax");
  const auto kf1 = s.Interprete<int>(parameters[1]);
  const auto min = s.Interprete<double>(parameters[2]);
  const auto max = s.Interprete<double>(parameters[3]);
  // Optional trailing arguments, as in Standard_Selector: minimum energy and
  // whether the angle is evaluated in the lab frame.
  double emin(0.);
  bool labframe(false);
  if (parameters.size() > 4) emin = s.Interprete<double>(parameters[4]);
  if (parameters.size() > 5) labframe = s.Interprete<int>(parameters[5]);
  // A negative kf code denotes the antiparticle; kf_code is unsigned, so it
  // has to be split off rather than passed through.
  Flavour flav1 = Flavour((kf_code)std::abs(kf1), kf1 < 0);
  Angle_Selector *sel = new Angle_Selector(key.p_proc);
  sel->SetRange(flav1, min, max, emin);
  sel->SetLabframe(labframe);
  return sel;
}

void ATOOLS::Getter<Selector_Base, Selector_Key, Angle_Selector>::PrintInfo(
    std::ostream &str, const size_t width) const {
  str << "[BhaBhaAngle, kf, angmin, angmax{, emin{, labframe}}]";
}

/*--------------------------------------------------------------------

  Bhabha 4 fermion Selector

  --------------------------------------------------------------------*/

FourFermion_Selector::FourFermion_Selector(Process_Base *const proc)
    : Selector_Base("BhaBha4f_Selector", proc), m_angmin(0.), m_angmax(10.0),
      m_emin(0.), m_labframe(false),
      m_flav1(Flavour(kf_none)), m_flav2(Flavour(kf_none)) {}

FourFermion_Selector::~FourFermion_Selector() {}

bool FourFermion_Selector::Trigger(Selector_List &sl) {
  if (!m_on)
    return true;

  // Check if at least one Bhabha pair exists with angle cuts satisfied.
  // Counted through the Selector_Log so this selector appears in the rejection
  // statistics and can raise the "phase space too tight" warning - it kept no
  // counts at all before, making it invisible to both.
  if (!m_sel_log->CountingIdentity(OneBhabha(sl)))
    return false;

  // Further cuts can be added here in the future
  return true;
}

void FourFermion_Selector::BuildCuts(Cut_Data *cuts) {
  // See Angle_Selector::BuildCuts: without this, YFS ISR sees smin~=0 and
  // never clamps vmax, so the t-channel pole at s'->0 goes unregulated.
  if (!m_on || m_emin <= 0.) return;
  if (m_isnlo) {
    cuts->smin = Max(cuts->smin, m_smin);
    return;
  }
  for (int i = m_nin; i < m_n; i++)
    if (m_flav1.Includes(p_fl[i]) || m_flav2.Includes(p_fl[i]))
      cuts->energymin[i] =
          Max(Max(m_emin, p_fl[i].SelMass()), cuts->energymin[i]);
}

bool FourFermion_Selector::PassesLeg(const ATOOLS::Vec4D &p,
                                     const ATOOLS::Vec4D &beam) const {
  Vec4D mom(p);
  // beam is a lab-frame momentum; see the m_labframe comment on the class.
  if (m_labframe) p_proc->Integrator()->Beam()->BoostBackLab(mom);
  const double ang = mom.Theta(beam);
  return (ang >= m_angmin) && (ang <= m_angmax) && (mom[0] >= m_emin);
}

bool FourFermion_Selector::OneBhabha(ATOOLS::Selector_List &sl) {
  // Find two DISTINCT final-state particles forming the configured pair: one
  // of flavour m_flav1 inside the acceptance around beam 0, one of flavour
  // m_flav2 inside the acceptance around beam 1.
  //
  // Two things were wrong before. The configured flavours were ignored
  // entirely - the legs were matched against the beam flavours instead, so the
  // Flavs setting had no effect on the selection. And the two legs were tested
  // independently, so in a four-fermion final state the leg of one pair could
  // satisfy the first condition while the leg of another satisfied the second,
  // accepting a combination that is not a pair at all. Requiring i != j fixes
  // that; the pair is what the "4f" in the name is about.
  const Vec4D beam1mom = p_proc->Integrator()->Beam()->GetBeam(0)->OutMomentum();
  const Vec4D beam2mom = p_proc->Integrator()->Beam()->GetBeam(1)->OutMomentum();
  for (size_t i = m_nin; i < sl.size(); i++) {
    if (!m_flav1.Includes(sl[i].Flavour())) continue;
    if (!PassesLeg(sl[i].Momentum(), beam1mom)) continue;
    for (size_t j = m_nin; j < sl.size(); j++) {
      if (j == i) continue;
      if (!m_flav2.Includes(sl[j].Flavour())) continue;
      if (!PassesLeg(sl[j].Momentum(), beam2mom)) continue;
      return true;
    }
  }
  return false;
}

void FourFermion_Selector::SetRange(Flavour flav1, Flavour flav2, double min,
                                    double max, double emin) {
  m_flav1 = flav1;
  m_flav2 = flav2;
  m_angmin = min;
  m_angmax = max;
  m_emin = emin;
  // Two legs are required above emin, so the smallest invariant mass this
  // admits is (2*emin)^2 for the back-to-back configuration the angular window
  // selects. This matches what Cut_Data::Complete derives from the energymin
  // array on the non-NLO path below (it sums the per-particle minima and then
  // squares), which sqr(emin) did not - the NLO clamp on vmax=1-smin/s was four
  // times looser than the LO one.
  m_smin = sqr(2. * emin);
  m_on = true;
}

DECLARE_GETTER(FourFermion_Selector, "BhaBha", Selector_Base, Selector_Key);

Selector_Base *
ATOOLS::Getter<Selector_Base, Selector_Key, FourFermion_Selector>::operator()(
    const Selector_Key &key) const {
  //   Scoped_Settings s{ key.m_settings };
  auto s = key.m_settings["BhaBha"]["BhBhaBha4f"];
//   Scoped_Settings s{ Settings::GetMainSettings()["BhaBha"]["BhBhaBha4f"] };
  s.DeclareVectorSettingsWithEmptyDefault({"Flavs"});
  s.DeclareVectorSettingsWithEmptyDefault({"Angles"});
  const auto bounds = s["Angles"].GetVector<double>();
  // Both bounds are indexed below, so an one-element list would read past the
  // end - checking for empty() alone was not enough.
  if (bounds.size() != 2)
    THROW(critical_error,
          "BhaBha4f selector needs exactly two \"Angles\" entries [min, max]");

  const auto flavs =
      s["Flavs"].SetSynonyms({"Flavours", "Flavors"}).GetVector<int>();
  if (flavs.size() != 2)
    THROW(critical_error, "Only two additional flavours allowed!");
  // A negative kf code denotes the antiparticle. kf_code is unsigned, so
  // passing the raw int through (as before) turned e.g. -11 into a huge code
  // that matched nothing, silently disabling that leg.
  const Flavour flav1((kf_code)std::abs(flavs[0]), flavs[0] < 0);
  const Flavour flav2((kf_code)std::abs(flavs[1]), flavs[1] < 0);
  const double emin = s["EMin"].SetDefault(0.).Get<double>();
  const bool labframe = s["LabFrame"].SetDefault(false).Get<bool>();

  FourFermion_Selector *sel = new FourFermion_Selector(key.p_proc);
  sel->SetRange(flav1, flav2, bounds[0], bounds[1], emin);
  sel->SetLabframe(labframe);
  return sel;
}

void ATOOLS::Getter<Selector_Base, Selector_Key,
                    FourFermion_Selector>::PrintInfo(std::ostream &str,
                                                     const size_t width) const {
  str << "BhaBha: {BhBhaBha4f: {Angles: [min, max], Flavours: [kf1, kf2]"
      << ", EMin: <GeV>, LabFrame: <bool>}}";
}
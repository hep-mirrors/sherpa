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
  double m_angmin, m_angmax;
  ATOOLS::Flavour m_flav;

public:
  Angle_Selector(Process_Base *const);
  ~Angle_Selector();
  void SetRange(ATOOLS::Flavour, double, double);
  bool Trigger(ATOOLS::Selector_List &);
  void BuildCuts(Cut_Data *);
};

class FourFermion_Selector : public Selector_Base {
  double m_angmin, m_angmax;
  ATOOLS::Flavour m_flav1, m_flav2;
  ATOOLS::Flavour m_beam1, m_beam2;

public:
  FourFermion_Selector(Process_Base *const);
  ~FourFermion_Selector();
  void SetRange(ATOOLS::Flavour, ATOOLS::Flavour, double, double);
  bool Trigger(ATOOLS::Selector_List &);
  void BuildCuts(Cut_Data *);
  bool OneBhabha(ATOOLS::Selector_List &);
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
      m_angmax(10.0), m_flav(Flavour(kf_none)) {}

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
      const Vec4D mom = sl[i].Momentum();
      if (m_flav == beamfl1)
        beamMom = p_proc->Integrator()->Beam()->GetBeam(0)->OutMomentum();
      else
        beamMom = p_proc->Integrator()->Beam()->GetBeam(1)->OutMomentum();
      //   if(m_labframe) p_proc->Integrator()->Beam()->BoostBackLab(mom);
      const double ang = mom.Theta(beamMom);
      //   if(!m_use_radians) ang *= 180./M_PI;
      if (m_sel_log->Hit(((ang < m_angmin) || (ang > m_angmax))))
        return false;
    }
  }
  return true;
}

void Angle_Selector::BuildCuts(Cut_Data *cuts) {}

void Angle_Selector::SetRange(Flavour flav, double min, double max) {
  m_flav = flav;
  m_angmin = min;
  m_angmax = max;
  m_on = true;
}

DECLARE_GETTER(Angle_Selector, "BhaBhaAngle", Selector_Base, Selector_Key);

Selector_Base *
ATOOLS::Getter<Selector_Base, Selector_Key, Angle_Selector>::operator()(
    const Selector_Key &key) const {
  Scoped_Settings s{key.m_settings};
  const auto parameters =
      s.SetDefault<std::string>({}).GetVector<std::string>();
  if (parameters.size() != 4)
    THROW(critical_error, "Invalid syntax");
  const auto kf1 = s.Interprete<int>(parameters[1]);
  const auto min = s.Interprete<double>(parameters[2]);
  const auto max = s.Interprete<double>(parameters[3]);
  Flavour flav1 = Flavour(kf1);
  Angle_Selector *sel = new Angle_Selector(key.p_proc);
  sel->SetRange(flav1, min, max);
  return sel;
}

void ATOOLS::Getter<Selector_Base, Selector_Key, Angle_Selector>::PrintInfo(
    std::ostream &str, const size_t width) const {
  str << "[BhaBhaAngle, kf1, kf2, min, max]";
}

/*--------------------------------------------------------------------

  Bhabha 4 fermion Selector

  --------------------------------------------------------------------*/

FourFermion_Selector::FourFermion_Selector(Process_Base *const proc)
    : Selector_Base("BhaBha4f_Selector", proc), m_angmin(0.), m_angmax(10.0),
      m_flav1(Flavour(kf_none)), m_flav2(Flavour(kf_none)) {}

FourFermion_Selector::~FourFermion_Selector() {}

bool FourFermion_Selector::Trigger(Selector_List &sl) {
  if (!m_on)
    return true;

  // Check if at least one Bhabha pair exists with angle cuts satisfied
  if (!OneBhabha(sl))
    return false;

  // Further cuts can be added here in the future
  return true;
}

void FourFermion_Selector::BuildCuts(Cut_Data *cuts) {}

bool FourFermion_Selector::OneBhabha(ATOOLS::Selector_List &sl) {
  // Get beam flavours
  const Flavour beamfl1 = p_proc->Integrator()->Beam()->GetBeam(0)->Beam();
  const Flavour beamfl2 = p_proc->Integrator()->Beam()->GetBeam(1)->Beam();

  bool found_beam1_pass = false;
  bool found_beam2_pass = false;

  for (size_t i = m_nin; i < sl.size(); i++) {
    const Flavour &particle_flav = sl[i].Flavour();
    const Vec4D mom = sl[i].Momentum();

    if (particle_flav == beamfl1) {
      Vec4D beamMom = p_proc->Integrator()->Beam()->GetBeam(0)->OutMomentum();
      const double ang = mom.Theta(beamMom);
      if ((ang >= m_angmin) && (ang <= m_angmax)) {
        found_beam1_pass = true;
      }
    }

    if (particle_flav == beamfl2) {
      Vec4D beamMom = p_proc->Integrator()->Beam()->GetBeam(1)->OutMomentum();
      const double ang = mom.Theta(beamMom);
      if ((ang >= m_angmin) && (ang <= m_angmax)) {
        found_beam2_pass = true;
      }
    }
    if (found_beam1_pass && found_beam2_pass) {
      return true;
    }
  }
  return found_beam1_pass && found_beam2_pass;
}

void FourFermion_Selector::SetRange(Flavour flav1, Flavour flav2, double min,
                                    double max) {
  m_flav1 = flav1;
  m_flav2 = flav2;
  m_angmin = min;
  m_angmax = max;
  m_on = true;
  // Store beam flavours for later use
  m_beam1 = p_proc->Integrator()->Beam()->GetBeam(0)->Beam();
  m_beam2 = p_proc->Integrator()->Beam()->GetBeam(1)->Beam();
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
  if (bounds.empty())
    THROW(critical_error,
          "Missing \"Angles\" specification in BhaBha4f selector");

  const auto flavs =
      s["Flavs"].SetSynonyms({"Flavours", "Flavors"}).GetVector<int>();
  if (flavs.size() != 2)
    THROW(critical_error, "Only two additional flavours allowed!");
  //   Flavour flav1 = Flavour((kf_code)abs(kf1),kf1<0);

  FourFermion_Selector *sel = new FourFermion_Selector(key.p_proc);
  sel->SetRange(flavs[0], flavs[1], bounds[0], bounds[1]);
  return sel;
}

void ATOOLS::Getter<Selector_Base, Selector_Key,
                    FourFermion_Selector>::PrintInfo(std::ostream &str,
                                                     const size_t width) const {
  str << "[BhaBhaAngle, kf1, kf2, min, max]";
}
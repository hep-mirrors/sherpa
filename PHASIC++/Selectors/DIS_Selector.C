#ifndef ATOOLS_Phys_Standard_Selector_DIS_H
#define ATOOLS_Phys_Standard_Selector_DIS_H
#include "PHASIC++/Selectors/Selector.H"

namespace PHASIC {
    class IPZIN_Selector : public Selector_Base {
    double m_pzmin, m_pzmax;
    ATOOLS::Flavour m_flav;
  public:
    IPZIN_Selector(Process_Base *const);
    ~IPZIN_Selector();
    void     SetRange(ATOOLS::Flavour,double,double);
    bool     Trigger(ATOOLS::Selector_List &);
    void     BuildCuts(Cut_Data *);
  };

  class IINEL_Selector : public Selector_Base {
    double m_ymin, m_ymax;
    ATOOLS::Flavour m_flav1,m_flav2;
  public:
    IINEL_Selector(Process_Base *const);
    ~IINEL_Selector();
    void     SetRange(ATOOLS::Flavour,ATOOLS::Flavour,double,double);
    bool     Trigger(ATOOLS::Selector_List &);
    void     BuildCuts(Cut_Data *);
  };


}

#endif

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"


using namespace PHASIC;
using namespace ATOOLS;

/*--------------------------------------------------------------------

  Inelesticity Selector

  --------------------------------------------------------------------*/

IINEL_Selector::IINEL_Selector(Process_Base *const proc):
  Selector_Base("INEL_Selector",proc), m_ymin(0.), m_ymax(1.0),
  m_flav1(Flavour(kf_none)), m_flav2(Flavour(kf_none))
{
}

IINEL_Selector::~IINEL_Selector() {
}

bool IINEL_Selector::Trigger(Selector_List &sl)
{
  if (!m_on) return true;
  for (int i=0;i<m_nin;i++) {
    for (int j=m_nin;j<sl.size();j++) {
      if ( (m_flav1.Includes(sl[i].Flavour()) &&
            m_flav2.Includes(sl[j].Flavour())) ||
           (m_flav1.Includes(sl[j].Flavour()) &&
            m_flav2.Includes(sl[i].Flavour())) ) {
        double yij = 1.0-(sl[j].Momentum()[0]/sl[i].Momentum()[0])*(1.0+sl[i].Momentum().CosTheta(sl[j].Momentum()))/2.0;
        if (m_sel_log->Hit( ((yij < m_ymin) || (yij > m_ymax)) )) return false;
      }
    }
  }
  return true;
}

void IINEL_Selector::BuildCuts(Cut_Data * cuts) {}

void IINEL_Selector::SetRange(Flavour flav1,Flavour flav2,double min,double max)
{
  m_flav1=flav1;
  m_flav2=flav2;
  m_ymin=min;
  m_ymax=max;
  m_on=true;
}

DECLARE_GETTER(IINEL_Selector,"INEL",Selector_Base,Selector_Key);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,IINEL_Selector>::
operator()(const Selector_Key &key) const
{
  Scoped_Settings s{ key.m_settings };
  const auto parameters = s.SetDefault<std::string>({}).GetVector<std::string>();
  if (parameters.size() != 5)
    THROW(critical_error, "Invalid syntax");
  const auto kf1 = s.Interprete<int>(parameters[1]);
  const auto kf2 = s.Interprete<int>(parameters[2]);
  const auto min = s.Interprete<double>(parameters[3]);
  const auto max = s.Interprete<double>(parameters[4]);
  Flavour flav1 = Flavour((kf_code)abs(kf1),kf1<0);
  Flavour flav2 = Flavour((kf_code)abs(kf2),kf2<0);
  IINEL_Selector *sel = new IINEL_Selector(key.p_proc);
  sel->SetRange(flav1,flav2,min,max);
  return sel;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,IINEL_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"[INEL, kf1, kf2, min, max]";
}

/*--------------------------------------------------------------------

  Pz Selector

  --------------------------------------------------------------------*/

IPZIN_Selector::IPZIN_Selector(Process_Base *const proc):
  Selector_Base("PZIN_Selector",proc), m_pzmin(0.), m_pzmax(0.),
  m_flav(Flavour(kf_none))
{
}

IPZIN_Selector::~IPZIN_Selector() {
}

bool IPZIN_Selector::Trigger(Selector_List &sl)
{
  if (!m_on) return true;
  for (int i=0;i<m_nin;i++) {
      if  (m_flav.Includes(sl[i].Flavour())) {
        double pz = std::abs(sl[i].Momentum()[3]);
        if (m_sel_log->Hit( ((pz < m_pzmin) || (pz > m_pzmax)) )) return false;
    }
  }
  return true;
}

void IPZIN_Selector::BuildCuts(Cut_Data * cuts)
{
}

void IPZIN_Selector::SetRange(Flavour flav,double min,double max)
{
  m_flav=flav;
  m_pzmin=min;
  m_pzmax=max;
  m_on=true;
}

DECLARE_GETTER(IPZIN_Selector,"PZIN",Selector_Base,Selector_Key);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,IPZIN_Selector>::
operator()(const Selector_Key &key) const
{
  Scoped_Settings s{ key.m_settings };
  const auto parameters = s.SetDefault<std::string>({}).GetVector<std::string>();
  if (parameters.size() != 4)
    THROW(critical_error, "Invalid syntax");
  const auto kf = s.Interprete<int>(parameters[1]);
  const auto min = s.Interprete<double>(parameters[2]);
  const auto max = s.Interprete<double>(parameters[3]);
  Flavour flav = Flavour((kf_code)abs(kf),kf<0);
  IPZIN_Selector *sel = new IPZIN_Selector(key.p_proc);
  sel->SetRange(flav,min,max);
  return sel;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,IPZIN_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"[PZIN, kf, min, max]";
}

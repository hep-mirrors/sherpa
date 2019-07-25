#include "PHASIC++/Selectors/Selector.H"

#include <cassert>

namespace PHASIC {

  // -----------------------------
  // single particle selectors
  // -----------------------------
  class PT_Selector : public Selector_Base {
    double  m_ptmin, m_ptmax;
    ATOOLS::Flavour m_flav;
  public:
    PT_Selector(Process_Base *const);
    ~PT_Selector();
    void     SetRange(ATOOLS::Flavour,double,double);
    bool     Trigger(ATOOLS::Selector_List &);
    void     BuildCuts(Cut_Data *);
  };

  class ET_Selector : public Selector_Base {
    double  m_etmin,m_etmax;
    ATOOLS::Flavour m_flav;
  public:
    ET_Selector(Process_Base *const);
    ~ET_Selector();
    void     SetRange(ATOOLS::Flavour,double,double);
    bool     Trigger(ATOOLS::Selector_List &);
    void     BuildCuts(Cut_Data *);
  };

  class Rapidity_Selector : public Selector_Base {
    double  m_ymin, m_ymax;
    ATOOLS::Flavour m_flav;
  public:
    Rapidity_Selector(Process_Base *const);
    ~Rapidity_Selector();
    void     SetRange(ATOOLS::Flavour,double,double);
    bool     Trigger(ATOOLS::Selector_List &);
    void     BuildCuts(Cut_Data *);
  };

  class PseudoRapidity_Selector : public Selector_Base {
    double  m_etamin, m_etamax;
    ATOOLS::Flavour m_flav;
  public:
    PseudoRapidity_Selector(Process_Base *const);
    ~PseudoRapidity_Selector();
    void     SetRange(ATOOLS::Flavour,double,double);
    bool     Trigger(ATOOLS::Selector_List &);
    void     BuildCuts(Cut_Data *);
  };

  // -----------------------------
  // two particle selectors
  // -----------------------------
  class IMass_Selector : public Selector_Base {
    double  m_massmin, m_massmax;
    ATOOLS::Flavour m_flav1,m_flav2;
  public:
    IMass_Selector(Process_Base *const);
    ~IMass_Selector();
    void     SetRange(ATOOLS::Flavour,ATOOLS::Flavour,double,double);
    bool     Trigger(ATOOLS::Selector_List &);
    void     BuildCuts(Cut_Data *);
  };

  class IQ2_Selector : public Selector_Base {
    double m_q2min, m_q2max;
    ATOOLS::Flavour m_flav1,m_flav2;
  public:
    IQ2_Selector(Process_Base *const);
    ~IQ2_Selector();
    void     SetRange(ATOOLS::Flavour,ATOOLS::Flavour,double,double);
    bool     Trigger(ATOOLS::Selector_List &);
    void     BuildCuts(Cut_Data *);
  };

  class PT2_Selector : public Selector_Base {
    double  m_ptmin, m_ptmax;
    ATOOLS::Flavour m_flav1,m_flav2;
  public:
    PT2_Selector(Process_Base *const);
    ~PT2_Selector();
    void     SetRange(ATOOLS::Flavour,ATOOLS::Flavour,double,double);
    bool     Trigger(ATOOLS::Selector_List &);
    void     BuildCuts(Cut_Data *);
  };

  class MT2_Selector : public Selector_Base {
    double  m_mtmin, m_mtmax;
    ATOOLS::Flavour m_flav1,m_flav2;
  public:
    MT2_Selector(Process_Base *const);
    ~MT2_Selector();
    void     SetRange(ATOOLS::Flavour,ATOOLS::Flavour,double,double);
    bool     Trigger(ATOOLS::Selector_List &);
    void     BuildCuts(Cut_Data *);
  };

  class MT2_v2_Selector : public Selector_Base {
    double  m_mtmin, m_mtmax;
    ATOOLS::Flavour m_flav1,m_flav2;
  public:
    MT2_v2_Selector(Process_Base *const);
    ~MT2_v2_Selector();
    void     SetRange(ATOOLS::Flavour,ATOOLS::Flavour,double,double);
    bool     Trigger(ATOOLS::Selector_List &);
    void     BuildCuts(Cut_Data *);
  };

  class DeltaY_Selector : public Selector_Base {
    double  m_dymin, m_dymax;
    ATOOLS::Flavour m_flav1,m_flav2;
  public:
    DeltaY_Selector(Process_Base *const);
    ~DeltaY_Selector();
    void     SetRange(ATOOLS::Flavour,ATOOLS::Flavour,double,double);
    bool     Trigger(ATOOLS::Selector_List &);
    void     BuildCuts(Cut_Data *);
  };

  class DeltaEta_Selector : public Selector_Base {
    double  m_detamin, m_detamax;
    ATOOLS::Flavour m_flav1,m_flav2;
  public:
    DeltaEta_Selector(Process_Base *const);
    ~DeltaEta_Selector();
    void     SetRange(ATOOLS::Flavour,ATOOLS::Flavour,double,double);
    bool     Trigger(ATOOLS::Selector_List &);
    void     BuildCuts(Cut_Data *);
  };

  class DeltaPhi_Selector : public Selector_Base {
    double  m_dphimin, m_dphimax;
    ATOOLS::Flavour m_flav1,m_flav2;
  public:
    DeltaPhi_Selector(Process_Base *const);
    ~DeltaPhi_Selector();
    void     SetRange(ATOOLS::Flavour,ATOOLS::Flavour,double,double);
    bool     Trigger(ATOOLS::Selector_List &);
    void     BuildCuts(Cut_Data *);
  };

  class DeltaR_Selector : public Selector_Base {
    double  m_dRmin, m_dRmax;
    ATOOLS::Flavour m_flav1,m_flav2;
  public:
    DeltaR_Selector(Process_Base *const);
    ~DeltaR_Selector();
    void     SetRange(ATOOLS::Flavour,ATOOLS::Flavour,double,double);
    bool     Trigger(ATOOLS::Selector_List &);
    void     BuildCuts(Cut_Data *);
  };

  class DeltaRy_Selector : public Selector_Base {
    double  m_dRymin, m_dRymax;
    ATOOLS::Flavour m_flav1,m_flav2;
  public:
    DeltaRy_Selector(Process_Base *const);
    ~DeltaRy_Selector();
    void     SetRange(ATOOLS::Flavour,ATOOLS::Flavour,double,double);
    bool     Trigger(ATOOLS::Selector_List &);
    void     BuildCuts(Cut_Data *);
  };

  class PhiStar_Selector : public Selector_Base {
    double  m_phistarmin, m_phistarmax;
    ATOOLS::Flavour m_flav1,m_flav2;
  public:
    PhiStar_Selector(Process_Base *const);
    ~PhiStar_Selector();
    void     SetRange(ATOOLS::Flavour,ATOOLS::Flavour,double,double);
    bool     Trigger(ATOOLS::Selector_List &);
    void     BuildCuts(Cut_Data *);
  };

  // -----------------------------
  // inclusive selectors
  // -----------------------------
  class Multiplicity_Selector : public Selector_Base {
    size_t  m_nmin, m_nmax;
    ATOOLS::Flavour m_flav;
  public:
    Multiplicity_Selector(Process_Base *const);
    ~Multiplicity_Selector();
    void     SetRange(ATOOLS::Flavour,size_t,size_t);
    bool     Trigger(ATOOLS::Selector_List &);
    void     BuildCuts(Cut_Data *);
  };

  class PTMIS_Selector : public Selector_Base {
    double  m_ptmismin, m_ptmismax;
    std::vector<ATOOLS::Flavour> m_flavs;
  public:
    PTMIS_Selector(Process_Base *const);
    ~PTMIS_Selector();
    void     SetRange(double,double);
    bool     Trigger(ATOOLS::Selector_List &);
    void     BuildCuts(Cut_Data *);
  };

  class ETMIS_Selector : public Selector_Base {
    double  m_etmismin, m_etmismax;
    std::vector<ATOOLS::Flavour> m_flavs;
  public:
    ETMIS_Selector(Process_Base *const);
    ~ETMIS_Selector();
    void     SetRange(double,double);
    bool     Trigger(ATOOLS::Selector_List &);
    void     BuildCuts(Cut_Data *);
  };

  class Isolation_Cut : public Selector_Base {
    double m_dR,m_exp,m_emax,m_massmax;
    std::vector<size_t> m_vf;
    ATOOLS::Flavour m_iflav;

    double Chi(double eg,double dr);
    double DR(const ATOOLS::Vec4D & p1,const ATOOLS::Vec4D & p2);
    double DEta12(const ATOOLS::Vec4D &,const ATOOLS::Vec4D &);
    double DPhi12(const ATOOLS::Vec4D &,const ATOOLS::Vec4D &);
  public:
    Isolation_Cut(Process_Base *const);

    void   SetRange(ATOOLS::Flavour,double,double,double,double);
    bool   Trigger(ATOOLS::Selector_List &);
    void   BuildCuts(Cut_Data *);
  };

  class NJettiness_Selector : public Selector_Base {
    size_t  m_N;
    double  m_njmin, m_njmax;
  public:
    NJettiness_Selector(Process_Base *const);
    ~NJettiness_Selector();
    void     SetRange(size_t,double,double);
    bool     Trigger(ATOOLS::Selector_List &);
    void     BuildCuts(Cut_Data *);
  };

}

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include <algorithm>

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

class Order_Y {
public:
  bool operator()(const Vec4D &a,const Vec4D &b)
  {
    return a.Y()>b.Y();
  }
};

/*--------------------------------------------------------------------

  Transverse Momentum Selector

  --------------------------------------------------------------------*/

PT_Selector::PT_Selector(Process_Base *const proc):
  Selector_Base("PT_Selector",proc), m_ptmin(0.), m_ptmax(0.),
  m_flav(Flavour(kf_none))
{
}

PT_Selector::~PT_Selector() {
}

bool PT_Selector::Trigger(Selector_List &sl)
{
  DEBUG_FUNC(m_on);
  if (!m_on) return true;
  for (size_t i=m_nin;i<sl.size();i++) {
    if (m_flav.Includes(sl[i].Flavour())) {
      double pti = sl[i].Momentum().PPerp();
      if (m_sel_log->Hit( ((pti<m_ptmin) || (pti>m_ptmax)) )) return false;
    }
  }
  return true;
}

void PT_Selector::BuildCuts(Cut_Data * cuts)
{
  if (m_isnlo) return;
  double sumM2=0.;
  for (int i=m_nin;i<m_n;i++) {
    sumM2+=sqr(p_fl[i].SelMass());
  }
  for (int i=m_nin;i<m_n;i++) {
    if (m_flav.Includes(p_fl[i])) {
      cuts->energymin[i] = Max(sqrt(sqr(m_ptmin)+sqr(p_fl[i].SelMass())),
                               cuts->energymin[i]);
      double Emax2 = sqr((m_smax+2.*sqr(p_fl[i].SelMass())-sumM2)
                       /(2.*sqrt(m_smax)));
      double cosmax = Min(cuts->cosmax[0][i],
                          sqrt(1.-sqr(m_ptmin)/(Emax2-sqr(p_fl[i].SelMass()))));
      cuts->cosmax[0][i] = cuts->cosmax[1][i] = cosmax;
      cuts->cosmax[i][0] = cuts->cosmax[i][1] = cosmax;
      cuts->etmin[i] = Max(sqrt(sqr(m_ptmin)+sqr(p_fl[i].SelMass())
                           *(1.-sqr(cuts->cosmax[0][i]))),cuts->etmin[i]);
    }
  }
}

void PT_Selector::SetRange(Flavour flav,double min,double max)
{
  m_flav=flav;
  m_ptmin=min;
  m_ptmax=max;

  for (size_t i=m_nin;i<m_n;i++) {
    if (m_flav.Includes(p_fl[i])) {
      m_on=true;
      // if more than one flav is found maybe increase m_smin
      m_smin = Max(m_smin,4.*sqr(m_ptmin));
    }
  }

  msg_Debugging()<<"flav="<<m_flav
                 <<", min="<<m_ptmin<<", max="<<m_ptmax
                 <<" -> smin="<<m_smin<<", on="<<m_on
                 <<std::endl;
}

DECLARE_ND_GETTER(PT_Selector,"PT",Selector_Base,Selector_Key,true);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,PT_Selector>::
operator()(const Selector_Key &key) const
{
  Scoped_Settings s{ key.m_settings };
  const auto parameters = s.SetDefault<std::string>({}).GetVector<std::string>();
  if (parameters[0] == "PTNLO")
    msg_Out()
      << "WARNING: Substituting PT selector for missing PTNLO selector\n";
  else
    assert(parameters[0] == "PT");
  if (parameters.size() != 4)
    THROW(critical_error, "Invalid syntax");
  const auto kf = s.Interprete<int>(parameters[1]);
  const auto min = s.Interprete<double>(parameters[2]);
  const auto max = s.Interprete<double>(parameters[3]);
  Flavour flav = Flavour((kf_code)abs(kf),kf<0);
  PT_Selector *sel = new PT_Selector(key.p_proc);
  sel->SetRange(flav,min,max);
  return sel;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,PT_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"transverse momentum selector";
}



/*--------------------------------------------------------------------

  Transverse Energy Selector

  --------------------------------------------------------------------*/

ET_Selector::ET_Selector(Process_Base *const proc):
  Selector_Base("ET_Selector",proc), m_etmin(0.), m_etmax(0.),
  m_flav(Flavour(kf_none))
{
}

ET_Selector::~ET_Selector() {
}

bool ET_Selector::Trigger(Selector_List &sl)
{
  DEBUG_FUNC(m_on);
  if (!m_on) return true;
  for (size_t i=m_nin;i<sl.size();i++) {
    if (m_flav.Includes(sl[i].Flavour())) {
      double eti = sl[i].Momentum().MPerp();
      if (m_sel_log->Hit( ((eti<m_etmin) || (eti>m_etmax)) )) return false;
    }
  }
  return true;
}

void ET_Selector::BuildCuts(Cut_Data * cuts)
{
  if (m_isnlo) return;
  double sumM2=0.;
  for (int i=m_nin;i<m_n;i++) {
    sumM2+=sqr(p_fl[i].SelMass());
  }
  for (int i=m_nin;i<m_n;i++) {
    if (m_flav.Includes(p_fl[i])) {
      cuts->energymin[i] = Max(sqrt(sqr(m_etmin)+sqr(p_fl[i].SelMass())),
                               cuts->energymin[i]);
      double Emax2 = sqr((m_smax+2.*sqr(p_fl[i].SelMass())-sumM2)
                       /(2.*sqrt(m_smax)));
      double cosmax = Min(cuts->cosmax[0][i],
                          sqrt(1.-sqr(m_etmin)/(Emax2-sqr(p_fl[i].SelMass()))));
      cuts->cosmax[0][i] = cuts->cosmax[1][i] = cosmax;
      cuts->cosmax[i][0] = cuts->cosmax[i][1] = cosmax;
      cuts->etmin[i] = Max(sqrt(sqr(m_etmin)+sqr(p_fl[i].SelMass())
                           *(1.-sqr(cuts->cosmax[0][i]))),cuts->etmin[i]);
    }
  }
}

void ET_Selector::SetRange(Flavour flav,double min,double max)
{
  m_flav=flav;
  m_etmin=min;
  m_etmax=max;

  for (size_t i=m_nin;i<m_n;i++) {
    if (m_flav.Includes(p_fl[i])) {
      m_on=true;
      // if more than one flav is found maybe increase m_smin
      m_smin = Max(m_smin,4.*sqr(m_etmin));
    }
  }

  msg_Debugging()<<"flav="<<m_flav
                 <<", min="<<m_etmin<<", max="<<m_etmax
                 <<" -> smin="<<m_smin<<", on="<<m_on
                 <<std::endl;
}

DECLARE_ND_GETTER(ET_Selector,"ET",Selector_Base,Selector_Key,true);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,ET_Selector>::
operator()(const Selector_Key &key) const
{
  Scoped_Settings s{ key.m_settings };
  const auto parameters = s.SetDefault<std::string>({}).GetVector<std::string>();
  assert(parameters[0] == "ET");
  if (parameters.size() != 4)
    THROW(critical_error, "Invalid syntax");
  const auto kf = s.Interprete<int>(parameters[1]);
  const auto min = s.Interprete<double>(parameters[2]);
  const auto max = s.Interprete<double>(parameters[3]);
  Flavour flav = Flavour((kf_code)abs(kf),kf<0);
  ET_Selector *sel = new ET_Selector(key.p_proc);
  sel->SetRange(flav,min,max);
  return sel;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,ET_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"transverse energy selector";
}



/*--------------------------------------------------------------------

  Rapidity Selector

  --------------------------------------------------------------------*/

Rapidity_Selector::Rapidity_Selector(Process_Base *const proc):
  Selector_Base("Rapidity_Selector",proc), m_ymin(0.), m_ymax(0.),
  m_flav(Flavour(kf_none))
{
}

Rapidity_Selector::~Rapidity_Selector() {
}

bool Rapidity_Selector::Trigger(Selector_List &sl)
{
  DEBUG_FUNC(m_on);
  if (!m_on) return true;
  for (size_t i=m_nin;i<sl.size();i++) {
    if (sl[i].Momentum()==Vec4D(0.,0.,0.,0.)) continue;
    if (m_flav.Includes(sl[i].Flavour())) {
      double yi = sl[i].Momentum().Y();
      if (m_sel_log->Hit( ((yi<m_ymin) || (yi>m_ymax)) )) return false;
    }
  }
  return true;
}

void Rapidity_Selector::BuildCuts(Cut_Data * cuts)
{
  if (m_isnlo) return;
  for (int i=m_nin;i<m_n;i++) {
    if (m_flav.Includes(p_fl[i])) {
      cuts->cosmax[0][i] = cuts->cosmax[i][0] =
        Min(cuts->cosmax[0][i],tanh(m_ymax)/sqrt(1.-sqr(p_fl[i].SelMass())
                                                 /sqr(cuts->energymin[i])));
      cuts->cosmax[1][i] = cuts->cosmax[i][1] =
        Min(cuts->cosmax[1][i],tanh(-m_ymin)/sqrt(1.-sqr(p_fl[i].SelMass())
                                                  /sqr(cuts->energymin[i])));
    }
  }
}

void Rapidity_Selector::SetRange(Flavour flav,double min,double max)
{
  m_flav=flav;
  m_ymin=min;
  m_ymax=max;

  for (size_t i=m_nin;i<m_n;i++) {
    if (m_flav.Includes(p_fl[i])) {
      m_on=true;
    }
  }

  msg_Debugging()<<"flav="<<m_flav
                 <<", min="<<m_ymin<<", max="<<m_ymax
                 <<" -> smin="<<m_smin<<", on="<<m_on
                 <<std::endl;
}

DECLARE_ND_GETTER(Rapidity_Selector,"Y",Selector_Base,Selector_Key,true);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,Rapidity_Selector>::
operator()(const Selector_Key &key) const
{
  Scoped_Settings s{ key.m_settings };
  const auto parameters = s.SetDefault<std::string>({}).GetVector<std::string>();
  assert(parameters[0] == "Y");
  if (parameters.size() != 4)
    THROW(critical_error, "Invalid syntax");
  const auto kf = s.Interprete<int>(parameters[1]);
  const auto min = s.Interprete<double>(parameters[2]);
  const auto max = s.Interprete<double>(parameters[3]);
  Flavour flav = Flavour((kf_code)abs(kf),kf<0);
  Rapidity_Selector *sel = new Rapidity_Selector(key.p_proc);
  sel->SetRange(flav,min,max);
  return sel;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,Rapidity_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"rapidity selector";
}



/*--------------------------------------------------------------------

  PseudoRapidity Selector

  --------------------------------------------------------------------*/

PseudoRapidity_Selector::PseudoRapidity_Selector(Process_Base *const proc):
  Selector_Base("PseudoRapidity_Selector",proc), m_etamin(0.), m_etamax(0.),
  m_flav(Flavour(kf_none))
{
}

PseudoRapidity_Selector::~PseudoRapidity_Selector() {
}

bool PseudoRapidity_Selector::Trigger(Selector_List &sl)
{
  DEBUG_FUNC(m_on);
  if (!m_on) return true;
  for (size_t i=m_nin;i<sl.size();i++) {
    if (sl[i].Momentum()==Vec4D(0.,0.,0.,0.)) continue;
    if (m_flav.Includes(sl[i].Flavour())) {
      double etai = sl[i].Momentum().Eta();
      if (m_sel_log->Hit( ((etai<m_etamin) || (etai>m_etamax)) )) return false;
    }
  }
  return true;
}

void PseudoRapidity_Selector::BuildCuts(Cut_Data * cuts)
{
  if (m_isnlo) return;
  for (int i=m_nin;i<m_n;i++) {
    if (m_flav.Includes(p_fl[i])) {
      cuts->cosmin[1][i] = cuts->cosmin[i][1]
        = Max(cuts->cosmin[1][i],tanh(-m_etamax));
      cuts->cosmin[0][i] = cuts->cosmin[i][0]
        = Max(cuts->cosmin[0][i],tanh(m_etamin));
      cuts->cosmax[0][i] = cuts->cosmax[i][0]
        = Min(cuts->cosmax[0][i],tanh(m_etamax));
      cuts->cosmax[1][i] = cuts->cosmax[i][1]
        = Min(cuts->cosmax[1][i],tanh(-m_etamin));
    }
  }
}

void PseudoRapidity_Selector::SetRange(Flavour flav,double min,double max)
{
  m_flav=flav;
  m_etamin=min;
  m_etamax=max;

  for (size_t i=m_nin;i<m_n;i++) {
    if (m_flav.Includes(p_fl[i])) {
      m_on=true;
    }
  }

  msg_Debugging()<<"flav="<<m_flav
                 <<", min="<<m_etamin<<", max="<<m_etamax
                 <<" -> smin="<<m_smin<<", on="<<m_on
                 <<std::endl;
}

DECLARE_ND_GETTER(PseudoRapidity_Selector,"Eta",Selector_Base,Selector_Key,true);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,PseudoRapidity_Selector>::
operator()(const Selector_Key &key) const
{
  Scoped_Settings s{ key.m_settings };
  const auto parameters = s.SetDefault<std::string>({}).GetVector<std::string>();
  assert(parameters[0] == "Eta");
  if (parameters.size() != 4)
    THROW(critical_error, "Invalid syntax");
  const auto kf = s.Interprete<int>(parameters[1]);
  const auto min = s.Interprete<double>(parameters[2]);
  const auto max = s.Interprete<double>(parameters[3]);
  Flavour flav = Flavour((kf_code)abs(kf),kf<0);
  PseudoRapidity_Selector *sel = new PseudoRapidity_Selector(key.p_proc);
  sel->SetRange(flav,min,max);
  return sel;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,PseudoRapidity_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"pseudorapidity selector";
}



/*--------------------------------------------------------------------

  Invariant Mass Selector

 --------------------------------------------------------------------*/

IMass_Selector::IMass_Selector(Process_Base *const proc):
  Selector_Base("Mass_Selector",proc), m_massmin(0.), m_massmax(0.),
  m_flav1(Flavour(kf_none)), m_flav2(Flavour(kf_none))
{
}

IMass_Selector::~IMass_Selector() {
}

bool IMass_Selector::Trigger(Selector_List &sl)
{
  DEBUG_FUNC(m_on);
  if (!m_on) return true;
  for (size_t i=m_nin;i<sl.size();i++) {
    for (size_t j=i+1;j<sl.size();j++) {
      if ( (m_flav1.Includes(sl[i].Flavour()) &&
            m_flav2.Includes(sl[j].Flavour()) ) ||
           (m_flav1.Includes(sl[j].Flavour()) &&
            m_flav2.Includes(sl[i].Flavour()) ) ) {
        double massij = sqrt((sl[i].Momentum()+sl[j].Momentum()).Abs2());
        if (m_sel_log->Hit( ((massij<m_massmin) ||
                             (massij>m_massmax)) )) return false;
      }
    }
  }
  return true;
}

void IMass_Selector::BuildCuts(Cut_Data * cuts)
{
  if (m_isnlo) return;
  for (size_t i=m_nin;i<m_n;i++) {
    for (size_t j=i+1;j<m_n;j++) {
      if ( (m_flav1.Includes(p_fl[i]) && m_flav2.Includes(p_fl[j])) ||
           (m_flav1.Includes(p_fl[j]) && m_flav2.Includes(p_fl[i])) ) {
        cuts->scut[i][j] = cuts->scut[j][i] = Max(cuts->scut[i][j],sqr(m_massmin));
      }
    }
  }
}

void IMass_Selector::SetRange(Flavour flav1,Flavour flav2,double min,double max)
{
  m_flav1=flav1;
  m_flav2=flav2;
  m_massmin=Max(min,m_flav1.SelMass()+m_flav2.SelMass());
  m_massmax=Min(max,(rpa->gen.PBeam(0)[0]+rpa->gen.PBeam(1)[0]));

  for (size_t i=m_nin;i<m_n;i++) {
    for (size_t j=i+1;j<m_n;j++) {
      if ( (m_flav1.Includes(p_fl[i]) && m_flav2.Includes(p_fl[j])) ||
           (m_flav1.Includes(p_fl[j]) && m_flav2.Includes(p_fl[i])) ) {
        m_on=true;
        // if more than one pair is found maybe increase m_smin
        if (sqr(m_massmin)>m_smin) m_smin = Max(sqr(m_massmin),m_smin);
      }
    }
  }
}

DECLARE_ND_GETTER(IMass_Selector,"Mass",Selector_Base,Selector_Key,true);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,IMass_Selector>::
operator()(const Selector_Key &key) const
{
  Scoped_Settings s{ key.m_settings };
  const auto parameters = s.SetDefault<std::string>({}).GetVector<std::string>();
  assert(parameters[0] == "Mass");
  if (parameters.size() != 5)
    THROW(critical_error, "Invalid syntax");
  const auto kf1 = s.Interprete<int>(parameters[1]);
  const auto kf2 = s.Interprete<int>(parameters[2]);
  const auto min = s.Interprete<double>(parameters[3]);
  const auto max = s.Interprete<double>(parameters[4]);
  Flavour flav1 = Flavour((kf_code)abs(kf1),kf1<0);
  Flavour flav2 = Flavour((kf_code)abs(kf2),kf2<0);
  IMass_Selector *sel = new IMass_Selector(key.p_proc);
  sel->SetRange(flav1,flav2,min,max);
  return sel;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,IMass_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"invariant mass selector";
}



/*--------------------------------------------------------------------

  Virtuality Selector

  --------------------------------------------------------------------*/

IQ2_Selector::IQ2_Selector(Process_Base *const proc):
  Selector_Base("Q2_Selector",proc), m_q2min(0.), m_q2max(0.),
  m_flav1(Flavour(kf_none)), m_flav2(Flavour(kf_none))
{
}

IQ2_Selector::~IQ2_Selector() {
}

bool IQ2_Selector::Trigger(Selector_List &sl)
{
  if (!m_on) return true;
  for (int i=0;i<m_nin;i++) {
    for (int j=m_nin;j<sl.size();j++) {
      if ( (m_flav1.Includes(sl[i].Flavour()) &&
            m_flav2.Includes(sl[j].Flavour())) ||
           (m_flav1.Includes(sl[j].Flavour()) &&
            m_flav2.Includes(sl[i].Flavour())) ) {
        double q2ij = -(sl[i].Momentum()-sl[j].Momentum()).Abs2();
        if (m_sel_log->Hit( ((q2ij < m_q2min) ||
                             (q2ij > m_q2max)) )) return false;
      }
    }
  }
  return true;
}

void IQ2_Selector::BuildCuts(Cut_Data * cuts)
{
  if (m_isnlo) return;
  for (int i=0;i<m_nin;i++) {
    for (int j=m_nin;j<m_n;j++) {
      if ( (m_flav1.Includes(p_fl[i]) && m_flav2.Includes(p_fl[j])) ||
           (m_flav1.Includes(p_fl[j]) && m_flav2.Includes(p_fl[i])) ) {
        cuts->scut[i][j] = Min(cuts->scut[i][j],-m_q2min);
      }
    }
  }
}

void IQ2_Selector::SetRange(Flavour flav1,Flavour flav2,double min,double max)
{
  m_flav1=flav1;
  m_flav2=flav2;
  m_q2min=min;
  m_q2max=max;
  m_on=true;
}

DECLARE_ND_GETTER(IQ2_Selector,"Q2",Selector_Base,Selector_Key,true);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,IQ2_Selector>::
operator()(const Selector_Key &key) const
{
  Scoped_Settings s{ key.m_settings };
  const auto parameters = s.SetDefault<std::string>({}).GetVector<std::string>();
  assert(parameters[0] == "Q2");
  if (parameters.size() != 5)
    THROW(critical_error, "Invalid syntax");
  const auto kf1 = s.Interprete<int>(parameters[1]);
  const auto kf2 = s.Interprete<int>(parameters[2]);
  const auto min = s.Interprete<double>(parameters[3]);
  const auto max = s.Interprete<double>(parameters[4]);
  Flavour flav1 = Flavour((kf_code)abs(kf1),kf1<0);
  Flavour flav2 = Flavour((kf_code)abs(kf2),kf2<0);
  IQ2_Selector *sel = new IQ2_Selector(key.p_proc);
  sel->SetRange(flav1,flav2,min,max);
  return sel;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,IQ2_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Q2 selector";
}



/*--------------------------------------------------------------------

  PT2 Selector

  --------------------------------------------------------------------*/

PT2_Selector::PT2_Selector(Process_Base *const proc):
  Selector_Base("PT2_Selector",proc), m_ptmin(0.), m_ptmax(0.),
  m_flav1(Flavour(kf_none)), m_flav2(Flavour(kf_none))
{
}

PT2_Selector::~PT2_Selector() {
}

bool PT2_Selector::Trigger(Selector_List &sl)
{
  DEBUG_FUNC(m_on);
  if (!m_on) return true;
  for (size_t i=m_nin;i<sl.size();i++) {
    for (size_t j=i+1;j<sl.size();j++) {
      if ( (m_flav1.Includes(sl[i].Flavour()) &&
            m_flav2.Includes(sl[j].Flavour())) ||
           (m_flav1.Includes(sl[j].Flavour()) &&
            m_flav2.Includes(sl[i].Flavour())) ) {
        double ptij = (sl[i].Momentum()+sl[j].Momentum()).PPerp();
        if (m_sel_log->Hit( ((ptij<m_ptmin) || (ptij>m_ptmax)) )) return false;
      }
    }
  }
  return true;
}

void PT2_Selector::BuildCuts(Cut_Data * cuts)
{
}

void PT2_Selector::SetRange(Flavour flav1,Flavour flav2,double min,double max)
{
  m_flav1=flav1;
  m_flav2=flav2;
  m_ptmin=min;
  m_ptmax=max;

  for (size_t i=m_nin;i<m_n;i++) {
    for (size_t j=i+1;j<m_n;j++) {
      if ( (m_flav1.Includes(p_fl[i]) && m_flav2.Includes(p_fl[j])) ||
           (m_flav1.Includes(p_fl[j]) && m_flav2.Includes(p_fl[i])) ) {
        m_on=true;
        // if more than one pair is found maybe increase m_smin
        m_smin = Max(m_smin,4.*sqr(m_ptmin));
      }
    }
  }
}

DECLARE_ND_GETTER(PT2_Selector,"PT2",Selector_Base,Selector_Key,true);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,PT2_Selector>::
operator()(const Selector_Key &key) const
{
  Scoped_Settings s{ key.m_settings };
  const auto parameters = s.SetDefault<std::string>({}).GetVector<std::string>();
  assert(parameters[0] == "PT2");
  if (parameters.size() != 5)
    THROW(critical_error, "Invalid syntax");
  const auto kf1 = s.Interprete<int>(parameters[1]);
  const auto kf2 = s.Interprete<int>(parameters[2]);
  const auto min = s.Interprete<double>(parameters[3]);
  const auto max = s.Interprete<double>(parameters[4]);
  Flavour flav1 = Flavour((kf_code)abs(kf1),kf1<0);
  Flavour flav2 = Flavour((kf_code)abs(kf2),kf2<0);
  PT2_Selector *sel = new PT2_Selector(key.p_proc);
  sel->SetRange(flav1,flav2,min,max);
  return sel;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,PT2_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"transverse momentum selector";
}



/*--------------------------------------------------------------------

  MT2 Selector

  --------------------------------------------------------------------*/

MT2_Selector::MT2_Selector(Process_Base *const proc):
  Selector_Base("MT2_Selector",proc), m_mtmin(0.), m_mtmax(0.),
  m_flav1(Flavour(kf_none)), m_flav2(Flavour(kf_none))
{
}

MT2_Selector::~MT2_Selector() {
}

bool MT2_Selector::Trigger(Selector_List &sl)
{
  DEBUG_FUNC(m_on);
  if (!m_on) return true;
  for (size_t i=m_nin;i<sl.size();i++) {
    for (size_t j=i+1;j<sl.size();j++) {
      if ( (m_flav1.Includes(sl[i].Flavour()) &&
            m_flav2.Includes(sl[j].Flavour())) ||
           (m_flav1.Includes(sl[j].Flavour()) &&
            m_flav2.Includes(sl[i].Flavour())) ) {
        double mtij = (sl[i].Momentum()+sl[j].Momentum()).MPerp();
        if (m_sel_log->Hit( ((mtij<m_mtmin) || (mtij>m_mtmax)) )) return false;
      }
    }
  }
  return true;
}

void MT2_Selector::BuildCuts(Cut_Data * cuts)
{
}

void MT2_Selector::SetRange(Flavour flav1,Flavour flav2,double min,double max)
{
  m_flav1=flav1;
  m_flav2=flav2;
  m_mtmin=min;
  m_mtmax=max;

  for (size_t i=m_nin;i<m_n;i++) {
    for (size_t j=i+1;j<m_n;j++) {
      if ( (m_flav1.Includes(p_fl[i]) && m_flav2.Includes(p_fl[j])) ||
           (m_flav1.Includes(p_fl[j]) && m_flav2.Includes(p_fl[i])) ) {
        m_on=true;
        // if more than one pair is found maybe increase m_smin
        m_smin = Max(m_smin,4.*sqr(m_mtmin));
      }
    }
  }
}

DECLARE_ND_GETTER(MT2_Selector,"MT2",Selector_Base,Selector_Key,true);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,MT2_Selector>::
operator()(const Selector_Key &key) const
{
  Scoped_Settings s{ key.m_settings };
  const auto parameters = s.SetDefault<std::string>({}).GetVector<std::string>();
  assert(parameters[0] == "MT2");
  if (parameters.size() != 5)
    THROW(critical_error, "Invalid syntax");
  const auto kf1 = s.Interprete<int>(parameters[1]);
  const auto kf2 = s.Interprete<int>(parameters[2]);
  const auto min = s.Interprete<double>(parameters[3]);
  const auto max = s.Interprete<double>(parameters[4]);
  Flavour flav1 = Flavour((kf_code)abs(kf1),kf1<0);
  Flavour flav2 = Flavour((kf_code)abs(kf2),kf2<0);
  MT2_Selector *sel = new MT2_Selector(key.p_proc);
  sel->SetRange(flav1,flav2,min,max);
  return sel;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,MT2_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"MT2 selector [ mT = sqrt( m(1,2)^2 + pT(1,2)^2 ) ]";
}



/*--------------------------------------------------------------------

  MT2_v2 Selector

  --------------------------------------------------------------------*/

MT2_v2_Selector::MT2_v2_Selector(Process_Base *const proc):
  Selector_Base("MT2_v2_Selector",proc), m_mtmin(0.), m_mtmax(0.),
  m_flav1(Flavour(kf_none)), m_flav2(Flavour(kf_none))
{
}

MT2_v2_Selector::~MT2_v2_Selector() {
}

bool MT2_v2_Selector::Trigger(Selector_List &sl)
{
  DEBUG_FUNC(m_on);
  if (!m_on) return true;
  for (size_t i=m_nin;i<sl.size();i++) {
    for (size_t j=i+1;j<sl.size();j++) {
      if ( (m_flav1.Includes(sl[i].Flavour()) &&
            m_flav2.Includes(sl[j].Flavour())) ||
           (m_flav1.Includes(sl[j].Flavour()) &&
            m_flav2.Includes(sl[i].Flavour())) ) {
        double mtij = sqrt(2.*sl[i].Momentum().PPerp()*sl[j].Momentum().PPerp()
                           *(1.-sl[i].Momentum().CosDPhi(sl[j].Momentum())));
        if (m_sel_log->Hit( ((mtij<m_mtmin) || (mtij>m_mtmax)) )) return false;
      }
    }
  }
  return true;
}

void MT2_v2_Selector::BuildCuts(Cut_Data * cuts)
{
}

void MT2_v2_Selector::SetRange(Flavour flav1,Flavour flav2,double min,double max)
{
  m_flav1=flav1;
  m_flav2=flav2;
  m_mtmin=min;
  m_mtmax=max;

  for (size_t i=m_nin;i<m_n;i++) {
    for (size_t j=i+1;j<m_n;j++) {
      if ( (m_flav1.Includes(p_fl[i]) && m_flav2.Includes(p_fl[j])) ||
           (m_flav1.Includes(p_fl[j]) && m_flav2.Includes(p_fl[i])) ) {
        m_on=true;
        // if more than one pair is found maybe increase m_smin
        m_smin = Max(m_smin,4.*sqr(m_mtmin));
      }
    }
  }
}

DECLARE_ND_GETTER(MT2_v2_Selector,"MT2_v2",Selector_Base,Selector_Key,true);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,MT2_v2_Selector>::
operator()(const Selector_Key &key) const
{
  Scoped_Settings s{ key.m_settings };
  const auto parameters = s.SetDefault<std::string>({}).GetVector<std::string>();
  assert(parameters[0] == "MT2_v2");
  if (parameters.size() != 5)
    THROW(critical_error, "Invalid syntax");
  const auto kf1 = s.Interprete<int>(parameters[1]);
  const auto kf2 = s.Interprete<int>(parameters[2]);
  const auto min = s.Interprete<double>(parameters[3]);
  const auto max = s.Interprete<double>(parameters[4]);
  Flavour flav1 = Flavour((kf_code)abs(kf1),kf1<0);
  Flavour flav2 = Flavour((kf_code)abs(kf2),kf2<0);
  MT2_v2_Selector *sel = new MT2_v2_Selector(key.p_proc);
  sel->SetRange(flav1,flav2,min,max);
  return sel;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,MT2_v2_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"MT2 selector [ mT = sqrt( 2 pT1 pT2 (1-cos(dphi(1,2)) ) ]";
}



/*--------------------------------------------------------------------

  DeltaY Selector

  --------------------------------------------------------------------*/

DeltaY_Selector::DeltaY_Selector(Process_Base *const proc):
  Selector_Base("DeltaY_Selector",proc), m_dymin(0.), m_dymax(0.),
  m_flav1(Flavour(kf_none)), m_flav2(Flavour(kf_none))
{
}

DeltaY_Selector::~DeltaY_Selector() {
}

bool DeltaY_Selector::Trigger(Selector_List &sl)
{
  DEBUG_FUNC(m_on);
  if (!m_on) return true;
  for (size_t i=m_nin;i<sl.size();i++) {
    for (size_t j=i+1;j<sl.size();j++) {
      if (sl[i].Momentum()==Vec4D(0.,0.,0.,0.)) continue;
      if ( (m_flav1.Includes(sl[i].Flavour()) &&
            m_flav2.Includes(sl[j].Flavour())) ||
           (m_flav1.Includes(sl[j].Flavour()) &&
            m_flav2.Includes(sl[i].Flavour())) ) {
        double dyij = dabs(sl[i].Momentum().DY(sl[j].Momentum()));
        if (m_sel_log->Hit( ((dyij<m_dymin) || (dyij>m_dymax)) )) return false;
      }
    }
  }
  return true;
}

void DeltaY_Selector::BuildCuts(Cut_Data * cuts)
{
}

void DeltaY_Selector::SetRange(Flavour fl1,Flavour fl2,double min,double max)
{
  m_flav1=fl1;
  m_flav2=fl2;
  m_dymin=min;
  m_dymax=max;

  for (size_t i=m_nin;i<m_n;i++) {
    for (size_t j=i+1;j<m_n;j++) {
      if ( (m_flav1.Includes(p_fl[i]) && m_flav2.Includes(p_fl[j])) ||
           (m_flav1.Includes(p_fl[j]) && m_flav2.Includes(p_fl[i])) ) {
        m_on=true;
      }
    }
  }
}

DECLARE_ND_GETTER(DeltaY_Selector,"DY",Selector_Base,Selector_Key,true);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,DeltaY_Selector>::
operator()(const Selector_Key &key) const
{
  Scoped_Settings s{ key.m_settings };
  const auto parameters = s.SetDefault<std::string>({}).GetVector<std::string>();
  assert(parameters[0] == "DY");
  if (parameters.size() != 5)
    THROW(critical_error, "Invalid syntax");
  const auto kf1 = s.Interprete<int>(parameters[1]);
  const auto kf2 = s.Interprete<int>(parameters[2]);
  const auto min = s.Interprete<double>(parameters[3]);
  const auto max = s.Interprete<double>(parameters[4]);
  Flavour flav1 = Flavour((kf_code)abs(kf1),kf1<0);
  Flavour flav2 = Flavour((kf_code)abs(kf2),kf2<0);
  DeltaY_Selector *sel = new DeltaY_Selector(key.p_proc);
  sel->SetRange(flav1,flav2,min,max);
  return sel;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,DeltaY_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"rapidity separation selector";
}



/*--------------------------------------------------------------------

  DeltaEta Selector

  --------------------------------------------------------------------*/

DeltaEta_Selector::DeltaEta_Selector(Process_Base *const proc):
  Selector_Base("DeltaEta_Selector",proc), m_detamin(0.), m_detamax(0.),
  m_flav1(Flavour(kf_none)), m_flav2(Flavour(kf_none))
{
}

DeltaEta_Selector::~DeltaEta_Selector() {
}

bool DeltaEta_Selector::Trigger(Selector_List &sl)
{
  DEBUG_FUNC(m_on);
  if (!m_on) return true;
  for (size_t i=m_nin;i<sl.size();i++) {
    for (size_t j=i+1;j<sl.size();j++) {
      if (sl[i].Momentum()==Vec4D(0.,0.,0.,0.)) continue;
      if ( (m_flav1.Includes(sl[i].Flavour()) &&
            m_flav2.Includes(sl[j].Flavour())) ||
           (m_flav1.Includes(sl[j].Flavour()) &&
            m_flav2.Includes(sl[i].Flavour())) ) {
        double detaij = dabs(sl[i].Momentum().DEta(sl[j].Momentum()));
        if (m_sel_log->Hit( ((detaij<m_detamin) || (detaij>m_detamax)) )) return false;
      }
    }
  }
  return true;
}

void DeltaEta_Selector::BuildCuts(Cut_Data * cuts)
{
}

void DeltaEta_Selector::SetRange(Flavour fl1,Flavour fl2,double min,double max)
{
  m_flav1=fl1;
  m_flav2=fl2;
  m_detamin=min;
  m_detamax=max;

  for (size_t i=m_nin;i<m_n;i++) {
    for (size_t j=i+1;j<m_n;j++) {
      if ( (m_flav1.Includes(p_fl[i]) && m_flav2.Includes(p_fl[j])) ||
           (m_flav1.Includes(p_fl[j]) && m_flav2.Includes(p_fl[i])) ) {
        m_on=true;
      }
    }
  }
}

DECLARE_ND_GETTER(DeltaEta_Selector,"DEta",Selector_Base,Selector_Key,true);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,DeltaEta_Selector>::
operator()(const Selector_Key &key) const
{
  Scoped_Settings s{ key.m_settings };
  const auto parameters = s.SetDefault<std::string>({}).GetVector<std::string>();
  assert(parameters[0] == "DEta");
  if (parameters.size() != 5)
    THROW(critical_error, "Invalid syntax");
  const auto kf1 = s.Interprete<int>(parameters[1]);
  const auto kf2 = s.Interprete<int>(parameters[2]);
  const auto min = s.Interprete<double>(parameters[3]);
  const auto max = s.Interprete<double>(parameters[4]);
  Flavour flav1 = Flavour((kf_code)abs(kf1),kf1<0);
  Flavour flav2 = Flavour((kf_code)abs(kf2),kf2<0);
  DeltaEta_Selector *sel = new DeltaEta_Selector(key.p_proc);
  sel->SetRange(flav1,flav2,min,max);
  return sel;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,DeltaEta_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"pseudorapidity separation selector";
}



/*--------------------------------------------------------------------

  DeltaPhi Selector

  --------------------------------------------------------------------*/

DeltaPhi_Selector::DeltaPhi_Selector(Process_Base *const proc):
  Selector_Base("DeltaPhi_Selector",proc), m_dphimin(0.), m_dphimax(0.),
  m_flav1(Flavour(kf_none)), m_flav2(Flavour(kf_none))
{
}

DeltaPhi_Selector::~DeltaPhi_Selector() {
}

bool DeltaPhi_Selector::Trigger(Selector_List &sl)
{
  DEBUG_FUNC(m_on);
  if (!m_on) return true;
  for (size_t i=m_nin;i<m_n;i++) {
    for (size_t j=i+1;j<m_n;j++) {
      if (sl[i].Momentum()==Vec4D(0.,0.,0.,0.)) continue;
      if ( (m_flav1.Includes(sl[i].Flavour()) &&
            m_flav2.Includes(sl[j].Flavour())) ||
           (m_flav1.Includes(sl[j].Flavour()) &&
            m_flav2.Includes(sl[i].Flavour())) ) {
        double dphiij = dabs(sl[i].Momentum().DPhi(sl[j].Momentum()));
        if (m_sel_log->Hit( ((dphiij<m_dphimin) || (dphiij>m_dphimax)) )) return false;
      }
    }
  }
  return true;
}

void DeltaPhi_Selector::BuildCuts(Cut_Data * cuts)
{
}

void DeltaPhi_Selector::SetRange(Flavour fl1,Flavour fl2,double min,double max)
{
  m_flav1=fl1;
  m_flav2=fl2;
  m_dphimin=min;
  m_dphimax=max;

  for (size_t i=m_nin;i<m_n;i++) {
    for (size_t j=i+1;j<m_n;j++) {
      if ( (m_flav1.Includes(p_fl[i]) && m_flav2.Includes(p_fl[j])) ||
           (m_flav1.Includes(p_fl[j]) && m_flav2.Includes(p_fl[i])) ) {
        m_on=true;
      }
    }
  }
}

DECLARE_ND_GETTER(DeltaPhi_Selector,"DPhi",Selector_Base,Selector_Key,true);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,DeltaPhi_Selector>::
operator()(const Selector_Key &key) const
{
  Scoped_Settings s{ key.m_settings };
  const auto parameters = s.SetDefault<std::string>({}).GetVector<std::string>();
  assert(parameters[0] == "DPhi");
  if (parameters.size() != 5)
    THROW(critical_error, "Invalid syntax");
  const auto kf1 = s.Interprete<int>(parameters[1]);
  const auto kf2 = s.Interprete<int>(parameters[2]);
  const auto min = s.Interprete<double>(parameters[3]);
  const auto max = s.Interprete<double>(parameters[4]);
  Flavour flav1 = Flavour((kf_code)abs(kf1),kf1<0);
  Flavour flav2 = Flavour((kf_code)abs(kf2),kf2<0);
  DeltaPhi_Selector *sel = new DeltaPhi_Selector(key.p_proc);
  sel->SetRange(flav1,flav2,min,max);
  return sel;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,DeltaPhi_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"azimuthal separation selector";
}



/*--------------------------------------------------------------------

  DeltaR Selector

  --------------------------------------------------------------------*/

DeltaR_Selector::DeltaR_Selector(Process_Base *const proc):
  Selector_Base("DeltaR_Selector",proc), m_dRmin(0.), m_dRmax(0.),
  m_flav1(Flavour(kf_none)), m_flav2(Flavour(kf_none))
{
}

DeltaR_Selector::~DeltaR_Selector() {
}

bool DeltaR_Selector::Trigger(Selector_List &sl)
{
  DEBUG_FUNC(m_on);
  if (!m_on) return true;
  for (size_t i=m_nin;i<sl.size();i++) {
    for (size_t j=i+1;j<sl.size();j++) {
      if (sl[i].Momentum()==Vec4D(0.,0.,0.,0.)) continue;
      if ( (m_flav1.Includes(sl[i].Flavour()) &&
            m_flav2.Includes(sl[j].Flavour())) ||
           (m_flav1.Includes(sl[j].Flavour()) &&
            m_flav2.Includes(sl[i].Flavour())) ) {
        double dRij = sl[i].Momentum().DR(sl[j].Momentum());
        if (m_sel_log->Hit( ((dRij<m_dRmin) || (dRij>m_dRmax)) )) return false;
      }
    }
  }
  return true;
}

void DeltaR_Selector::BuildCuts(Cut_Data * cuts)
{
}

void DeltaR_Selector::SetRange(Flavour fl1,Flavour fl2,double min,double max)
{
  m_flav1=fl1;
  m_flav2=fl2;
  m_dRmin=min;
  m_dRmax=max;

  for (size_t i=m_nin;i<m_n;i++) {
    for (size_t j=i+1;j<m_n;j++) {
      if ( (m_flav1.Includes(p_fl[i]) && m_flav2.Includes(p_fl[j])) ||
           (m_flav1.Includes(p_fl[j]) && m_flav2.Includes(p_fl[i])) ) {
        m_on=true;
      }
    }
  }
}

DECLARE_ND_GETTER(DeltaR_Selector,"DR",Selector_Base,Selector_Key,true);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,DeltaR_Selector>::
operator()(const Selector_Key &key) const
{
  Scoped_Settings s{ key.m_settings };
  const auto parameters = s.SetDefault<std::string>({}).GetVector<std::string>();
  assert(parameters[0] == "DR");
  if (parameters.size() != 5)
    THROW(critical_error, "Invalid syntax");
  const auto kf1 = s.Interprete<int>(parameters[1]);
  const auto kf2 = s.Interprete<int>(parameters[2]);
  const auto min = s.Interprete<double>(parameters[3]);
  const auto max = s.Interprete<double>(parameters[4]);
  Flavour flav1 = Flavour((kf_code)abs(kf1),kf1<0);
  Flavour flav2 = Flavour((kf_code)abs(kf2),kf2<0);
  DeltaR_Selector *sel = new DeltaR_Selector(key.p_proc);
  sel->SetRange(flav1,flav2,min,max);
  return sel;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,DeltaR_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"angular separation (using phi and eta) selector";
}



/*--------------------------------------------------------------------

  DeltaRy Selector

  --------------------------------------------------------------------*/

DeltaRy_Selector::DeltaRy_Selector(Process_Base *const proc):
  Selector_Base("DeltaRy_Selector",proc), m_dRymin(0.), m_dRymax(0.),
  m_flav1(Flavour(kf_none)), m_flav2(Flavour(kf_none))
{
}

DeltaRy_Selector::~DeltaRy_Selector() {
}

bool DeltaRy_Selector::Trigger(Selector_List &sl)
{
  DEBUG_FUNC(m_on);
  if (!m_on) return true;
  for (size_t i=m_nin;i<sl.size();i++) {
    for (size_t j=i+1;j<sl.size();j++) {
      if (sl[i].Momentum()==Vec4D(0.,0.,0.,0.)) continue;
      if ( (m_flav1.Includes(sl[i].Flavour()) &&
            m_flav2.Includes(sl[j].Flavour())) ||
           (m_flav1.Includes(sl[j].Flavour()) &&
            m_flav2.Includes(sl[i].Flavour())) ) {
        double dRyij = sl[i].Momentum().DRy(sl[j].Momentum());
        if (m_sel_log->Hit( ((dRyij<m_dRymin) || (dRyij>m_dRymax)) )) return false;
      }
    }
  }
  return true;
}

void DeltaRy_Selector::BuildCuts(Cut_Data * cuts)
{
}

void DeltaRy_Selector::SetRange(Flavour fl1,Flavour fl2,double min,double max)
{
  m_flav1=fl1;
  m_flav2=fl2;
  m_dRymin=min;
  m_dRymax=max;

  for (size_t i=m_nin;i<m_n;i++) {
    for (size_t j=i+1;j<m_n;j++) {
      if ( (m_flav1.Includes(p_fl[i]) && m_flav2.Includes(p_fl[j])) ||
           (m_flav1.Includes(p_fl[j]) && m_flav2.Includes(p_fl[i])) ) {
        m_on=true;
      }
    }
  }
}

DECLARE_ND_GETTER(DeltaRy_Selector,"DR(y)",Selector_Base,Selector_Key,true);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,DeltaRy_Selector>::
operator()(const Selector_Key &key) const
{
  Scoped_Settings s{ key.m_settings };
  const auto parameters = s.SetDefault<std::string>({}).GetVector<std::string>();
  assert(parameters[0] == "DR(y)");
  if (parameters.size() != 5)
    THROW(critical_error, "Invalid syntax");
  const auto kf1 = s.Interprete<int>(parameters[1]);
  const auto kf2 = s.Interprete<int>(parameters[2]);
  const auto min = s.Interprete<double>(parameters[3]);
  const auto max = s.Interprete<double>(parameters[4]);
  Flavour flav1 = Flavour((kf_code)abs(kf1),kf1<0);
  Flavour flav2 = Flavour((kf_code)abs(kf2),kf2<0);
  DeltaRy_Selector *sel = new DeltaRy_Selector(key.p_proc);
  sel->SetRange(flav1,flav2,min,max);
  return sel;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,DeltaRy_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"angular separation (using phi and y) selector";
}



/*--------------------------------------------------------------------

  PhiStar Selector

  --------------------------------------------------------------------*/

PhiStar_Selector::PhiStar_Selector(Process_Base *const proc):
  Selector_Base("PhiStar_Selector",proc), m_phistarmin(0.), m_phistarmax(0.),
  m_flav1(Flavour(kf_none)), m_flav2(Flavour(kf_none))
{
}

PhiStar_Selector::~PhiStar_Selector() {
}

bool PhiStar_Selector::Trigger(Selector_List &sl)
{
  DEBUG_FUNC(m_on);
  if (!m_on) return true;
  for (size_t i=m_nin;i<sl.size();i++) {
    for (size_t j=i+1;j<sl.size();j++) {
      if (sl[i].Momentum()==Vec4D(0.,0.,0.,0.)) continue;
      if ( (m_flav1.Includes(sl[i].Flavour()) &&
            m_flav2.Includes(sl[j].Flavour())) ||
           (m_flav1.Includes(sl[j].Flavour()) &&
            m_flav2.Includes(sl[i].Flavour())) ) {
        // phi* = tan((pi-dphi)/2) * sin(theta*)
        //      = tan((pi-dphi)/2) * sqrt(1-cos^2(theta*))
        // with
        // cos(theta*) = tanh(deta/2)
        double phistarij = tan(0.5*(M_PI-sl[i].Momentum()
                                              .DPhi(sl[j].Momentum())))
                           *sqrt(1.-sqr(tanh(0.5*sl[i].Momentum()
                                                      .DEta(sl[j].Momentum()))));
        if (m_sel_log->Hit( ((phistarij<m_phistarmin) ||
                             (phistarij>m_phistarmax)) )) return false;
      }
    }
  }
  return true;
}

void PhiStar_Selector::BuildCuts(Cut_Data * cuts)
{
}

void PhiStar_Selector::SetRange(Flavour fl1,Flavour fl2,double min,double max)
{
  m_flav1=fl1;
  m_flav2=fl2;
  m_phistarmin=min;
  m_phistarmax=max;

  for (size_t i=m_nin;i<m_n;i++) {
    for (size_t j=i+1;j<m_n;j++) {
      if ( (m_flav1.Includes(p_fl[i]) && m_flav2.Includes(p_fl[j])) ||
           (m_flav1.Includes(p_fl[j]) && m_flav2.Includes(p_fl[i])) ) {
        m_on=true;
      }
    }
  }
}

DECLARE_ND_GETTER(PhiStar_Selector,"PhiStar",Selector_Base,Selector_Key,true);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,PhiStar_Selector>::
operator()(const Selector_Key &key) const
{
  Scoped_Settings s{ key.m_settings };
  const auto parameters = s.SetDefault<std::string>({}).GetVector<std::string>();
  assert(parameters[0] == "PhiStar");
  if (parameters.size() != 5)
    THROW(critical_error, "Invalid syntax");
  const auto kf1 = s.Interprete<int>(parameters[1]);
  const auto kf2 = s.Interprete<int>(parameters[2]);
  const auto min = s.Interprete<double>(parameters[3]);
  const auto max = s.Interprete<double>(parameters[4]);
  Flavour flav1 = Flavour((kf_code)abs(kf1),kf1<0);
  Flavour flav2 = Flavour((kf_code)abs(kf2),kf2<0);
  PhiStar_Selector *sel = new PhiStar_Selector(key.p_proc);
  sel->SetRange(flav1,flav2,min,max);
  return sel;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,PhiStar_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"\\phi^* selector";
}



/*--------------------------------------------------------------------

 Multiplicity Selector

 --------------------------------------------------------------------*/

Multiplicity_Selector::Multiplicity_Selector(Process_Base *const proc):
  Selector_Base("Multiplicity_Selector",proc), m_nmin(0), m_nmax(0),
  m_flav(Flavour(kf_none))
{
}

Multiplicity_Selector::~Multiplicity_Selector() {
}

bool Multiplicity_Selector::Trigger(Selector_List &sl)
{
  DEBUG_FUNC(m_on);
  if (!m_on) return true;
  double cnt(0);
  for (size_t i=m_nin;i<sl.size();i++) {
    if (m_flav.Includes(sl[i].Flavour())) {
      cnt++;
    }
  }
  if (m_sel_log->Hit( ((cnt<m_nmin) || (cnt>m_nmax)) ))
    return false;
  return true;
}

void Multiplicity_Selector::BuildCuts(Cut_Data * cuts)
{
}

void Multiplicity_Selector::SetRange(Flavour fl,size_t min,size_t max)
{
  m_flav=fl;
  m_nmin=min;
  m_nmax=max;
  // on only if at least one in proc
//  for (size_t i=m_nin;i<m_n;i++) {
//    if (m_flav1.Includes(p_fl[i])) {
//      m_on=true;
//    }
//  }
  // always on
  m_on=true;
}

DECLARE_ND_GETTER(Multiplicity_Selector,"N",Selector_Base,Selector_Key,true);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,Multiplicity_Selector>::
operator()(const Selector_Key &key) const
{
  Scoped_Settings s{ key.m_settings };
  const auto parameters = s.SetDefault<std::string>({}).GetVector<std::string>();
  assert(parameters[0] == "N");
  if (parameters.size() != 4)
    THROW(critical_error, "Invalid syntax");
  const auto kf = s.Interprete<int>(parameters[1]);
  const auto min = s.Interprete<double>(parameters[2]);
  const auto max = s.Interprete<double>(parameters[3]);
  Flavour flav = Flavour((kf_code)abs(kf),kf<0);
  Multiplicity_Selector *sel = new Multiplicity_Selector(key.p_proc);
  sel->SetRange(flav,min,max);
  return sel;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,Multiplicity_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"multiplicity selector";
}



/*--------------------------------------------------------------------

 PTMIS Selector

 --------------------------------------------------------------------*/

PTMIS_Selector::PTMIS_Selector(Process_Base *const proc):
 Selector_Base("PTMIS_Selector",proc), m_ptmismin(0.), m_ptmismax(0.)
{
  // maybe extend to other invisible flavours depending on model
  m_flavs=Flavour_Vector(1,Flavour(kf_neutrino));
}

PTMIS_Selector::~PTMIS_Selector() {
}

bool PTMIS_Selector::Trigger(Selector_List &sl)
{
  DEBUG_FUNC(m_on);
  if (!m_on) return true;
  Vec4D mismom(0.,0.,0.,0.);
  for (size_t k=0;k<m_flavs.size();k++) {
    for (size_t i=m_nin;i<sl.size();i++) {
      if (m_flavs[k].Includes(sl[i].Flavour())) {
        mismom+=sl[i].Momentum();
      }
    }
  }
  double ptmis(mismom.PPerp());
  if (m_sel_log->Hit( ((ptmis<m_ptmismin) || (ptmis>m_ptmismax)) ))
    return false;
  return true;
}

void PTMIS_Selector::BuildCuts(Cut_Data * cuts)
{
}

void PTMIS_Selector::SetRange(double min,double max)
{
  m_ptmismin=min;
  m_ptmismax=max;
  m_smin = Max(m_smin,m_ptmismin*m_ptmismin);
  // on only if at least one invisible in proc
//  for (size_t i=m_nin;i<m_n;i++) {
//    for (size_t j(0);j<m_flavs.size();++i) {
//      if (m_flavs[j].Includes(p_fl[i])) {
//        m_on=true;
//      }
//    }
//  }
  // always on
  m_on=true;
}

DECLARE_ND_GETTER(PTMIS_Selector,"PTmis",Selector_Base,Selector_Key,true);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,PTMIS_Selector>::
operator()(const Selector_Key &key) const
{
  Scoped_Settings s{ key.m_settings };
  const auto parameters = s.SetDefault<std::string>({}).GetVector<std::string>();
  assert(parameters[0] == "PTmis");
  if (parameters.size() != 3)
    THROW(critical_error, "Invalid syntax");
  const auto min = s.Interprete<double>(parameters[2]);
  const auto max = s.Interprete<double>(parameters[3]);
  PTMIS_Selector *sel = new PTMIS_Selector(key.p_proc);
  sel->SetRange(min,max);
  return sel;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,PTMIS_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"missing transverse momentum selector";
}



/*--------------------------------------------------------------------

 ETMIS Selector

 --------------------------------------------------------------------*/

ETMIS_Selector::ETMIS_Selector(Process_Base *const proc):
 Selector_Base("ETMIS_Selector",proc), m_etmismin(0.), m_etmismax(0.)
{
  // maybe extend to other invisible flavours depending on model
  m_flavs=Flavour_Vector(1,Flavour(kf_neutrino));
}

ETMIS_Selector::~ETMIS_Selector() {
}

bool ETMIS_Selector::Trigger(Selector_List &sl)
{
  DEBUG_FUNC(m_on);
  if (!m_on) return true;
  Vec4D mismom(0.,0.,0.,0.);
  for (size_t k=0;k<m_flavs.size();k++) {
    for (size_t i=m_nin;i<sl.size();i++) {
      if (m_flavs[k].Includes(sl[i].Flavour())) {
        mismom+=sl[i].Momentum();
      }
    }
  }
  double etmis(mismom.MPerp());
  if (m_sel_log->Hit( ((etmis<m_etmismin) || (etmis>m_etmismax)) ))
    return false;
  return true;
}

void ETMIS_Selector::BuildCuts(Cut_Data * cuts)
{
}

void ETMIS_Selector::SetRange(double min,double max)
{
  m_etmismin=min;
  m_etmismax=max;
  m_smin = Max(m_smin,m_etmismin*m_etmismin);
  // on only if at least one invisible in proc
//  for (size_t i=m_nin;i<m_n;i++) {
//    for (size_t j(0);j<m_flavs.size();++i) {
//      if (m_flavs[j].Includes(p_fl[i])) {
//        m_on=true;
//      }
//    }
//  }
  // always on
  m_on=true;
}

DECLARE_ND_GETTER(ETMIS_Selector,"ETmis",Selector_Base,Selector_Key,true);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,ETMIS_Selector>::
operator()(const Selector_Key &key) const
{
  Scoped_Settings s{ key.m_settings };
  const auto parameters = s.SetDefault<std::string>({}).GetVector<std::string>();
  assert(parameters[0] == "ETmis");
  if (parameters.size() != 3)
    THROW(critical_error, "Invalid syntax");
  const auto min = s.Interprete<double>(parameters[2]);
  const auto max = s.Interprete<double>(parameters[3]);
  ETMIS_Selector *sel = new ETMIS_Selector(key.p_proc);
  sel->SetRange(min,max);
  return sel;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,ETMIS_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"missing transverse energy selector";
}



/*--------------------------------------------------------------------

  photon isolation cut: hep-ph/9801442

  --------------------------------------------------------------------*/

Isolation_Cut::Isolation_Cut(Process_Base *const proc) :
  Selector_Base("IsolationCut",proc), m_dR(0.), m_exp(0.), m_emax(0.),
  m_massmax(0.), m_iflav(Flavour(kf_none))
{
}

void Isolation_Cut::SetRange(Flavour flav,double dR,double exp,double emax,
                                          double massmax)
{
  m_dR    = dR;
  m_exp   = exp;
  m_emax  = emax;
  m_massmax = massmax;
  m_iflav = flav;
  m_vf.clear();
  for (int i=m_nin;i<m_nin+m_nout;i++) {
    if (m_iflav.Includes(p_fl[i])) {
      m_on=true;
      m_vf.push_back(i);
    }
  }
}

class edrt {
public:
  double E;
  double dr;
  edrt(double _e,double _dr) : E(_e), dr(_dr) {}
};
class Order_edrt {
public:
  int operator()(const edrt a, const edrt b);
};
int Order_edrt::operator()(const edrt a, const edrt b) {
  if (a.dr<b.dr) return 1;
  return 0;
}

bool Isolation_Cut::Trigger(Selector_List &sl)
{
  DEBUG_FUNC(m_on);
  if (!m_on) return true;
  vector<size_t> vfsub;
  for (size_t i=m_nin;i<m_n;i++)
    if (m_iflav.Includes(sl[i].Flavour())) vfsub.push_back(i);
  const vector<size_t> *const vf(&vfsub);
  for (size_t k=0;k<vf->size();k++) {
    double egamma=sl[(*vf)[k]].Momentum().PPerp();
    vector<edrt> edrlist;
    for (size_t i=m_nin;i<m_n;i++) {
      if (Flavour(kf_jet).Includes(sl[i].Flavour()) ||
          (sl[i].Flavour().Strong() && sl[i].Flavour().Mass()<m_massmax)) {
        double dr=DR(sl[(*vf)[k]].Momentum(),sl[i].Momentum());
        if (dr<m_dR) edrlist.push_back(edrt(sl[i].Momentum().PPerp(),dr));
      }
    }
    if (edrlist.size()>0) {
      stable_sort(edrlist.begin(),edrlist.end(),Order_edrt());
      double etot=0.;
      for (size_t i=0;i<edrlist.size();i++) {
        etot+=edrlist[i].E;
        if (m_sel_log->Hit(etot>Chi(egamma,edrlist[i].dr))) return false;
      }
      edrlist.clear();
    }
  }
  return true;
}


void Isolation_Cut::BuildCuts(Cut_Data * cuts)
{
}

double Isolation_Cut::Chi(double eg,double dr)
{
  if      (m_exp<0.)  return 0.;//rpa->gen.Ecms();
  else if (m_exp==0.) return m_emax;
  else if (m_exp==1.) return eg*m_emax*(1.-cos(dr))/(1.-cos(m_dR));
  else if (m_exp==2.) return eg*m_emax*sqr((1.-cos(dr))/(1.-cos(m_dR)));
  else                return eg*m_emax*pow((1.-cos(dr))/(1.-cos(m_dR)),m_exp);
}

double Isolation_Cut::DR(const Vec4D & p1,const Vec4D & p2)
{
  return  sqrt(sqr(DEta12(p1,p2)) + sqr(DPhi12(p1,p2)));
}
double Isolation_Cut::DEta12(const Vec4D & p1,const Vec4D & p2)
{
  // eta1,2 = -log(tan(theta_1,2)/2)
  double c1=p1[3]/Vec3D(p1).Abs();
  double c2=p2[3]/Vec3D(p2).Abs();
  return  0.5 *log( (1 + c1)*(1 - c2)/((1-c1)*(1+c2)));
}

double Isolation_Cut::DPhi12(const Vec4D & p1,const Vec4D & p2)
{
  double pt1=sqrt(p1[1]*p1[1]+p1[2]*p1[2]);
  double pt2=sqrt(p2[1]*p2[1]+p2[2]*p2[2]);
  return acos((p1[1]*p2[1]+p1[2]*p2[2])/(pt1*pt2));
}

DECLARE_ND_GETTER(Isolation_Cut,"IsolationCut",Selector_Base,Selector_Key,true);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,Isolation_Cut>::
operator()(const Selector_Key &key) const
{
  Scoped_Settings s{ key.m_settings };
  const auto parameters = s.SetDefault<std::string>({}).GetVector<std::string>();
  assert(parameters[0] == "IsolationCut");
  if (parameters.size() < 5)
    THROW(critical_error, "Invalid syntax");
  const auto kf = s.Interprete<int>(parameters[1]);
  const auto dR = s.Interprete<double>(parameters[2]);
  const auto exp = s.Interprete<int>(parameters[3]);
  const auto emax = s.Interprete<double>(parameters[4]);
  auto massmax = 0.0;
  if (parameters.size() > 5)
    massmax = s.Interprete<double>(parameters[5]);
  Flavour flav = Flavour((kf_code)abs(kf),kf<0);
  Isolation_Cut *sel = new Isolation_Cut(key.p_proc);
  sel->SetRange(flav,dR,exp,emax,massmax);
  return sel;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,Isolation_Cut>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Isolation_Cut selector [hep-ph/9801442]";
}



/*--------------------------------------------------------------------

 NJettiness arXiv:1004.2489

 --------------------------------------------------------------------*/

NJettiness_Selector::NJettiness_Selector(Process_Base *const proc):
 Selector_Base("NJettiness_Selector",proc), m_N(0), m_njmin(0.), m_njmax(0.)
{
  THROW(not_implemented,"Not fully implemented yet.");
}

NJettiness_Selector::~NJettiness_Selector() {
}

bool NJettiness_Selector::Trigger(Selector_List &sl)
{
  DEBUG_FUNC(m_on);
  if (!m_on) return true;
  // find all reference vectors n_i
  Vec4D_Vector ni(m_N+2,Vec4D(0.,0.,0.,0));
  ni[0]=Vec4D(1.,0.,0.,1.);
  ni[1]=Vec4D(1.,0.,0.,-1.);
  // construct all q_i and non hadronic q
  Vec4D_Vector qi(m_N+2,Vec4D(0.,0.,0.,0));
  Vec4D q(0.,0.,0.,0.);
  // find x_a, x_b, Q2
  Vec4D sumq(q);
  for (size_t i(0);i<qi.size();++i) sumq+=qi[i];
  return true;
}

void NJettiness_Selector::BuildCuts(Cut_Data * cuts)
{
}

void NJettiness_Selector::SetRange(size_t N,double min,double max)
{
  m_N=N;
  m_njmin=min;
  m_njmax=max;
}

DECLARE_ND_GETTER(NJettiness_Selector,"NJ",Selector_Base,Selector_Key,true);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,NJettiness_Selector>::
operator()(const Selector_Key &key) const
{
  Scoped_Settings s{ key.m_settings };
  const auto parameters = s.SetDefault<std::string>({}).GetVector<std::string>();
  assert(parameters[0] == "NJ");
  if (parameters.size() != 4)
    THROW(critical_error, "Invalid syntax");
  const auto n = s.Interprete<size_t>(parameters[1]);
  const auto min = s.Interprete<double>(parameters[2]);
  const auto max = s.Interprete<double>(parameters[3]);
  NJettiness_Selector *sel = new NJettiness_Selector(key.p_proc);
  sel->SetRange(n,min,max);
  return sel;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,NJettiness_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"NJettiness selector [arXiv:1004.2489]";
}




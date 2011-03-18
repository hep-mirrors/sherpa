#ifndef ATOOLS_Phys_Standard_Selector_H
#define ATOOLS_Phys_Standard_Selector_H

#include "PHASIC++/Selectors/Selector.H"
#include "ATOOLS/Org/Data_Reader.H"

namespace PHASIC {

  class Energy_Selector : public Selector_Base {
    double * emin, * emax, * value;
    bool     m_strong;
  public:
    Energy_Selector(int,int,ATOOLS::Flavour *);
    ~Energy_Selector();
    void     SetRange(ATOOLS::Flavour_Vector,double,double);
    bool     Trigger(const ATOOLS::Vec4D_Vector & );
    bool     JetTrigger(const ATOOLS::Vec4D_Vector &,ATOOLS::NLO_subevtlist *const);
    bool     NoJetTrigger(const ATOOLS::Vec4D_Vector &);
    void     BuildCuts(Cut_Data *);
  };

  class ET_Selector : public Selector_Base {
    double * etmin, * etmax, * value;
    bool     m_strong;
  public:
    ET_Selector(int,int,ATOOLS::Flavour *);
    ~ET_Selector();
    void     SetRange(ATOOLS::Flavour_Vector,double,double);
    bool     Trigger(const ATOOLS::Vec4D_Vector & );
    bool     JetTrigger(const ATOOLS::Vec4D_Vector &,ATOOLS::NLO_subevtlist *const);
    bool     NoJetTrigger(const ATOOLS::Vec4D_Vector &);
    void     BuildCuts(Cut_Data *);
  };

  class PT_Selector : public Selector_Base {
    double * ptmin, * ptmax, * value;
    bool     m_strong;
  public:
    PT_Selector(int,int,ATOOLS::Flavour *);
    ~PT_Selector();
    void     SetRange(ATOOLS::Flavour_Vector,double,double);
    bool     Trigger(const ATOOLS::Vec4D_Vector & );
    bool     JetTrigger(const ATOOLS::Vec4D_Vector &,ATOOLS::NLO_subevtlist *const);
    bool     NoJetTrigger(const ATOOLS::Vec4D_Vector &);
    void     BuildCuts(Cut_Data *);
  };

  class Rapidity_Selector : public Selector_Base {
    double * ymin, * ymax, * value;
    bool     m_strong;
  public:
    Rapidity_Selector(int,int,ATOOLS::Flavour *);
    ~Rapidity_Selector();
    void     SetRange(ATOOLS::Flavour_Vector,double,double);
    bool     Trigger(const ATOOLS::Vec4D_Vector & );
    bool     JetTrigger(const ATOOLS::Vec4D_Vector &,ATOOLS::NLO_subevtlist *const);
    bool     NoJetTrigger(const ATOOLS::Vec4D_Vector &);
    void     BuildCuts(Cut_Data *);
    int      IsConditional() { return 1; }
  };

  class PseudoRapidity_Selector : public Selector_Base {
    double * etamin, * etamax, * value;
    bool     m_strong;
  public:
    PseudoRapidity_Selector(int,int,ATOOLS::Flavour *);
    ~PseudoRapidity_Selector();
    void     SetRange(ATOOLS::Flavour_Vector,double,double);
    bool     Trigger(const ATOOLS::Vec4D_Vector & );
    bool     JetTrigger(const ATOOLS::Vec4D_Vector &,ATOOLS::NLO_subevtlist *const);
    bool     NoJetTrigger(const ATOOLS::Vec4D_Vector &);
    void     BuildCuts(Cut_Data *);
  };

  class Angle_Selector : public Selector_Base {
    double ** cosmin, ** cosmax, * value;
    bool     m_strong;
  public:
    Angle_Selector(int,int,ATOOLS::Flavour *);
    ~Angle_Selector();
    void     SetRange(ATOOLS::Flavour_Vector,double,double);
    void     SetRange(ATOOLS::Flavour_Vector,int,double,double);
    bool     Trigger(const ATOOLS::Vec4D_Vector & );
    bool     JetTrigger(const ATOOLS::Vec4D_Vector &,ATOOLS::NLO_subevtlist *const);
    bool     NoJetTrigger(const ATOOLS::Vec4D_Vector &);
    void     BuildCuts(Cut_Data *);
  };

  class PT2_Selector : public Selector_Base {
    double ** pt2min, ** pt2max, * value;
    bool     m_strong;
  public:
    PT2_Selector(int,int,ATOOLS::Flavour *);
    ~PT2_Selector();
    void     SetRange(ATOOLS::Flavour_Vector,double,double);
    bool     Trigger(const ATOOLS::Vec4D_Vector & );
    bool     JetTrigger(const ATOOLS::Vec4D_Vector &,ATOOLS::NLO_subevtlist *const);
    bool     NoJetTrigger(const ATOOLS::Vec4D_Vector &);
    void     BuildCuts(Cut_Data *);
  };

  class IMass_Selector : public Selector_Base {
    double ** massmin, ** massmax, * value;
    bool     m_strong;
  public:
    IMass_Selector(int,int,ATOOLS::Flavour *);
    ~IMass_Selector();
    void     SetRange(ATOOLS::Flavour_Vector,double,double);
    bool     Trigger(const ATOOLS::Vec4D_Vector & );
    bool     JetTrigger(const ATOOLS::Vec4D_Vector &,ATOOLS::NLO_subevtlist *const);
    bool     NoJetTrigger(const ATOOLS::Vec4D_Vector &);
    void     BuildCuts(Cut_Data *);
  };

  class IQ2_Selector : public Selector_Base {
    double ** massmin, ** massmax, * value;
    bool     m_strong;
  public:
    IQ2_Selector(int,int,ATOOLS::Flavour *);
    ~IQ2_Selector();
    void     SetRange(ATOOLS::Flavour_Vector,double,double);
    bool     Trigger(const ATOOLS::Vec4D_Vector & );
    bool     JetTrigger(const ATOOLS::Vec4D_Vector &,ATOOLS::NLO_subevtlist *const);
    bool     NoJetTrigger(const ATOOLS::Vec4D_Vector &);
    void     BuildCuts(Cut_Data *);
  };

  class Delta_Eta_Selector : public Selector_Base {
    double ** detamin, ** detamax, * value;
  public:
    Delta_Eta_Selector(int,int,ATOOLS::Flavour *);
    ~Delta_Eta_Selector();
    void     SetRange(ATOOLS::Flavour_Vector,double,double);
    bool     Trigger(const ATOOLS::Vec4D_Vector & );
    void     BuildCuts(Cut_Data *);
  };

  class Delta_Phi_Selector : public Selector_Base {
    double ** dphimin, ** dphimax, * value;
  public:
    Delta_Phi_Selector(int,int,ATOOLS::Flavour *);
    ~Delta_Phi_Selector();
    void     SetRange(ATOOLS::Flavour_Vector,double,double);
    bool     Trigger(const ATOOLS::Vec4D_Vector & );
    void     BuildCuts(Cut_Data *);
  };

  class Delta_R_Selector : public Selector_Base {
    double ** drmin, ** drmax, * value;
  public:
    Delta_R_Selector(int,int,ATOOLS::Flavour *);
    ~Delta_R_Selector();
    void     SetRange(ATOOLS::Flavour_Vector,double,double);
    bool     Trigger(const ATOOLS::Vec4D_Vector & );
    void     BuildCuts(Cut_Data *);
  };

}

#endif

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

  Energy Selector

  --------------------------------------------------------------------*/

Energy_Selector::Energy_Selector(int _nin,int _nout, Flavour * _fl):
  Selector_Base("Energy_Selector")
{
  m_nin  = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl   = _fl;
  m_smin = 0.;
  m_smax = rpa.gen.Ecms()*rpa.gen.Ecms();
  m_strong = 0;
  if (m_nin==2) if (m_fl[0].Strong()&&m_fl[1].Strong()) m_strong = 1;
  
  double Emax = rpa.gen.PBeam(0)[0]+rpa.gen.PBeam(1)[0];
  emin  = new double[m_n];
  emax  = new double[m_n];
  value  = new double[m_n];
  for (int i=0;i<m_n;i++) { emin[i] = 0.; emax[i] = Emax; }
  m_sel_log = new Selector_Log(m_name);
}

Energy_Selector::~Energy_Selector() 
{
  delete [] emin;
  delete [] emax;
  delete [] value;
}


bool Energy_Selector::Trigger(const Vec4D_Vector & mom) 
{
  double ei;
  for (int i=m_nin;i<m_n;i++) {
    ei = value[i] = mom[i][0];
    if (m_sel_log->Hit( ((ei<emin[i]) || (ei>emax[i])) )) return 0;
  }
  return 1;
}

bool Energy_Selector::JetTrigger(const Vec4D_Vector &,NLO_subevtlist *const subs)
{
  if (!m_strong) return 1;
  msg_Error()<<"Energy_Selector::JetTrigger: IR unsave cut"<<std::endl;
  return 0;
}

bool Energy_Selector::NoJetTrigger(const Vec4D_Vector & mom)
{
  if (!m_strong) return Trigger(mom);
  msg_Error()<<"Energy_Selector::NoJetTrigger: IR unsave cut"<<std::endl;
  return 0;
}


void Energy_Selector::BuildCuts(Cut_Data * cuts) 
{
  for (int i=0;i<m_n;i++) {
    cuts->energymin[i] = Max(emin[i],cuts->energymin[i]);
  }
}

void Energy_Selector::SetRange(std::vector<Flavour> crit,double _min,
			       double _max)
{
  if (crit.size() != 1) {
    msg_Error()<<"Wrong number of arguments in Energy_Selector::SetRange : "
			  <<crit.size()<<endl;
    return;
  }


  double MaxEmin = 0.;
  for (int i=m_nin;i<m_n;i++) {
    if (crit[0].Includes(m_fl[i])) {
      emin[i] = Max(_min,m_fl[i].SelMass()); 
      emax[i] = Min(_max,(rpa.gen.PBeam(0)[0]+rpa.gen.PBeam(1)[0]));
      if (emin[i]>MaxEmin ) MaxEmin = emin[i];
      if (m_fl[i].Strong()) m_strong = 1;
    }
  }
  m_smin = Max(MaxEmin*MaxEmin,m_smin);
}

namespace PHASIC{

DECLARE_ND_GETTER(Energy_Selector_Getter,"Energy",Selector_Base,Selector_Key,true);

Selector_Base *Energy_Selector_Getter::operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<3) THROW(critical_error,"Invalid syntax");
  int crit1=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][0]));
  double min=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][1]));
  double max=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][2]));
  Flavour flav = Flavour((kf_code)abs(crit1));
  if (crit1<0) flav = flav.Bar();
  Flavour_Vector critflavs(1,flav);
  Energy_Selector *sel = new Energy_Selector
    (key.p_proc->NIn(),key.p_proc->NOut(),
     (Flavour*)&key.p_proc->Process()->Flavours().front());
  sel->SetRange(critflavs,min,max);
  return sel;
}

void Energy_Selector_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"energy selector"; 
}

}

/*--------------------------------------------------------------------

  Transverse Energy Selector

  --------------------------------------------------------------------*/

ET_Selector::ET_Selector(int _nin,int _nout, Flavour * _fl):
  Selector_Base("ET_Selector")
{
  m_nin  = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl   = _fl;
  m_smin = 0.;
  m_smax = rpa.gen.Ecms()*rpa.gen.Ecms();
  m_strong = 0;
  if (m_nin==2) if (m_fl[0].Strong()&&m_fl[1].Strong()) m_strong = 1;
  
  double Emax = rpa.gen.PBeam(0)[0]+rpa.gen.PBeam(1)[0];
  etmin  = new double[m_n];
  etmax  = new double[m_n];
  value  = new double[m_n];
  for (int i=0;i<m_n;i++) { etmin[i] = 0.; etmax[i] = Emax; }
  m_sel_log = new Selector_Log(m_name);
}

ET_Selector::~ET_Selector() 
{
  delete [] etmin;
  delete [] etmax;
  delete [] value;
}

bool ET_Selector::Trigger(const Vec4D_Vector & mom) 
{
  double eti;
  for (int i=m_nin;i<m_n;i++) {   
    eti = value[i] = mom[i][0]*mom[i].PPerp()/mom[i].P();
    if (m_sel_log->Hit( ((eti<etmin[i]) || (eti>etmax[i])) )) return 0;
  }
  return 1;
}

bool ET_Selector::JetTrigger(const Vec4D_Vector &,NLO_subevtlist *const subs)
{
  if (!m_strong) return 1;
  msg_Error()<<"ET_Selector::JetTrigger: IR unsave cut"<<std::endl;
  return 0;
}

bool ET_Selector::NoJetTrigger(const Vec4D_Vector & mom)
{
  if (!m_strong) return Trigger(mom);
  msg_Error()<<"ET_Selector::NoJetTrigger: IR unsave cut"<<std::endl;
  return 0;
}


void ET_Selector::BuildCuts(Cut_Data * cuts) 
{
  for (int i=m_nin;i<m_n;i++) {
    cuts->energymin[i] = Max(etmin[i],cuts->energymin[i]);
    cuts->cosmax[0][i] = cuts->cosmax[1][i] = cuts->cosmax[i][0] = cuts->cosmax[i][1] =  
      Min(cuts->cosmax[0][i],sqrt(1.-4.*sqr(etmin[i])/m_smax));
    cuts->etmin[i] = Max(etmin[i],cuts->etmin[i]);
  }
}

void ET_Selector::SetRange(std::vector<Flavour> crit,double _min,double _max)
{
  if (crit.size() != 1) {
    msg_Error()<<"Wrong number of arguments in ET_Selector::SetRange : "
			  <<crit.size()<<endl;
    return;
  }

  double MaxEtmin = 0.;
  for (int i=m_nin;i<m_n;i++) {
    if (crit[0].Includes(m_fl[i])) {
      etmin[i] = _min; 
      etmax[i] = Min(_max,(rpa.gen.PBeam(0)[0]+rpa.gen.PBeam(1)[0]));
      if (etmin[i] > MaxEtmin) MaxEtmin = etmin[i];
      if (m_fl[i].Strong()) m_strong = 1;
    }
  }
  m_smin = Max(MaxEtmin*MaxEtmin,m_smin);
}

namespace PHASIC {

DECLARE_ND_GETTER(ET_Selector_Getter,"ET",Selector_Base,Selector_Key,true);

Selector_Base *ET_Selector_Getter::operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<3) THROW(critical_error,"Invalid syntax");
  int crit1=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][0]));
  double min=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][1]));
  double max=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][2]));
  Flavour flav = Flavour((kf_code)abs(crit1));
  if (crit1<0) flav = flav.Bar();
  Flavour_Vector critflavs(1,flav);
  ET_Selector *sel = new ET_Selector
    (key.p_proc->NIn(),key.p_proc->NOut(),
     (Flavour*)&key.p_proc->Process()->Flavours().front());
  sel->SetRange(critflavs,min,max);
  return sel;
}

void ET_Selector_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"transverse energy selector"; 
}

}

/*--------------------------------------------------------------------

  PT Selector

  --------------------------------------------------------------------*/

PT_Selector::PT_Selector(int _nin,int _nout, Flavour * _fl) :
  Selector_Base("PT_Selector")
{
  m_nin  = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl   = _fl;
  m_smin = 0.;
  m_smax = rpa.gen.Ecms()*rpa.gen.Ecms();
  m_strong = 0;
  if (m_nin==2) if (m_fl[0].Strong()&&m_fl[1].Strong()) m_strong = 1;
  
  double Emax = rpa.gen.PBeam(0)[0]+rpa.gen.PBeam(1)[0];
  ptmin  = new double[m_n];
  ptmax  = new double[m_n];
  value = new double[m_n];
  for (int i=0;i<m_n;i++) { ptmin[i] = 0.; ptmax[i] = Emax; }
  m_sel_log = new Selector_Log(m_name);
}

PT_Selector::~PT_Selector() {
  delete [] ptmin;
  delete [] ptmax;
  delete [] value;
}


bool PT_Selector::Trigger(const Vec4D_Vector & mom) 
{
   double sumM2=0.;
  for (int i=m_nin;i<m_n;i++) {
    sumM2+=sqr(m_fl[i].SelMass());
  }
  double pti;
  for (int i=m_nin;i<m_n;i++) {
    pti = value[i] = sqrt(sqr(mom[i][1]) + sqr(mom[i][2]));
    if (m_sel_log->Hit( ((pti<ptmin[i]) || (pti>ptmax[i])) )) return 0;
  }
  return 1;
}

bool PT_Selector::JetTrigger(const Vec4D_Vector &,NLO_subevtlist *const subs)
{
  if (!m_strong) return 1;
  msg_Error()<<"PT_Selector::JetTrigger: IR unsave cut"<<std::endl;
  return 0;
}

bool PT_Selector::NoJetTrigger(const Vec4D_Vector &mom)
{
  if (!m_strong) return Trigger(mom);
  msg_Error()<<"PT_Selector::NoJetTrigger: IR unsave cut"<<std::endl;
  return 0;
}


void PT_Selector::BuildCuts(Cut_Data * cuts) 
{
  double sumM2=0.;
  for (int i=m_nin;i<m_n;i++) {
    sumM2+=sqr(m_fl[i].SelMass());
  }
  for (int i=m_nin;i<m_n;i++) {
    cuts->energymin[i] = Max(sqrt(sqr(ptmin[i])+sqr(m_fl[i].SelMass())),cuts->energymin[i]);
    double Emax2 = sqr((m_smax+2.*sqr(m_fl[i].SelMass())-sumM2)/(2.*sqrt(m_smax)));
    cuts->cosmax[0][i] = cuts->cosmax[1][i] = cuts->cosmax[i][0] = cuts->cosmax[i][1] =  
      Min(cuts->cosmax[0][i],sqrt(1.-sqr(ptmin[i])/(Emax2-sqr(m_fl[i].SelMass()))));
    cuts->etmin[i] = Max(sqrt(sqr(ptmin[i])+sqr(m_fl[i].SelMass())*(1.-sqr(cuts->cosmax[0][i]))),cuts->etmin[i]);
  }
}

void PT_Selector::SetRange(std::vector<Flavour> crit,double _min,double _max)
{
  if (crit.size() != 1) {
    msg_Error()<<"Wrong number of arguments in PT_Selector::SetRange : "
	       <<crit.size()<<endl;
    return;
  }

  double MaxPTmin = 0.;
  for (int i=m_nin;i<m_n;i++) {
    if (crit[0].Includes(m_fl[i])) {
      ptmin[i] = _min; 
      ptmax[i] = Min(_max,(rpa.gen.PBeam(0)[0]+rpa.gen.PBeam(1)[0]));
      if (ptmin[i]>MaxPTmin) MaxPTmin = ptmin[i];
      if (m_fl[i].Strong()) m_strong = 1;
    }
  }
  m_smin = Max(m_smin,4.*MaxPTmin*MaxPTmin);
}

namespace PHASIC{

DECLARE_ND_GETTER(PT_Selector_Getter,"PT",Selector_Base,Selector_Key,true);

Selector_Base *PT_Selector_Getter::operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<3) THROW(critical_error,"Invalid syntax");
  int crit1=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][0]));
  double min=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][1]));
  double max=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][2]));
  Flavour flav = Flavour((kf_code)abs(crit1));
  if (crit1<0) flav = flav.Bar();
  Flavour_Vector critflavs(1,flav);
  PT_Selector *sel = new PT_Selector
    (key.p_proc->NIn(),key.p_proc->NOut(),
     (Flavour*)&key.p_proc->Process()->Flavours().front());
  sel->SetRange(critflavs,min,max);
  return sel;
}

void PT_Selector_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"transverse momentum selector"; 
}

}

/*--------------------------------------------------------------------

  Rapidity Selector

  --------------------------------------------------------------------*/


Rapidity_Selector::Rapidity_Selector(int _nin,int _nout, Flavour * _fl) :
  Selector_Base("Rapidity_Selector")
{
  m_nin  = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl   = _fl;
  m_smin = 0.;
  m_smax = rpa.gen.Ecms()*rpa.gen.Ecms();
  m_strong = 0;
  if (m_nin==2) if (m_fl[0].Strong()&&m_fl[1].Strong()) m_strong = 1;
  
  double Emax = rpa.gen.PBeam(0)[0]+rpa.gen.PBeam(1)[0];
  double pl;

  ymin  = new double[m_n];
  ymax  = new double[m_n];
  value = new double[m_n];
  for (int i=0;i<m_n;i++) {
    pl      = sqrt(Emax*Emax-sqr(_fl[i].SelMass()));
    ymax[i] = log( (Emax+pl)/(Emax-pl) );
    ymin[i] = -ymax[i];
    if (_fl[i].SelMass()==0.) {
      ymax[i] = 100.; ymin[i] = -ymax[i];  
    }
  }
  m_sel_log = new Selector_Log(m_name);
}

Rapidity_Selector::~Rapidity_Selector() {
  delete [] ymin;
  delete [] ymax;
  delete [] value;
}


bool Rapidity_Selector::Trigger(const Vec4D_Vector & mom) 
{
  double yi;
  for (int i=m_nin;i<m_n;i++) {
    yi = value[i] = 0.5 * log( (mom[i][0]+mom[i][3])/(mom[i][0]-mom[i][3]) );
    if (m_sel_log->Hit( ((yi<ymin[i]) || (yi>ymax[i])) )) return 0;
  }
  return 1;
}

bool Rapidity_Selector::JetTrigger(const Vec4D_Vector &,NLO_subevtlist *const subs)
{
  if (!m_strong) return 1;
  msg_Error()<<"Rapidity_Selector::JetTrigger: IR unsave cut"<<std::endl;
  return 0;
}

bool Rapidity_Selector::NoJetTrigger(const Vec4D_Vector &mom)
{
  if (!m_strong) return Trigger(mom);
  msg_Error()<<"Rapidity_Selector::NoJetTrigger: IR unsave cut"<<std::endl;
  return 0;
}


void Rapidity_Selector::BuildCuts(Cut_Data * cuts) 
{
  for (int i=m_nin;i<m_n;i++) {
    cuts->cosmax[0][i] = cuts->cosmax[i][0] =  
      Min(cuts->cosmax[0][i],1./sqrt(1.-sqr(m_fl[i].SelMass())/sqr(cuts->energymin[i]))*tanh(ymax[i]));
    cuts->cosmax[1][i] = cuts->cosmax[i][1] = 
      Min(cuts->cosmax[0][i],1./sqrt(1.-sqr(m_fl[i].SelMass())/sqr(cuts->energymin[i]))*tanh(-ymin[i]));
  }
}


void Rapidity_Selector::SetRange(std::vector<Flavour> crit,double _min, 
				 double _max)
{
  if (crit.size() != 1) {
    msg_Error()<<"Wrong number of arguments in Rapidity_Selector::SetRange : "
			  <<crit.size()<<endl;
    return;
  }
  
  double E = rpa.gen.PBeam(0)[0]+rpa.gen.PBeam(1)[0];
  double pl,y;
  
  for (int i=m_nin;i<m_n;i++) {
    if (crit[0].Includes(m_fl[i])) {
      pl      = sqrt(E*E-sqr(m_fl[i].SelMass())); 
      y       = log((E+pl)/(E-pl));
      std::cout<<" determined y yields : "<<y<<std::endl;
      ymin[i] = Max(_min,-y);
      ymax[i] = Min(_max,y);
      if (m_fl[i].Strong()) m_strong = 1;
    }
  }
}

namespace PHASIC{

DECLARE_ND_GETTER(Rapidity_Selector_Getter,"Rapidity",Selector_Base,Selector_Key,true);

Selector_Base *Rapidity_Selector_Getter::operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<3) THROW(critical_error,"Invalid syntax");
  int crit1=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][0]));
  double min=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][1]));
  double max=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][2]));
  Flavour flav = Flavour((kf_code)abs(crit1));
  if (crit1<0) flav = flav.Bar();
  Flavour_Vector critflavs(1,flav);
  Rapidity_Selector *sel = new Rapidity_Selector
    (key.p_proc->NIn(),key.p_proc->NOut(),
     (Flavour*)&key.p_proc->Process()->Flavours().front());
  sel->SetRange(critflavs,min,max);
  return sel;
}

void Rapidity_Selector_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"rapidity selector"; 
}

}

/*--------------------------------------------------------------------

  PseudoRapidity Selector

  --------------------------------------------------------------------*/


PseudoRapidity_Selector::PseudoRapidity_Selector(int _nin,int _nout, Flavour * _fl) :
  Selector_Base("PseudoRapidity_Selector")
{
  m_nin  = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl   = _fl;
  m_smin = 0.;
  m_smax = rpa.gen.Ecms()*rpa.gen.Ecms();
  m_strong = 0;
  if (m_nin==2) if (m_fl[0].Strong()&&m_fl[1].Strong()) m_strong = 1;
 
  etamin  = new double[m_n];
  etamax  = new double[m_n];
  value = new double[m_n];
  for (int i=0;i<m_n;i++) {
      etamax[i] = 100.;
      etamin[i] = -etamax[i];
  }
  m_sel_log = new Selector_Log(m_name);
}

PseudoRapidity_Selector::~PseudoRapidity_Selector() {
  delete [] etamin;
  delete [] etamax;
  delete [] value;
}


bool PseudoRapidity_Selector::Trigger(const Vec4D_Vector & mom) 
{
  double etai,theta;
  
  for (int i=m_nin;i<m_n;i++) {
    theta = acos(Vec3D(mom[i])*Vec3D(mom[0])/(Vec3D(mom[i]).Abs()*Vec3D(mom[0]).Abs())); 
    etai  = value[i] = -log(tan(theta/2.));
    if (m_sel_log->Hit( ((etai<etamin[i]) || (etai>etamax[i])) )) return 0;
  }
  return 1;
}

bool PseudoRapidity_Selector::JetTrigger(const Vec4D_Vector &,NLO_subevtlist *const subs)
{
  if (!m_strong) return 1;
  msg_Error()<<"PseudoRapidity_Selector::JetTrigger: IR unsave cut"<<std::endl;
  return 0;
}

bool PseudoRapidity_Selector::NoJetTrigger(const Vec4D_Vector &mom)
{
  if (!m_strong) return Trigger(mom);
  msg_Error()<<"PseudoRapidity_Selector::NoJetTrigger: IR unsave cut"<<std::endl;
  return 0;
}


void PseudoRapidity_Selector::BuildCuts(Cut_Data * cuts) 
{
  for (int i=m_nin;i<m_n;i++) {
    cuts->cosmin[1][i] = cuts->cosmin[i][1] = Max(cuts->cosmin[1][i],tanh(-etamax[i]));
    cuts->cosmin[0][i] = cuts->cosmin[i][0] = Max(cuts->cosmin[0][i],tanh(etamin[i])); 
    cuts->cosmax[0][i] = cuts->cosmax[i][0] = Min(cuts->cosmax[0][i],tanh(etamax[i]));
    cuts->cosmax[1][i] = cuts->cosmax[i][1] = Min(cuts->cosmax[1][i],tanh(-etamin[i]));
  }
}


void PseudoRapidity_Selector::SetRange(std::vector<Flavour> crit,double _min, 
			       double _max)
{
  if (crit.size() != 1) {
    msg_Error()<<"Wrong number of arguments in PseudoRapidity_Selector::SetRange : "
			  <<crit.size()<<endl;
    return;
  }
  
  for (int i=m_nin;i<m_n;i++) {
    if (crit[0].Includes(m_fl[i])) {
      etamin[i] = _min;
      etamax[i] = _max;
      if (m_fl[i].Strong()) m_strong = 1;
    }
  }
}

namespace PHASIC{

DECLARE_ND_GETTER(PseudoRapidity_Selector_Getter,"PseudoRapidity",Selector_Base,Selector_Key,true);

Selector_Base *PseudoRapidity_Selector_Getter::operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<3) THROW(critical_error,"Invalid syntax");
  int crit1=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][0]));
  double min=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][1]));
  double max=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][2]));
  Flavour flav = Flavour((kf_code)abs(crit1));
  if (crit1<0) flav = flav.Bar();
  Flavour_Vector critflavs(1,flav);
  PseudoRapidity_Selector *sel = new PseudoRapidity_Selector
    (key.p_proc->NIn(),key.p_proc->NOut(),
     (Flavour*)&key.p_proc->Process()->Flavours().front());
  sel->SetRange(critflavs,min,max);
  return sel;
}

void PseudoRapidity_Selector_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"pseudorapidity selector"; 
}

}

/*--------------------------------------------------------------------

  Angle Selector

  --------------------------------------------------------------------*/

Angle_Selector::Angle_Selector(int _nin,int _nout, Flavour * _fl):
  Selector_Base("Angle_Selector")
{
  m_nin  = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl   = _fl;
  m_smin = 0.;
  m_smax = rpa.gen.Ecms()*rpa.gen.Ecms();
  m_strong = 0;
  if (m_nin==2) if (m_fl[0].Strong()&&m_fl[1].Strong()) m_strong = 1;

  cosmin = new double*[m_n];
  cosmax = new double*[m_n];
  value  = new double[m_n*m_n];
  for (int i=0;i<m_n;i++) { cosmin[i] = new double[m_n]; cosmax[i] = new double[m_n]; }
  for (int i=0;i<m_n;i++) {
    for (int j=i+1;j<m_n;j++) {
      //for numerical reason min and max should be slightly larger than -/+1
      cosmin[i][j] = cosmin[j][i] = -1.1; 
      cosmax[i][j] = cosmax[j][i] =  1.1; 
    }
  }
    
  m_sel_log = new Selector_Log(m_name);
}

Angle_Selector::~Angle_Selector() {
  for (int i=0;i<m_n;i++) {
    delete [] cosmin[i];
    delete [] cosmax[i];
  }
  delete [] cosmin;
  delete [] cosmax;
  delete [] value;
}

bool Angle_Selector::Trigger(const Vec4D_Vector & mom) 
{
  double cosij;
  for (int i=0;i<m_n;i++) {
    for (int j=i+1;j<m_n;j++) {
      cosij = value[m_n*i+j] = 
	Vec3D(mom[i])*Vec3D(mom[j])/(Vec3D(mom[i]).Abs()*Vec3D(mom[j]).Abs());
      if (m_sel_log->Hit( ((cosij < cosmin[i][j]) || 
			 (cosij > cosmax[i][j])) )) return 0;
    }
  }
  return 1;
}


bool Angle_Selector::JetTrigger(const Vec4D_Vector &,NLO_subevtlist *const subs)
{
  if (!m_strong) return 1;
  msg_Error()<<"Angle_Selector::JetTrigger: IR unsave cut"<<std::endl;
  return 0;
}

bool Angle_Selector::NoJetTrigger(const Vec4D_Vector &mom)
{
  if (!m_strong) return Trigger(mom);
  msg_Error()<<"Angle_Selector::NoJetTrigger: IR unsave cut"<<std::endl;
  return 0;
}

void Angle_Selector::BuildCuts(Cut_Data * cuts) 
{
  for (int i=0;i<m_n-1;i++) {
    for (int j=i+1;j<m_n;j++) {
//       cuts->cosmin[i][j] = cuts->cosmin[j][i] = 
// 	Max(cosmin[i][j],cuts->cosmin[i][j]);
      cuts->cosmax[i][j] = cuts->cosmax[j][i] = 
	Min(cosmax[i][j],cuts->cosmax[i][j]);
    }
    if (i==1) {
      for (int j=i+1;j<m_n;j++) {
	cuts->cosmin[i][j] = cuts->cosmin[j][i] = 
	  Max(cuts->cosmin[i][j],-cuts->cosmax[0][j]);
	cuts->cosmax[i][j] = cuts->cosmax[j][i] = 
	  Min(cuts->cosmax[i][j],-cuts->cosmin[0][j]);
      }
    }
  }
}

void Angle_Selector::SetRange(std::vector<Flavour> crit,
			      double _min, double _max)
{
  if (crit.size() != 2) {
    msg_Error()<<"Wrong number of arguments in Angle_Selector::SetRange : "
			  <<crit.size()<<endl;
    return;
  }

  //for numerical reasons exact +1 or -1 may be prolematic as borders
  if (IsEqual(_min,-1.)) _min = -1.1;
  if (IsEqual(_max,1.))  _max = 1.1;


  for (int i=m_nin;i<m_n;i++) {
    for (int j=i+1;j<m_n;j++) {
      if ( ((crit[0].Includes(m_fl[i])) && (crit[1].Includes(m_fl[j])) ) || 
	   ((crit[0].Includes(m_fl[j])) && (crit[1].Includes(m_fl[i])) ) ) {
	cosmin[i][j] = cosmin[j][i] = _min; 
	cosmax[i][j] = cosmax[j][i] = _max; 
      if (m_fl[i].Strong()||m_fl[j].Strong()) m_strong = 1;
      }
    }
  }
}

void Angle_Selector::SetRange(std::vector<Flavour> crit,int beam,
			      double _min, double _max)
{
  if (crit.size() != 1) {
    msg_Error()<<"Wrong number of arguments in Angle_Selector::SetRange : "
			  <<crit.size()<<endl;
    return;
  }
  for (int i=m_nin;i<m_n;i++) {
    if ( (crit[0].Includes(m_fl[i])) || ((crit[0].Bar()).Includes(m_fl[i]) ) ) {
      cosmin[i][beam] = cosmin[beam][i] = Max(_min,-1.1); 
      cosmax[i][beam] = cosmax[beam][i] = Min(_max, 1.1); 
      if (m_fl[i].Strong()) m_strong = 1;
    }
  }
}

namespace PHASIC{

DECLARE_ND_GETTER(Angle_Selector_Getter,"Angle",Selector_Base,Selector_Key,true);

Selector_Base *Angle_Selector_Getter::operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<4) THROW(critical_error,"Invalid syntax");
  int crit1=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][0]));
  int crit2=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][1]));
  double min=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][2]));
  double max=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][3]));
  Flavour flav = Flavour((kf_code)abs(crit1));
  if (crit1<0) flav = flav.Bar();
  Flavour_Vector critflavs(1,flav);
  flav = Flavour((kf_code)abs(crit2));
  if (crit2<0) flav = flav.Bar();
  critflavs.push_back(flav);
  Angle_Selector *sel = new Angle_Selector
    (key.p_proc->NIn(),key.p_proc->NOut(),
     (Flavour*)&key.p_proc->Process()->Flavours().front());
  sel->SetRange(critflavs,min,max);
  return sel;
}

void Angle_Selector_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"angle selector"; 
}

DECLARE_ND_GETTER(BeamAngle_Selector_Getter,"BeamAngle",Selector_Base,Selector_Key,true);

Selector_Base *BeamAngle_Selector_Getter::operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<4) THROW(critical_error,"Invalid syntax");
  int crit1=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][0]));
  int crit2=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][1]));
  double min=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][2]));
  double max=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][3]));
  Flavour flav = Flavour((kf_code)abs(crit1));
  if (crit1<0) flav = flav.Bar();
  Flavour_Vector critflavs(1,flav);
  flav = Flavour((kf_code)abs(crit2));
  if (crit2<0) flav = flav.Bar();
  Angle_Selector *sel = new Angle_Selector
    (key.p_proc->NIn(),key.p_proc->NOut(),
     (Flavour*)&key.p_proc->Process()->Flavours().front());
  sel->SetRange(critflavs,flav,min,max);
  return sel;
}

void BeamAngle_Selector_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"beam angle selector"; 
}

}
 /*--------------------------------------------------------------------
 
  PT2 Selector

  --------------------------------------------------------------------*/

PT2_Selector::PT2_Selector(int _nin,int _nout, Flavour * _fl):
  Selector_Base("PT2_Selector")
{
  m_nin  = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl   = _fl;
  m_smin = 0.;
  m_smax = rpa.gen.Ecms()*rpa.gen.Ecms();
  m_strong = 0;
  if (m_nin==2) if (m_fl[0].Strong()&&m_fl[1].Strong()) m_strong = 1;

  pt2min = new double*[m_n];
  pt2max = new double*[m_n];
  value  = new double[m_n*m_n];
  for (int i=0;i<m_n;i++) { pt2min[i] = new double[m_n]; pt2max[i] = new double[m_n]; }
  for (int i=0;i<m_n;i++) {
    for (int j=i+1;j<m_n;j++) {
      //for numerical reason min and max should be slightly larger than -/+1
      pt2min[i][j] = pt2min[j][i] = 0.; 
      pt2max[i][j] = pt2max[j][i] = 2.*(rpa.gen.PBeam(0)[0]+rpa.gen.PBeam(1)[0]);
    }
  }
    
  m_sel_log = new Selector_Log(m_name);
}

PT2_Selector::~PT2_Selector() {
  for (int i=0;i<m_n;i++) {
    delete [] pt2min[i];
    delete [] pt2max[i];
  }
  delete [] pt2min;
  delete [] pt2max;
  delete [] value;
}

bool PT2_Selector::Trigger(const Vec4D_Vector &mom) 
{
  double ptij;
  for (int i=0;i<m_n;i++) {
    for (int j=i+1;j<m_n;j++) {
      ptij = value[m_n*i+j] = (mom[i]+mom[j]).PPerp();
      if (m_sel_log->Hit( ((ptij < pt2min[i][j]) || 
			   (ptij > pt2max[i][j])) )) return 0;
    }
  }
  return 1;
}


bool PT2_Selector::JetTrigger(const Vec4D_Vector &,NLO_subevtlist *const subs)
{
  if (!m_strong) return 1;
  msg_Error()<<"PT2_Selector::JetTrigger: IR unsave cut"<<std::endl;
  return 0;
}

bool PT2_Selector::NoJetTrigger(const Vec4D_Vector &mom)
{
  if (!m_strong) return Trigger(mom);
  msg_Error()<<"PT2_Selector::NoJetTrigger: IR unsave cut"<<std::endl;
  return 0;
}

void PT2_Selector::BuildCuts(Cut_Data * cuts) 
{
}

void PT2_Selector::SetRange(std::vector<Flavour> crit,
			      double _min, double _max)
{
  if (crit.size() != 2) {
    msg_Error()<<"Wrong number of arguments in PT2_Selector::SetRange : "
			  <<crit.size()<<endl;
    return;
  }

  for (int i=m_nin;i<m_n;i++) {
    for (int j=i+1;j<m_n;j++) {
      if ( ((crit[0].Includes(m_fl[i])) && (crit[1].Includes(m_fl[j])) ) || 
	   ((crit[0].Includes(m_fl[j])) && (crit[1].Includes(m_fl[i])) ) ) {
	pt2min[i][j] = pt2min[j][i] = _min; 
	pt2max[i][j] = pt2max[j][i] = _max; 
      if (m_fl[i].Strong()||m_fl[j].Strong()) m_strong = 1;
      }
    }
  }
}

namespace PHASIC{

DECLARE_ND_GETTER(PT2_Selector_Getter,"PT2",Selector_Base,Selector_Key,true);

Selector_Base *PT2_Selector_Getter::operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<4) THROW(critical_error,"Invalid syntax");
  int crit1=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][0]));
  int crit2=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][1]));
  double min=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][2]));
  double max=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][3]));
  Flavour flav = Flavour((kf_code)abs(crit1));
  if (crit1<0) flav = flav.Bar();
  Flavour_Vector critflavs(1,crit1);
  flav = Flavour((kf_code)abs(crit2));
  if (crit2<0) flav = flav.Bar();
  critflavs.push_back(flav);
  PT2_Selector *sel = new PT2_Selector
    (key.p_proc->NIn(),key.p_proc->NOut(),
     (Flavour*)&key.p_proc->Process()->Flavours().front());
  sel->SetRange(critflavs,min,max);
  return sel;
}

void PT2_Selector_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"PT2 selector"; 
}

}

/*--------------------------------------------------------------------

  Invariant Mass Selector

  --------------------------------------------------------------------*/

IMass_Selector::IMass_Selector(int _nin,int _nout, Flavour * _fl):
  Selector_Base("Mass_Selector")
{
  m_nin  = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl   = _fl;
  m_smin = 0.;
  m_smax = rpa.gen.Ecms()*rpa.gen.Ecms();
  m_strong = 0;
  
  massmin = new double*[m_n];
  massmax = new double*[m_n];
  value   = new double[m_n*m_n];

  for (int i=0;i<m_n;i++) { 
    massmin[i] = new double[m_n]; 
    massmax[i] = new double[m_n];
  }

  for (int i=m_nin;i<m_n;i++) {
    for (int j=i+1;j<m_n;j++) {
      massmin[i][j] = massmin[j][i] = 0.; 
      massmax[i][j] = massmax[j][i] = 2.*(rpa.gen.PBeam(0)[0]+rpa.gen.PBeam(1)[0]); 
    }
  }
  m_sel_log = new Selector_Log(m_name);
}

IMass_Selector::~IMass_Selector() {
  for (int i=0;i<m_n;i++) {
    delete [] massmin[i];
    delete [] massmax[i];
  }
  delete [] massmin;
  delete [] massmax;
  delete [] value;
}

bool IMass_Selector::Trigger(const Vec4D_Vector & mom) 
{
  double massij;
  for (int i=m_nin;i<m_n;i++) {
    for (int j=i+1;j<m_n;j++) {
      massij = value[i*m_n+j] = sqrt((mom[i]+mom[j]).Abs2());
      if (m_sel_log->Hit( ((massij < massmin[i][j]) || 
			 (massij > massmax[i][j])) )) return 0;
    }
  }
  return 1;
}

bool IMass_Selector::JetTrigger(const Vec4D_Vector &,NLO_subevtlist *const subs)
{
  if (!m_strong) return 1;
  msg_Error()<<"Mass_Selector::JetTrigger: IR unsave cut"<<std::endl;
  return 0;
}

bool IMass_Selector::NoJetTrigger(const Vec4D_Vector &mom)
{
  if (!m_strong) return Trigger(mom);
  msg_Error()<<"Mass_Selector::NoJetTrigger: IR unsave cut"<<std::endl;
  return 0;
}

void IMass_Selector::BuildCuts(Cut_Data * cuts) 
{
  for (int i=m_nin;i<m_n-1;i++) {
    for (int j=i+1;j<m_n;j++) {
      cuts->scut[i][j] = cuts->scut[j][i] = Max(cuts->scut[i][j],sqr(massmin[i][j]));
    }
  }
}

void IMass_Selector::SetRange(std::vector<Flavour> crit,double _min, double _max)
{
  if (crit.size() != 2) {
    msg_Error()<<"Wrong number of arguments in Mass_Selector::SetRange : "
			  <<crit.size()<<endl;
    return;
  }

  for (int i=m_nin;i<m_n;i++) {
    for (int j=m_nin+1;j<m_n;j++) {
      if ( ((crit[0].Includes(m_fl[i])) && (crit[1].Includes(m_fl[j])) ) || 
	   ((crit[0].Includes(m_fl[j])) && (crit[1].Includes(m_fl[i])) ) ) {
	massmin[i][j] = massmin[j][i] = Max(_min,m_fl[i].SelMass()+m_fl[j].SelMass()); 
	massmax[i][j] = massmax[j][i] = _max;
	if (sqr(massmin[i][j])>m_smin) m_smin = Max(sqr(massmin[i][j]),m_smin);
	if (m_fl[i].Strong()||m_fl[j].Strong()) m_strong=true;
      }
    }
  }
}

namespace PHASIC{

DECLARE_ND_GETTER(IMass_Selector_Getter,"Mass",Selector_Base,Selector_Key,true);

Selector_Base *IMass_Selector_Getter::operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<4) THROW(critical_error,"Invalid syntax");
  int crit1=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][0]));
  int crit2=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][1]));
  double min=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][2]));
  double max=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][3]));
  Flavour flav = Flavour((kf_code)abs(crit1));
  if (crit1<0) flav = flav.Bar();
  Flavour_Vector critflavs(1,flav);
  flav = Flavour((kf_code)abs(crit2));
  if (crit2<0) flav = flav.Bar();
  critflavs.push_back(flav);
  IMass_Selector *sel = new IMass_Selector
    (key.p_proc->NIn(),key.p_proc->NOut(),
     (Flavour*)&key.p_proc->Process()->Flavours().front());
  sel->SetRange(critflavs,min,max);
  return sel;
}

void IMass_Selector_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"mass selector"; 
}

}

/*--------------------------------------------------------------------

  Virtuality Selector

  --------------------------------------------------------------------*/

IQ2_Selector::IQ2_Selector(int _nin,int _nout, Flavour * _fl):
  Selector_Base("Q2_Selector")
{
  m_nin  = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl   = _fl;
  m_smin = 0.;
  m_smax = rpa.gen.Ecms()*rpa.gen.Ecms();
  m_strong = 0;
  
  massmin = new double*[m_n];
  massmax = new double*[m_n];
  value   = new double[m_n*m_n];

  for (int i=0;i<m_n;i++) { 
    massmin[i] = new double[m_n]; 
    massmax[i] = new double[m_n];
  }

  for (int i=0;i<m_nin;i++) {
    for (int j=m_nin;j<m_n;j++) {
      massmin[i][j] = massmin[j][i] = 0.; 
      massmax[i][j] = massmax[j][i] = m_smax; 
    }
  }
  m_sel_log = new Selector_Log(m_name);
}

IQ2_Selector::~IQ2_Selector() {
  for (int i=0;i<m_n;i++) {
    delete [] massmin[i];
    delete [] massmax[i];
  }
  delete [] massmin;
  delete [] massmax;
  delete [] value;
}

bool IQ2_Selector::Trigger(const Vec4D_Vector & mom) 
{
  double massij;
  for (int i=0;i<m_nin;i++) {
    for (int j=m_nin;j<m_n;j++) {
      massij = value[i*m_n+j] = -(mom[i]-mom[j]).Abs2();
      if (m_sel_log->Hit( ((massij < massmin[i][j]) || 
			 (massij > massmax[i][j])) )) return 0;
    }
  }
  return 1;
}

bool IQ2_Selector::JetTrigger(const Vec4D_Vector &,NLO_subevtlist *const subs)
{
  if (!m_strong) return 1;
  msg_Error()<<"Q2_Selector::JetTrigger: IR unsafe cut"<<std::endl;
  return 0;
}

bool IQ2_Selector::NoJetTrigger(const Vec4D_Vector &mom)
{
  if (!m_strong) return Trigger(mom);
  msg_Error()<<"Q2_Selector::NoJetTrigger: IR unsafe cut"<<std::endl;
  return 0;
}

void IQ2_Selector::BuildCuts(Cut_Data * cuts) 
{
  for (int i=0;i<m_nin;i++) {
    for (int j=m_nin;j<m_n;j++) {
      cuts->scut[i][j] = Min(cuts->scut[i][j],-massmin[i][j]);
    }
  }
}

void IQ2_Selector::SetRange(std::vector<Flavour> crit,double _min, double _max)
{
  if (crit.size() != 2) {
    msg_Error()<<"Wrong number of arguments in Mass_Selector::SetRange : "
			  <<crit.size()<<endl;
    return;
  }

  for (int i=0;i<m_nin;i++) {
    for (int j=m_nin;j<m_n;j++) {
      if ( ((crit[0].Includes(m_fl[i])) && (crit[1].Includes(m_fl[j])) ) || 
	   ((crit[0].Includes(m_fl[j])) && (crit[1].Includes(m_fl[i])) ) ) {
	massmin[i][j] = massmin[j][i] = _min;
	massmax[i][j] = massmax[j][i] = _max;
	if (m_fl[i].Strong()||m_fl[j].Strong()) m_strong=true;
      }
    }
  }
}

namespace PHASIC{

DECLARE_ND_GETTER(IQ2_Selector_Getter,"Q2",Selector_Base,Selector_Key,true);

Selector_Base *IQ2_Selector_Getter::operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<4) THROW(critical_error,"Invalid syntax");
  int crit1=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][0]));
  int crit2=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][1]));
  double min=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][2]));
  double max=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][3]));
  Flavour flav = Flavour((kf_code)abs(crit1));
  if (crit1<0) flav = flav.Bar();
  Flavour_Vector critflavs(1,flav);
  flav = Flavour((kf_code)abs(crit2));
  if (crit2<0) flav = flav.Bar();
  critflavs.push_back(flav);
  IQ2_Selector *sel = new IQ2_Selector
    (key.p_proc->NIn(),key.p_proc->NOut(),
     (Flavour*)&key.p_proc->Process()->Flavours().front());
  sel->SetRange(critflavs,min,max);
  return sel;
}

void IQ2_Selector_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"mass selector"; 
}

}

/*--------------------------------------------------------------------

  Delta Eta Selector

  --------------------------------------------------------------------*/

Delta_Eta_Selector::Delta_Eta_Selector(int _nin,int _nout, Flavour * _fl) :
  Selector_Base("Delta_Eta_Selector")
{
  m_nin  = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl   = _fl;
  m_smin = 0.;
  m_smax = rpa.gen.Ecms()*rpa.gen.Ecms();
  
  detamin = new double*[m_n];
  detamax = new double*[m_n];
  value = new double[m_n*m_n];

  for (int i=0;i<m_n;i++) { 
    detamin[i] = new double[m_n]; 
    detamax[i] = new double[m_n];
  }

  for (int i=m_nin;i<m_n;i++) {
    for (int j=i+1;j<m_n;j++) {
      detamin[i][j] = detamin[j][i] = 0.; 
      detamax[i][j] = detamax[j][i] = 200.;
    }
  }
  m_sel_log = new Selector_Log(m_name);
}

Delta_Eta_Selector::~Delta_Eta_Selector() {
  for (int i=0;i<m_n;i++) {
    delete [] detamin[i];
    delete [] detamax[i];
  }
  delete [] detamin;
  delete [] detamax;
  delete [] value;
}

bool Delta_Eta_Selector::Trigger(const Vec4D_Vector & mom) 
{
  double detaij;
  for (int i=m_nin;i<m_n;i++) {
    for (int j=i+1;j<m_n;j++) {
      detaij = value[i*m_n+j] = mom[i].DEta(mom[j]);
//       PRINT_INFO("("<<m_fl[i]<<" "<<m_fl[j]<<") : "<<detaij
// 		 <<" in {"<<detamin[i][j]<<", "<<detamax[i][j]<<"}");
      if (m_sel_log->Hit( ((detaij < detamin[i][j]) || 
			   (detaij > detamax[i][j])) )) return 0;
    }
  }
  return 1;
}

void Delta_Eta_Selector::BuildCuts(Cut_Data * cuts) {}

void Delta_Eta_Selector::SetRange(std::vector<Flavour> crit,double _min, double _max)
{
  if (crit.size() != 2) {
    msg_Error()<<"Wrong number of arguments in Delta_Eta_Selector::SetRange : "
			  <<crit.size()<<endl;
    return;
  }

  for (int i=m_nin;i<m_n;i++) {
    for (int j=m_nin+1;j<m_n;j++) {
      if ( ((crit[0].Includes(m_fl[i])) && (crit[1].Includes(m_fl[j])) ) || 
	   ((crit[0].Includes(m_fl[j])) && (crit[1].Includes(m_fl[i])) ) ) {
	detamin[i][j] = detamin[j][i] = _min;
	detamax[i][j] = detamax[j][i] = _max;
      }
    }
  }
}

namespace PHASIC{

DECLARE_ND_GETTER(Delta_Eta_Selector_Getter,"DeltaEta",Selector_Base,Selector_Key,true);

Selector_Base *Delta_Eta_Selector_Getter::operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<4) THROW(critical_error,"Invalid syntax");
  int crit1=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][0]));
  int crit2=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][1]));
  double min=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][2]));
  double max=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][3]));
  Flavour flav = Flavour((kf_code)abs(crit1));
  if (crit1<0) flav = flav.Bar();
  Flavour_Vector critflavs(1,flav);
  flav = Flavour((kf_code)abs(crit2));
  if (crit2<0) flav = flav.Bar();
  critflavs.push_back(flav);
  Delta_Eta_Selector *sel = new Delta_Eta_Selector
    (key.p_proc->NIn(),key.p_proc->NOut(),
     (Flavour*)&key.p_proc->Process()->Flavours().front());
  sel->SetRange(critflavs,min,max);
  return sel;
}

void Delta_Eta_Selector_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"\\Delta\\eta selector"; 
}

}

/*--------------------------------------------------------------------

  Delta Phi Selector

  --------------------------------------------------------------------*/

Delta_Phi_Selector::Delta_Phi_Selector(int _nin,int _nout, Flavour * _fl):
  Selector_Base("Delta_Phi_Selector")
{
  m_nin  = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl   = _fl;
  m_smin = 0.;
  m_smax = rpa.gen.Ecms()*rpa.gen.Ecms();
  
  dphimin = new double*[m_n];
  dphimax = new double*[m_n];
  value = new double[m_n*m_n];

  for (int i=0;i<m_n;i++) { 
    dphimin[i] = new double[m_n]; 
    dphimax[i] = new double[m_n];
  }

  for (int i=m_nin;i<m_n;i++) {
    for (int j=i+1;j<m_n;j++) {
      dphimin[i][j] = dphimin[j][i] = 0.; 
      dphimax[i][j] = dphimax[j][i] = 200.;
    }
  }
  m_sel_log = new Selector_Log(m_name);
}

Delta_Phi_Selector::~Delta_Phi_Selector() {
  for (int i=0;i<m_n;i++) {
    delete [] dphimin[i];
    delete [] dphimax[i];
  }
  delete [] dphimin;
  delete [] dphimax;
  delete [] value;
}

bool Delta_Phi_Selector::Trigger(const Vec4D_Vector & mom) 
{
  double dphiij;
  for (int i=m_nin;i<m_n;i++) {
    for (int j=i+1;j<m_n;j++) {
      dphiij = value[i*m_n+j] = mom[i].DPhi(mom[j]);
//       PRINT_INFO("("<<m_fl[i]<<" "<<m_fl[j]<<") : "<<dphiij
// 		 <<" in {"<<dphimin[i][j]<<", "<<dphimax[i][j]<<"}");
      if (m_sel_log->Hit( ((dphiij < dphimin[i][j]) || 
			   (dphiij > dphimax[i][j])) )) return 0;
    }
  }
  return 1;
}

void Delta_Phi_Selector::BuildCuts(Cut_Data * cuts) {}

void Delta_Phi_Selector::SetRange(std::vector<Flavour> crit,double _min, double _max)
{
  if (crit.size() != 2) {
    msg_Error()<<"Wrong number of arguments in Delta_Phi_Selector::SetRange : "
			  <<crit.size()<<endl;
    return;
  }

  for (int i=m_nin;i<m_n;i++) {
    for (int j=m_nin+1;j<m_n;j++) {
      if ( ((crit[0].Includes(m_fl[i])) && (crit[1].Includes(m_fl[j])) ) || 
	   ((crit[0].Includes(m_fl[j])) && (crit[1].Includes(m_fl[i])) ) ) {
	dphimin[i][j] = dphimin[j][i] = _min;
	dphimax[i][j] = dphimax[j][i] = _max;
      }
    }
  }
}

namespace PHASIC{

DECLARE_ND_GETTER(Delta_Phi_Selector_Getter,"DeltaPhi",Selector_Base,Selector_Key,true);

Selector_Base *Delta_Phi_Selector_Getter::operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<4) THROW(critical_error,"Invalid syntax");
  int crit1=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][0]));
  int crit2=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][1]));
  double min=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][2]));
  double max=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][3]));
  Flavour flav = Flavour((kf_code)abs(crit1));
  if (crit1<0) flav = flav.Bar();
  Flavour_Vector critflavs(1,flav);
  flav = Flavour((kf_code)abs(crit2));
  if (crit2<0) flav = flav.Bar();
  critflavs.push_back(flav);
  Delta_Phi_Selector *sel = new Delta_Phi_Selector
    (key.p_proc->NIn(),key.p_proc->NOut(),
     (Flavour*)&key.p_proc->Process()->Flavours().front());
  sel->SetRange(critflavs,min,max);
  return sel;
}

void Delta_Phi_Selector_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"\\Delta\\phi selector"; 
}

}

/*--------------------------------------------------------------------

  Delta R Selector

  --------------------------------------------------------------------*/

Delta_R_Selector::Delta_R_Selector(int _nin,int _nout, Flavour * _fl):
  Selector_Base("Delta_R_Selector")
{
  m_nin  = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl   = _fl;
  m_smin = 0.;
  m_smax = rpa.gen.Ecms()*rpa.gen.Ecms();
  
  drmin = new double*[m_n];
  drmax = new double*[m_n];
  value = new double[m_n*m_n];

  for (int i=0;i<m_n;i++) { 
    drmin[i] = new double[m_n]; 
    drmax[i] = new double[m_n];
  }

  for (int i=m_nin;i<m_n;i++) {
    for (int j=i+1;j<m_n;j++) {
      drmin[i][j] = drmin[j][i] = 0.; 
      drmax[i][j] = drmax[j][i] = 200.;
    }
  }
  m_sel_log = new Selector_Log(m_name);
}

Delta_R_Selector::~Delta_R_Selector() {
  for (int i=0;i<m_n;i++) {
    delete [] drmin[i];
    delete [] drmax[i];
  }
  delete [] drmin;
  delete [] drmax;
  delete [] value;
}

bool Delta_R_Selector::Trigger(const Vec4D_Vector & mom) 
{
  double drij;
  for (int i=m_nin;i<m_n;i++) {
    for (int j=i+1;j<m_n;j++) {
      drij = value[i*m_n+j] = mom[i].DR(mom[j]);
//       PRINT_INFO("("<<m_fl[i]<<" "<<m_fl[j]<<") : "<<drij
// 		 <<" in {"<<drmin[i][j]<<", "<<drmax[i][j]<<"}");
      if (m_sel_log->Hit( ((drij < drmin[i][j]) || 
			   (drij > drmax[i][j])) )) return 0;
    }
  }
  return 1;
}

void Delta_R_Selector::BuildCuts(Cut_Data * cuts) {}

void Delta_R_Selector::SetRange(std::vector<Flavour> crit,double _min, double _max)
{
  if (crit.size() != 2) {
    msg_Error()<<"Wrong number of arguments in Delta_R_Selector::SetRange : "
			  <<crit.size()<<endl;
    return;
  }

  for (int i=m_nin;i<m_n;i++) {
    for (int j=m_nin+1;j<m_n;j++) {
      if ( ((crit[0].Includes(m_fl[i])) && (crit[1].Includes(m_fl[j])) ) || 
	   ((crit[0].Includes(m_fl[j])) && (crit[1].Includes(m_fl[i])) ) ) {
	drmin[i][j] = drmin[j][i] = _min;
	drmax[i][j] = drmax[j][i] = _max;
      }
    }
  }
}

namespace PHASIC{

DECLARE_ND_GETTER(Delta_R_Selector_Getter,"DeltaR",Selector_Base,Selector_Key,true);

Selector_Base *Delta_R_Selector_Getter::operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<3) THROW(critical_error,"Invalid syntax");
  int crit1=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][0]));
  int crit2=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][1]));
  double min=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][2]));
  double max=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][3]));
  Flavour flav = Flavour((kf_code)abs(crit1));
  if (crit1<0) flav = flav.Bar();
  Flavour_Vector critflavs(1,flav);
  flav = Flavour((kf_code)abs(crit2));
  if (crit2<0) flav = flav.Bar();
  critflavs.push_back(flav);
  Delta_R_Selector *sel = new Delta_R_Selector
    (key.p_proc->NIn(),key.p_proc->NOut(),
     (Flavour*)&key.p_proc->Process()->Flavours().front());
  sel->SetRange(critflavs,min,max);
  return sel;
}

void Delta_R_Selector_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"\\Delta R selector"; 
}

}

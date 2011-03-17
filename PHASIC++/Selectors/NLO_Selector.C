#ifndef ATOOLS_Phys_Standard_Selector_H
#define ATOOLS_Phys_Standard_Selector_H

#include "PHASIC++/Selectors/Selector.H"
#include "ATOOLS/Org/Data_Reader.H"

namespace PHASIC {

  class PTNLO_Selector : public Selector_Base {
    std::vector<double>  ptmin, ptmax;
    std::vector<ATOOLS::Flavour> flav;
    int     m_strong;
  public:
    PTNLO_Selector(int,int,ATOOLS::Flavour *);
    ~PTNLO_Selector();
    void     SetRange(ATOOLS::Flavour_Vector,double,double);
    bool     Trigger(const ATOOLS::Vec4D_Vector & );
    bool     JetTrigger(const ATOOLS::Vec4D_Vector &,ATOOLS::NLO_subevtlist *const);
    bool     NoJetTrigger(const ATOOLS::Vec4D_Vector &);
    void     BuildCuts(Cut_Data *);
  };

  class RapidityNLO_Selector : public Selector_Base {
    std::vector<double>  ymin, ymax;
    std::vector<ATOOLS::Flavour> flav;
    int     m_strong;
  public:
    RapidityNLO_Selector(int,int,ATOOLS::Flavour *);
    ~RapidityNLO_Selector();
    void     SetRange(ATOOLS::Flavour_Vector,double,double);
    bool     Trigger(const ATOOLS::Vec4D_Vector & );
    bool     JetTrigger(const ATOOLS::Vec4D_Vector &,ATOOLS::NLO_subevtlist *const);
    bool     NoJetTrigger(const ATOOLS::Vec4D_Vector &);
    void     BuildCuts(Cut_Data *);
  };

  class PseudoRapidityNLO_Selector : public Selector_Base {
    std::vector<double>  etamin, etamax;
    std::vector<ATOOLS::Flavour> flav;
    int     m_strong;
  public:
    PseudoRapidityNLO_Selector(int,int,ATOOLS::Flavour *);
    ~PseudoRapidityNLO_Selector();
    void     SetRange(ATOOLS::Flavour_Vector,double,double);
    bool     Trigger(const ATOOLS::Vec4D_Vector & );
    bool     JetTrigger(const ATOOLS::Vec4D_Vector &,ATOOLS::NLO_subevtlist *const);
    bool     NoJetTrigger(const ATOOLS::Vec4D_Vector &);
    void     BuildCuts(Cut_Data *);
  };

  class PT2NLO_Selector : public Selector_Base {
    std::vector<double>  ptmin, ptmax;
    std::vector<ATOOLS::Flavour> flav1,flav2;
    int     m_strong;
  public:
    PT2NLO_Selector(int,int,ATOOLS::Flavour *);
    ~PT2NLO_Selector();
    void     SetRange(ATOOLS::Flavour_Vector,double,double);
    bool     Trigger(const ATOOLS::Vec4D_Vector & );
    bool     JetTrigger(const ATOOLS::Vec4D_Vector &,ATOOLS::NLO_subevtlist *const);
    bool     NoJetTrigger(const ATOOLS::Vec4D_Vector &);
    void     BuildCuts(Cut_Data *);
  };

  class Isolation_Cut : public Selector_Base {
    int    m_mode;
    double m_d0;
    double m_emax;
    std::vector<int> m_if;
    ATOOLS::Flavour m_iflav;
    
    double Chi(double eg,double dr);
    double DR(const ATOOLS::Vec4D & p1,const ATOOLS::Vec4D & p2);
    double DEta12(const ATOOLS::Vec4D &,const ATOOLS::Vec4D &);
    double DPhi12(const ATOOLS::Vec4D &,const ATOOLS::Vec4D &);
  public:
    Isolation_Cut(int,int,ATOOLS::Flavour*,int);
    
    void   SetRange(ATOOLS::Flavour_Vector,double,double);
    bool   Trigger(const ATOOLS::Vec4D_Vector &);
    bool   JetTrigger(const ATOOLS::Vec4D_Vector &,ATOOLS::NLO_subevtlist *const);
    bool   NoJetTrigger(const ATOOLS::Vec4D_Vector & p) {return 1;}
    void   BuildCuts(Cut_Data *);
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


 /*--------------------------------------------------------------------
 
  PTNLO Selector

  --------------------------------------------------------------------*/

PTNLO_Selector::PTNLO_Selector(int _nin,int _nout, Flavour * _fl):
  Selector_Base("PTNLO_Selector")
{
  m_nin  = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl   = _fl;
  m_smin = 0.;
  m_smax = rpa.gen.Ecms()*rpa.gen.Ecms();
  m_strong = 0;
  if (m_nin==2) if (m_fl[0].Strong()&&m_fl[1].Strong()) m_strong = -1;
  
  m_sel_log = new Selector_Log(m_name);
}

PTNLO_Selector::~PTNLO_Selector() {
}


bool PTNLO_Selector::Trigger(const Vec4D_Vector & mom) 
{
  double pti;
  for (size_t k=0;k<flav.size();k++) {
    for (int i=m_nin;i<m_n;i++) {
      if (flav[k].Includes(m_fl[i])) {
	pti = sqrt(sqr(mom[i][1]) + sqr(mom[i][2]));
	if (m_sel_log->Hit( ((pti<ptmin[k]) || (pti>ptmax[k])) )) return 0;
      }
    }
  }
  return 1;
}

bool PTNLO_Selector::JetTrigger(const Vec4D_Vector &mom,ATOOLS::NLO_subevtlist *const subs)
{
  if (m_strong==0) return 1;
  if (m_strong==-1) {
    double pti;
    for (size_t k=0;k<flav.size();k++) {
      for (size_t i=m_nin;i<subs->back()->m_n;i++) {
	if (flav[k].Includes(subs->back()->p_fl[i])) {
	  pti = sqrt(sqr(mom[i][1]) + sqr(mom[i][2]));
	  if (m_sel_log->Hit( ((pti<ptmin[k]) || (pti>ptmax[k])) )) return 0;
	} 
      }
    }
    return 1;
  }
  msg_Error()<<"PTNLO_Selector::JetTrigger: IR unsave cut"<<std::endl;
  return 0;
}

bool PTNLO_Selector::NoJetTrigger(const Vec4D_Vector & mom)
{
  if (m_strong==0) return Trigger(mom);
  return 1;
}


void PTNLO_Selector::BuildCuts(Cut_Data * cuts) 
{
}

void PTNLO_Selector::SetRange(std::vector<Flavour> crit,double _min, 
			      double _max)
{
  if (crit.size() != 1) {
    msg_Error()<<"Wrong number of arguments in PTNLO_Selector::SetRange : "
	       <<crit.size()<<endl;
    return;
  }

  double MaxPTmin = 0.;
  flav.push_back(crit[0]);
  ptmin.push_back(_min);
  ptmax.push_back(Min(_max,(rpa.gen.PBeam(0)[0]+rpa.gen.PBeam(1)[0])));
  bool used=0;
  for (int i=m_nin;i<m_n;i++) {
    if (crit[0].Includes(m_fl[i])) {
      used=1;
      if (m_fl[i].Strong()) m_strong = 1;
      MaxPTmin=_min;
      break;
    }
  }
  if (!used) {
    flav.pop_back();
    ptmin.pop_back();
    ptmax.pop_back();
  }
  m_smin = Max(m_smin,4.*MaxPTmin*MaxPTmin);
}


namespace PHASIC{

DECLARE_ND_GETTER(PTNLO_Selector_Getter,"PTNLO",Selector_Base,Selector_Key,true);

Selector_Base *PTNLO_Selector_Getter::operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<3) THROW(critical_error,"Invalid syntax");
  int crit1=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][0]));
  double min=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][1]));
  double max=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][2]));
  Flavour flav = Flavour((kf_code)abs(crit1));
  if (crit1<0) flav = flav.Bar();
  Flavour_Vector critflavs(1,flav);
  PTNLO_Selector *sel = new PTNLO_Selector
    (key.p_proc->NIn(),key.p_proc->NOut(),
     (Flavour*)&key.p_proc->Process()->Flavours().front());
  sel->SetRange(critflavs,min,max);
  return sel;
}

void PTNLO_Selector_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"transverse momentum selector for NLO"; 
}

}


/*--------------------------------------------------------------------

  RapidityNLO Selector

  --------------------------------------------------------------------*/

RapidityNLO_Selector::RapidityNLO_Selector(int _nin,int _nout, Flavour * _fl):
  Selector_Base("RapidityNLO_Selector")
{
  m_nin  = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl   = _fl;
  m_smin = 0.;
  m_smax = rpa.gen.Ecms()*rpa.gen.Ecms();
  m_strong = 0;
  if (m_nin==2) if (m_fl[0].Strong()&&m_fl[1].Strong()) m_strong = -1;
  
  m_sel_log = new Selector_Log(m_name);
}

RapidityNLO_Selector::~RapidityNLO_Selector() {
}


bool RapidityNLO_Selector::Trigger(const Vec4D_Vector & mom) 
{
  double yi;
  for (size_t k=0;k<flav.size();k++) {
    for (int i=m_nin;i<m_n;i++) {
      if (flav[k].Includes(m_fl[i])) {
	yi = mom[i].Y();
	if (m_sel_log->Hit( ((yi<ymin[k]) || (yi>ymax[k])) )) return 0;
      }
    }
  }
  return 1;
}

bool RapidityNLO_Selector::JetTrigger(const Vec4D_Vector &mom,ATOOLS::NLO_subevtlist *const subs)
{
  if (m_strong==0) return 1;
  if (m_strong==-1) {
    double yi;
    for (size_t k=0;k<flav.size();k++) {
      for (size_t i=m_nin;i<subs->back()->m_n;i++) {
	if (flav[k].Includes(subs->back()->p_fl[i])) {
	  yi = mom[i].Y();
	  if (m_sel_log->Hit( ((yi<ymin[k]) || (yi>ymax[k])) )) return 0;
	} 
      }
    }
    return 1;
  }
  msg_Error()<<"RapidityNLO_Selector::JetTrigger: IR unsave cut"<<std::endl;
  return 0;
}

bool RapidityNLO_Selector::NoJetTrigger(const Vec4D_Vector & mom)
{
  if (m_strong==0) return Trigger(mom);
  return 1;
}


void RapidityNLO_Selector::BuildCuts(Cut_Data * cuts) 
{
}

void RapidityNLO_Selector::SetRange(std::vector<Flavour> crit,double _min, 
			       double _max)
{
  if (crit.size() != 1) {
    msg_Error()<<"Wrong number of arguments in RapidityNLO_Selector::SetRange : "
	       <<crit.size()<<endl;
    return;
  }
  if (_min != -_max) {
    msg_Error()<<"Asymetric cuts not allowed in RapidityNLO_Selector::SetRange : "<<endl;
    return;
  }

  flav.push_back(crit[0]);
  ymin.push_back(_min);
  ymax.push_back(_max);
  bool used=0;
  for (int i=m_nin;i<m_n;i++) {
    if (crit[0].Includes(m_fl[i])) {
      used=1;
      if (m_fl[i].Strong()) m_strong = 1;
      break;
    }
  }
  if (!used) {
    flav.pop_back();
    ymin.pop_back();
    ymax.pop_back();
  }
}

namespace PHASIC{

DECLARE_ND_GETTER(RapidityNLO_Selector_Getter,"RapidityNLO",Selector_Base,Selector_Key,true);

Selector_Base *RapidityNLO_Selector_Getter::operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<3) THROW(critical_error,"Invalid syntax");
  int crit1=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][0]));
  double min=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][1]));
  double max=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][2]));
  Flavour flav = Flavour((kf_code)abs(crit1));
  if (crit1<0) flav = flav.Bar();
  Flavour_Vector critflavs(1,flav);
  RapidityNLO_Selector *sel = new RapidityNLO_Selector
    (key.p_proc->NIn(),key.p_proc->NOut(),
     (Flavour*)&key.p_proc->Process()->Flavours().front());
  sel->SetRange(critflavs,min,max);
  return sel;
}

void RapidityNLO_Selector_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"rapidity selector for NLO"; 
}

}
/*--------------------------------------------------------------------

  PseudoRapidityNLO Selector

  --------------------------------------------------------------------*/

PseudoRapidityNLO_Selector::PseudoRapidityNLO_Selector(int _nin,int _nout, Flavour * _fl):
  Selector_Base("PseudoRapidityNLO_Selector")
 {
  m_nin  = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl   = _fl;
  m_smin = 0.;
  m_smax = rpa.gen.Ecms()*rpa.gen.Ecms();
  m_strong = 0;
  if (m_nin==2) if (m_fl[0].Strong()&&m_fl[1].Strong()) m_strong = -1;
  
  m_sel_log = new Selector_Log(m_name);
}

PseudoRapidityNLO_Selector::~PseudoRapidityNLO_Selector() {
}


bool PseudoRapidityNLO_Selector::Trigger(const Vec4D_Vector & mom) 
{
  double etai;
  for (size_t k=0;k<flav.size();k++) {
    for (int i=m_nin;i<m_n;i++) {
      if (flav[k].Includes(m_fl[i])) {
	etai = mom[i].Eta();
	if (m_sel_log->Hit( ((etai<etamin[k]) || (etai>etamax[k])) )) return 0;
      }
    }
  }
  return 1;
}

bool PseudoRapidityNLO_Selector::JetTrigger(const Vec4D_Vector &mom,ATOOLS::NLO_subevtlist *const subs)
{
  if (m_strong==0) return 1;
  if (m_strong==-1) {
    double etai;
    for (size_t k=0;k<flav.size();k++) {
      for (size_t i=m_nin;i<subs->back()->m_n;i++) {
	if (flav[k].Includes(subs->back()->p_fl[i])) {
	  etai = mom[i].Eta();
	  if (m_sel_log->Hit( ((etai<etamin[k]) || (etai>etamax[k])) )) return 0;
	} 
      }
    }
    return 1;
  }
  msg_Error()<<"PseudoRapidityNLO_Selector::JetTrigger: IR unsave cut"<<std::endl;
  return 0;
}

bool PseudoRapidityNLO_Selector::NoJetTrigger(const Vec4D_Vector & mom)
{
  if (m_strong==0) return Trigger(mom);
  return 1;
}


void PseudoRapidityNLO_Selector::BuildCuts(Cut_Data * cuts) 
{
}

void PseudoRapidityNLO_Selector::SetRange(std::vector<Flavour> crit,double _min, 
			       double _max)
{
  if (crit.size() != 1) {
    msg_Error()<<"Wrong number of arguments in PseudoRapidityNLO_Selector::SetRange : "
	       <<crit.size()<<endl;
    return;
  }
  if (_min != -_max) {
    msg_Error()<<"Asymetric cuts not allowed in PseudoRapidityNLO_Selector::SetRange : "<<endl;
    return;
  }

  flav.push_back(crit[0]);
  etamin.push_back(_min);
  etamax.push_back(_max);
  bool used=0;
  for (int i=m_nin;i<m_n;i++) {
    if (crit[0].Includes(m_fl[i])) {
      used=1;
      if (m_fl[i].Strong()) m_strong = 1;
      break;
    }
  }
  if (!used) {
    flav.pop_back();
    etamin.pop_back();
    etamax.pop_back();
  }
}

namespace PHASIC{

DECLARE_ND_GETTER(PseudoRapidityNLO_Selector_Getter,"PseudoRapidityNLO",Selector_Base,Selector_Key,true);

Selector_Base *PseudoRapidityNLO_Selector_Getter::operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<3) THROW(critical_error,"Invalid syntax");
  int crit1=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][0]));
  double min=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][1]));
  double max=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][2]));
  Flavour flav = Flavour((kf_code)abs(crit1));
  if (crit1<0) flav = flav.Bar();
  Flavour_Vector critflavs(1,flav);
  PseudoRapidityNLO_Selector *sel = new PseudoRapidityNLO_Selector
    (key.p_proc->NIn(),key.p_proc->NOut(),
     (Flavour*)&key.p_proc->Process()->Flavours().front());
  sel->SetRange(critflavs,min,max);
  return sel;
}

void PseudoRapidityNLO_Selector_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"pseudorapidity selector for NLO"; 
}

}
/*--------------------------------------------------------------------

  PT2NLO Selector

  --------------------------------------------------------------------*/

PT2NLO_Selector::PT2NLO_Selector(int _nin,int _nout, Flavour * _fl):
  Selector_Base("PT2NLO_Selector")
{
  m_nin  = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl   = _fl;
  m_smin = 0.;
  m_smax = rpa.gen.Ecms()*rpa.gen.Ecms();
  m_strong = 0;
  if (m_nin==2) if (m_fl[0].Strong()&&m_fl[1].Strong()) m_strong = -1;
  
  m_sel_log = new Selector_Log(m_name);
}

PT2NLO_Selector::~PT2NLO_Selector() {
}


bool PT2NLO_Selector::Trigger(const Vec4D_Vector & mom) 
{
  double ptij;
  for (size_t k=0;k<flav1.size();k++) {
    for (int i=m_nin;i<m_n;i++) {
      for (int j=i+1;j<m_n;j++) {
	if ( ((flav1[k].Includes(m_fl[i])) && (flav2[k].Includes(m_fl[j])) ) || 
	     ((flav1[k].Includes(m_fl[j])) && (flav2[k].Includes(m_fl[i])) ) ) {
	  ptij = (mom[i]+mom[j]).PPerp();
	  if (m_sel_log->Hit( ((ptij<ptmin[k]) || (ptij>ptmax[k])) )) return 0;
	}
      }
    }
  }
  return 1;
}

bool PT2NLO_Selector::JetTrigger(const Vec4D_Vector &mom,ATOOLS::NLO_subevtlist *const subs)
{
  if (m_strong==0) return 1;
  if (m_strong==-1) {
    double ptij;
    for (size_t k=0;k<flav1.size();k++) {
      for (size_t i=m_nin;i<subs->back()->m_n;i++) {
	for (size_t j=i+1;j<subs->back()->m_n;j++) {
	  if ( ((flav1[k].Includes(subs->back()->p_fl[i])) &&
		(flav2[k].Includes(subs->back()->p_fl[j])) ) || 
	       ((flav1[k].Includes(subs->back()->p_fl[j])) &&
		(flav2[k].Includes(subs->back()->p_fl[i])) ) ) {
	    ptij = (mom[i]+mom[j]).PPerp();
	    if (m_sel_log->Hit( ((ptij<ptmin[k]) || (ptij>ptmax[k])) )) return 0;
	  }
	}
      }
    }
    return 1;
  }
  msg_Error()<<"PTNLO_Selector::JetTrigger: IR unsave cut"<<std::endl;
  return 0;
}

bool PT2NLO_Selector::NoJetTrigger(const Vec4D_Vector &mom)
{
  if (m_strong==0) return Trigger(mom);
  return 1;
}


void PT2NLO_Selector::BuildCuts(Cut_Data * cuts) 
{
}

void PT2NLO_Selector::SetRange(std::vector<Flavour> crit,double _min, 
			       double _max)
{
  if (crit.size() != 2) {
    msg_Error()<<"Wrong number of arguments in PTNLO_Selector::SetRange : "
	       <<crit.size()<<endl;
    return;
  }

  double MaxPTmin = 0.;
  flav1.push_back(crit[0]);
  flav2.push_back(crit[1]);
  ptmin.push_back(_min);
  ptmax.push_back(Min(_max,(rpa.gen.PBeam(0)[0]+rpa.gen.PBeam(1)[0])));
  bool used=0;
  for (int i=m_nin;i<m_n;i++) {
    for (int j=i+1;j<m_n;j++) {
      if ( ((crit[0].Includes(m_fl[i])) && (crit[1].Includes(m_fl[j])) ) || 
	   ((crit[0].Includes(m_fl[j])) && (crit[1].Includes(m_fl[i])) ) ) {
	used=1;
	MaxPTmin = _min;
	if (m_fl[i].Strong()||m_fl[j].Strong()) m_strong = 1;
	break;
      }
    }
  }
  if (!used) {
    flav1.pop_back();
    flav2.pop_back();
    ptmin.pop_back();
    ptmax.pop_back();
  }
  m_smin = Max(m_smin,4.*MaxPTmin*MaxPTmin);
}


namespace PHASIC{

DECLARE_ND_GETTER(PT2NLO_Selector_Getter,"PT2NLO",Selector_Base,Selector_Key,true);

Selector_Base *PT2NLO_Selector_Getter::operator()(const Selector_Key &key) const
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
  PT2NLO_Selector *sel = new PT2NLO_Selector
    (key.p_proc->NIn(),key.p_proc->NOut(),
     (Flavour*)&key.p_proc->Process()->Flavours().front());
  sel->SetRange(critflavs,min,max);
  return sel;
}

void PT2NLO_Selector_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"PT2NLO selector"; 
}

}

/*--------------------------------------------------------------------

  photon isolation cut: hep-ph/9801442
 
  --------------------------------------------------------------------*/

Isolation_Cut::Isolation_Cut(int nin,int nout,Flavour * _fl,int mode) : 
  Selector_Base("IsolationCut"), m_mode(mode)
{
  m_nin  = nin; 
  m_nout = nout; 
  m_smin = 0.;
  m_fl   = _fl;
  
  m_sel_log = new Selector_Log(m_name);
}

void Isolation_Cut::SetRange(std::vector<Flavour> crit,double _min, 
			     double _max)
{
  if (crit.size() != 1) {
    msg_Error()<<"Wrong number of arguments in Isolation_Cut::SetRange : "
			  <<crit.size()<<endl;
    return;
  }
  m_d0    = _min;
  m_emax  = _max;
  m_iflav = crit[0];
  m_if.clear();
  for (int i=m_nin;i<m_nin+m_nout;i++) {
    if (crit[0].Includes(m_fl[i])) {
      m_if.push_back(i);
    }
  }
}

class edr {
public:
  double E;
  double dr;
  edr(double _e,double _dr) : E(_e), dr(_dr) {}
};
class Order_edr {
public:
  int operator()(const edr a, const edr b);
};
int Order_edr::operator()(const edr a, const edr b) {
  if (a.dr<b.dr) return 1;
  return 0;
}

bool Isolation_Cut::Trigger(const Vec4D_Vector & p)
{
  for (size_t k=0;k<m_if.size();k++) {
    double egamma=p[m_if[k]].PPerp();
    vector<edr> edrlist;
    for (int i=m_nin;i<m_nin+m_nout;i++) if (m_fl[i].Strong()) {
      double dr=DR(p[m_if[k]],p[i]);
      if (dr<m_d0) edrlist.push_back(edr(p[i].PPerp(),dr));
    }
    if (edrlist.size()>0) {
      stable_sort(edrlist.begin(),edrlist.end(),Order_edr());
      double etot=0.;
      for (size_t i=0;i<edrlist.size();i++) {
	etot+=edrlist[i].E;
 	if (m_sel_log->Hit(etot>Chi(egamma,edrlist[i].dr))) return 0;
      }
      edrlist.clear();
    }
  }
  return 1;
}

bool Isolation_Cut::JetTrigger(const Vec4D_Vector &p,ATOOLS::NLO_subevtlist *const subs)
{
  vector<int> vf;
  for (size_t i=m_nin;i<subs->back()->m_n;i++) {
    if (m_iflav.Includes(subs->back()->p_fl[i])) {
      vf.push_back(i);
    }
  }
  for (size_t k=0;k<vf.size();k++) {
    double egamma=p[vf[k]].PPerp();
    vector<edr> edrlist;
    for (size_t i=m_nin;i<subs->back()->m_n;i++)
      if (subs->back()->p_fl[i].Strong()) {
      double dr=DR(p[vf[k]],p[i]);
      if (dr<m_d0) edrlist.push_back(edr(p[i].PPerp(),dr));
    }
    if (edrlist.size()>0) {
      stable_sort(edrlist.begin(),edrlist.end(),Order_edr());
      double etot=0.;
      for (size_t i=0;i<edrlist.size();i++) {
	etot+=edrlist[i].E;
	if (m_sel_log->Hit(etot>Chi(egamma,edrlist[i].dr))) return 0;
      }
      edrlist.clear();
    }
  }
  return 1;
}


void Isolation_Cut::BuildCuts(Cut_Data * cuts) 
{
}

double Isolation_Cut::Chi(double eg,double dr)
{
  if (m_mode==0) return m_emax;
  if (m_mode<0) return 0.;//rpa.gen.Ecms();
  return eg*m_emax*pow((1.-cos(dr))/(1.-cos(m_d0)),m_mode);
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

namespace PHASIC{

DECLARE_ND_GETTER(Isolation_Cut_Getter,"IsolationCut",Selector_Base,Selector_Key,true);

Selector_Base *Isolation_Cut_Getter::operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<3) THROW(critical_error,"Invalid syntax");
  int crit1=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][0]));
  double min=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][1]));
  int mode=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][2]));
  double emax=1.;
  if (key.front().size()>3)
    emax=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][3]));
  Flavour flav = Flavour((kf_code)abs(crit1));
  if (crit1<0) flav = flav.Bar();
  Flavour_Vector critflavs(1,flav);
  Isolation_Cut *sel = new Isolation_Cut
    (key.p_proc->NIn(),key.p_proc->NOut(),
     (Flavour*)&key.p_proc->Process()->Flavours().front(),mode);
  sel->SetRange(critflavs,min,emax);
  return sel;
}

void Isolation_Cut_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Isolation_Cut selector: hep-ph/9801442"; 
}

}


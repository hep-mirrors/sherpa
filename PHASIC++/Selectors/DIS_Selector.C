#ifndef ATOOLS_Phys_Standard_Selector_DIS_H
#define ATOOLS_Phys_Standard_Selector_DIS_H
#include "PHASIC++/Selectors/Selector.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Phys/Flavour.H"
namespace PHASIC {  
  class IPZIN_Selector : public Selector_Base {
    double * pzinmin, * pzinmax;
    bool     m_strong;
  public:
    IPZIN_Selector(int,int,ATOOLS::Flavour *);
    ~IPZIN_Selector();
    void     SetRange(ATOOLS::Flavour_Vector,double,double);
    bool     Trigger(const ATOOLS::Vec4D_Vector & );
    void     BuildCuts(Cut_Data *);
  };

  class IINEL_Selector : public Selector_Base {
    double ** ymin, ** ymax, * value;
    bool     m_strong;
  public:
    IINEL_Selector(int,int,ATOOLS::Flavour *);
    ~IINEL_Selector();
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

/*--------------------------------------------------------------------

  Inelesticity Selector

  --------------------------------------------------------------------*/

IINEL_Selector::IINEL_Selector(int _nin,int _nout, Flavour * _fl):
  Selector_Base("INEL_Selector")
{
  m_nin  = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl   = _fl;
  m_smin = 0.;
  m_smax = 1;//rpa->gen.Ecms()*rpa->gen.Ecms();
  m_strong = 0;
  
  ymin = new double*[m_n];
  ymax = new double*[m_n];
  value   = new double[m_n*m_n];

  for (int i=0;i<m_n;i++) { 
    ymin[i] = new double[m_n]; 
    ymax[i] = new double[m_n];
  }

  for (int i=0;i<m_nin;i++) {
    for (int j=m_nin;j<m_n;j++) {
      ymin[i][j] = ymin[j][i] = 0.; 
      ymax[i][j] = ymax[j][i] = m_smax; 
    }
  }
  m_sel_log = new Selector_Log(m_name);
}

IINEL_Selector::~IINEL_Selector() {
  for (int i=0;i<m_n;i++) {
    delete [] ymin[i];
    delete [] ymax[i];
  }
  delete [] ymin;
  delete [] ymax;
  delete [] value;
}


bool IINEL_Selector::Trigger(const ATOOLS::Vec4D_Vector & mom )
{
  double yij;
  for (int i=0;i<m_nin;i++) {
    for (int j=m_nin;j<m_n;j++) {
      //massij = value[i*m_n+j] = -(mom[i]-mom[j]).Abs2();
      yij = value[i*m_n+j]  = 1.0-(mom[j][0]/mom[i][0])*(1.0+mom[i].CosTheta(mom[j]))/2.0;//
      if (m_sel_log->Hit( ((yij < ymin[i][j]) || 
			   (yij > ymax[i][j])) )) return 0;
    }
  }
  return 1;
}

void IINEL_Selector::BuildCuts(Cut_Data * cuts) {}

void IINEL_Selector::SetRange(std::vector<Flavour> crit,double _min, double _max)
{
  if (crit.size() != 2) {
    msg_Error()<<"Wrong number of arguments in Mass_Selector::SetRange : "
			  <<crit.size()<<std::endl;
    return;
  }

  for (int i=0;i<m_nin;i++) {
    for (int j=m_nin;j<m_n;j++) {
      if ( ((crit[0].Includes(m_fl[i])) && (crit[1].Includes(m_fl[j])) ) || 
	   ((crit[0].Includes(m_fl[j])) && (crit[1].Includes(m_fl[i])) ) ) {
	ymin[i][j] = ymin[j][i] = _min;
	ymax[i][j] = ymax[j][i] = _max;
	if (m_fl[i].Strong()||m_fl[j].Strong()) m_strong=true;
      }
    }
  }
}

DECLARE_ND_GETTER(IINEL_Selector,"INEL",Selector_Base,Selector_Key,true);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,IINEL_Selector>::
operator()(const Selector_Key &key) const
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
  IINEL_Selector *sel = new IINEL_Selector
    (key.p_proc->NIn(),key.p_proc->NOut(),
     (Flavour*)&key.p_proc->Process()->Flavours().front());
  sel->SetRange(critflavs,min,max);
  return sel;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,IINEL_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"INEL selector";
}


/*--------------------------------------------------------------------

  Pz Selector

  --------------------------------------------------------------------*/

IPZIN_Selector::IPZIN_Selector(int _nin,int _nout, Flavour * _fl) :
  Selector_Base("IPZIN_Selector")
{
  m_nin  = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl   = _fl;
  m_smin = 0.;
  m_smax = rpa->gen.Ecms()*rpa->gen.Ecms();
  m_strong = 0;
  if (m_nin==2) if (m_fl[0].Strong()&&m_fl[1].Strong()) m_strong = 1;
  
  double Emax = rpa->gen.PBeam(0)[0]+rpa->gen.PBeam(1)[0];
  pzinmin  = new double[m_nin];
  pzinmax  = new double[m_nin];
  for (int i=0;i<m_nin;i++) { pzinmin[i] = 0.; pzinmax[i] = Emax; }
  m_sel_log = new Selector_Log(m_name);
}

IPZIN_Selector::~IPZIN_Selector() {
  delete [] pzinmin;
  delete [] pzinmax;
}


bool IPZIN_Selector::Trigger(const Vec4D_Vector & mom) 
{

  double pzini;
  for (int i=0;i<m_nin;i++) {
    pzini  = std::abs(mom[i][3]);
    if (m_sel_log->Hit( ((pzini<pzinmin[i]) || (pzini>pzinmax[i])) )) return 0;
  }
  return 1;
}

void IPZIN_Selector::BuildCuts(Cut_Data * cuts) 
{
}

void IPZIN_Selector::SetRange(std::vector<Flavour> crit,double _min,double _max)
{
  if (crit.size() != 1) {
    msg_Error()<<"Wrong number of arguments in IPZIN_Selector::SetRange : "
	       <<crit.size()<<std::endl;
    return;
  }
  for (int i=0;i<m_nin;i++) {
    if (crit[0].Includes(m_fl[i])) {
      pzinmin[i] = _min; 
      pzinmax[i] = Min(_max,(rpa->gen.PBeam(0)[0]+rpa->gen.PBeam(1)[0]));
      if (m_fl[i].Strong()) m_strong = 1;
    }    
  }
  //m_smin = Max(m_smin,4.*MaxPZmax*MaxPZmax);
}

DECLARE_ND_GETTER(IPZIN_Selector,"PZIN",Selector_Base,Selector_Key,true);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,IPZIN_Selector>::
operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<3) THROW(critical_error,"Invalid syntax");
  int crit1=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][0]));
  double min=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][1]));
  double max=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][2]));
  Flavour flav = Flavour((kf_code)abs(crit1));
  if (crit1<0) flav = flav.Bar();
  Flavour_Vector critflavs(1,flav);
  IPZIN_Selector *sel = new IPZIN_Selector
    (key.p_proc->NIn(),key.p_proc->NOut(),
     (Flavour*)&key.p_proc->Process()->Flavours().front());
  sel->SetRange(critflavs,min,max);
  return sel;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,IPZIN_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"PZIN selector"; 
}

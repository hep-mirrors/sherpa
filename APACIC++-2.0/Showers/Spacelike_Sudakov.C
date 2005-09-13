#include "Spacelike_Sudakov.H"
#include "Spacelike_Kinematics.H"
#include "Timelike_Sudakov.H"
#include "Sudakov_Tools.H"
#include "Knot.H"

#include "QCD_Splitting_Functions.H"
#include "QED_Splitting_Functions.H"
#include "Run_Parameter.H"
#include "Remnant_Base.H"
#include "Beam_Base.H"

#include <iomanip>

#ifdef PROFILE__all
#include "prof.hh"
#else 
#define PROFILE_HERE {}
#define PROFILE_LOCAL(LOCALNAME) {}
#endif


using namespace APACIC;
using namespace PDF;
using namespace ATOOLS;


Spacelike_Sudakov::Spacelike_Sudakov(PDF_Base * pdf,Sudakov_Tools * tools,Spacelike_Kinematics * kin,
				     double pt2min,ATOOLS::Data_Read * dataread,const size_t beam) : 
  Backward_Splitting_Group(0,0), p_tools(tools), p_kin(kin), 
  p_pdfa(pdf->GetBasicPDF()->GetCopy()),
  m_pt2min(dabs(pt2min)), 
  m_s_hadron(dataread->GetValue<double>("IS_MAX_SCALE",sqr(rpa.gen.Ecms()))),
  m_t0(-m_pt2min),
  m_emin(dataread->GetValue<double>("IS_MINIMAL_E",0.5)),
  m_cpl_scheme(dataread->GetValue<int>("IS_COUPLINGS",1)),
  m_pdf_scheme(dataread->GetValue<int>("IS_PDF_SCALE",1)),
  m_ordering_scheme(dataread->GetValue<int>("IS_ORDERING",0)),
  m_facscale(rpa.gen.Ycut()*m_s_hadron),
  m_pdf_fac(5.),
  m_pdf_scale_fac(dataread->GetValue<double>("IS_PDF_SCALE_FACTOR",1.)),
  m_cpl_scale_fac(dataread->GetValue<double>("IS_CPL_SCALE_FACTOR",0.25))
{
  s_kfactorscheme = dataread->GetValue<int>("S_KFACTOR_SCHEME",1);        
  p_pdf           = pdf->GetBasicPDF()->GetCopy();
  p_tools->CalculateMaxCouplings(m_cpl_scheme,m_pt2min,m_s_hadron);

  if (p_pdf->Bunch().IsHadron()) {
    for (int i=1;i<6;++i) {
      Flavour fl = Flavour(kf::code(i));
      Add(new q_qg(fl,p_tools));
      Add(new q_qg(fl.Bar(),p_tools));
      Add(new q_gq(fl,p_tools));
      Add(new q_gq(fl.Bar(),p_tools));
      Add(new g_qq(fl,p_tools));
      Add(new g_qq(fl.Bar(),p_tools));
    }
    Add(new g_gg(p_tools));
    Add(new g_gg(p_tools));
  }

  if (p_pdf->Bunch().IsLepton()) {
    Add(new f_fp(Flavour(kf::e),p_tools));
    Add(new f_fp(Flavour(kf::e).Bar(),p_tools));
    Add(new f_pf(Flavour(kf::e),p_tools));
    Add(new f_pf(Flavour(kf::e).Bar(),p_tools));
    Add(new p_ff(Flavour(kf::e),p_tools));
    Add(new p_ff(Flavour(kf::e).Bar(),p_tools));
  }

  PrintStat();
}


bool Spacelike_Sudakov::Dice(Knot * mo,double sprime,bool jetveto,int & extra_pdf) 
{
  PROFILE_HERE;
  mo->tmax = mo->t;  // last start t
  m_inflav = mo->part->Flav(); 
  m_t      = mo->t;
  m_x      = mo->x;

  if (!((m_t-m_t0)<rpa.gen.Accu())) {
    if (mo->prev) {
      mo->t    = mo->tout;
      mo->stat = 0;
      return false; 
    }
  }
  
  double xe(2.*m_emin*sqrt(sprime)/sqr(rpa.gen.Ecms()));  
  m_zmin = m_x/(1.-xe);
  m_zmax = m_x/(m_x+xe);
  if (m_zmin>m_zmax) {
    mo->t    = mo->tout;
    mo->stat = 0;
    return false; 
  }
  
  p_pdf->Calculate(m_x,(-m_t));
  
  while (m_t<m_t0) {
    CrudeInt(m_zmin,m_zmax);
    ProduceT();
    if (m_t>m_t0) {
      mo->t    = mo->tout;
      mo->stat = 0;
      return false;     
    }
    SelectOne();
    m_z   = GetZ();   
    m_pt2 = -(1.-m_z)*m_t;

    double uhat(-m_t - sprime* (1.-m_z)/m_z);

    if (uhat<0. && !Veto(mo,jetveto,extra_pdf)) {
      if (!RemnantVeto(mo)) {
	UniformPhi();
	mo->z      = m_z;
	mo->t      = m_t;
	mo->phi    = m_phi;
	mo->pt2lcm = m_pt2;
	mo->thcrit = 4.*m_z*m_z*m_t/(4.*m_z*m_z*m_t-(1.-m_z)*m_x*m_x*m_s_hadron);
	return true;
      }
      else break;
    }
  }
  mo->t    = mo->tout;
  mo->stat = 0;
  return false; 
}

void Spacelike_Sudakov::ProduceT() {
  if (m_lastint <0.) m_t = +1.;            
  else {
    double rn(ran.Get());
    if (IsZero( 2.*M_PI*log(rn) / (m_lastint*m_pdf_fac) )) m_t=m_t0;
    else m_t *= exp( 2.*M_PI*log(rn) / (m_lastint*m_pdf_fac) );
  }
  return;
}

bool Spacelike_Sudakov::Veto(Knot * mo,bool jetveto,int & extra_pdf) 
{  
  PROFILE_HERE;
  if ((1.-m_z)*m_t>m_t0)         return true;

  // 1. masses, z-range and splitting function
  if (!extra_pdf && MassVeto(0)) return true;

  // 2. alphaS
  if (CplVeto())                 return true;

  // 3. ordering
  if (OrderingVeto(mo))          return true;

  // 4. extra pdf weight for first branch below  Q_jet
  if (extra_pdf) {
    extra_pdf=0;
    if (MassVeto(1))             return true;
  }

  return false;
}

bool Spacelike_Sudakov::MassVeto(int extra_pdf) 
{
  PROFILE_HERE;
  double weight(p_pdf->GetXPDF(GetFlB())/(p_pdf->GetXPDF(GetFlA())*m_pdf_fac)); 
  double q2(-m_t), firstq2(m_facscale);
  double wb_jet(0.);

  switch (m_pdf_scheme) {
  case 0:
    firstq2/=(1.-m_z);
    break;
  default:
    q2 = m_pt2;
  }

  q2 *= m_pdf_scale_fac;


  if (extra_pdf) {
    p_pdf->Calculate(m_x,0.,0.,firstq2);
    wb_jet   = p_pdf->GetXPDF(GetFlB());
  }
  p_pdf->Calculate(m_x,q2);
  if (!extra_pdf) {
    wb_jet   = p_pdf->GetXPDF(GetFlB());
  }
  p_pdfa->Calculate(m_x/m_z,q2);
  double test = p_pdfa->GetXPDF(GetFlA());
  if (IsZero(test)) return true;
  weight     *= test/wb_jet * GetWeight(m_z,-m_t,0);

  if (ran.Get() > weight) return true;
  return false;
}

bool Spacelike_Sudakov::CplVeto() 
{
  switch (m_cpl_scheme) {
  case 0 :  return false;
  case 2 :  return (GetCoupling(m_cpl_scale_fac*m_t)/GetCoupling()   > ran.Get()) ? false : true;   
  default : return (GetCoupling(m_cpl_scale_fac*m_pt2)/GetCoupling() > ran.Get()) ? false : true;   
  }
  return true;
}

bool Spacelike_Sudakov::OrderingVeto(Knot * mo) 
{
  if (!m_inflav.Strong()) {
    mo->maxpt2 = m_pt2;
    return false;
  }
  double th(4.*m_z*m_z*m_t/(4.*m_z*m_z*m_t-(1.-m_z)*m_x*m_x*m_s_hadron));
  switch (m_ordering_scheme) {
  case 0 :
    mo->maxpt2 = m_pt2;
    return false;
  case 2 : if (th > mo->thcrit) return true;
    mo->maxpt2 = m_pt2;
    return false;
  default :
    if (m_pt2 > mo->maxpt2)     return true;
    mo->maxpt2 = m_pt2;
    return false;
  }
  return false;
}

bool Spacelike_Sudakov::RemnantVeto(Knot * mo) 
{
  double E(p_remnant->GetBeam()->Energy()*mo->x/m_z);
  Particle part(1,GetFlA(),Vec4D(E,0.0,0.0,E));
  if (!p_remnant->TestExtract(&part)) return true;
  return false;
}

void Spacelike_Sudakov::Add(Splitting_Function * spl) 
{
  for (SplFunIter iter(m_group);iter();++iter) {
    if (iter()->GetFlB()==spl->GetFlB()) {
      iter()->Add(spl);
      return ;
    }
  }
  m_group.Append(new Backward_Splitting_Group(spl,p_pdf));
  p_selected=spl;
}


double Spacelike_Sudakov::CrudeInt(double _zmin, double _zmax) 
{
  PROFILE_HERE;
  SplFunIter iter(m_group);
  for (;iter();++iter)
    if (iter()->GetFlB()==m_inflav) { p_selected=iter(); break; }
  if (!iter()) { p_selected=NULL; return m_lastint = -1.; }
  return m_lastint = p_selected->CrudeInt(_zmin,_zmax);
}     

#include "Spacelike_Sudakov.H"
#include "Spacelike_Kinematics.H"
#include "Timelike_Sudakov.H"
#include "Sudakov_Tools.H"
#include "Knot.H"

#include "QCD_Splitting_Functions.H"
#include "QED_Splitting_Functions.H"
#include "Run_Parameter.H"

#include <iomanip>


using namespace APACIC;
using namespace PDF;
using namespace ATOOLS;


Spacelike_Sudakov::Spacelike_Sudakov(PDF_Base * _pdf,Sudakov_Tools * _tools,Spacelike_Kinematics * _kin,
				     double _pt2min,ATOOLS::Data_Read * _dataread) : 
  Backward_Splitting_Group(0,0), p_tools(_tools), p_kin(_kin), m_pt2min(dabs(_pt2min)), 
  m_miveto_scheme(0), m_last_veto(0)
{
  p_pdf             = _pdf->GetBasicPDF(); 
  p_pdfa            = p_pdf->GetBasicPDF()->GetCopy();
  m_ordering_scheme = _dataread->GetValue<int>("IS ORDERING",0);  /* Switch for ordering due to coherence:  
                                                                     0 = none, 1 = pt^2, 2 = pt^2/E^2     */
  m_cpl_scheme      = _dataread->GetValue<int>("IS COUPLINGS",3); /*  (0=fix, 1=pt^2, 2=t/4)              */ 
  m_pdf_scheme      = _dataread->GetValue<int>("IS PDF SCALE",1); /*  0 = -Q^2, 1 = -(1-z)*Q^2 */
  m_pdf_scale_fac   = _dataread->GetValue<double>("IS PDF SCALE FACTOR",1.);
  m_jetveto_scheme  = _dataread->GetValue<int>("IS JETVETOSCHEME",2);

  m_emin            = .5;
  m_pt2max          = sqr(rpa.gen.Ecms());
  m_qjet            = rpa.gen.Ycut()*sqr(rpa.gen.Ecms());
  m_facscale        = m_qjet;
  m_pdf_fac         = 5.; 
  m_lambda2         = p_tools->GetLambda2(); 
  m_b               = p_tools->GetBnorm();   

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


bool Spacelike_Sudakov::Dice(Knot * mo,double sprime,bool jetveto,int & extra_pdf) {
  mo->tmax = mo->t;  // last start t
  msg.Tracking()<<" tstart="<<mo->tmax<<" ("<<mo->kn_no<<")"<<std::endl;
  m_last_veto = 0;
  m_inflav = mo->part->Flav(); 
  m_t      = mo->t;
  m_x      = mo->x;
  m_t0     = - m_pt2min;

  if (!((m_t-m_t0)<rpa.gen.Accu())) {
    if (mo->prev) {
      mo->t    = mo->tout;
      mo->stat = 0;
      return 0; 
    }
  }
  
  m_xe   = 2.*m_emin*sqrt(sprime)/sqr(rpa.gen.Ecms());  
  m_zmin = m_x/(1.-m_xe);
  m_zmax = m_x/(m_x+m_xe);
  if (m_zmin>m_zmax) {
    mo->t    = mo->tout;
    mo->stat = 0;
    return 0; 
  }
  
  p_pdf->Calculate(m_x,0.,0.,sqrt(-m_t));
  
  while (m_t<m_t0) {
    CrudeInt(m_zmin,m_zmax);
    ProduceT();
    if (m_t>m_t0) {
      mo->t    = mo->tout;
      mo->stat = 0;
      return 0;     
    }
    SelectOne();
    m_z   = GetZ();   
    m_pt2 = -(1.-m_z)*m_t;

    double uhat = -m_t - sprime* (1.-m_z)/m_z;
    if (uhat>=0.) m_last_veto=9;

    if (uhat<0. && !Veto(mo,jetveto,extra_pdf)) {
      UniformPhi();
      mo->z      = m_z;
      mo->t      = m_t;
      mo->phi    = m_phi;
      if (m_ordering_scheme==2) {
	double th = 4.*m_z*m_z*m_t/(4.*m_z*m_z*m_t-(1.-m_z)*m_x*m_x*m_pt2max);
	mo->thcrit = th;
      }
      return 1;
    }
  }
  mo->t    = mo->tout;
  mo->stat = 0;
  return 0; 
}

void Spacelike_Sudakov::ProduceT() {
  if (m_lastint <0.) m_t = +1.;            
  else m_t *= exp( 2.*M_PI*log(ran.Get()) / (m_lastint*m_pdf_fac) );
  return;
}

static double rem_ran, rem_weight;

bool Spacelike_Sudakov::Veto(Knot * mo,bool jetveto,int & extra_pdf) 
{  
  m_last_veto=0;

  if ((1.-m_z)*m_t>m_t0) {
    m_last_veto=1;
    return 1;
  }
  // 1. masses, z-range and splitting function
  if (!extra_pdf) {
    if (MassVeto(0)) {
      m_last_veto=2;
      return 1;
    }
  }
  // 2. alphaS
  if (CplVeto()) {
    m_last_veto=3;
    return 1;
  }
  // 3. angular ordering
  if (PTVeto(mo)) {
    m_last_veto=4;
    return 1;
  }
  // 4. jet veto
  if (jetveto) {      // test only!!!
    if (JetVeto(mo)) {
      m_last_veto=5;
      return 1;
    }
   }

  // 5. extra pdf weight for first branch below  Q_jet
  if (extra_pdf) {
    if (MassVeto(1)) {
      m_last_veto=6;
      extra_pdf=0;
      return 1;
    }
    extra_pdf=0;
  }
  // 6. energy conservation for multiple iteractions
  if (m_miveto_scheme) {
    if (MIVeto(mo->part->Flav(),mo->x,mo->t)) {
      // std::cout<<"mi veto caught ("<<this<<") "<<rem_weight<<" "<<rem_ran<<std::endl;
      m_last_veto=7;
      return 1;
    }
  }
  return 0;
}

bool Spacelike_Sudakov::MIVeto(const ATOOLS::Flavour &flavour,
			       const double x,const double scale) 
{
  p_pdf->Calculate(x,0,0,scale);
  // std::cout<<"mi veto check ("<<this<<") "<<p_pdf->GetXPDF(flavour)<<std::endl;
  if (p_pdf->GetXPDF(flavour)==0.) return true;
  return false;
}

bool Spacelike_Sudakov::MassVeto(int extra_pdf) 
{
  double weight  = p_pdf->GetXPDF(GetFlB())/(p_pdf->GetXPDF(GetFlA())*m_pdf_fac); 

  double q = sqrt(-m_t);
  double firstq = sqrt(m_facscale);
  double wb_jet;
  switch (m_pdf_scheme) {
  case 0:
    firstq/=sqrt(1.-m_z);
    break;
  default:
    q = sqrt(m_pt2);
  }

  q *= m_pdf_scale_fac;


  if (extra_pdf) {
    firstq*=m_pdf_scale_fac;
    p_pdf->Calculate(m_x,0.,0.,firstq);
    wb_jet   = p_pdf->GetXPDF(GetFlB());
  }
  p_pdf->Calculate(m_x,0.,0.,q);
  if (!extra_pdf) {
    wb_jet   = p_pdf->GetXPDF(GetFlB());
  }
  if (m_x/m_z>=1.) {
    std::cout<<"x="<<m_x<<"    z="<<m_z<<"  x/z="<<m_x/m_z<<std::endl;
  }
  p_pdfa->Calculate(m_x/m_z,0.,0.,q);
  weight        *= p_pdfa->GetXPDF(GetFlA())/wb_jet;
  weight        *= GetWeight(m_z,-m_t,0);

  if (ran.Get() > weight) return 1;
  return 0;
}

bool Spacelike_Sudakov::CplVeto() 
{
  switch (m_cpl_scheme) {
  case 0 : 
    return 0;
  case 2 : 
    return (GetCoupling(0.25*m_t)/GetCoupling() > ran.Get()) ? 0 : 1;   
  case 3 : {
    double a3 = GetCoupling();
    double b3 = GetCoupling(0.25*m_pt2);
//     std::cout<<" coupling "<<a3<<"  bigger? than "<<b3;
//     if (b3>a3) std::cout<<"IS!!!"<<std::endl;
//     else  std::cout<<std::endl;
    double r3 = ran.Get();
    double w3 = b3/a3;
    return (w3 > r3) ? 0 : 1;
  }
  default : 
    double a = GetCoupling();
    double b = GetCoupling(m_pt2);
    double r = ran.Get();
    double w = b/a;
    return (w > r) ? 0 : 1;
    return (GetCoupling(m_pt2)/GetCoupling() > ran.Get()) ? 0 : 1;   
  }
}

bool Spacelike_Sudakov::PTVeto(Knot * mo) 
{
  double th = 4.*m_z*m_z*m_t/(4.*m_z*m_z*m_t-(1.-m_z)*m_x*m_x*m_pt2max);
  if (!m_inflav.Strong()) {
    mo->maxpt2 = m_pt2;
    return 0;
  }

  switch (m_ordering_scheme) {
  case 0 : 
    mo->maxpt2 = m_pt2;
    return 0;
  case 2 : 
    if (th > mo->thcrit) return 1;
    mo->maxpt2 = m_pt2;
    return 0;
  default :
    if (m_pt2 > mo->maxpt2) return 1;
    mo->thcrit = th;
    mo->maxpt2 = m_pt2;
    return 0;
  }
}

bool Spacelike_Sudakov::JetVeto(Knot * mo) 
{
  std::cout<<" Spacelike_Sudakov::JetVeto("<<mo->kn_no<<") called "<<m_pt2<<" vs. "<<m_qjet;
  if (m_pt2>m_qjet) {
    std::cout<<"would be veto "<<std::endl;
    return 0;  // test only!!!!

    return 1;
  }
  std::cout<<std::endl;
  return 0;
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
  SplFunIter iter(m_group);
  for (;iter();++iter)
    if (iter()->GetFlB()==m_inflav) { p_selected=iter(); break; }
  if (!iter()) return m_lastint = -1.;
  return m_lastint = p_selected->CrudeInt(_zmin,_zmax);
}     

void Spacelike_Sudakov::SetJetvetoPt2(const double pt2) 
{ 
  p_kin->SetJetvetoPt2(pt2); 
}

void Spacelike_Sudakov::SetFactorisationScale(const double scale)
{
  m_facscale=scale;
}

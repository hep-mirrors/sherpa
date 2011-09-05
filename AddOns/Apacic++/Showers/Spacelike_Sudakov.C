#include "AddOns/Apacic++/Showers/Spacelike_Sudakov.H"
#include "AddOns/Apacic++/Showers/Spacelike_Kinematics.H"
#include "AddOns/Apacic++/Showers/Timelike_Sudakov.H"
#include "AddOns/Apacic++/Showers/Sudakov_Tools.H"
#include "AddOns/Apacic++/Main/Knot.H"

#include "AddOns/Apacic++/Showers/QCD_Splitting_Functions.H"
#include "AddOns/Apacic++/Showers/QED_Splitting_Functions.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "PDF/Remnant/Remnant_Base.H"
#include "BEAM/Main/Beam_Base.H"

#include <iomanip>
#include <limits>

#ifdef PROFILE__all
#include "prof.hh"
#else 
#define PROFILE_HERE {}
#define PROFILE_LOCAL(LOCALNAME) {}
#endif


using namespace APACIC;
using namespace PDF;
using namespace ATOOLS;


size_t Spacelike_Sudakov::s_nem(0);
size_t Spacelike_Sudakov::s_maxem(std::numeric_limits<size_t>::max());

Spacelike_Sudakov::
Spacelike_Sudakov(PDF_Base * pdf,Sudakov_Tools * tools,
		  Spacelike_Kinematics * kin,
		  double pt2min,ATOOLS::Data_Reader * dataread,
		  const size_t beam) : 
  Backward_Splitting_Group(p_rms,0,0), p_rms(NULL), p_tools(tools), p_kin(kin), 
  p_pdfa(pdf->GetBasicPDF()->GetCopy()),
  m_pt2min(dabs(pt2min)), 
  m_t0(-m_pt2min),
  m_pdf_fac(5.), m_rbmax(1.0)
{
  p_pdf = pdf->GetBasicPDF()->GetCopy();
}

Spacelike_Sudakov::~Spacelike_Sudakov() 
{
  delete p_pdf;
}

bool Spacelike_Sudakov::Initialize()
{
  p_tools->CalculateMaxCouplings
    (m_cpl_scheme,m_pt2min*m_cpl_scale_fac,m_s_hadron*m_cpl_scale_fac);

  if (p_pdf->Bunch().IsHadron()) {
    for (int i=1;i<6;++i) {
      Flavour fl = Flavour((kf_code)(i));
      Add(new q_qg(p_rms,fl,p_tools));
      Add(new q_qg(p_rms,fl.Bar(),p_tools));
      Add(new q_gq(p_rms,fl,p_tools));
      Add(new q_gq(p_rms,fl.Bar(),p_tools));
      Add(new g_qq(p_rms,fl,p_tools));
      Add(new g_qq(p_rms,fl.Bar(),p_tools));
    }
    Add(new g_gg(p_rms,p_tools));
    Add(new g_gg(p_rms,p_tools));
  }

  if (p_pdf->Bunch().IsLepton()) {
    Add(new f_fp(p_rms,Flavour(kf_e),p_tools));
    Add(new f_fp(p_rms,Flavour(kf_e).Bar(),p_tools));
    Add(new f_pf(p_rms,Flavour(kf_e),p_tools));
    Add(new f_pf(p_rms,Flavour(kf_e).Bar(),p_tools));
    Add(new p_ff(p_rms,Flavour(kf_e),p_tools));
    Add(new p_ff(p_rms,Flavour(kf_e).Bar(),p_tools));
  }

  PrintStat();
  return true;
}

void Spacelike_Sudakov::AcceptBranch(const Knot *const mo) 
{
  msg_Debugging()<<METHOD<<"("<<mo->kn_no<<"):"<<std::endl;
  if (mo->prev==NULL) return;
  Knot *mother(mo->prev), *sister(mo->prev->left);
  if (!mother->part->Flav().Strong()) {
    mother->maxpt2 = sister->maxpt2 = mo->maxpt2;
    mother->thcrit = sister->thcrit = mo->thcrit;
    return;
  }
  mother->maxpt2 = sister->maxpt2 = Min(mo->smaxpt2,mo->maxpt2);
  mother->thcrit = sister->thcrit = mo->sthcrit;
  msg_Debugging()<<"  accept mother = "<<mother->kn_no
		 <<", set maxpt = "<<sqrt(mother->maxpt2)
		 <<", thcrit = "<<mother->thcrit<<"\n";
  msg_Debugging()<<"  accept sister = "<<sister->kn_no
		 <<", set maxpt = "<<sqrt(sister->maxpt2)
		 <<", thcrit = "<<sister->thcrit<<"\n";
}

bool Spacelike_Sudakov::Dice(Knot * mo,double sprime,
			     double m2p,int & extra_pdf) 
{
  if (s_nem>=s_maxem) {
    mo->t    = mo->tout;
    mo->stat = 0;
    return false; 
  }
  PROFILE_HERE;
  mo->tmax = mo->t;  // last start t
  m_inflav = mo->part->Flav(); 
  m_t      = mo->t;
  m_x      = mo->x;

  if (!((m_t-m_t0)<rpa->gen.Accu())) {
    if (mo->prev) {
      mo->t    = mo->tout;
      mo->stat = 0;
      return false; 
    }
  }
  
  double xe(2.*m_emin*sqrt(sprime)/sqr(rpa->gen.Ecms()));  
  m_zmin = m_x/(1.-xe);
  m_zmax = m_x/(m_x+xe);
  if (m_zmin>m_zmax) {
    mo->t    = mo->tout;
    mo->stat = 0;
    return false; 
  }
  
  double q2pdf(mo->pt2lcm);
  if (!extra_pdf) {
    switch (m_pdf_scheme) {
    case 0: q2pdf=Max(-m_t,-m_t0); break;
    default: q2pdf=Max(mo->pt2lcm+mo->tout,-m_t0);
    }
    q2pdf*=m_pdf_scale_fac;
  }
  msg_Debugging()<<METHOD<<"(): s' = "<<sprime<<", Q = "
		 <<sqrt(q2pdf)<<", x = "<<m_x<<", "
		 <<m_zmin<<" < z < "<<m_zmax<<"\n";
  if (q2pdf<p_pdf->Q2Min()) {
    mo->t    = mo->tout;
    mo->stat = 0;
    return false;
  }
  p_pdf->Calculate(m_x,q2pdf);

  CrudeInt(m_zmin,m_zmax);

  while (m_t<m_t0) {
    ProduceT();
    if (mo->right && m_t<mo->right->t) continue;
    if (m_t>m_t0) {
      mo->t    = mo->tout;
      mo->stat = 0;
      return false;     
    }
    SelectOne();
    m_z   = GetZ();
    m_ta  = sqr(p_ms->Mass(GetFlA())); 
    m_tc  = sqr(p_ms->Mass(GetFlC()));
    m_pt2 = m_z*(1.-m_z)*m_ta-(1.-m_z)*m_t-m_z*m_tc;
    if (m_pt2+m_ta<m_pt2min) continue;

    double uhat(m_ta+m_tc+m2p-m_t-sprime*(1.-m_z)/m_z);

    if (uhat<0. && !Veto(mo)) {
      if (!RemnantVeto(mo)) {
	UniformPhi();
	mo->z      = m_z;
	mo->t      = m_t;
	mo->phi    = m_phi;
	return true;
      }
    }
  }
  mo->t    = mo->tout;
  mo->stat = 0;
  return false; 
}

void Spacelike_Sudakov::ProduceT() {
  if (m_lastint <0.) m_t = +1.;            
  else {
    double rn(ran.Get()), ne(2.*M_PI*log(rn)/(m_lastint*m_pdf_fac));
    if (IsZero(ne)) m_t=m_t0;
    else m_t*=exp(ne);
  }
  return;
}

bool Spacelike_Sudakov::Veto(Knot * mo) 
{  
  PROFILE_HERE;
  // 1. masses, z-range and splitting function
  if (MassVeto()) return true;

  // 2. alphaS
  if (CplVeto())                 return true;

  // 3. ordering
  if (OrderingVeto(mo))          return true;

  return false;
}

bool Spacelike_Sudakov::MassVeto() 
{
  PROFILE_HERE;
  // cancel weight from CrudeInt
  double weight(p_pdf->GetXPDF(GetFlB())/
		(p_pdf->GetXPDF(GetFlA())*m_pdf_fac)); 
  double q2(0.0);
  switch (m_pdf_scheme) {
  case 0: q2=-m_t; break;
  default: q2=m_pt2+m_ta;
  }
  q2*=m_pdf_scale_fac;
  if (q2<p_pdf->Q2Min() || m_x>m_z) return true;
  p_pdfa->Calculate(m_x/m_z,q2);
  double test(p_pdfa->GetXPDF(GetFlA()));
  if (IsZero(test)) return true;
  p_pdfa->Calculate(m_x,q2);
  // weight: P(z)*(x/z*f(x/z,t))/(x*f(x,t));
  weight*=test/p_pdfa->GetXPDF(GetFlB())*GetWeight(m_z,-m_t,0);
  if (ran.Get()>weight) return true;
  return false;
}

bool Spacelike_Sudakov::CplVeto() 
{
  switch (m_cpl_scheme) {
  case 0 :  return false;
  case 2 :  return GetCoupling(m_cpl_scale_fac*m_t)/GetCoupling()<ran.Get();   
  default : return GetCoupling(m_cpl_scale_fac*m_pt2)/GetCoupling()<ran.Get();   
  }
  return true;
}

bool Spacelike_Sudakov::OrderingVeto(Knot * mo) 
{
  msg_Debugging()<<METHOD<<"("<<mo->kn_no<<")\n";
  msg_Indent();
  double th(4.*m_z*m_z*m_t/(4.*m_z*m_z*m_t-(1.-m_z)*m_x*m_x*m_s_hadron));
  mo->sthcrit = th; //asin(sqrt(th));
  if (mo->sthcrit<0.0) mo->sthcrit=M_PI-mo->sthcrit;
  msg_Debugging()<<"ss: thcrit = "<<mo->thcrit<<", th = "<<th<<std::endl;
  msg_Debugging()<<"ss: maxpt2 = "<<mo->maxpt2<<", pt2 = "<<m_pt2<<std::endl;
  if (!m_inflav.Strong()) return false;
  switch (m_ordering_scheme) {
  case 0 :
    return false;
  case 2 : 
    if (th > mo->thcrit) {
      msg_Debugging()<<"ss: th veto\n";
      return true;
    }
    return false;
  default :
    return false;
  }
  return false;
}

bool Spacelike_Sudakov::PTVeto(Knot *knot)
{
  if (!knot->part->Flav().Strong()) {
    knot->smaxpt2=knot->right->maxpt2;
    msg_Debugging()<<"ss: ew knot "<<knot->kn_no<<", maxpt = "
		   <<sqrt(knot->right->maxpt2)<<std::endl;
    return false;
  }
  double E2(sqr(sqrt(knot->left->E2)+sqrt(knot->right->E2)));
  double z(p_kin->Kinematics()->
	   LightConeZ(knot->right->z,E2,knot->t,
		      knot->right->t,knot->left->t));
  double pt2(z*(1.0-z)*knot->t-(1.0-z)*knot->right->t-z*knot->left->t);
  if (knot->part->Flav().IsGluon() && knot->right->part->Flav().IsGluon()) 
    if (ran.Get()<0.5) pt2=knot->right->maxpt2;
  knot->smaxpt2=Max(0.0,pt2);
  msg_Debugging()<<"ss: knot "<<knot->kn_no<<", maxpt = "
		 <<sqrt(knot->maxpt2)<<", pt = "<<sqrt(dabs(pt2))
		 <<"("<<pt2<<")"<<std::endl;
  if (m_ordering_scheme==1 &&
      knot->part->Info()!='H' && 
      (pt2<0.0 || pt2>knot->maxpt2)) {
    msg_Debugging()<<"ss: pt veto\n";
    return true;
  }
  return false;
}

bool Spacelike_Sudakov::RemnantVeto(Knot * mo) 
{
  double E(p_remnant->GetBeam()->Energy()*mo->x/m_z);
  return !p_remnant->TestExtract
    (GetFlA(),Vec4D(E,0.0,0.0,sqrt(E*E-sqr(p_ms->Mass(GetFlA())))));
}

void Spacelike_Sudakov::Add(Splitting_Function * spl) 
{
  for (Splitting_Vector::iterator sit(m_splittings.begin());
       sit!=m_splittings.end();++sit) {
    if ((*sit)->GetFlB()==spl->GetFlB()) {
      (*sit)->Add(spl);
      return;
    }
  }
  m_splittings.push_back(new Backward_Splitting_Group(p_rms,spl,p_pdf));
  p_selected=spl;
}


double Spacelike_Sudakov::CrudeInt(double zmin,double zmax) 
{
  PROFILE_HERE;
  Splitting_Vector::iterator sit(m_splittings.begin());
  for (;sit!=m_splittings.end();++sit)
    if ((*sit)->GetFlB()==m_inflav) { 
      p_selected=*sit; 
      break; 
    }
  if (sit==m_splittings.end()) { 
    p_selected=NULL; 
    return m_lastint=-1.0; 
  }
  return m_lastint = p_selected->CrudeInt(zmin,zmax)*m_rbmax;
}     

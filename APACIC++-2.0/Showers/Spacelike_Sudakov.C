#include "Spacelike_Sudakov.H"
#include "Timelike_Sudakov.H"
#include "Sudakov_Tools.H"
#include "Knot.H"

#include "QCD_Splitting_Functions.H"
#include "QED_Splitting_Functions.H"
#include "Run_Parameter.H"


using namespace APACIC;
using namespace PDF;
using namespace AMATOOLS;
using namespace APHYTOOLS;
using namespace AORGTOOLS;

Spacelike_Sudakov::Spacelike_Sudakov(PDF_Base * _pdf,Sudakov_Tools * _tools) : 
  p_tools(_tools), Backward_Splitting_Group(0,0) 
{
  p_pdf    = _pdf;
  p_pdfa   = p_pdf->GetCopy();

  m_cpl_scheme      = 1;  /*  (0=fix, 1=pt^2 2=t/4)                             */   
  m_ordering_scheme = 0;  /* Switch for ordering due to coherence:  
                           0 = none, 1 = pt^2, 2 = pt^2/E^2                   */

  m_emin        = .5;

  m_pt2max      = sqr(rpa.gen.Ecms());
  m_pt2min      = rpa.pshower.InitialQ02();
  m_lambda2     = p_tools->GetLambda2(); // from as(pt2max) = 1/(beta_0 * log(pt2max/lambda^2))
  m_b           = p_tools->GetBnorm();   //      2*pi*beta_0*scalefactor

  msg.Events()<<"Init QCD splitting functions ..."<<std::endl
		 <<"   ("<<rpa.gen.Beam1()<<","<<rpa.gen.Beam1().IsHadron()<<","
		 <<rpa.gen.Beam2()<<","<<rpa.gen.Beam2().IsHadron()<<")"<<std::endl;
  if ( (rpa.gen.Beam1().IsHadron()) || (rpa.gen.Beam2().IsHadron())) {
    // -- initialise QCD Splittingfunctions --
    for (int i=1;i<6;++i) {
      msg.Debugging()<<"   ... : "<<Flavour(i)<<" -> "<<Flavour(i)<<std::endl;
      // add gluon quark & quark gluon (loop over active Flavours)
      Add(new q_qg(Flavour(kf::code(i)),p_tools));
      Add(new q_qg(Flavour(i).Bar(),p_tools));
      Add(new q_gq(Flavour(i),p_tools));
      Add(new q_gq(Flavour(i).Bar(),p_tools));
      // add q qbar & qbar q (loop over active Flavours)
      Add(new g_qq(Flavour(i),p_tools));
      Add(new g_qq(Flavour(i).Bar(),p_tools));
    }
    // add gluon gluon twice!
    Add(new g_gg(p_tools));
    Add(new g_gg(p_tools));
  }

  /*
  if ( (rpa.gen.Beam1().IsLepton()) || (rpa.gen.Beam2().IsLepton())) {
    for (int i=7;i<13;i+=2) {
      // add gluon quark & quark gluon (loop over active Flavours)
      Add(new f_fp(Flavour(i),p_tools));
      Add(new f_fp(Flavour(i).Bar(),p_tools));
      Add(new f_pf(Flavour(i),p_tools));
      Add(new f_pf(Flavour(i).Bar(),p_tools));
      // add q qbar & qbar q (loop over active Flavours)
      Add(new p_ff(Flavour(i),p_tools));
      Add(new p_ff(Flavour(i).Bar(),p_tools));
    }
  }
  */

  if ( (rpa.gen.Beam1().IsLepton()) || (rpa.gen.Beam2().IsLepton())) {
    Add(new f_fp(Flavour(kf::e),p_tools));                
    Add(new f_fp(Flavour(kf::e).Bar(),p_tools));
    Add(new f_pf(Flavour(kf::e),p_tools));
    Add(new f_pf(Flavour(kf::e).Bar(),p_tools));
    Add(new p_ff(Flavour(kf::e),p_tools));
    Add(new p_ff(Flavour(kf::e).Bar(),p_tools));
  }

  PrintStat();
    
}


bool Spacelike_Sudakov::Dice(Knot * mo,double sprime) {
  m_inflav = mo->part->Flav(); 
  m_t      = mo->t;
  m_x      = mo->x;
  m_t0     = m_pt2min;
  
  msg.Debugging()<<"Spacelike_Sudakov::Dice (t,x): "<<m_t<<" / "<<m_x<<" / for ("<<mo->kn_no
		 <<"), "<<m_inflav<<std::endl;
  
  if (!((m_t-m_t0)<rpa.gen.Accu())) {
    msg.Debugging()<<"Spacelike_Sudakov::Dice : mo can't branch (t > t_0) : "
		   <<m_t<<" < "<<m_t0<<std::endl
		   <<"      mo = "<<m_inflav<<", mass = "<<m_inflav.Mass()<<", "
		   <<"status = "<<mo->stat<<std::endl;
    if (mo->prev) {
      msg.Debugging()<<"      prev = "<<mo->prev->part->Flav()
		     <<", stat = "<<mo->prev->stat<<", t = "<<mo->prev->t<<std::endl;
      mo->t    = mo->tout;
      mo->stat = 0;
      return 0; 
    }
  }
  
  // regulator of parton splitting functions
  if (m_inflav.Strong()) m_xe   = 2.*m_emin*sprime/sqr(2.*rpa.gen.Ecms());
  else                 m_xe   = 0.0001;
  m_zmin = m_x/(1.-m_xe);
  m_zmax = m_x/(m_x+m_xe);
  msg.Debugging()<<"Spacelike_Sudakov::Dice : zrange "<<m_zmin<<" < "<<m_zmax<<std::endl;
  // *AS*
  /*
    xe   = 2.*emin*sqrt(sprime)/sqr(rpa.gen.Ecms());  
    zmin = x/(1.-xe);
    zmax = x/(x+xe);
    msg.Out()<<"Spacelike_Sudakov::Dice : zrange "<<zmin<<" < "<<zmax<<std::endl;
  */

  if (m_zmin>m_zmax) {
    msg.Debugging()<<"Spacelike_Sudakov::Dice : mother can't branch (zmax<zmin) : "
		   <<m_x<<" < "<<m_t0<<std::endl
		   <<"      mother = "<<m_inflav<<", mass = "<<m_inflav.Mass()<<", "
		   <<"status = "<<mo->stat<<std::endl;
    if (mo->prev) 
      msg.Debugging()<<"      prev = "<<mo->prev->part->Flav()
		     <<", stat = "<<mo->prev->stat<<", t = "<<mo->prev->t<<std::endl;
    
    mo->t    = mo->tout;
    mo->stat = 0;
    return 0; 
  }
  
  
  p_pdf->Calculate(m_x,sqrt(-m_t));
  
  while (m_t<m_t0) {
    CrudeInt(m_zmin,m_zmax);  // using above pdf !!! and inflav
    ProduceT();
    if (m_t>m_t0) {
      msg.Debugging()<<"Spacelike_Sudakov::No Branch for ("<<mo->kn_no<<"), "<<m_inflav
		     <<", "<<m_t<<", set on t="<<m_t0<<std::endl;
      mo->t    = mo->tout;
      mo->stat = 0;
      return 0;      // no further branching
    }
    SelectOne();
    m_z   = GetZ();   
    m_pt2 = -(1.-m_z)*m_t;
    if (!Veto(mo)) {
      msg.Debugging()<<"Spacelike_Sudakov::Dice Branch with t="
		     <<m_t<<", z="<<m_z<<", "<<m_inflav<<" for "<<m_lastint<<std::endl;
      UniformPhi();
      mo->z      = m_z;
      mo->t      = m_t;
      mo->phi    = m_phi;
      return 1;
    }    
  }
  msg.Debugging()<<"Spacelike_Sudakov::Banged out of Dice !"<<std::endl;
  mo->t    = mo->tout;
  mo->stat = 0;
  return 0; 
}

void Spacelike_Sudakov::ProduceT() {
  if (m_lastint <rpa.gen.Accu()) m_t = m_t0;
  m_t *= exp( 2.*M_PI*log(ran.Get()) / m_lastint );
  return;
}

bool Spacelike_Sudakov::Veto(Knot * mo) 
{  
  
  // "lower" cutoff reached / still spacelike
  if ((1.-m_z)*m_t>m_t0)   return 1;
  // 1. masses, z-range and splitting function
  if (MassVeto())    return 1;
  // 2. alphaS
  if (CplVeto())  return 1;
  // 3. angular ordering
  if (PTVeto(mo)) return 1;
  // passed vetos
  return 0;
}



bool Spacelike_Sudakov::MassVeto() 
{
  // same pt2 of actual branch but with the different values x and x/z.
  double weight  = p_pdf->GetXPDF(GetFlB())/p_pdf->GetXPDF(GetFlA()); 
  p_pdf->Calculate(m_x,sqrt(-m_t));
  p_pdfa->Calculate(m_x/m_z,sqrt(-m_t));
  weight        *= GetWeight(m_z,-m_t,0) * p_pdfa->GetXPDF(GetFlA())/p_pdf->GetXPDF(GetFlB()); 
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
  default : 
    return (GetCoupling(m_pt2)/GetCoupling() > ran.Get()) ? 0 : 1;   
  }
}

bool Spacelike_Sudakov::PTVeto(Knot * mo) 
{
  // approximately pt^2/pl^2 with virtual masses neglected. Check this !

  double th = 4.*m_z*m_z*m_t/(4.*m_z*m_z*m_t-(1.-m_z)*m_x*m_x*m_pt2max);
  if (!m_inflav.Strong()) {
    mo->thcrit = th;
    mo->maxpt2 = m_pt2;
    return 0;
  }

  switch (m_ordering_scheme) {
  case 0 : 
    mo->thcrit = th;
    mo->maxpt2 = m_pt2;
    return 0;
  case 2 : 
    if (th > mo->thcrit) return 1;
    mo->thcrit = th;
    mo->maxpt2 = m_pt2;
    return 0;
  default :
    if (m_pt2 > mo->maxpt2) return 1;
    mo->thcrit = th;
    mo->maxpt2 = m_pt2;
    return 0;
  }
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
  if (!iter()) std::cout<<"warning "<<m_inflav<<" not implemented in Splitting_Group"<<std::endl;
  return m_lastint = p_selected->CrudeInt(_zmin,_zmax);
}     

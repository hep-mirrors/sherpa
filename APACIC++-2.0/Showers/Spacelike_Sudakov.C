#include "Spacelike_Sudakov.H"
#include "QCD_Splitting_Functions.H"
#include "QED_Splitting_Functions.H"
#include "Run_Parameter.H"

using namespace APACIC;
using namespace PDF;
using namespace AMATOOLS;
using namespace APHYTOOLS;
using namespace AORGTOOLS;

Spacelike_Sudakov::Spacelike_Sudakov(PDF_Base * _pdf,Sudakov_Tools * _tools) : 
  tools(_tools), Backward_Splitting_Group(0,0) {
  pdf    = _pdf;
  pdfa   = pdf->GetCopy();

  cpl_scheme      = 1;  /*  (0=fix, 1=pt^2 2=t/4)                             */   
  ordering_scheme = 0;  /* Switch for ordering due to coherence:  
                           0 = none, 1 = pt^2, 2 = pt^2/E^2                   */

  emin        = .5;

  pt2max      = sqr(rpa.gen.Ecms());
  pt2min      = rpa.pshower.InitialQ02();
  lambda2     = tools->GetLambda2(); // from as(pt2max) = 1/(beta_0 * log(pt2max/lambda^2))
  b           = tools->GetBnorm();   //      2*pi*beta_0*scalefactor

  msg.Debugging()<<"Init QCD splitting functions ..."<<std::endl
		 <<"   ("<<rpa.gen.Beam1()<<","<<rpa.gen.Beam1().ishadron()<<","
		 <<rpa.gen.Beam2()<<","<<rpa.gen.Beam2().ishadron()<<")"<<std::endl;
  if ( (rpa.gen.Beam1().ishadron()) || (rpa.gen.Beam2().ishadron())) {
    // -- initialise QCD Splittingfunctions --
    for (int i=1;i<6;++i) {
      msg.Debugging()<<"   ... : "<<Flavour(i)<<" -> "<<Flavour(i)<<std::endl;
      // add gluon quark & quark gluon (loop over active Flavours)
      Add(new q_qg(Flavour(kf::code(i)),tools));
      Add(new q_qg(Flavour(i).bar(),tools));
      Add(new q_gq(Flavour(i),tools));
      Add(new q_gq(Flavour(i).bar(),tools));
      // add q qbar & qbar q (loop over active Flavours)
      Add(new g_qq(Flavour(i),tools));
      Add(new g_qq(Flavour(i).bar(),tools));
    }
    // add gluon gluon twice!
    Add(new g_gg(tools));
    Add(new g_gg(tools));
  }

  /*
  if ( (rpa.gen.Beam1().islepton()) || (rpa.gen.Beam2().islepton())) {
    for (int i=7;i<13;i+=2) {
      // add gluon quark & quark gluon (loop over active Flavours)
      Add(new f_fp(Flavour(i),tools));
      Add(new f_fp(Flavour(i).bar(),tools));
      Add(new f_pf(Flavour(i),tools));
      Add(new f_pf(Flavour(i).bar(),tools));
      // add q qbar & qbar q (loop over active Flavours)
      Add(new p_ff(Flavour(i),tools));
      Add(new p_ff(Flavour(i).bar(),tools));
    }
  }
  */

  if ( (rpa.gen.Beam1().islepton()) || (rpa.gen.Beam2().islepton())) {
    Add(new f_fp(Flavour(kf::e),tools));                
    Add(new f_fp(Flavour(kf::e).bar(),tools));
    Add(new f_pf(Flavour(kf::e),tools));
    Add(new f_pf(Flavour(kf::e).bar(),tools));
    Add(new p_ff(Flavour(kf::e),tools));
    Add(new p_ff(Flavour(kf::e).bar(),tools));
  }

  PrintStat();
    
};


bool Spacelike_Sudakov::Dice(Knot * mother,double sprime) {
  inflav = mother->part->flav(); 
  t      = mother->t;
  x      = mother->x;
  t0     = pt2min;
  
  msg.Debugging()<<"Spacelike_Sudakov::Dice (t,x): "<<t<<" / "<<x<<" / for ("<<mother->kn_no
		 <<"), "<<inflav<<std::endl;
  
  if (!((t-t0)<rpa.gen.Accu())) {
    msg.Debugging()<<"Spacelike_Sudakov::Dice : mother can't branch (t > t_0) : "
		   <<t<<" < "<<t0<<std::endl
		   <<"      mother = "<<inflav<<", mass = "<<inflav.mass()<<", "
		   <<"status = "<<mother->stat<<std::endl;
    if (mother->prev) {
      msg.Debugging()<<"      prev = "<<mother->prev->part->flav()
		     <<", stat = "<<mother->prev->stat<<", t = "<<mother->prev->t<<std::endl;
      mother->t    = mother->tout;
      mother->stat = 0;
      return 0; 
    }
  }
  
  // regulator of parton splitting functions
  if (inflav.strong()) xe   = 2.*emin*sprime/sqr(2.*rpa.gen.Ecms());
  else                 xe   = 0.0001;
  zmin = x/(1.-xe);
  zmax = x/(x+xe);
  msg.Debugging()<<"Spacelike_Sudakov::Dice : zrange "<<zmin<<" < "<<zmax<<std::endl;
  
  if (zmin>zmax) {
    msg.Debugging()<<"Spacelike_Sudakov::Dice : mother can't branch (zmax<zmin) : "
		   <<x<<" < "<<t0<<std::endl
		   <<"      mother = "<<inflav<<", mass = "<<inflav.mass()<<", "
		   <<"status = "<<mother->stat<<std::endl;
    if (mother->prev) 
      msg.Debugging()<<"      prev = "<<mother->prev->part->flav()
		     <<", stat = "<<mother->prev->stat<<", t = "<<mother->prev->t<<std::endl;
    
    mother->t    = mother->tout;
    mother->stat = 0;
    return 0; 
  }
  
  
  pdf->Calculate(x,sqrt(-t));
  
  while (t<t0) {
    CrudeInt(zmin,zmax);  // using above pdf !!! and inflav
    ProduceT();
    if (t>t0) {
      msg.Debugging()<<"Spacelike_Sudakov::No Branch for ("<<mother->kn_no<<"), "<<inflav
		     <<", "<<t<<", set on t="<<t0<<std::endl;
      mother->t    = mother->tout;
      mother->stat = 0;
      return 0;      // no further branching
    }
    SelectOne();
    z   = GetZ();   
    pt2 = -(1.-z)*t;
    if (!Veto(mother)) {
      msg.Debugging()<<"Spacelike_Sudakov::Dice Branch with t="
		     <<t<<", z="<<z<<", "<<inflav<<" for "<<lastint<<std::endl;
      UniformPhi();
      mother->z      = z;
      mother->t      = t;
      mother->phi    = phi;
      return 1;
    }    
  }
  msg.Debugging()<<"Spacelike_Sudakov::Banged out of Dice !"<<std::endl;
  mother->t    = mother->tout;
  mother->stat = 0;
  return 0; 
}

void Spacelike_Sudakov::ProduceT() {
  if (lastint <rpa.gen.Accu()) t = t0;
  t *= exp( 2.*M_PI*log(Ran.get()) / lastint );
  return;
}

bool Spacelike_Sudakov::Veto(Knot * mo) 
{  
  msg.Debugging()<<"      Enter the vetos with t, z, pt2 = "
		 <<t<<", "<<z<<", "<<-(1.-z)*t<<std::endl;
  
  // "lower" cutoff reached / still spacelike
  if ((1.-z)*t>t0)   return 1;

  // 1. masses, z-range and splitting function
  if (MassVeto())    return 1;
  // 2. alphaS
  if (CplVeto())  return 1;
  // 3. angular ordering
  if (PTVeto(mo)) return 1;
  msg.Debugging()<<"            Passed Vetos "<<std::endl;
  return 0;
}



bool Spacelike_Sudakov::MassVeto() 
{
  // same pt2 of actual branch but with the different values x and x/z.
  double weight  = pdf->GetXPDF(GetFlB())/pdf->GetXPDF(GetFlA()); 
  pdf->Calculate(x,sqrt(-t));
  pdfa->Calculate(x/z,sqrt(-t));
  weight        *= GetWeight(z,-t,0) * pdfa->GetXPDF(GetFlA())/pdf->GetXPDF(GetFlB()); 
  msg.Debugging()<<"            MassVeto "<<weight<<":"<<x<<","<<x/z<<","<<t<<std::endl;
  if (Ran.get() > weight) return 1;
  return 0;
}

bool Spacelike_Sudakov::CplVeto() {
  msg.Debugging()<<"            CplVeto "<<std::endl;
  switch (cpl_scheme) {
  case 0 : 
    return 0;
  case 2 : 
    return (GetCoupling(0.25*t)/GetCoupling() > Ran.get()) ? 0 : 1;   
  default : 
    msg.Debugging()<<"            AlphaSVeto "<<GetCoupling(pt2)<<" / "<<GetCoupling()<<std::endl;
    return (GetCoupling(pt2)/GetCoupling() > Ran.get()) ? 0 : 1;   
  }
}

bool Spacelike_Sudakov::PTVeto(Knot * mo) {
  msg.Debugging()<<"            PTVeto "<<pt2<<std::endl;

  // approximately pt^2/pl^2 with virtual masses neglected. Check this !
  double th = 4.*z*z*t/(4.*z*z*t-(1.-z)*x*x*pt2max);
  if (!inflav.strong()) {
    mo->thcrit = th;
    mo->maxpt2 = pt2;
    return 0;
  }

  switch (ordering_scheme) {
  case 0 : 
    mo->thcrit = th;
    mo->maxpt2 = pt2;
    return 0;
  case 2 : 
    if (th > mo->thcrit) return 1;
    mo->thcrit = th;
    mo->maxpt2 = pt2;
    return 0;
  default :
    if (pt2 > mo->maxpt2) return 1;
    mo->thcrit = th;
    mo->maxpt2 = pt2;
    return 0;
  }
}


#include "Spacelike_Sudakov.H"
#include "Timelike_Sudakov.H"
#include "Sudakov_Tools.H"

#include "QCD_Splitting_Functions.H"
#include "QED_Splitting_Functions.H"
#include "Run_Parameter.H"


using namespace APACIC;
using namespace PDF;
using namespace AMATOOLS;
using namespace APHYTOOLS;
using namespace AORGTOOLS;

Spacelike_Sudakov::Spacelike_Sudakov(PDF_Base * _pdf,Sudakov_Tools * _tools) : 
  tools(_tools), Backward_Splitting_Group(0,0) 
{
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

  msg.Events()<<"Init QCD splitting functions ..."<<std::endl
		 <<"   ("<<rpa.gen.Beam1()<<","<<rpa.gen.Beam1().IsHadron()<<","
		 <<rpa.gen.Beam2()<<","<<rpa.gen.Beam2().IsHadron()<<")"<<std::endl;
  if ( (rpa.gen.Beam1().IsHadron()) || (rpa.gen.Beam2().IsHadron())) {
    // -- initialise QCD Splittingfunctions --
    for (int i=1;i<6;++i) {
      msg.Debugging()<<"   ... : "<<Flavour(i)<<" -> "<<Flavour(i)<<std::endl;
      // add gluon quark & quark gluon (loop over active Flavours)
      Add(new q_qg(Flavour(kf::code(i)),tools));
      Add(new q_qg(Flavour(i).Bar(),tools));
      Add(new q_gq(Flavour(i),tools));
      Add(new q_gq(Flavour(i).Bar(),tools));
      // add q qbar & qbar q (loop over active Flavours)
      Add(new g_qq(Flavour(i),tools));
      Add(new g_qq(Flavour(i).Bar(),tools));
    }
    // add gluon gluon twice!
    Add(new g_gg(tools));
    Add(new g_gg(tools));
  }

  /*
  if ( (rpa.gen.Beam1().IsLepton()) || (rpa.gen.Beam2().IsLepton())) {
    for (int i=7;i<13;i+=2) {
      // add gluon quark & quark gluon (loop over active Flavours)
      Add(new f_fp(Flavour(i),tools));
      Add(new f_fp(Flavour(i).Bar(),tools));
      Add(new f_pf(Flavour(i),tools));
      Add(new f_pf(Flavour(i).Bar(),tools));
      // add q qbar & qbar q (loop over active Flavours)
      Add(new p_ff(Flavour(i),tools));
      Add(new p_ff(Flavour(i).Bar(),tools));
    }
  }
  */

  if ( (rpa.gen.Beam1().IsLepton()) || (rpa.gen.Beam2().IsLepton())) {
    Add(new f_fp(Flavour(kf::e),tools));                
    Add(new f_fp(Flavour(kf::e).Bar(),tools));
    Add(new f_pf(Flavour(kf::e),tools));
    Add(new f_pf(Flavour(kf::e).Bar(),tools));
    Add(new p_ff(Flavour(kf::e),tools));
    Add(new p_ff(Flavour(kf::e).Bar(),tools));
  }

  PrintStat();
    
}


bool Spacelike_Sudakov::Dice(Knot * mo,double sprime) {
  inflav = mo->part->Flav(); 
  t      = mo->t;
  x      = mo->x;
  t0     = pt2min;
  
  msg.Debugging()<<"Spacelike_Sudakov::Dice (t,x): "<<t<<" / "<<x<<" / for ("<<mo->kn_no
		 <<"), "<<inflav<<std::endl;
  
  if (!((t-t0)<rpa.gen.Accu())) {
    msg.Debugging()<<"Spacelike_Sudakov::Dice : mo can't branch (t > t_0) : "
		   <<t<<" < "<<t0<<std::endl
		   <<"      mo = "<<inflav<<", mass = "<<inflav.Mass()<<", "
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
  if (inflav.Strong()) xe   = 2.*emin*sprime/sqr(2.*rpa.gen.Ecms());
  else                 xe   = 0.0001;
  zmin = x/(1.-xe);
  zmax = x/(x+xe);
  msg.Debugging()<<"Spacelike_Sudakov::Dice : zrange "<<zmin<<" < "<<zmax<<std::endl;
  // *AS*
  /*
    xe   = 2.*emin*sqrt(sprime)/sqr(rpa.gen.Ecms());  
    zmin = x/(1.-xe);
    zmax = x/(x+xe);
    msg.Out()<<"Spacelike_Sudakov::Dice : zrange "<<zmin<<" < "<<zmax<<std::endl;
  */

  if (zmin>zmax) {
    msg.Debugging()<<"Spacelike_Sudakov::Dice : mother can't branch (zmax<zmin) : "
		   <<x<<" < "<<t0<<std::endl
		   <<"      mother = "<<inflav<<", mass = "<<inflav.Mass()<<", "
		   <<"status = "<<mo->stat<<std::endl;
    if (mo->prev) 
      msg.Debugging()<<"      prev = "<<mo->prev->part->Flav()
		     <<", stat = "<<mo->prev->stat<<", t = "<<mo->prev->t<<std::endl;
    
    mo->t    = mo->tout;
    mo->stat = 0;
    return 0; 
  }
  
  
  pdf->Calculate(x,sqrt(-t));
  
  while (t<t0) {
    CrudeInt(zmin,zmax);  // using above pdf !!! and inflav
    ProduceT();
    if (t>t0) {
      msg.Debugging()<<"Spacelike_Sudakov::No Branch for ("<<mo->kn_no<<"), "<<inflav
		     <<", "<<t<<", set on t="<<t0<<std::endl;
      mo->t    = mo->tout;
      mo->stat = 0;
      return 0;      // no further branching
    }
    SelectOne();
    z   = GetZ();   
    pt2 = -(1.-z)*t;
    if (!Veto(mo)) {
      msg.Debugging()<<"Spacelike_Sudakov::Dice Branch with t="
		     <<t<<", z="<<z<<", "<<inflav<<" for "<<lastint<<std::endl;
      UniformPhi();
      mo->z      = z;
      mo->t      = t;
      mo->phi    = phi;
      return 1;
    }    
  }
  msg.Debugging()<<"Spacelike_Sudakov::Banged out of Dice !"<<std::endl;
  mo->t    = mo->tout;
  mo->stat = 0;
  return 0; 
}

void Spacelike_Sudakov::ProduceT() {
  if (lastint <rpa.gen.Accu()) t = t0;
  t *= exp( 2.*M_PI*log(ran.Get()) / lastint );
  return;
}

bool Spacelike_Sudakov::Veto(Knot * mo) 
{  
  
  // "lower" cutoff reached / still spacelike
  if ((1.-z)*t>t0)   return 1;
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
  double weight  = pdf->GetXPDF(GetFlB())/pdf->GetXPDF(GetFlA()); 
  pdf->Calculate(x,sqrt(-t));
  pdfa->Calculate(x/z,sqrt(-t));
  weight        *= GetWeight(z,-t,0) * pdfa->GetXPDF(GetFlA())/pdf->GetXPDF(GetFlB()); 
  if (ran.Get() > weight) return 1;
  return 0;
}

bool Spacelike_Sudakov::CplVeto() 
{
  switch (cpl_scheme) {
  case 0 : 
    return 0;
  case 2 : 
    return (GetCoupling(0.25*t)/GetCoupling() > ran.Get()) ? 0 : 1;   
  default : 
    return (GetCoupling(pt2)/GetCoupling() > ran.Get()) ? 0 : 1;   
  }
}

bool Spacelike_Sudakov::PTVeto(Knot * mo) 
{
  // approximately pt^2/pl^2 with virtual masses neglected. Check this !

  double th = 4.*z*z*t/(4.*z*z*t-(1.-z)*x*x*pt2max);
  if (!inflav.Strong()) {
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



void Spacelike_Sudakov::Add(Splitting_Function * spl) 
{
  for (SplFunIter iter(group);iter();++iter) {
    if (iter()->GetFlB()==spl->GetFlB()) {
      iter()->Add(spl);
      return ;
    }
  }
  group.Append(new Backward_Splitting_Group(spl,pdf));
  selected=spl;
}


double Spacelike_Sudakov::CrudeInt(double _zmin, double _zmax) 
{
  SplFunIter iter(group);
  for (;iter();++iter)
    if (iter()->GetFlB()==inflav) { selected=iter(); break; }
  if (!iter()) std::cout<<"warning "<<inflav<<" not implemented in Splitting_Group"<<std::endl;
  return lastint = selected->CrudeInt(_zmin,_zmax);
}     

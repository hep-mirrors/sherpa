#include "CSSHOWER++/Showers/Kinematics_Base.H"
#include "CSSHOWER++/Tools/Singlet.H"
#include "CSSHOWER++/Showers/Sudakov.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Math/Poincare.H"

using namespace CSSHOWER;
using namespace PHASIC;
using namespace ATOOLS;

double Kinematics_FF::GetY(const double &Q2,const double &kt2,const double &z,
			   const double &mi2,const double &mj2,const double &mk2) const
{
  if (z<=0.0 || z>=1.0 || Q2<=mi2+mj2+mk2) return -1.0;
  return (kt2/(z*(1.0-z))+(1.0-z)/z*mi2+z/(1.0-z)*mj2)/(Q2-mi2-mj2-mk2);
}

int Kinematics_FF::MakeKinematics
(Parton *const split,const ATOOLS::Flavour & fli,
 const ATOOLS::Flavour & flj,Parton *&pc,const int mode)
{
  Parton * spect = split->GetSpect();
  Vec4D p1 = split->Momentum(), p2 = spect->Momentum();

  double mi2=sqr(p_ms->Mass(fli)),mj2=sqr(p_ms->Mass(flj));
  double mij2 = sqr(p_ms->Mass(split->GetFlavour())); 
  double mk2 = sqr(p_ms->Mass(spect->GetFlavour()));

  double y=GetY((p1+p2).Abs2(),split->KtTest(),split->ZTest(),mi2,mj2,mk2);
  Kin_Args ff(y,split->ZTest(),split->Phi());
  if (ConstructFFDipole(mi2,mj2,mij2,mk2,p1,p2,ff)<0) return -1;
  if (mode==0 && (ff.m_pi+ff.m_pj).Abs2()<split->TMin()) return -1;

  split->SetMomentum(ff.m_pi);
  spect->SetMomentum(ff.m_pk);
  if (pc==NULL) pc = new Parton(flj,ff.m_pj,pst::FS);
  else pc->SetMomentum(ff.m_pj);

  return 1;
}

double Kinematics_FI::GetY(const double &Q2,const double &kt2,const double &z,
			   const double &mi2,const double &mj2,const double &ma2) const
{
  if (z<=0.0 || z>=1.0 || Q2>=mi2+mj2+ma2) return -1.0;
  return 1.0/(1.0-(kt2/(z*(1.0-z))+mi2*(1.0-z)/z+mj2*z/(1.0-z))/(Q2-ma2-mi2-mj2));
}

int Kinematics_FI::MakeKinematics
(Parton *const split,const ATOOLS::Flavour & fli,
 const ATOOLS::Flavour & flj,Parton *&pc,const int mode)
{ 
  Parton * spect = split->GetSpect();
  Vec4D p1 = split->Momentum(), p2 = spect->Momentum(), rp1 = p1;

  double ma2 = p_ms->Mass2(spect->GetFlavour());
  double mi2=sqr(p_ms->Mass(fli)), mj2=sqr(p_ms->Mass(flj));
  double mij2 = sqr(p_ms->Mass(split->GetFlavour())); 
  
  double y=1.0-GetY((p1-p2).Abs2(),split->KtTest(),split->ZTest(),mi2,mj2,ma2);
  Kin_Args fi(y,split->ZTest(),split->Phi());
  if (ConstructFIDipole(mi2,mj2,mij2,ma2,p1,p2,fi)<0) return -1;
  if (mode==0 && (fi.m_pi+fi.m_pj).Abs2()<split->TMin()) return -1;

  split->SetMomentum(fi.m_pi);
  spect->SetMomentum(fi.m_pk);
  if (pc==NULL) pc = new Parton(flj,fi.m_pj,pst::FS);
  else pc->SetMomentum(fi.m_pj);
  
  return 1;
}

double Kinematics_IF::GetY(const double &Q2,const double &kt2,const double &z,
			   const double &ma2,const double &mi2,const double &mk2) const
{
  if (z<=0.0 || z>=1.0 || Q2>=ma2+mi2+mk2) return -1.0;
  return -z/(Q2-ma2-mi2-mk2)*((kt2+mi2)/(1.0-z)+(1.0-z)*ma2);
}

int Kinematics_IF::MakeKinematics
(Parton *const split,const ATOOLS::Flavour & fla,
 const ATOOLS::Flavour & fli,Parton *&pc,const int mode)
{
  Parton *b(NULL);
  for (PLiter pit(split->GetSing()->begin());pit!=split->GetSing()->end();++pit)
    if ((*pit)->GetType()==pst::IS && *pit!=split) {
      b=*pit;
      break;
    }
  if (b==NULL) THROW(fatal_error,"Corrupted singlet");
  double mb2(p_ms->Mass2(b->GetFlavour()));

  Parton * spect = split->GetSpect();
  Vec4D p1 = split->Momentum(), p2 = spect->Momentum();

  double mk2 = p_ms->Mass2(spect->GetFlavour());
  double mi2 = p_ms->Mass2(fli), ma2 = p_ms->Mass2(fla);
  double mai2 = p_ms->Mass2(split->GetFlavour()); 
  
  double y=GetY((p2-p1).Abs2(),split->KtTest(),split->ZTest(),ma2,mi2,mk2);
  Kin_Args ifp(y,split->ZTest(),split->Phi(),split->Kin());
  if (dabs(y-split->ZTest())<Kin_Args::s_uxeps) ifp.m_mode=1;
  if (ConstructIFDipole(ma2,mi2,mai2,mk2,mb2,p1,p2,b->Momentum(),ifp)<0) return -1;
  if (mode==0 && -(ifp.m_pi-ifp.m_pj).Abs2()<split->TMin()) return -1;

  split->SetLT(ifp.m_lam);
  split->SetMomentum(ifp.m_pi);
  spect->SetMomentum(ifp.m_pk);
  if (pc==NULL) pc = new Parton(fli,ifp.m_pj,pst::FS);
  else pc->SetMomentum(ifp.m_pj);
  
  return 1;
}

double Kinematics_II::GetY(const double &Q2,const double &kt2,const double &z,
			   const double &ma2,const double &mi2,const double &mb2) const
{
  if (z<=0.0 || z>=1.0 || Q2<=ma2+mi2+mb2) return -1.0;
  return z/(Q2-ma2-mb2-mi2)*((kt2+mi2)/(1.0-z)+(1.0-z)*ma2);
}

int Kinematics_II::MakeKinematics
(Parton *const split,const ATOOLS::Flavour & fli,
 const ATOOLS::Flavour & newfl,Parton *&pc,const int mode)
{
  Parton * spect = split->GetSpect();
  Vec4D p1 = split->Momentum(), p2 = spect->Momentum();
  
  double ma2 = sqr(p_ms->Mass(fli)), mi2 = sqr(p_ms->Mass(newfl));
  double mai2 = sqr(p_ms->Mass(split->GetFlavour()));
  double mb2 = sqr(p_ms->Mass(spect->GetFlavour()));

  double y=GetY((p1+p2).Abs2(),split->KtTest(),split->ZTest(),ma2,mi2,mb2);
  Kin_Args ii(y,split->ZTest(),split->Phi());
  if (ConstructIIDipole(ma2,mi2,mai2,mb2,p1,p2,ii)<0) return -1;
  if (mode==0 && -(ii.m_pi-ii.m_pj).Abs2()<split->TMin()) return -1;

  split->SetLT(ii.m_lam);
  split->SetMomentum(ii.m_pi);
  spect->SetMomentum(ii.m_pk);
  if (pc==NULL) pc = new Parton(newfl,ii.m_pj,pst::FS);
  else pc->SetMomentum(ii.m_pj);

  return 1;
}

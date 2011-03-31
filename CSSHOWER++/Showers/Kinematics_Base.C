#include "CSSHOWER++/Showers/Kinematics_Base.H"
#include "CSSHOWER++/Tools/Singlet.H"
#include "CSSHOWER++/Showers/Sudakov.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Math/Poincare.H"

using namespace CSSHOWER;
using namespace PHASIC;
using namespace ATOOLS;

void Kinematics_Base::SetFixVec(Parton *const p,Vec4D mom,
				const Kin_Args &lt,const int mode) const
{
  if (p->FixSpec()==Vec4D()) return;
  Vec4D oldp(p->OldMomentum()), ref(p->FixSpec());
  if (mode==3 || (mode==1 && lt.m_mode==0)) {
    Poincare_Sequence lam(lt.m_lam);
    lam.Invert();
    mom=lam*mom;
  }
  Poincare oldcms(oldp), newcms(mom);
  oldcms.Boost(ref);
  newcms.BoostBack(ref);
//   ref*=mom[0]/ref[0];
//   ref[0]=ref.PSpat();
  p->SetFixSpec(ref);
  p->SetOldMomentum(mom);
}

double Kinematics_FF::GetY(const double &Q2,const double &kt2,const double &z,
			   const double &mi2,const double &mj2,const double &mk2) const
{
  if (z<=0.0 || z>=1.0 || Q2<=mi2+mj2+mk2) return -1.0;
  return (kt2/(z*(1.0-z))+(1.0-z)/z*mi2+z/(1.0-z)*mj2)/(Q2-mi2-mj2-mk2);
}

int Kinematics_FF::MakeKinematics
(Parton *const split,const double & mi2,const double & mj2,
 const ATOOLS::Flavour &flj,Parton *&pc,const int mode)
{
  Parton * spect = split->GetSpect();
  Vec4D p1 = split->Momentum(), p2 = spect->Momentum();

  double mij2 = split->Mass2(), mk2 = spect->Mass2();
  if (mk2 && !spect->GetFlavour().Strong()) mk2=p2.Abs2();
  bool nospec(false);
  if (p_ms->Mass2(split->GetFlavour())>10.0 &&
      !split->GetFlavour().Strong()) {
    if (split->GetPrev()==NULL)
      THROW(fatal_error,"Missing splitting information");
    mij2=p1.Abs2();
    p2=split->GetPrev()->FixSpec();
    mk2=0.0;
    nospec=true;
  }

  double y=GetY((p1+p2).Abs2(),split->KtTest(),split->ZTest(),mi2,mj2,mk2);
  Kin_Args ff(y,split->ZTest(),split->Phi());
  if (ConstructFFDipole(mi2,mj2,mij2,mk2,p1,p2,ff)<0) return -1;
  if (mode==0 && (ff.m_pi+ff.m_pj).Abs2()<split->TMin()) return -1;

  split->SetMomentum(ff.m_pi);
  if (mi2) SetFixVec(split,ff.m_pi,ff,0);
  if (!nospec) {
    if (mk2) SetFixVec(spect,ff.m_pk,ff,0);
    spect->SetMomentum(ff.m_pk);
  }
  else if (!IsEqual(ff.m_pk,p2,1.0e-3))
    msg_Error()<<METHOD<<"(): Error in EW splitting.\n  Shifted p_k = "
	       <<p2<<" -> "<<ff.m_pk<<std::endl;
  if (pc==NULL) {
    pc = new Parton(flj,ff.m_pj,pst::FS);
    pc->SetMass2(p_ms->Mass2(flj));
  }
  else {
    if (mj2) SetFixVec(pc,ff.m_pj,ff,0);
    pc->SetMomentum(ff.m_pj);
  }

  return 1;
}

double Kinematics_FI::GetY(const double &Q2,const double &kt2,const double &z,
			   const double &mi2,const double &mj2,const double &ma2) const
{
  if (z<=0.0 || z>=1.0 || Q2>=mi2+mj2+ma2) return -1.0;
  return 1.0/(1.0-(kt2/(z*(1.0-z))+mi2*(1.0-z)/z+mj2*z/(1.0-z))/(Q2-ma2-mi2-mj2));
}

int Kinematics_FI::MakeKinematics
(Parton *const split,const double & mi2,const double & mj2,
 const ATOOLS::Flavour &flj,Parton *&pc,const int mode)
{ 
  Parton * spect = split->GetSpect();
  Vec4D p1 = split->Momentum(), p2 = spect->Momentum(), rp1 = p1;

  double ma2 = spect->Mass2(), mij2 = split->Mass2(); 
  bool nospec(false);
  if (p_ms->Mass2(split->GetFlavour())>10.0 &&
      !split->GetFlavour().Strong()) {
    if (split->GetPrev()==NULL)
      THROW(fatal_error,"Missing splitting information");
    mij2=p1.Abs2();
    p2=split->GetPrev()->FixSpec();
    ma2=0.0;
    nospec=true;
  }
  
  double y=1.0-GetY((p1-p2).Abs2(),split->KtTest(),split->ZTest(),mi2,mj2,ma2);
  Kin_Args fi(y,split->ZTest(),split->Phi());
  if (ConstructFIDipole(mi2,mj2,mij2,ma2,p1,p2,fi)<0) return -1;
  if (mode==0 && (fi.m_pi+fi.m_pj).Abs2()<split->TMin()) return -1;

  split->SetMomentum(fi.m_pi);
  if (mi2) SetFixVec(split,fi.m_pi,fi,2);
  if (!nospec) spect->SetMomentum(fi.m_pk);
  else if (!IsEqual(fi.m_pk,p2,1.0e-3))
    msg_Error()<<METHOD<<"(): Error in EW splitting.\n  Shifted p_k = "
	       <<p2<<" -> "<<fi.m_pk<<std::endl;
  if (pc==NULL) {
    pc = new Parton(flj,fi.m_pj,pst::FS);
    pc->SetMass2(p_ms->Mass2(flj));
  }
  else {
    if (mj2) SetFixVec(pc,fi.m_pj,fi,2);
    pc->SetMomentum(fi.m_pj);
  }
  
  return 1;
}

double Kinematics_IF::GetY(const double &Q2,const double &kt2,const double &z,
			   const double &ma2,const double &mi2,const double &mk2) const
{
  if (z<=0.0 || z>=1.0 || Q2>=ma2+mi2+mk2) return -1.0;
  return -z/(Q2-ma2-mi2-mk2)*((kt2+mi2)/(1.0-z)+(1.0-z)*ma2);
}

int Kinematics_IF::MakeKinematics
(Parton *const split,const double & ma2,const double & mi2,
 const ATOOLS::Flavour &fli,Parton *&pc,const int mode)
{
  Parton *b(NULL);
  for (PLiter pit(split->GetSing()->begin());pit!=split->GetSing()->end();++pit)
    if ((*pit)->GetType()==pst::IS && *pit!=split) {
      b=*pit;
      break;
    }
  if (b==NULL) THROW(fatal_error,"Corrupted singlet");
  double mb2(b->Mass2());

  Parton * spect = split->GetSpect();
  Vec4D p1 = split->Momentum(), p2 = spect->Momentum();

  double mk2 = spect->Mass2(), mai2 = split->Mass2(); 
  if (mk2 && !spect->GetFlavour().Strong()) mk2=p2.Abs2();
  
  double y=GetY((p2-p1).Abs2(),split->KtTest(),split->ZTest(),ma2,mi2,mk2);
  Kin_Args ifp(y,split->ZTest(),split->Phi(),split->Kin());
  if (dabs(y-split->ZTest())<Kin_Args::s_uxeps) ifp.m_mode=1;
  if (ConstructIFDipole(ma2,mi2,mai2,mk2,mb2,p1,p2,b->Momentum(),ifp)<0) return -1;
  if (mode==0 && -(ifp.m_pi-ifp.m_pj).Abs2()<split->TMin()) return -1;

  split->SetLT(ifp.m_lam);
  split->SetMomentum(ifp.m_pi);
  spect->SetMomentum(ifp.m_pk);
  if (mk2) SetFixVec(spect,ifp.m_pk,ifp,1);
  if (pc==NULL) {
    pc = new Parton(fli,ifp.m_pj,pst::FS);
    pc->SetMass2(p_ms->Mass2(fli));
  }
  else {
    if (mi2) SetFixVec(pc,ifp.m_pj,ifp,1);
    pc->SetMomentum(ifp.m_pj);
  }
  
  return 1;
}

double Kinematics_II::GetY(const double &Q2,const double &kt2,const double &z,
			   const double &ma2,const double &mi2,const double &mb2) const
{
  if (z<=0.0 || z>=1.0 || Q2<=ma2+mi2+mb2) return -1.0;
  return z/(Q2-ma2-mb2-mi2)*((kt2+mi2)/(1.0-z)+(1.0-z)*ma2);
}

int Kinematics_II::MakeKinematics
(Parton *const split,const double & ma2,const double & mi2,
 const ATOOLS::Flavour &newfl,Parton *&pc,const int mode)
{
  Parton * spect = split->GetSpect();
  Vec4D p1 = split->Momentum(), p2 = spect->Momentum();
  
  double mai2 = split->Mass2(), mb2 = spect->Mass2();

  double y=GetY((p1+p2).Abs2(),split->KtTest(),split->ZTest(),ma2,mi2,mb2);
  Kin_Args ii(y,split->ZTest(),split->Phi());
  if (ConstructIIDipole(ma2,mi2,mai2,mb2,p1,p2,ii)<0) return -1;
  if (mode==0 && -(ii.m_pi-ii.m_pj).Abs2()<split->TMin()) return -1;

  split->SetLT(ii.m_lam);
  split->SetMomentum(ii.m_pi);
  spect->SetMomentum(ii.m_pk);
  if (pc==NULL) {
    pc = new Parton(newfl,ii.m_pj,pst::FS);
    pc->SetMass2(p_ms->Mass2(newfl));
  }
  else {
    if (mi2) SetFixVec(pc,ii.m_pj,ii,3);
    pc->SetMomentum(ii.m_pj);
  }

  return 1;
}

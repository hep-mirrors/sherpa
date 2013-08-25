#include "MCATNLO/Showers/Kinematics_Base.H"
#include "MCATNLO/Tools/Singlet.H"
#include "MCATNLO/Showers/Sudakov.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Math/Poincare.H"

using namespace MCATNLO;
using namespace PHASIC;
using namespace ATOOLS;

double Kinematics_FF::GetY(const double &Q2,const double &kt2,const double &z,
			   const double &mi2,const double &mj2,const double &mk2,
			   const ATOOLS::Flavour &fla,const ATOOLS::Flavour &flc,
			   const bool force) const
{
  if (!force && (z<=0.0 || z>=1.0 || Q2<=mi2+mj2+mk2)) return -1.0;
  if (fla.IsFermion()) {
    if (flc.IsFermion()) return kt2/z/(Q2-mi2-mj2-mk2);
    return kt2/(1.0-z)/(Q2-mi2-mj2-mk2);
  }
  if (flc.IsFermion()) return kt2/(Q2-mi2-mj2-mk2);
  return kt2/(z*(1.0-z))/(Q2-mi2-mj2-mk2);
}

int Kinematics_FF::MakeKinematics
(Parton *const split,const ATOOLS::Flavour &fli,
 const ATOOLS::Flavour &flj,Parton *&pc)
{
  Parton * spect = split->GetSpect();
  Vec4D p1 = split->Momentum(), p2 = spect->Momentum();

  double mi2 = p_ms->Mass2(fli), mj2 = p_ms->Mass2(flj);
  double mij2 = p_ms->Mass2(split->GetFlavour()), mk2 = p_ms->Mass2(spect->GetFlavour());
  if (mk2 && !spect->GetFlavour().Strong()) mk2=p2.Abs2();

  double y=GetY((p1+p2).Abs2(),split->KtTest(),split->ZTest(),mi2,mj2,mk2,
		split->GetFlavour(),flj,1);
  Kin_Args ff(y,split->ZTest(),split->Phi());
  if (ConstructFFDipole(mi2,mj2,mij2,mk2,p1,p2,ff)<0) return -1;

  split->SetMomentum(ff.m_pi);
  spect->SetMomentum(ff.m_pk);
  if (pc==NULL) {
    pc = new Parton(flj,ff.m_pj,pst::FS);
  }
  else {
    pc->SetMomentum(ff.m_pj);
  }

  return 1;
}

double Kinematics_FI::GetY(const double &Q2,const double &kt2,const double &z,
			   const double &mi2,const double &mj2,const double &ma2,
			   const ATOOLS::Flavour &fla,const ATOOLS::Flavour &flc,
			   const bool force) const
{
  if (!force && (z<=0.0 || z>=1.0 || Q2>=mi2+mj2+ma2)) return -1.0;
  if (fla.IsFermion()) {
    if (flc.IsFermion()) return 1.0/(1.0-kt2/z/(Q2-ma2-mi2-mj2));
    return 1.0/(1.0-kt2/(1.0-z)/(Q2-ma2-mi2-mj2));
  }
  if (flc.IsFermion()) return 1.0/(1.0-kt2/(Q2-ma2-mi2-mj2));
  return 1.0/(1.0-kt2/(z*(1.0-z))/(Q2-ma2-mi2-mj2));
}

int Kinematics_FI::MakeKinematics
(Parton *const split,const ATOOLS::Flavour &fli,
 const ATOOLS::Flavour &flj,Parton *&pc)
{ 
  Parton * spect = split->GetSpect();
  Vec4D p1 = split->Momentum(), p2 = spect->Momentum(), rp1 = p1;

  double mi2 = p_ms->Mass2(fli), mj2 = p_ms->Mass2(flj);
  double ma2 = p_ms->Mass2(spect->GetFlavour()), mij2 = p_ms->Mass2(split->GetFlavour()); 
  
  double Q2((p1-p2).Abs2());
  double y=GetY(Q2,split->KtTest(),split->ZTest(),mi2,mj2,ma2,
		split->GetFlavour(),flj,1);
  y=1.0-y*(Q2-mij2-ma2)/(Q2-mi2-mj2-ma2);
  Kin_Args fi(y,split->ZTest(),split->Phi());
  if (ConstructFIDipole(mi2,mj2,mij2,ma2,p1,p2,fi)<0) return -1;

  split->SetMomentum(fi.m_pi);
  spect->SetMomentum(fi.m_pk);
  if (pc==NULL) {
    pc = new Parton(flj,fi.m_pj,pst::FS);
  }
  else {
    pc->SetMomentum(fi.m_pj);
  }
  
  return 1;
}

double Kinematics_IF::GetY(const double &Q2,const double &kt2,const double &z,
			   const double &ma2,const double &mi2,const double &mk2,
			   const ATOOLS::Flavour &flb,const ATOOLS::Flavour &flc,
			   const bool force) const
{
  if (!force && (z<=0.0 || z>=1.0 || Q2>=ma2+mi2+mk2)) return -1.0;
  if (flc.IsFermion()) return -z/(Q2-ma2-mi2-mk2)*kt2;
  return -z/(Q2-ma2-mi2-mk2)*kt2/(1.0-z);
}

int Kinematics_IF::MakeKinematics
(Parton *const split,const ATOOLS::Flavour &fla,
 const ATOOLS::Flavour &fli,Parton *&pc)
{
  Parton *b(NULL);
  for (PLiter pit(split->GetSing()->begin());pit!=split->GetSing()->end();++pit)
    if ((*pit)->GetType()==pst::IS && *pit!=split) {
      b=*pit;
      break;
    }
  if (b==NULL) THROW(fatal_error,"Corrupted singlet");
  double ma2(p_ms->Mass2(fla)), mi2(p_ms->Mass2(fli));
  double mb2(p_ms->Mass2(b->GetFlavour()));

  Parton * spect = split->GetSpect();
  Vec4D p1 = split->Momentum(), p2 = spect->Momentum();

  double mk2 = p_ms->Mass2(spect->GetFlavour()), mai2 = p_ms->Mass2(split->GetFlavour()); 
  if (mk2 && !spect->GetFlavour().Strong()) mk2=p2.Abs2();
  
  double y=GetY((p2-p1).Abs2(),split->KtTest(),split->ZTest(),ma2,mi2,mk2,
		split->GetFlavour(),fli,1);
  Kin_Args ifp(y,split->ZTest(),split->Phi(),split->Kin());
  if (dabs(y-split->ZTest())<Kin_Args::s_uxeps) ifp.m_mode=1;
  if (ConstructIFDipole(ma2,mi2,mai2,mk2,mb2,p1,p2,b->Momentum(),ifp)<0) return -1;

  split->SetLT(ifp.m_lam);
  split->SetMomentum(ifp.m_pi);
  spect->SetMomentum(ifp.m_pk);
  if (pc==NULL) {
    pc = new Parton(fli,ifp.m_pj,pst::FS);
  }
  else {
    pc->SetMomentum(ifp.m_pj);
  }
  
  return 1;
}

double Kinematics_II::GetY(const double &Q2,const double &kt2,const double &z,
			   const double &ma2,const double &mi2,const double &mb2,
			   const ATOOLS::Flavour &flb,const ATOOLS::Flavour &flc,
			   const bool force) const
{
  if (!force && (z<=0.0 || z>=1.0 || Q2<=ma2+mi2+mb2)) return -1.0;
  if (flc.IsFermion()) return z/(Q2-ma2-mb2-mi2)*kt2;
  return z/(Q2-ma2-mb2-mi2)*kt2/(1.0-z);
}

int Kinematics_II::MakeKinematics
(Parton *const split,const ATOOLS::Flavour &fla,
 const ATOOLS::Flavour &newfl,Parton *&pc)
{
  Parton * spect = split->GetSpect();
  Vec4D p1 = split->Momentum(), p2 = spect->Momentum();
  
  double ma2 = p_ms->Mass2(fla), mi2 = p_ms->Mass2(newfl);
  double mai2 = p_ms->Mass2(split->GetFlavour()), mb2 = p_ms->Mass2(spect->GetFlavour());

  double y=GetY((p1+p2).Abs2(),split->KtTest(),split->ZTest(),ma2,mi2,mb2,
		split->GetFlavour(),newfl,1);
  Kin_Args ii(y,split->ZTest(),split->Phi(),split->Kin());
  if (ConstructIIDipole(ma2,mi2,mai2,mb2,p1,p2,ii)<0) return -1;

  split->SetLT(ii.m_lam);
  split->SetMomentum(ii.m_pi);
  spect->SetMomentum(ii.m_pk);
  if (pc==NULL) {
    pc = new Parton(newfl,ii.m_pj,pst::FS);
  }
  else {
    pc->SetMomentum(ii.m_pj);
  }

  return 1;
}

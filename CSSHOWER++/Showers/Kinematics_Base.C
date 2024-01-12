#include "CSSHOWER++/Showers/Kinematics_Base.H"
#include "CSSHOWER++/Tools/Singlet.H"
#include "CSSHOWER++/Showers/Sudakov.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace CSSHOWER;
using namespace PHASIC;
using namespace ATOOLS;

double Kinematics_FF::GetKT2(const double &Q2,const double &y,const double &z,
			     const double &mi2,const double &mj2,const double &mk2,
			     const ATOOLS::Flavour &fla,const ATOOLS::Flavour &flc) const
{
  double pipj=(Q2-mi2-mj2-mk2)*y;
  if (m_evolscheme==0 || m_evolscheme==2 || m_evolscheme == 20) {
    double kt2=pipj*z*(1.0-z)-sqr(1.0-z)*mi2-sqr(z)*mj2;
    if (m_evolscheme==0) return kt2;
    if (m_evolscheme==2) return kt2+mi2+mj2;
    // like scheme 2, but only applied to gluon splitters
    if (m_evolscheme==20)
      return fla.IsGluon()?(kt2+mi2+mj2):kt2;
  }
  else if (m_evolscheme==1 || m_evolscheme==3 || m_evolscheme == 30) {
    double kt2=pipj*z*(1.0-z);
    if (fla.IsFermion()) kt2=pipj*(flc.IsVector()?(1.0-z):z);
    else if (flc.IsFermion()) kt2=pipj;
    if (m_evolscheme==1) return kt2;
    if (m_evolscheme==3) return kt2+mi2+mj2;
        // like scheme 3, but only applied to gluon splitters
    if (m_evolscheme==30)
      return fla.IsGluon()?(kt2+mi2+mj2):kt2;
  }
  else THROW(fatal_error, "Not implemented");
  return 0.0;
}

double Kinematics_FF::GetY(const double &Q2,const double &_kt2,const double &z,
			   const double &mi2,const double &mj2,const double &mk2,
			   const ATOOLS::Flavour &fla,const ATOOLS::Flavour &flc,
			   const bool force) const
{
  if (!force && (z<=0.0 || z>=1.0 || Q2<=mi2+mj2+mk2)) return -1.0;
  double kt2=_kt2;
  if (m_evolscheme==2  || m_evolscheme==3)  kt2=kt2-mi2-mj2;
  if (m_evolscheme==20 || m_evolscheme==30)
    kt2=(fla.IsGluon())?(kt2-mi2-mj2):kt2;

  if (m_evolscheme==0 || m_evolscheme==2 || m_evolscheme == 20) {
    return (kt2/(z*(1.0-z))+(1.0-z)/z*mi2+z/(1.0-z)*mj2)/(Q2-mi2-mj2-mk2);
  }
  else if (m_evolscheme==1 || m_evolscheme==3 || m_evolscheme == 30) {
    if (fla.IsFermion()) {
      if (flc.IsFermion()) return kt2/z/(Q2-mi2-mj2-mk2);
      return kt2/(1.0-z)/(Q2-mi2-mj2-mk2);
    }
    if (flc.IsFermion()) return kt2/(Q2-mi2-mj2-mk2);
    return kt2/(z*(1.0-z))/(Q2-mi2-mj2-mk2);
  }
  else THROW(fatal_error, "Not implemented");
  return 0.0;
}

int Kinematics_FF::MakeKinematics
(Parton *const split,const double & mi2,const double & mj2,
 const ATOOLS::Flavour &flj,Parton *&pc,const int mode)
{
  Parton * spect = split->GetSpect();
  Vec4D p1 = split->Momentum(), p2 = spect->Momentum();

  double mij2 = split->Mass2(), mk2 = spect->Mass2();

  double y=GetY((p1+p2).Abs2(),split->KtTest(),split->ZTest(),mi2,mj2,mk2,
		split->GetFlavour(),flj,1);
  Kin_Args ff(y,split->ZTest(),split->Phi());
  if (ConstructFFDipole(mi2,mj2,mij2,mk2,p1,p2,ff)<0 ||
      !ValidateDipoleKinematics(mi2, mj2, mk2, ff)) return -1;

  split->SetMomentum(ff.m_pi);
  spect->SetMomentum(ff.m_pk);
  if (pc==NULL) {
    pc = new Parton(flj,ff.m_pj,pst::FS);
    pc->SetMass2(p_ms->Mass2(flj));
  }
  else {
    pc->SetMomentum(ff.m_pj);
  }

  return 1;
}

double Kinematics_FI::GetKT2(const double &Q2,const double &y,const double &z,
			     const double &mi2,const double &mj2,const double &ma2,
			     const ATOOLS::Flavour &fla,const ATOOLS::Flavour &flc) const
{
  double pipj=-(Q2-ma2-mi2-mj2)*(1.0-y)/y;
  if (m_evolscheme==0 || m_evolscheme==2 || m_evolscheme == 20) {
    double kt2=pipj*z*(1.0-z)-sqr(1.0-z)*mi2-sqr(z)*mj2;
    if (m_evolscheme==0) return kt2;
    if (m_evolscheme==2) return kt2+mi2+mj2;
    if (m_evolscheme==20)
      return fla.IsGluon()?(kt2+mi2+mj2):kt2;
  }
  else if (m_evolscheme==1 || m_evolscheme==3 || m_evolscheme == 30) {
    double kt2=pipj*z*(1.0-z);
    if (fla.IsFermion()) kt2=pipj*(flc.IsVector()?(1.0-z):z);
    else if (flc.IsFermion()) kt2=pipj;
    if (m_evolscheme==1) return kt2;
    if (m_evolscheme==3) return kt2+mi2+mj2;
    if (m_evolscheme==30)
      return fla.IsGluon()?(kt2+mi2+mj2):kt2;
  }
  else THROW(fatal_error, "Not implemented");
  return 0.0;

}

double Kinematics_FI::GetY(const double &Q2,const double &_kt2,const double &z,
			   const double &mi2,const double &mj2,const double &ma2,
			   const ATOOLS::Flavour &fla,const ATOOLS::Flavour &flc,
			   const bool force) const
{
  if (!force && (z<=0.0 || z>=1.0 || Q2>=mi2+mj2+ma2)) return -1.0;
  double kt2=_kt2;
  if (m_evolscheme==2 || m_evolscheme==3) kt2=kt2-mi2-mj2;
  if (m_evolscheme==20 || m_evolscheme==30)
    kt2=(fla.IsGluon())?(kt2-mi2-mj2):kt2;

  if (m_evolscheme==0 || m_evolscheme==2 || m_evolscheme == 20) {
    return 1.0/(1.0-(kt2/(z*(1.0-z))+mi2*(1.0-z)/z+mj2*z/(1.0-z))/(Q2-ma2-mi2-mj2));
  }
  else if (m_evolscheme==1 || m_evolscheme==3 || m_evolscheme == 30) {
    if (fla.IsFermion()) {
      if (flc.IsFermion()) return 1.0/(1.0-kt2/z/(Q2-ma2-mi2-mj2));
      return 1.0/(1.0-kt2/(1.0-z)/(Q2-ma2-mi2-mj2));
    }
    if (flc.IsFermion()) return 1.0/(1.0-kt2/(Q2-ma2-mi2-mj2));
    return 1.0/(1.0-kt2/(z*(1.0-z))/(Q2-ma2-mi2-mj2));
  }
  else THROW(fatal_error, "Not implemented");
  return 0.0;
}

int Kinematics_FI::MakeKinematics
(Parton *const split,const double & mi2,const double & mj2,
 const ATOOLS::Flavour &flj,Parton *&pc,const int mode)
{
  Parton * spect = split->GetSpect();
  Vec4D p1 = split->Momentum(), p2 = spect->Momentum(), rp1 = p1;

  double ma2 = spect->Mass2(), mij2 = split->Mass2();

  double Q2((p1-p2).Abs2());
  double y=GetY(Q2,split->KtTest(),split->ZTest(),mi2,mj2,ma2,
		split->GetFlavour(),flj,1);
  Kin_Args fi(1.0-y,split->ZTest(),split->Phi(),8);
  if (ConstructFIDipole(mi2,mj2,mij2,ma2,p1,p2,fi)<0 ||
      !ValidateDipoleKinematics(mi2, mj2, ma2, fi)) return -1;

  split->SetMomentum(fi.m_pi);
  spect->SetMomentum(fi.m_pk);
  if (pc==NULL) {
    pc = new Parton(flj,fi.m_pj,pst::FS);
    pc->SetMass2(p_ms->Mass2(flj));
  }
  else {
    pc->SetMomentum(fi.m_pj);
  }

  return 1;
}

double Kinematics_IF::GetKT2(const double &Q2,const double &y,const double &z,
			     const double &ma2,const double &mi2,const double &mk2,
			     const ATOOLS::Flavour &flb,const ATOOLS::Flavour &flc) const
{
  /// for is splitting we pass only the flavour of the parton
  /// that enters the ME (b) and the final state one (c).
  /// in a g -> q q~ splitting fl(b) = - fl(c)!
  double pipj=(Q2-ma2-mi2-mk2)*y/z;
  const bool isgluonsplitting((flb.Kfcode() == flc.Kfcode())) ;
  if (m_evolscheme==0 || m_evolscheme==2 || m_evolscheme == 20) {
    double kt2=-pipj*(1.0-z)-mi2-sqr(1.0-z)*ma2;
    if (m_evolscheme==0) return kt2;
    if (m_evolscheme==2) return kt2+mi2+ma2;
    if (m_evolscheme==20)
      return (isgluonsplitting) ? (kt2+mi2+ma2):kt2 ;
  }
  else if (m_evolscheme==1 || m_evolscheme==3 || m_evolscheme == 30) {
    double kt2=-pipj*(1.0-z);
    if (flc.IsFermion()) kt2=-pipj;
    if (m_evolscheme==1) return kt2;
    if (m_evolscheme==3) return kt2+mi2+ma2;
    if (m_evolscheme==30)
      return (isgluonsplitting) ? (kt2+mi2+ma2):kt2 ;
  }
  else THROW(fatal_error, "Not implemented");
  return 0.0;
}

double Kinematics_IF::GetY(const double &Q2,const double &_kt2,const double &z,
			   const double &ma2,const double &mi2,const double &mk2,
			   const ATOOLS::Flavour &flb,const ATOOLS::Flavour &flc,
			   const bool force) const
{
  if (!force && (z<=0.0 || z>=1.0 || Q2>=ma2+mi2+mk2)) return -1.0;
  double kt2=_kt2;
  const bool isgluonsplitting((flb.Kfcode() == flc.Kfcode())) ;
  if (m_evolscheme==2 || m_evolscheme==3) kt2=kt2-mi2-ma2;
  if (m_evolscheme==20 || m_evolscheme==30)
    kt2=(isgluonsplitting)?(kt2-mi2-ma2):kt2;
  if (m_evolscheme==0 || m_evolscheme==2 || m_evolscheme == 20) {
    return -z/(Q2-ma2-mi2-mk2)*((kt2+mi2)/(1.0-z)+(1.0-z)*ma2);
  }
  else if (m_evolscheme==1 || m_evolscheme==3 || m_evolscheme == 30) {
    if (flc.IsFermion()) return -z/(Q2-ma2-mi2-mk2)*kt2;
    return -z/(Q2-ma2-mi2-mk2)*kt2/(1.0-z);
  }
  else THROW(fatal_error, "Not implemented");
  return 0.0;
}

int Kinematics_IF::MakeKinematics
(Parton *const split,const double & ma2,const double & mi2,
 const ATOOLS::Flavour &fli,Parton *&pc,const int mode)
{
  if (split->ForcedDecay()) {
    Kin_Args args(0,split->ZTest(),split->Phi(),split->Kin());
    return MakeForcedDecayKinematics(split,pc,args,mode);
  }
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

  double y=GetY((p2-p1).Abs2(),split->KtTest(),split->ZTest(),ma2,mi2,mk2,
		split->GetFlavour(),fli,1);
  Kin_Args ifp(y,split->ZTest(),split->Phi(),split->Kin());
  if (dabs(y-split->ZTest())<Kin_Args::s_uxeps) ifp.m_mode=1;
  if (ConstructIFDipole(ma2,mi2,mai2,mk2,mb2,p1,p2,b->Momentum(),ifp)<0 ||
      !ValidateDipoleKinematics(ma2, mi2, ifp.m_mk2>=0.0?ifp.m_mk2:mk2, ifp)) return -1;

  split->SetLT(ifp.m_lam);
  ifp.m_lam.Invert();
  split->SetMomentum(ifp.m_lam*ifp.m_pi);
  spect->SetMomentum(ifp.m_lam*ifp.m_pk);
  if (p1[3]*ifp.m_pi[3]<0.) return -1;
  if (pc==NULL) {
    pc = new Parton(fli,ifp.m_lam*ifp.m_pj,pst::FS);
    pc->SetMass2(p_ms->Mass2(fli));
  }
  else {
    pc->SetMomentum(ifp.m_lam*ifp.m_pj);
  }
  return 1;
}

double Kinematics_II::GetKT2(const double &Q2,const double &y,const double &z,
			     const double &ma2,const double &mi2,const double &mb2,
			     const ATOOLS::Flavour &flb,const ATOOLS::Flavour &flc) const
{
  /// for is splitting we pass only the flavour of the parton
  /// that enters the ME (b) and the final state one (c).
  /// in a g -> q q~ splitting fl(b) = - fl(c)!
  double pipj=(Q2-ma2-mi2-mb2)*y/z;
  const bool isgluonsplitting((flb.Kfcode() == flc.Kfcode())) ;
  if (m_evolscheme==0 || m_evolscheme==2 || m_evolscheme == 20) {
    double kt2=pipj*(1.0-z)-mi2-sqr(1.0-z)*ma2;
    if (m_evolscheme==0) return kt2;
    if (m_evolscheme==2) return kt2+mi2+ma2;
    if (m_evolscheme==20)
      return (isgluonsplitting) ? (kt2+mi2+ma2):kt2 ;
  }
  else if (m_evolscheme==1 || m_evolscheme==3 || m_evolscheme == 30) {
    double kt2=pipj*(1.0-z);
    if (flc.IsFermion()) kt2=pipj;
    if (m_evolscheme==1) return kt2;
    if (m_evolscheme==3) return kt2+mi2+ma2;
    if (m_evolscheme==30)
      return (isgluonsplitting) ? (kt2+mi2+ma2):kt2 ;
  }
  else THROW(fatal_error, "Not implemented");
  return 0.0;
}

double Kinematics_II::GetY(const double &Q2,const double &_kt2,const double &z,
			   const double &ma2,const double &mi2,const double &mb2,
			   const ATOOLS::Flavour &flb,const ATOOLS::Flavour &flc,
			   const bool force) const
{
  if (!force && (z<=0.0 || z>=1.0 || Q2<=ma2+mi2+mb2)) return -1.0;
  double kt2=_kt2;
  const bool isgluonsplitting((flb.Kfcode() == flc.Kfcode()));
  if (m_evolscheme==2 || m_evolscheme==3) kt2=kt2-mi2-ma2;
  if (m_evolscheme==20 || m_evolscheme==30)
    kt2=(isgluonsplitting)?(kt2-mi2-ma2):kt2;

  if (m_evolscheme==0 || m_evolscheme==2 || m_evolscheme == 20) {
    return z/(Q2-ma2-mb2-mi2)*((kt2+mi2)/(1.0-z)+(1.0-z)*ma2);
  }
  else if (m_evolscheme==1 || m_evolscheme==3 || m_evolscheme == 30) {
    if (flc.IsFermion()) return z/(Q2-ma2-mb2-mi2)*kt2;
    return z/(Q2-ma2-mb2-mi2)*kt2/(1.0-z);
  }
  else THROW(fatal_error, "Not implemented");
  return 0.0;
}

int Kinematics_II::MakeKinematics
(Parton *const split,const double & ma2,const double & mi2,
 const ATOOLS::Flavour &newfl,Parton *&pc,const int mode)
{
  if (split->ForcedDecay()) {
    Kin_Args args(0,split->ZTest(),split->Phi(),split->Kin());
    return MakeForcedDecayKinematics(split,pc,args,mode);
  }
  Parton * spect = split->GetSpect();
  Vec4D p1 = split->Momentum(), p2 = spect->Momentum();

  double mai2 = split->Mass2(), mb2 = spect->Mass2();

  double y=GetY((p1+p2).Abs2(),split->KtTest(),split->ZTest(),ma2,mi2,mb2,
		split->GetFlavour(),newfl,1);
  Kin_Args ii(y,split->ZTest(),split->Phi(),split->Kin());
  if (ConstructIIDipole(ma2,mi2,mai2,mb2,p1,p2,ii)<0 ||
      !ValidateDipoleKinematics(ma2, mi2, mb2, ii)) return -1;

  split->SetLT(ii.m_lam);
  ii.m_lam.Invert();
  split->SetMomentum(ii.m_lam*ii.m_pi);
  spect->SetMomentum(ii.m_lam*ii.m_pk);
  if (p1[3]*ii.m_pi[3]<0.) return -1;
  if (pc==NULL) {
    pc = new Parton(newfl,ii.m_lam*ii.m_pj,pst::FS);
    pc->SetMass2(p_ms->Mass2(newfl));
  }
  else {
    pc->SetMomentum(ii.m_lam*ii.m_pj);
  }
  return 1;
}

int Kinematics_Base::
MakeForcedDecayKinematics(Parton *const split,Parton *&pnew,Kin_Args & args,const int mode) const {
  Vec4D p1_t      = split->Momentum(), p2_t = Vec4D(0.,0.,0.,0.);
  for (PLiter plit=split->GetSing()->begin();plit!=split->GetSing()->end();plit++) {
    if ((*plit)->GetType()==pst::FS) p2_t += (*plit)->Momentum();
  }
  double s12_t    = (p1_t+p2_t).Abs2(),    E12 = sqrt(s12_t), z = split->ZTest();
  Vec4D  pplus    = E12/2. * (p1_t[3]>0. ? Vec4D(1.,0,0,1.) : Vec4D(1.,0,0,-1.));
  Vec4D  pminus   = E12/2. * (p1_t[3]>0. ? Vec4D(1.,0,0,-1.) : Vec4D(1.,0,0,1.));
  double m12      = Max(0.,p1_t.Abs2()), m22 = Max(0.,p2_t.Abs2());
  if (m22<1.e-8) m22 = 0.;
  double alpha1_t = 2.*p1_t*pminus/s12_t, beta1_t = 2.*p1_t*pplus/s12_t;
  double alpha2_t = 2.*p2_t*pminus/s12_t, beta2_t = 2.*p2_t*pplus/s12_t;
  Vec4D  P1       = alpha1_t*pplus + beta1_t*pminus;
  Vec4D  P2       = alpha2_t*pplus + beta2_t*pminus;
  Vec4D  kperp    = p2_t - P2;
  double kt2      = dabs(kperp.Abs2())>1.e-12 ? dabs(kperp.Abs2()) : 0.;
  double M2       = dabs(P2.Abs2())   >1.e-12 ? dabs(P2.Abs2())    : 0.;
  double MT2      = M2+kt2;
  double At       = alpha1_t*(1.-z)/z + alpha2_t;
  double c        = MT2 * alpha1_t - m12 * alpha2_t;
  double d        = m12 * alpha1_t * alpha2_t;
  double g        = MT2 * alpha1_t * alpha2_t;
  double zmax     = alpha1_t * c / ( sqr( sqrt(d)+sqrt(g) ) + (alpha1_t - alpha2_t) * c);
  if (z>zmax) {
    split->SetZTest(z = zmax);
    At = alpha1_t*(1.-z)/z + alpha2_t;
  }
  double alphaj   = 1./(2.*c) * (At*c + (d-g) + sqrt( Max(0., sqr(At*c+d-g) - 4.*At*c*d)));
  double betaj    = m12/(s12_t*alphaj);
  Vec4D  pj       = alphaj*pplus + betaj*pminus;
  double alpha1   = alpha1_t/z;
  double beta1    = 0.;
  Vec4D  p1       = alpha1*pplus + beta1*pminus;
  double alpha2   = alpha1_t*(1.-z)/z+alpha2_t-alphaj;
  double beta2    = MT2/(s12_t*alpha2);
  Vec4D  p2long   = alpha2*pplus + beta2*pminus;
  Vec4D  p2       = p2long + kperp;
  Vec4D  bb2      = Vec4D(p2long[0],0.,0.,-p2long[3]);
  args.m_lam.push_back(Poincare(bb2));
  args.m_lam.push_back(P2);
  Vec4D test, tmp;
  for (PLiter plit=split->GetSing()->begin();plit!=split->GetSing()->end();plit++) {
    if ((*plit)->GetType()==pst::FS) {
      tmp   = (*plit)->Momentum();
      test += args.m_lam*tmp;
    }
  }
  split->SetMomentum(p1);
  split->SetLT(args.m_lam);
  pnew = new Parton(split->GetFlavour().Bar(),pj,pst::FS);
  pnew->SetMass2(p_ms->Mass2(split->GetFlavour()));
  return 1;
}


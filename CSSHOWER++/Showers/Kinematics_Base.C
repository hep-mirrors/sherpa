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
  Vec4D np(0.0,cross(Vec3D(oldp),Vec3D(mom)));
  if (np.PSpat2()<=1.0e-6) {
    msg_Debugging()<<"Set fixed n_perp\n";
    np=Vec4D(0.0,cross(Vec3D(oldp),Vec3D(1.0,1.0,0.0)));
  }
  np*=1.0/np.PSpat();
  Vec4D lp(0.0,cross(Vec3D(oldp),Vec3D(np)));
  lp*=1.0/lp.PSpat();
  Vec4D pl(0.0,(Vec3D(ref)*Vec3D(lp))*lp);
  Vec4D pn(0.0,(Vec3D(ref)*Vec3D(np))*np);
  double D(oldp*ref), T(oldp.PSpat()), F(ref[0]);
  double Q(mom[0]), P(mom.PSpat()), S(mom.Abs2());
  Poincare rot(oldp,mom);
  if (oldp.Abs2()>1.0e-3 && mom.Abs2()>1.0e-3) {
    Poincare oldcms(oldp), newcms(mom);
    oldcms.Boost(ref);
    rot.Rotate(ref);
    newcms.BoostBack(ref);
  }
  else {
    double E((Q*D+P/T*(F*S-Q*D))/S);
    ref=Vec4D(E,Vec3D(mom)*(Q*E-D)/(P*P));
    ref+=pn+rot*pl;
  }
  p->SetFixSpec(ref);
  p->SetOldMomentum(mom);
}

double Kinematics_FF::GetKT2(const double &Q2,const double &y,const double &z,
			     const double &mi2,const double &mj2,const double &mk2,
			     const ATOOLS::Flavour &fla,const ATOOLS::Flavour &flc) const
{
  double pipj=(Q2-mi2-mj2-mk2)*y;
  if (m_evolscheme==0) {
    return pipj*z*(1.0-z)-sqr(1.0-z)*mi2-sqr(z)*mj2;
  }
  double kt2=pipj*z*(1.0-z);
  if (fla.IsFermion()) kt2=pipj*(flc.IsVector()?(1.0-z):z);
  else if (flc.IsFermion()) kt2=pipj;
  return kt2;
}

double Kinematics_FF::GetY(const double &Q2,const double &kt2,const double &z,
			   const double &mi2,const double &mj2,const double &mk2,
			   const ATOOLS::Flavour &fla,const ATOOLS::Flavour &flc,
			   const bool force) const
{
  if (!force && (z<=0.0 || z>=1.0 || Q2<=mi2+mj2+mk2)) return -1.0;
  if (m_evolscheme==0) {
    return (kt2/(z*(1.0-z))+(1.0-z)/z*mi2+z/(1.0-z)*mj2)/(Q2-mi2-mj2-mk2);
  }
  if (fla.IsFermion()) {
    if (flc.IsFermion()) return kt2/z/(Q2-mi2-mj2-mk2);
    return kt2/(1.0-z)/(Q2-mi2-mj2-mk2);
  }
  if (flc.IsFermion()) return kt2/(Q2-mi2-mj2-mk2);
  return kt2/(z*(1.0-z))/(Q2-mi2-mj2-mk2);
}

int Kinematics_FF::MakeKinematics
(Parton *const split,const double & mi2,const double & mj2,
 const ATOOLS::Flavour &flj,Parton *&pc,const int mode)
{
  Parton * spect = split->GetSpect();
  Vec4D p1 = split->Momentum(), p2 = spect->Momentum();

  double mij2 = split->Mass2(), mk2 = spect->Mass2();
  if (spect->KScheme()) mk2=p2.Abs2();
  bool nospec(false);
  if (mode && split->GetPrev() && split->GetPrev()->KScheme()) {
    mij2=p1.Abs2();
    p2=split->GetPrev()->FixSpec();
    mk2=0.0;
    nospec=true;
  }

  double y=GetY((p1+p2).Abs2(),split->KtTest(),split->ZTest(),mi2,mj2,mk2,
		split->GetFlavour(),flj,1);
  Kin_Args ff(y,split->ZTest(),split->Phi());
  if (ConstructFFDipole(mi2,mj2,mij2,mk2,p1,p2,ff)<0) return -1;

  split->SetMomentum(ff.m_pi);
  if (mi2) SetFixVec(split,ff.m_pi,ff,0);
  if (!nospec) {
    if (mk2) SetFixVec(spect,ff.m_pk,ff,0);
    spect->SetMomentum(ff.m_pk);
  }
  else if (!IsEqual(ff.m_pk,p2,1.0e-3))
    msg_Error()<<METHOD<<"(): Error in EW splitting ( y = "<<y
	       <<" ).\n  Shifted p_k = "
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

double Kinematics_FI::GetKT2(const double &Q2,const double &y,const double &z,
			     const double &mi2,const double &mj2,const double &ma2,
			     const ATOOLS::Flavour &fla,const ATOOLS::Flavour &flc) const
{
  double pipj=-(Q2-ma2-mi2-mj2)*(1.0-y)/y;
  if (m_evolscheme==0) {
    return pipj*z*(1.0-z)-sqr(1.0-z)*mi2-sqr(z)*mj2;
  }
  double kt2=pipj*z*(1.0-z);
  if (fla.IsFermion()) kt2=pipj*(flc.IsVector()?(1.0-z):z);
  else if (flc.IsFermion()) kt2=pipj;
  return kt2;
}

double Kinematics_FI::GetY(const double &Q2,const double &kt2,const double &z,
			   const double &mi2,const double &mj2,const double &ma2,
			   const ATOOLS::Flavour &fla,const ATOOLS::Flavour &flc,
			   const bool force) const
{
  if (!force && (z<=0.0 || z>=1.0 || Q2>=mi2+mj2+ma2)) return -1.0;
  if (m_evolscheme==0) {
    return 1.0/(1.0-(kt2/(z*(1.0-z))+mi2*(1.0-z)/z+mj2*z/(1.0-z))/(Q2-ma2-mi2-mj2));
  }
  if (fla.IsFermion()) {
    if (flc.IsFermion()) return 1.0/(1.0-kt2/z/(Q2-ma2-mi2-mj2));
    return 1.0/(1.0-kt2/(1.0-z)/(Q2-ma2-mi2-mj2));
  }
  if (flc.IsFermion()) return 1.0/(1.0-kt2/(Q2-ma2-mi2-mj2));
  return 1.0/(1.0-kt2/(z*(1.0-z))/(Q2-ma2-mi2-mj2));
}

int Kinematics_FI::MakeKinematics
(Parton *const split,const double & mi2,const double & mj2,
 const ATOOLS::Flavour &flj,Parton *&pc,const int mode)
{ 
  Parton * spect = split->GetSpect();
  Vec4D p1 = split->Momentum(), p2 = spect->Momentum(), rp1 = p1;

  double ma2 = spect->Mass2(), mij2 = split->Mass2(); 
  bool nospec(false);
  if (mode && split->GetPrev() && split->GetPrev()->KScheme()) {
    mij2=p1.Abs2();
    p2=split->GetPrev()->FixSpec();
    ma2=0.0;
    nospec=true;
  }
  
  double Q2((p1-p2).Abs2());
  double y=GetY(Q2,split->KtTest(),split->ZTest(),mi2,mj2,ma2,
		split->GetFlavour(),flj,1);
  Kin_Args fi(1.0-y,split->ZTest(),split->Phi(),8);
  if (ConstructFIDipole(mi2,mj2,mij2,ma2,p1,p2,fi)<0) return -1;

  split->SetMomentum(fi.m_pi);
  if (mi2) SetFixVec(split,fi.m_pi,fi,2);
  if (!nospec) {
    if (ma2) SetFixVec(spect,fi.m_pk,fi,2);
    spect->SetMomentum(fi.m_pk);
  }
  else if (!IsEqual(fi.m_pk,p2,1.0e-3))
    msg_Error()<<METHOD<<"(): Error in EW splitting ( y = "<<y
	       <<" ).\n  Shifted p_k = "
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

double Kinematics_IF::GetKT2(const double &Q2,const double &y,const double &z,
			     const double &ma2,const double &mi2,const double &mk2,
			     const ATOOLS::Flavour &flb,const ATOOLS::Flavour &flc) const
{
  double pipj=(Q2-ma2-mi2-mk2)*y/z;
  if (m_evolscheme==0) {
    return -pipj*(1.0-z)-mi2-sqr(1.0-z)*ma2;
  }
  double kt2=-pipj*(1.0-z);
  if (flc.IsFermion()) kt2=-pipj;
  return kt2;
}

double Kinematics_IF::GetY(const double &Q2,const double &kt2,const double &z,
			   const double &ma2,const double &mi2,const double &mk2,
			   const ATOOLS::Flavour &flb,const ATOOLS::Flavour &flc,
			   const bool force) const
{
  if (!force && (z<=0.0 || z>=1.0 || Q2>=ma2+mi2+mk2)) return -1.0;
  if (m_evolscheme==0) {
    return -z/(Q2-ma2-mi2-mk2)*((kt2+mi2)/(1.0-z)+(1.0-z)*ma2);
  }
  if (flc.IsFermion()) return -z/(Q2-ma2-mi2-mk2)*kt2;
  return -z/(Q2-ma2-mi2-mk2)*kt2/(1.0-z);
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
  if (spect->KScheme()) mk2=p2.Abs2();
  bool nospec(false);
  if (mode && split->GetPrev() && split->GetPrev()->KScheme()) {
    mai2=p1.Abs2();
    p2=split->GetPrev()->FixSpec();
    mk2=0.0;
    nospec=true;
  }

  double y=GetY((p2-p1).Abs2(),split->KtTest(),split->ZTest(),ma2,mi2,mk2,
		split->GetFlavour(),fli,1);
  Kin_Args ifp(y,split->ZTest(),split->Phi(),split->Kin());
  if (dabs(y-split->ZTest())<Kin_Args::s_uxeps) ifp.m_mode=1;
  if (ConstructIFDipole(ma2,mi2,mai2,mk2,mb2,p1,p2,b->Momentum(),ifp)<0) return -1;

  split->SetLT(ifp.m_lam);
  split->SetMomentum(ifp.m_pi);
  if (ma2) SetFixVec(split,ifp.m_pi,ifp,1);
  if (nospec) {
    Vec4D ps(ifp.m_lam*spect->Momentum());
    if (mk2) SetFixVec(spect,ps,ifp,1);
    spect->SetMomentum(ps);
  }
  else {
    if (mk2) SetFixVec(spect,ifp.m_pk,ifp,1);
    spect->SetMomentum(ifp.m_pk);
  }
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

double Kinematics_II::GetKT2(const double &Q2,const double &y,const double &z,
			     const double &ma2,const double &mi2,const double &mb2,
			     const ATOOLS::Flavour &flb,const ATOOLS::Flavour &flc) const
{
  double pipj=(Q2-ma2-mi2-mb2)*y/z;
  if (m_evolscheme==0) {
    return pipj*(1.0-z)-mi2-sqr(1.0-z)*ma2;
  }
  double kt2=pipj*(1.0-z);
  if (flc.IsFermion()) kt2=pipj;
  return kt2;
}

double Kinematics_II::GetY(const double &Q2,const double &kt2,const double &z,
			   const double &ma2,const double &mi2,const double &mb2,
			   const ATOOLS::Flavour &flb,const ATOOLS::Flavour &flc,
			   const bool force) const
{
  if (!force && (z<=0.0 || z>=1.0 || Q2<=ma2+mi2+mb2)) return -1.0;
  if (m_evolscheme==0) {
    return z/(Q2-ma2-mb2-mi2)*((kt2+mi2)/(1.0-z)+(1.0-z)*ma2);
  }
  if (flc.IsFermion()) return z/(Q2-ma2-mb2-mi2)*kt2;
  return z/(Q2-ma2-mb2-mi2)*kt2/(1.0-z);
}

int Kinematics_II::MakeKinematics
(Parton *const split,const double & ma2,const double & mi2,
 const ATOOLS::Flavour &newfl,Parton *&pc,const int mode)
{
  Parton * spect = split->GetSpect();
  Vec4D p1 = split->Momentum(), p2 = spect->Momentum();
  
  double mai2 = split->Mass2(), mb2 = spect->Mass2();
  bool nospec(false);
  if (mode && split->GetPrev() && split->GetPrev()->KScheme()) {
    mai2=p1.Abs2();
    // temporary !!!
    p2[1]=p2[2]=0.0;
    p2[0]=dabs(p2[3]);
    // end temporary
    mb2=0.0;
    nospec=true;
  }

  double y=GetY((p1+p2).Abs2(),split->KtTest(),split->ZTest(),ma2,mi2,mb2,
		split->GetFlavour(),newfl,1);
  Kin_Args ii(y,split->ZTest(),split->Phi(),split->Kin());
  if (ConstructIIDipole(ma2,mi2,mai2,mb2,p1,p2,ii)<0) return -1;

  split->SetLT(ii.m_lam);
  split->SetMomentum(ii.m_pi);
  if (ma2) SetFixVec(split,ii.m_pi,ii,3);
  if (nospec) {
    Vec4D ps(ii.m_lam*spect->Momentum());
    if (mb2) SetFixVec(spect,ps,ii,3);
    spect->SetMomentum(ps);
  }
  else {
    if (mb2) SetFixVec(spect,ii.m_pk,ii,3);
    spect->SetMomentum(ii.m_pk);
  }
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

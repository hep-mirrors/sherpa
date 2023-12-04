#include "PHASIC++/Channels/Antenna_Kinematics.H"

#include "PHASIC++/Channels/Transverse_Kinematics.H"
#include "ATOOLS/Math/ZAlign.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace PHASIC;
using namespace ATOOLS;

int PHASIC::ClusterAntenna
(Ant_Args &ffp,const size_t &i,const size_t &j,
 const size_t &k,const double &mij2)
{
  Vec4D pi(ffp.m_p[i]), pj(ffp.m_p[j]), pk(ffp.m_p[k]), K;
  for (size_t l(0);l<ffp.m_p.size();++l)
    if (ffp.m_b[l]&2) {
      if (l!=j) K+=ffp.m_p[l];
      if (l==i) ffp.m_mode|=2;
      if (l==k) ffp.m_mode|=4;
    }
  if (ffp.m_mode&2) { K=-K; }
  Vec4D n(K+pj);
  double z((pi*n)/((pi+pj)*n));
  ffp.m_pijt=pi/z;
  ffp.m_Kt=n-(1.0-z)*ffp.m_pijt;
  ffp.m_z=(ffp.m_mode&1)?1.0/z:z;
  ffp.m_y=dabs((pi*pj)/(pi*n));
  ffp.m_phi=ComputePhi(pi,n,pj);
  ffp.m_pk=pk;
  Vec4D Kt(K-ffp.m_pijt+pi+pj);
  Poincare oldcm, newcm;
  if (ffp.m_mode&2) {
    oldcm=Poincare(-Kt);
    newcm=Poincare(-K);
    oldcm.Boost(ffp.m_Kt);
    newcm.BoostBack(ffp.m_Kt);
    oldcm.Boost(ffp.m_pijt);
    newcm.BoostBack(ffp.m_pijt);
  }
  else {
    newcm=Poincare(Kt);
    oldcm=Poincare(K);
  }
  for (size_t l(0);l<ffp.m_p.size();++l)
    if (ffp.m_b[l]&1) {
      oldcm.Boost(ffp.m_p[l]);
      newcm.BoostBack(ffp.m_p[l]);
    }
  if (ffp.m_mode&4) {
    oldcm.Boost(ffp.m_pk);
    newcm.BoostBack(ffp.m_pk);
  }
  ffp.m_p[i]=ffp.m_pijt;
  ffp.m_p.erase(ffp.m_p.begin()+j);
  ffp.m_stat=1;
  return 1;
}

int PHASIC::ConstructAntenna
(Ant_Args &ffp,const size_t &ij,const size_t &k,
 const double &mi2,const double &mj2)
{
  Vec4D pij(ffp.m_p[ij]), pk(ffp.m_p[k]), Kt;
  for (size_t i(0);i<ffp.m_p.size();++i)
    if (ffp.m_b[i]&2) {
      Kt+=ffp.m_p[i];
      if (i==ij) ffp.m_mode|=2;
      if (i==k) ffp.m_mode|=4;
    }
  if (ffp.m_mode&2) { Kt=-Kt; }
  ffp.m_Kt=Kt;
  ffp.m_pijt=pij;
  double v(Kt[0]>0.0?ffp.m_y:-ffp.m_y);
  double z((ffp.m_mode&1)?1.0/ffp.m_z:ffp.m_z);
  double Q2(2.*pij*Kt), Kt2(Kt.Abs2());
  double kap(Kt.Abs2()/Q2), mit2(pij.Abs2());
  double n2((1.0-z)*Q2+Kt2+mit2-mi2);
  double vo(sqrt(1.0-mit2*Kt2/sqr(Q2/2.0)));
  double gam(Q2*(1+vo)/2.0);
  double muit2(mit2/Q2), mui2(mi2/Q2), muj2(mj2/Q2);
  double muitb2(mit2/gam), muib2(mi2/gam);
  double mujb2(mj2/gam), kapb(Kt2/gam);
  double zb((z+2.0*mui2)/(1.0+vo+2.0*muit2));
  if (zb*zb<muib2/(1.0+muitb2)*(1.0+kapb)) {
    msg_IODebugging()<<METHOD<<"(): Invalid kinematics."<<std::endl;
    return -1;
  }
  zb+=sqrt(zb*zb-muib2/(1.0+muitb2)*(1.0+kapb));
  ffp.m_pk=pk;
  ffp.m_pi=zb*pij+(mui2/zb-muit2*zb)/vo*(Kt-kapb*pij);
  Vec4D n(Kt+pij-ffp.m_pi), n_perp(ffp.m_pk);
  LN_Pair ln(GetLN(pij,Kt));
  n_perp-=((pk*ln.m_n)*ln.m_l+(pk*ln.m_l)*ln.m_n)/(gam/2.0);
  if (n_perp.PSpat2()<=rpa->gen.SqrtAccu()) {
    msg_IODebugging()<<METHOD<<"(): Set fixed n_perp\n";
    n_perp=LT(pij,Kt,Vec4D::ZVEC);
    if (n_perp.PSpat2()<=rpa->gen.SqrtAccu())
      n_perp=LT(pij,Kt,Vec4D::YVEC);
    if (n_perp.PSpat2()<=rpa->gen.SqrtAccu())
      n_perp=LT(pij,Kt,Vec4D::XVEC);
  }
  n_perp*=1.0/sqrt(dabs(n_perp.Abs2()));
  Vec4D l_perp(LT(n,ffp.m_pi,n_perp));
  l_perp*=1.0/sqrt(dabs(l_perp.Abs2()));
  double zeta(2/(1+vo)*(1.0-z+kap+muit2-mui2)/(1.0-zb+kapb));
  double vb(z*v*2/(1+vo)-muib2/zb*(1.0-zb+(mujb2-kapb*(1-zeta))/zeta));
  vb/=zb-muib2/zb*(1.0-zb+kapb)/zeta;
  double ktt(vb/zeta*((zeta-vb)*(1.0-zb)-(vb+(1-zeta))*kapb+mujb2)*gam-mj2);
  if (ktt<0.0) {
    msg_IODebugging()<<METHOD<<"(): Invalid kinematics."<<std::endl;
    return -1;
  }
  ktt=sqrt(ktt);
  ffp.m_pj=ktt*(cos(ffp.m_phi)*n_perp+sin(ffp.m_phi)*l_perp);
  ffp.m_pj+=(1.0-z+muit2-mui2+muj2)/(vo*zeta)*(pij-muitb2*Kt);
  ffp.m_pj+=vb*(1.0+vo)/(2.0*vo)*
    ((Kt-kapb*pij)-(1.0-zb+kapb)/zeta*(pij-muitb2*Kt));
#ifdef DEBUG__Kinematics
  {
    DEBUG_FUNC("pi-nb frame");
    DEBUG_VAR((ln.m_l-ln.m_n).Abs2());
    DEBUG_VAR(z<<" "<<v<<" "<<kap<<" "<<-v*(1.+kap));
    Poincare cms(ln.m_l-ln.m_n);
    Vec4D pi(ffp.m_pi), pj(ffp.m_pj), pk(ffp.m_pk);
    cms.Boost(pi);
    cms.Boost(pj);
    cms.Boost(pk);
    Poincare zax(pi,Vec4D::ZVEC);
    zax.Rotate(pi);
    zax.Rotate(pj);
    zax.Rotate(pk);
    double ei(sqrt(-z*pij*Kt/2.));
    DEBUG_VAR(pi[0]<<" "<<ei<<" "<<ei/pi[0]-1.);
    double ct(1.+2.*v*z/(1.-z-v*(1.+kap)));
    DEBUG_VAR(pj.CosTheta(pi)<<" "<<ct<<" "<<ct/pj.CosTheta(pi)-1.);
    double ej(-2.*v*ei/(1.-ct));
    DEBUG_VAR(pj[0]<<" "<<ej<<" "<<ej/pj[0]-1.);
  }
#endif
#ifdef DEBUG__Kinematics
  {
    DEBUG_FUNC("n frame");
    double n2((1-z+kap)*Q2);
    DEBUG_VAR(n.Abs2()<<" "<<n2<<" "<<n2/n.Abs2()-1.);
    Poincare cms(-n);
    Vec4D pi(ffp.m_pi), pj(ffp.m_pj), pk(ffp.m_pk);
    Vec4D K(Kt+pij-ffp.m_pi-ffp.m_pj), nn(n);
    Vec4D KKt(K+Kt);
    double kkt2((4.*kap+v*(1.-z))*Q2);
    DEBUG_VAR(KKt.Abs2()<<" "<<kkt2<<" "<<kkt2/KKt.Abs2()-1.);
    double pjKt(Q2/2.*(1.-z)*(1.-v));
    DEBUG_VAR(pj*Kt<<" "<<pjKt<<" "<<pjKt/(pj*Kt)-1.);
    double x((pij*n)/(Kt*n)), y(1./(1.-z+2.*kap));
    DEBUG_VAR(x<<" "<<y<<" "<<x/y-1.);
    cms.Boost(pi);
    cms.Boost(pj);
    cms.Boost(pk);
    cms.Boost(K);
    cms.Boost(nn);
    Poincare zax(pi,Vec4D::ZVEC);
    zax.Rotate(pi);
    zax.Rotate(pj);
    zax.Rotate(pk);
    zax.Rotate(K);
    double KtKt(Q2*(v*(1.-z)+2.*kap));
    DEBUG_VAR(K*Kt<<" "<<KtKt<<" "<<KtKt/(K*Kt)-1.);
    double eK(Q2/2.*(1.-z+2.*kap)/sqrt(Q2*(1.-z+kap)));
    DEBUG_VAR(K[0]<<" "<<eK<<" "<<eK/K[0]-1.);
    double ei(z/sqrt((1.-z+kap)/Q2)/2.);
    DEBUG_VAR(pi[0]<<" "<<ei<<" "<<ei/pi[0]-1.);
    double ej((1.-z)/sqrt((1.-z+kap)/Q2)/2.);
    DEBUG_VAR(pj[0]<<" "<<ej<<" "<<ej/pj[0]-1.);
    double ctj(1.-2.*v*(1.-z+kap)/(1.-z));
    DEBUG_VAR(pj.CosTheta(pi)<<" "<<ctj<<" "<<ctj/pj.CosTheta(pi)-1.);
    double ejoeK(-(1.-z)/(1.-z+2.*kap));
    DEBUG_VAR(pj[0]/-K[0]<<" "<<ejoeK<<" "<<ejoeK/(pj[0]/-K[0])-1.);
    double eioeK(-z/(1.-z+2.*kap));
    DEBUG_VAR(pi[0]/-K[0]<<" "<<eioeK<<" "<<eioeK/(pi[0]/-K[0])-1.);
    double cp(CosPhi(pk,pj,nn,pi));
    DEBUG_VAR(pk.CosDPhi(pj)<<" "<<cp<<" "<<cp/pk.CosDPhi(pj)-1.);
    double ctk(pk.CosTheta(pi));
    double Aikji((pi[0]/pk[0]*(1.-ctj)+1.-ctj*ctk)/(1.-ctk));
    double Bikji(sqrt((1.-ctj*ctj)*(1.-ctk*ctk))/(1.-ctk));
    double Wikji(1./(Aikji-Bikji*cp));
    double tWikji(pj[0]/pi[0]*(pi*pk)/((pi+pk)*pj));
    DEBUG_VAR(tWikji<<" "<<Wikji<<" "<<Wikji/tWikji-1.);
    double twikji(1./(pi*pj)*(pi*pk)/((pi+pk)*pj));
    double wikji(1./(Q2*v)*2.*Wikji/(1.-z));
    DEBUG_VAR(wikji<<" "<<twikji<<" "<<twikji/wikji-1.);
    DEBUG_VAR(ej*ej*(1.-ctj));
  }
#endif
#ifdef DEBUG__Kinematics
  {
    DEBUG_FUNC("Kt frame "<<Kt);
    Poincare cms(-Kt);
    Vec4D pi(ffp.m_pi), pj(ffp.m_pj), pk(ffp.m_pk), nn(n);
    cms.Boost(pi);
    cms.Boost(pj);
    cms.Boost(pk);
    cms.Boost(nn);
    Poincare zax(pi,Vec4D::ZVEC);
    zax.Rotate(pi);
    zax.Rotate(pj);
    zax.Rotate(pk);
    zax.Rotate(nn);
    double ei(z/sqrt(kap/Q2)/2.);
    DEBUG_VAR(pi[0]<<" "<<ei<<" "<<ei/pi[0]-1.);
    double ej((1.-z)*(1.-v)/sqrt(kap/Q2)/2.);
    DEBUG_VAR(pj[0]<<" "<<ej<<" "<<ej/pj[0]-1.);
    double ctj(1.-2.*v/(1.-v)*kap/(1.-z));
    DEBUG_VAR(pj.CosTheta(pi)<<" "<<ctj<<" "<<ctj/pj.CosTheta(pi)-1.);
    double cp(CosPhi(pk,pj,nn,pi));
    DEBUG_VAR(pk.CosDPhi(pj)<<" "<<cp<<" "<<cp/pk.CosDPhi(pj)-1.);
    double ctk(pk.CosTheta(pi));
    double Aikji((pi[0]/pk[0]*(1.-ctj)+1.-ctj*ctk)/(1.-ctk));
    double Bikji(sqrt((1.-ctj*ctj)*(1.-ctk*ctk))/(1.-ctk));
    double Wikji(1./(Aikji-Bikji*cp));
    double tWikji(pj[0]/pi[0]*(pi*pk)/((pi+pk)*pj));
    DEBUG_VAR(tWikji<<" "<<Wikji<<" "<<Wikji/tWikji-1.);
    double twikji(1./(pi*pj)*(pi*pk)/((pi+pk)*pj));
    double wikji(1./(Q2*v*(1.-v))*2.*Wikji/(1.-z));
    DEBUG_VAR(wikji<<" "<<twikji<<" "<<twikji/wikji-1.);
  }
#endif
  Vec4D K(Kt+pij-ffp.m_pi-ffp.m_pj);
  if (ffp.m_mode&2) {
    Poincare oldcm(-Kt), newcm(-K);
    for (size_t i(0);i<ffp.m_p.size();++i)
      if (ffp.m_b[i]&1) {
	newcm.Boost(ffp.m_p[i]);
	oldcm.BoostBack(ffp.m_p[i]);
      }
    newcm.Boost(ffp.m_pi);
    oldcm.BoostBack(ffp.m_pi);
    newcm.Boost(ffp.m_pj);
    oldcm.BoostBack(ffp.m_pj);
    if (ffp.m_mode&4) {
      newcm.Boost(ffp.m_pk);
      oldcm.BoostBack(ffp.m_pk);
    }
    newcm.Boost(ffp.m_Kt);
    oldcm.BoostBack(ffp.m_Kt);
    newcm.Boost(ffp.m_pijt);
    oldcm.BoostBack(ffp.m_pijt);
  }
  else {
    Poincare oldcm(Kt), newcm(K);
    for (size_t i(0);i<ffp.m_p.size();++i)
      if (ffp.m_b[i]&1) {
	oldcm.Boost(ffp.m_p[i]);
	newcm.BoostBack(ffp.m_p[i]);
      }
    if (ffp.m_mode&4) {
      oldcm.Boost(ffp.m_pk);
      newcm.BoostBack(ffp.m_pk);
    }
  }
#ifdef DEBUG__Kinematics
  {
    DEBUG_FUNC("K frame (Kt = "<<Kt<<", kap = "<<kap<<", v = "<<v<<")");
    Vec4D pi(ffp.m_pi), pj(ffp.m_pj), pk(ffp.m_pk);
    double eji((1.-z)/z/(1.-v));
    DEBUG_VAR(pj[0]/pi[0]<<" "<<eji<<" "<<eji/(pj[0]/pi[0])-1.);
    double ei(-z*pij*Kt/sqrt(Kt.Abs2())*(1.-v));
    DEBUG_VAR(pi[0]<<" "<<ei<<" "<<ei/pi[0]-1.);
    double eioeK(-z*(1.-v)/kap/2.);
    DEBUG_VAR(pi[0]/-Kt[0]<<" "<<eioeK<<" "<<eioeK/(pi[0]/-Kt[0])-1.);
    double ej(-pij*Kt/sqrt(Kt.Abs2())*(1.-z));
    DEBUG_VAR(pj[0]<<" "<<ej<<" "<<ej/pj[0]-1.);
    double ejoeK(-(1.-z)/kap/2.);
    DEBUG_VAR(pj[0]/-Kt[0]<<" "<<ejoeK<<" "<<ejoeK/(pj[0]/-Kt[0])-1.);
    double ct(1.-2.*kap/(1.-z)*v/(1.-v));
    DEBUG_VAR(pj.CosTheta(pi)<<" "<<ct<<" "<<ct/pj.CosTheta(pi)-1.);
    double t(pij*Kt*(1.-z)*v/(1.-v));
    DEBUG_VAR(pj[0]*pj[0]*(1.-pj.CosTheta(pi))<<" "<<t<<" "<<t/(pj[0]*pj[0]*(1.-pj.CosTheta(pi)))-1.);
    DEBUG_VAR(ffp.m_pk*K<<" "<<pk*Kt<<" "<<(ffp.m_pk*K)/(pk*Kt)-1.);
    DEBUG_VAR(ej*ej*(1.-ct));
  }
#endif
  return 1;
}

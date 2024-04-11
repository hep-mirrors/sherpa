#include "PHASIC++/Channels/CSS_Kinematics.H"

#include "PHASIC++/Channels/Transverse_Kinematics.H"
#include "ATOOLS/Math/ZAlign.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace PHASIC;
using namespace ATOOLS;

double Kin_Args::s_uxeps=1.0e-3;

Kin_Args PHASIC::ClusterFFDipole
(const double &mi2,const double &mj2,const double &mij2,
 const double &mk2,const Vec4D &pi,const Vec4D &pj,
 const Vec4D &pk,const int mode)
{
  double pipj(pi*pj), pipk(pi*pk), pjpk(pj*pk);
  double yijk(pipj/(pipj+pipk+pjpk)), zi(pipk/(pipk+pjpk));
  double sij((pi+pj).Abs2()), Q2((pi+pj+pk).Abs2());
  double po(sqr(Q2-mij2-mk2)-4.0*mij2*mk2);
  double pn(sqr(Q2-sij-mk2)-4.0*sij*mk2);
  if (pn<0.0 || po<0.0) {
    msg_IODebugging()<<METHOD<<"(): Invalid kinematics."<<std::endl;
    return Kin_Args();
  }
  Vec4D Q(pi+pj+pk);
  Kin_Args res(yijk,zi,0.0,(mode&4)?1:0,1);
  res.m_pk=sqrt(po/pn)*(pk-(Q*pk)/Q2*Q)+(Q2+mk2-mij2)/(2.0*Q2)*Q;
  res.m_pi=Q-res.m_pk;
  if (res.m_pi[0]<0.0 || res.m_pk[0]*pk[0]<0.0) {
    msg_IODebugging()<<METHOD<<"(): Invalid kinematics."<<std::endl;
    return Kin_Args();
  }
  if (mode&1) res.m_phi=ComputePhi(res.m_pi,res.m_pk,pi);
  pn=Sign(Q2-sij-mk2)*sqrt(pn);
  res.m_nb=((Q2-sij+mk2+pn)/2.*pk-mk2*Q)/pn;
  return res;
}

int PHASIC::ConstructFFDipole
(const double &mi2,const double &mj2,const double &mij2,
 const double &mk2,const Vec4D &pij,const Vec4D &pk,Kin_Args &ffp)
{
  Vec4D n_perp(0.0,cross(Vec3D(pij),Vec3D(pk)));
  if (n_perp.PSpat2()<=rpa->gen.SqrtAccu()) {
    msg_IODebugging()<<"Set fixed n_perp\n";
    n_perp=Vec4D(0.0,1.0,1.0,0.0);
    Poincare zrot(pij,Vec4D::ZVEC);
    zrot.RotateBack(n_perp);
  }
  n_perp*=1.0/n_perp.PSpat();
  Vec4D l_perp(LT(pij,pk,n_perp));
  l_perp*=1.0/sqrt(dabs(l_perp.Abs2()));
  Vec4D Q(pij+pk);
  double Q2(Q.Abs2()), sij(ffp.m_y*(Q2-mk2)+(1.0-ffp.m_y)*(mi2+mj2));
  double po(sqr(Q2-mij2-mk2)-4.0*mij2*mk2);
  double pn(sqr(Q2-sij-mk2)-4.0*sij*mk2);
  if (po<0.0 || pn<0.0) {
    msg_IODebugging()<<METHOD<<"(): Kinematics does not fit."<<std::endl;
    return -1;
  }
  double ecm(Q2-sij-mk2);
  po=Sign(ecm)*sqrt(po);
  pn=Sign(ecm)*sqrt(pn);
  double gam(0.5*(ecm+pn));
  double zt(ecm/pn*(ffp.m_z-mk2/gam*(sij+mi2-mj2)/ecm));
  double ktt(sij*zt*(1.0-zt)-(1.0-zt)*mi2-zt*mj2);
  if (ktt<0.0) {
    msg_IODebugging()<<METHOD<<"(): Invalid kinematics."<<std::endl;
    return -1;
  }
  ktt=sqrt(ktt);
  ffp.m_pk=pn/po*(pk-(Q2-mij2+mk2)/(2.0*Q2)*Q)+(Q2-sij+mk2)/(2.0*Q2)*Q;
  ffp.m_pj=Q-ffp.m_pk;
  ffp.m_pi=ktt*sin(ffp.m_phi)*l_perp;
  ffp.m_pi+=ktt*cos(ffp.m_phi)*n_perp+zt/pn*(gam*ffp.m_pj-sij*ffp.m_pk)+
    (mi2+ktt*ktt)/zt/pn*(ffp.m_pk-mk2/gam*ffp.m_pj);
  ffp.m_pj=Q-ffp.m_pk-ffp.m_pi;
  ffp.m_nb=(((ecm+pn)/2.+mk2)*ffp.m_pk-mk2*Q)/pn;
  return 1;
}

Kin_Args PHASIC::ClusterFIDipole
(const double &mi2,const double &mj2,const double &mij2,
 const double &ma2,const Vec4D &pi,const Vec4D &pj,
 const Vec4D &pa,const int mode)
{
  Vec4D Q(pa-pi-pj), Ql(Vec4D(Q[0],0.0,0.0,Q[3]));
  double pipj(pi*pj), pipa(pi*pa), pjpa(pj*pa);
  double xija((pipa+pjpa-pipj)/(pipa+pjpa)), zi(pipa/(pipa+pjpa));
  double sij((pi+pj).Abs2()), Q2(Q.Abs2()), kt2(Q.PPerp2());
  double po(sqr(Q2-mij2-ma2)-4.0*ma2*(mij2+kt2));
  double pn(sqr(Q2-sij-ma2)-4.0*ma2*(sij+kt2));
  if (!(mode&8)) xija/=(Q2-mi2-mj2-ma2)/(Q2-mij2-ma2);
  if (pn<0.0 || po<0.0) {
    msg_IODebugging()<<METHOD<<"(): Invalid kinematics."<<std::endl;
    return Kin_Args();
  }
  Kin_Args res(1.0-xija,zi,0.0,(mode&4)?1:0,1);
  res.m_pk=sqrt(po/pn)*(pa-(Q*pa)/(Q2+kt2)*Ql)
    +(Q2+ma2-mij2)/(2.0*(Q2+kt2))*Ql;
  res.m_pi=res.m_pk-Q;
  if (res.m_pi[0]<0.0 || res.m_pk[0]<0.0) {
    msg_IODebugging()<<METHOD<<"(): Invalid kinematics."<<std::endl;
    return Kin_Args();
  }
  if (mode&1) res.m_phi=ComputePhi(pi+pj,pa,pi);
  return res;
}

int PHASIC::ConstructFIDipole
(const double &mi2,const double &mj2,const double &mij2,
 const double &ma2,const Vec4D &pij,const Vec4D &pa,Kin_Args &fip)
{
  Vec4D Q(pa-pij), Ql(Vec4D(Q[0],0.0,0.0,Q[3]));
  double Q2(Q.Abs2()), kt2(Q.PPerp2()), yt(1.0-fip.m_y);
  if (fip.m_mode&8) yt=(1.0-yt)/yt;
  else yt=((Q2-mij2-ma2)/(Q2-mi2-mj2-ma2)-yt)/yt;
  double sij(-yt*(Q2-ma2)+(1.0+yt)*(mi2+mj2));
  double po(sqr(Q2-mij2-ma2)-4.0*ma2*(mij2+kt2));
  double ecm(Q2-sij-ma2), pn(sqr(ecm)-4.0*ma2*(sij+kt2));
  if (pn<0.0 || po<0.0) {
    msg_IODebugging()<<METHOD<<"(): Invalid kinematics."<<std::endl;
    return -1;
  }
  pn=sqrt(pn);
  po=sqrt(po);
  fip.m_pk=pn/po*(pa-(Q*pa)/(Q2+kt2)*Ql)
    +(Q2+ma2-sij)/(2.0*(Q2+kt2))*Ql;
  fip.m_pi=fip.m_pj=fip.m_pk-Q;
  Vec4D x(fip.m_pi), y(fip.m_pk);
  Vec4D n_perp(0.0,cross(Vec3D(fip.m_pi),Vec3D(fip.m_pk)));
  if (n_perp.PSpat2()<=rpa->gen.SqrtAccu()) {
    msg_IODebugging()<<"Set fixed n_perp\n";
    n_perp=Vec4D(0.0,1.0,1.0,0.0);
    Poincare zrot(fip.m_pi,Vec4D::ZVEC);
    zrot.RotateBack(n_perp);
  }
  n_perp*=1.0/n_perp.PSpat();
  Vec4D l_perp(LT(fip.m_pi,fip.m_pk,n_perp));
  l_perp*=1.0/sqrt(dabs(l_perp.Abs2()));
  double pnn(Sign(ecm)*sqrt(sqr(ecm)-4.0*sij*ma2)), gam(0.5*(ecm+pnn));
  double zt(ecm/pnn*(fip.m_z-ma2/gam*(sij+mi2-mj2)/ecm));
  double ktt(sij*zt*(1.0-zt)-(1.0-zt)*mi2-zt*mj2);
  if (ktt<0.0 || gam==0.0) {
    msg_IODebugging()<<METHOD<<"(): Invalid kinematics."<<std::endl;
    return -1;
  }
  ktt=sqrt(ktt);
  fip.m_pi=ktt*sin(fip.m_phi)*l_perp;
  fip.m_pi+=ktt*cos(fip.m_phi)*n_perp+zt/pnn*(gam*fip.m_pj+sij*fip.m_pk)-
    (mi2+ktt*ktt)/zt/pnn*(fip.m_pk+ma2/gam*fip.m_pj);
  fip.m_pj=fip.m_pk-Q-fip.m_pi;
  return 1;
}

Kin_Args PHASIC::ClusterIFDipole
(const double &ma2,const double &mj2,const double &maj2,
 const double &mk2,const double &mb2,const Vec4D &pa,
 const Vec4D &pj,const Vec4D &pk,const Vec4D &pb,const int mode)
{
  Vec4D Q(pa-pj-pk), Ql(Vec4D(Q[0],0.0,0.0,Q[3]));
  double papj(pa*pj), papk(pa*pk), pjpk(pj*pk), Q2(Q.Abs2());
  double xjka((papj+papk-pjpk)/(papj+papk)), uj(papj/(papj+papk));
  Kin_Args res(uj,xjka,0.0,(mode&4)?1:0,1);
  if (dabs(xjka-uj)<Kin_Args::s_uxeps ||
      dabs(Q2)<Kin_Args::s_uxeps) res.m_mode=1;
  if (res.m_mode==1) {
    double sjk((pj+pk).Abs2()), kt2(Q.PPerp2());
    double po(sqr(Q2-mk2-maj2)-4.0*maj2*(mk2+kt2));
    double pn(sqr(Q2-sjk-ma2)-4.0*ma2*(sjk+kt2));
    if (pn<0.0 || po<0.0) {
      msg_IODebugging()<<METHOD<<"(): Invalid kinematics."<<std::endl;
      return Kin_Args();
    }
    res.m_pi=sqrt(po/pn)*(pa-(Q2+ma2-sjk)/(2.0*(Q2+kt2))*Ql)
      +(Q2+maj2-mk2)/(2.0*(Q2+kt2))*Ql;
    res.m_pk=res.m_pi-Q;
    if (res.m_pi[0]<0.0 || res.m_pk[0]<0.0) {
      msg_IODebugging()<<METHOD<<"(): Invalid kinematics."<<std::endl;
      return Kin_Args();
    }
    if (mode&1) res.m_phi=ComputePhi(pj+pk,pa,pj);
  }
  else {
    double saj((pa-pj).Abs2());
    double po(sqr(Q2-maj2-mk2)-4.0*maj2*mk2);
    double pn(sqr(Q2-saj-mk2)-4.0*saj*mk2);
    if (pn<0.0 || po<0.0) {
      msg_IODebugging()<<METHOD<<"(): Invalid kinematics."<<std::endl;
      return -1;
    }
    res.m_pk=sqrt(po/pn)*(pk-(Q*pk)/Q2*Q)-(Q2+mk2-maj2)/(2.0*Q2)*Q;
    res.m_pi=Q+res.m_pk;
    if (dabs(res.m_pk.Abs2()-mk2)>rpa->gen.SqrtAccu() ||
	dabs(res.m_pi.Abs2()-maj2)>rpa->gen.SqrtAccu()) {
      msg_IODebugging()<<METHOD<<"(): Unstable point. Using scheme 1.\n";
      return ClusterIFDipole(ma2,mj2,maj2,mk2,mb2,pa,pj,pk,pb,mode|4);
    }
    ZAlign lam(res.m_pi,pb,maj2,mb2);
    if (lam.Status()<0) {
      msg_IODebugging()<<METHOD<<"(): Invalid kinematics."<<std::endl;
      return Kin_Args();
    }
    res.m_lam=lam;
    res.m_pi=res.m_lam*res.m_pi;
    res.m_pk=res.m_lam*res.m_pk;
    if (mode&1) res.m_phi=ComputePhi(res.m_pi,res.m_pk,res.m_lam*pj);
  }
  return res;
}

int PHASIC::ConstructIFDipole
(const double &ma2,const double &mj2,const double &maj2,
 const double &mkt2,const double &mb2,const Vec4D &paj,
 const Vec4D &pk,const Vec4D &pb,Kin_Args &ifp)
{
  if (ifp.m_mode==1) {
    Vec4D Q(paj-pk), Ql(Vec4D(Q[0],0.0,0.0,Q[3]));
    double Q2(Q.Abs2()), kt2(Q.PPerp2()), yt((1.0-ifp.m_z)/ifp.m_z);
    double mk2(ifp.m_mk2>=0.0?ifp.m_mk2:mkt2);
    double sjk(-yt*(Q2-ma2)+(1.0+yt)*(mj2+mk2));
    double po(sqr(Q2-mkt2-maj2)-4.0*maj2*(mkt2+kt2));
    double ecm(Q2-sjk-ma2), pn(sqr(ecm)-4.0*ma2*(sjk+kt2));
    if (pn<0.0 || po<0.0) {
      msg_IODebugging()<<METHOD<<"(): Invalid kinematics."<<std::endl;
      return -1;
    }
    pn=sqrt(pn);
    po=sqrt(po);
    Vec4D pa(pn/po*(paj-(Q2+maj2-mkt2)/(2.0*(Q2+kt2))*Ql)
	     +(Q2+ma2-sjk)/(2.0*(Q2+kt2))*Ql);
    ifp.m_pk=ifp.m_pj=(ifp.m_pi=pa)-Q;
    Vec4D n_perp(0.0,cross(Vec3D(ifp.m_pj),Vec3D(ifp.m_pi)));
    if (n_perp.PSpat2()<=rpa->gen.SqrtAccu()) {
      msg_IODebugging()<<"Set fixed n_perp\n";
      n_perp=Vec4D(0.0,1.0,1.0,0.0);
      Poincare zrot(ifp.m_pj,Vec4D::ZVEC);
      zrot.RotateBack(n_perp);
    }
    n_perp*=1.0/n_perp.PSpat();
    Vec4D l_perp(LT(ifp.m_pj,ifp.m_pi,n_perp));
    l_perp*=1.0/sqrt(dabs(l_perp.Abs2()));
    double pnn(Sign(ecm)*sqrt(sqr(ecm)-4.0*sjk*ma2)), gam(0.5*(ecm+pnn));
    double zt(ecm/pnn*(ifp.m_y-ma2/gam*(sjk+mj2-mk2)/ecm));
    double ktt(sjk*zt*(1.0-zt)-(1.0-zt)*mj2-zt*mk2);
    if (ktt<0.0) {
      msg_IODebugging()<<METHOD<<"(): Invalid kinematics."<<std::endl;
      return -1;
    }
    ktt=sqrt(ktt);
    ifp.m_pj=ktt*sin(ifp.m_phi)*l_perp;
    ifp.m_pj+=ktt*cos(ifp.m_phi)*n_perp+zt/pnn*(gam*ifp.m_pk+sjk*ifp.m_pi)-
      (sjk*(1.0-zt)+mj2-mk2)/pnn*(ifp.m_pi+ma2/gam*ifp.m_pk);
    ifp.m_pk=ifp.m_pi-Q-ifp.m_pj;
  }
  else {
    Vec4D n_perp(0.,1.,0.,0.), l_perp(0.,0.,1.,0.);
    if (cross(Vec3D(paj),Vec3D(pk)).Sqr()>0.) {
      n_perp = Vec4D(0.0,cross(Vec3D(paj),Vec3D(pk)));
      n_perp*=1.0/n_perp.PSpat();
      l_perp = Vec4D(LT(paj,pk,n_perp));
      l_perp*=1.0/sqrt(dabs(l_perp.Abs2()));
    }
    Vec4D Q(paj-pk);
    double mk2(ifp.m_mk2>=0.0?ifp.m_mk2:mkt2);
    double Q2(Q.Abs2()), po(sqr(Q2-maj2-mkt2)-4.0*maj2*mkt2);
    double yt(ifp.m_y/ifp.m_z), saj(yt*(Q2-mk2)+(1.0-yt)*(ma2+mj2));
    double ecm(Q2-saj-mk2), pn(sqr(ecm)-4.0*saj*mk2);
    if (pn<0.0 || po<0.0) {
      msg_IODebugging()<<METHOD<<"(): Invalid kinematics."<<std::endl;
      return -1;
    }
    pn=Sign(ecm)*sqrt(pn);
    po=Sign(ecm)*sqrt(po);
    ifp.m_pk=pn/po*(pk+(Q2+mkt2-maj2)/(2.0*Q2)*Q)-(Q2+mk2-saj)/(2.0*Q2)*Q;
    ifp.m_pj=ifp.m_pi=Q+ifp.m_pk;
    double gam(0.5*(ecm+pn)), rf(ifp.m_z-ifp.m_y);
    double zt(ecm/pn*((ifp.m_z-1.0)/rf-mk2/gam*(saj+mj2-ma2)/ecm));
    double ktt(saj*zt*(1.0-zt)-(1.0-zt)*mj2-zt*ma2);
    if (ktt<0.0 || gam==0.0) {
      msg_IODebugging()<<METHOD<<"(): Invalid kinematics."<<std::endl;
      return -1;
    }
    rf=dabs(rf);
    ktt=rf*sqrt(ktt);
    ifp.m_pj=ktt*sin(ifp.m_phi)*l_perp;
    ifp.m_pj+=ktt*cos(ifp.m_phi)*n_perp-
      zt*rf/pn*(gam*ifp.m_pi+saj*ifp.m_pk)+
      (mj2*rf*rf+ktt*ktt)/pn/(zt*rf)*(ifp.m_pk+mk2/gam*ifp.m_pi);
    Vec4D testpj = ifp.m_pj;
    ifp.m_pj[0]=sqrt(mj2*rf*rf+ifp.m_pj.PSpat2());
    ifp.m_pj*=1.0/rf;
    ifp.m_pi=Q+ifp.m_pj+ifp.m_pk;
    ZAlign lam(ifp.m_pi,pb,ma2,mb2);
    if (lam.Status()<0) {
      msg_IODebugging()<<METHOD<<"(): Invalid kinematics."<<std::endl;
      return -1;
    }
    ifp.m_lam=lam;
    ifp.m_pi=ifp.m_lam*ifp.m_pi;
    ifp.m_pj=ifp.m_lam*ifp.m_pj;
    ifp.m_pk=ifp.m_lam*ifp.m_pk;
  }
  return 1;
}

Kin_Args PHASIC::ClusterIIDipole
(const double &ma2,const double &mj2,const double &maj2,
 const double &mb2,const Vec4D &pa,const Vec4D &pj,
 const Vec4D &pb,const int mode)
{
  double papb(pa*pb), pjpa(pj*pa), pjpb(pj*pb);
  double xjab(1.0-(pjpa+pjpb)/papb), vj(pjpa/papb);
  if (xjab<0.0 || vj<0.0 || 1.0<xjab+vj) {
    msg_IODebugging()<<METHOD<<"(): Invalid kinematics."<<std::endl;
    return Kin_Args();
  }
  double sab((pa+pb).Abs2()), Q2((pa-pj+pb).Abs2());
  double po(sqr(Q2-maj2-mb2)-4.0*maj2*mb2);
  double pn(sqr(sab-ma2-mb2)-4.0*ma2*mb2);
  if (po<0.0 || pn<0.0) {
    msg_IODebugging()<<METHOD<<"(): Invalid kinematics."<<std::endl;
    return Kin_Args();
  }
  pn=sqrt(pn);
  po=sqrt(po);
  Kin_Args res(vj,xjab,0.0,(mode&4)?1:0,1);
  if (mode&2) {
    if (res.m_mode==0) {
      Vec4D pbh(pb);
      res.m_lam.push_back(Poincare(pa-pj+pb));
      res.m_lam.back().Boost(pbh);
      res.m_lam.push_back(Poincare(pbh,pb));
    }
  }
  res.m_pi=po/pn*(pa-2.0*ma2/(sab-ma2-mb2+pn)*pb)+
    2.0*maj2/(Q2-maj2-mb2+po)*pb;
  res.m_pk=pb;
  if (res.m_pi[0]<0.0 || IsEqual(Q2+po,maj2+mb2)) {
    msg_IODebugging()<<METHOD<<"(): Invalid kinematics."<<std::endl;
    return Kin_Args();
  }
  if (mode&2) {
    if (res.m_mode==0) {
      res.m_lam.push_back(Poincare(res.m_pi+res.m_pk));
      res.m_lam.back().Invert();
    }
    else {
      res.m_lam.push_back(Poincare(pa-pj+pb,res.m_pi+res.m_pk,1));
    }
  }
  if (mode&1) res.m_phi=ComputePhi(pa,pb,pj);
  return res;
}

int PHASIC::ConstructIIDipole
(const double &ma2,const double &mj2,const double &maj2,
 const double &mb2,const Vec4D &paj,const Vec4D &pb,Kin_Args &iip)
{
  double Q2((paj+pb).Abs2()), xjab(iip.m_z), vj(iip.m_y);
  double sab((Q2-mj2)/xjab-(ma2+mb2)*(1-xjab)/xjab);
  double po(sqr(Q2-maj2-mb2)-4.0*maj2*mb2);
  double pn(sqr(sab-ma2-mb2)-4.0*ma2*mb2);
  if (po<0.0 || pn<0.0) {
    msg_IODebugging()<<METHOD<<"(): Invalid kinematics."<<std::endl;
    return -1;
  }
  pn=sqrt(pn);
  po=sqrt(po);
  double ecm(sab-ma2-mb2), gam(0.5*(ecm+pn));
  iip.m_pi=pn/po*(paj-2.0*maj2/(Q2-maj2-mb2+po)*pb)+ma2/gam*pb;
  iip.m_pk=pb;
  double saj(-vj*(sab-ma2-mb2)+ma2+mj2);
  double x(iip.m_x?iip.m_x:xjab+vj);
  double zt(ecm/pn*(x-mb2/gam*(saj+ma2-mj2)/ecm));
  double ktt(vj*(sab-ma2-mb2)*(1.0-zt)-sqr(1.0-zt)*ma2-mj2);
  if (ktt<0.0) {
    msg_IODebugging()<<METHOD<<"(): Invalid kinematics."<<std::endl;
    return -1;
  }
  ktt=sqrt(ktt);
  msg_IODebugging()<<"Set fixed n_perp\n";
  Vec4D n_perp(0.0,1.0,1.0,0.0);
  Poincare zrot(iip.m_pi,Vec4D::ZVEC);
  zrot.RotateBack(n_perp);
  n_perp*=1.0/n_perp.PSpat();
  Vec4D l_perp(LT(iip.m_pi,pb,n_perp));
  l_perp*=1.0/sqrt(dabs(l_perp.Abs2()));
  iip.m_pj=ktt*(sin(iip.m_phi)*l_perp+cos(iip.m_phi)*n_perp);
  iip.m_pj+=(1.0-zt)/pn*(gam*iip.m_pi-ma2*iip.m_pk)+
    (mj2+ktt*ktt)/(1.0-zt)/pn*(iip.m_pk-mb2/gam*iip.m_pi);
  if (iip.m_mode==0) {
    Vec4D pbh(pb);
    iip.m_lam.push_back(iip.m_pi-iip.m_pj+pb);
    iip.m_lam.back().Boost(pbh);
    iip.m_lam.push_back(Poincare(pbh,pb));
    iip.m_lam.push_back(Poincare(paj+pb));
    iip.m_lam.back().Invert();
    iip.m_lam.Invert();
  }
  else {
    iip.m_lam.push_back(Poincare(paj+pb,iip.m_pi-iip.m_pj+pb,1));
  }
  return 1;
}

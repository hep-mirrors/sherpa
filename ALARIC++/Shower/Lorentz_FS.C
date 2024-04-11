#include "ALARIC++/Shower/Lorentz_FS.H"

#include "ALARIC++/Shower/Shower.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace ALARIC;
using namespace PHASIC;
using namespace ATOOLS;

Lorentz_FS_Rad::Lorentz_FS_Rad(const Kernel_Key &k):
  Lorentz(k,0)
{
}

double Lorentz_FS_Rad::Jacobian(const Splitting &s) const
{
  return s.m_J;
}

int Lorentz_FS_Rad::Construct(Splitting &s,const int mode) const
{
  if (mode&1) return Update(s,mode);
  size_t ij, k;
  Vec4D qa(s.p_c->Mom()), Kt(s.m_Kt);
  const Amplitude &a(*s.p_c->Ampl());
  for (size_t i(0);i<a.size();++i) {
    if (a[i]==s.p_c) ij=i;
    if (a[i]==s.p_s) k=i;
  }
  double gam(2.*qa*Kt), t(s.m_t/gam);
  if (p_sk->PS()->EvolScheme(0)==0) {
    s.m_y=t/(1.-s.m_z);
    s.m_J=1.;
  }
  else if (p_sk->PS()->EvolScheme(0)==1) {
    s.m_y=t/(1.-s.m_z+t);
    s.m_J=(1.-s.m_z)/(1.-s.m_z+t);
  }
  else {
    THROW(not_implemented,"Unknown evolution scheme");
  }
  s.m_x=s.m_z;
  Ant_Args ff(s.m_y,s.m_x,s.m_phi);
  ff.m_b=s.m_rcl;
  ff.m_p.reserve(a.size());
  for (size_t i(0);i<a.size();++i) {
    if (a[i]==s.p_c) ff.m_ij=ff.m_p.size();
    else if (a[i]==s.p_s) ff.m_k=ff.m_p.size();
    ff.m_p.push_back(a[i]->Mom());
  }
  if (ConstructAntenna(ff,ij,k,s.m_mi2,s.m_mj2)<0) return -1;
  if (p_ms->Mass2(s.p_c->Flav())) {
    double mu2(-ff.m_pijt.Abs2()/gam), kap(-ff.m_Kt.Abs2()/gam);
    double siij(ff.m_pi*ff.m_pijt/(-mu2*gam));
    s.m_J=1.0-2.0*mu2*(2.0*(kap-mu2)-s.m_z+siij*(1.0+2.0*mu2));
    s.m_J/=pow(sqr(1.0+2.0*mu2*(1.0-siij))-4.0*mu2*(1.0-s.m_z+kap),1.5);
  }
  s.m_pi=ff.m_pi;
  s.m_pj=ff.m_pj;
  s.m_pk=ff.m_pk;
  s.m_kap=Kt.Abs2()/gam;
  s.m_zi=s.m_x;
  if (p_sk->PS()->KernelScheme()&1)
    s.m_zi=(s.m_pi*ff.m_nb)/((s.m_pi+s.m_pj)*ff.m_nb);
  s.m_p=ff.m_p;
  s.m_kt2=s.m_t;
  s.m_K=Vec4D();
  for (size_t i(0);i<a.size();++i)
    if (s.m_rcl[i]&2 && a[i]!=s.p_c) s.m_K+=s.m_p[i];
  if (s.m_iink) s.m_K=-s.m_K-s.m_pi-s.m_pj;
  if (s.m_kt2<s.m_t0) return -1;
  return 1;
}

bool Lorentz_FS_Rad::Cluster
(Splitting &s,const int mode,
 const int i, const int j,const int k,
 const ATOOLS::Cluster_Amplitude *a) const
{
  Ant_Args ff;
  ff.m_b=s.m_rcl;
  s.m_K=Vec4D();
  for (size_t l(0);l<a->Legs().size();++l) {
    if (s.m_rcl[l]&2) {
      if (l!=j) s.m_K+=a->Leg(l)->Mom();
      if (l==i) s.m_iink=1;
    }
    ff.m_p.push_back(a->Leg(l)->Mom());
  }
  if (s.m_iink) s.m_K=-s.m_K-s.m_pj;
  ClusterAntenna(ff,i,j,k,s.m_mij2);
  if (ff.m_stat<0) return false;
  SetParams(s);
  double gam(2.*ff.m_pijt*ff.m_Kt);
  s.m_x=s.m_z=ff.m_z;
  if (p_sk->PS()->EvolScheme(0)==0) {
    Vec4D n(s.m_K+s.m_pj);
    s.m_t=2.*(s.m_pi*s.m_pj)*(s.m_pj*n)/(s.m_pi*n);
    s.m_J=1.;
  }
  else if (p_sk->PS()->EvolScheme(0)==1) {
    s.m_t=2.*(s.m_pi*s.m_pj)*(s.m_pj*s.m_K)/(s.m_pi*s.m_K);
    s.m_J=(1.-s.m_z)/(1.-s.m_z+s.m_t/gam);
  }
  else {
    THROW(not_implemented,"Unknown evolution scheme");
  }
  s.m_y=ff.m_y;
  s.m_phi=ff.m_phi;
  s.m_kap=s.m_K.Abs2()/gam;
  s.m_p=ff.m_p;
  s.m_Kt=ff.m_Kt;
  s.m_kt2=s.m_t;
  s.m_zi=s.m_x;
  if (p_sk->PS()->KernelScheme()&1)
    s.m_zi=(s.m_x*s.m_p[i]*ff.m_nb)/
      ((s.m_x*s.m_p[i]+s.m_p[j])*ff.m_nb);
  s.m_p.erase(s.m_p.begin()+j);
  if (s.m_kt2<s.m_t0) return false;
  return true;
}

Lorentz_FS_Split::Lorentz_FS_Split(const Kernel_Key &k):
  Lorentz(k,0)
{
}

double Lorentz_FS_Split::Jacobian(const Splitting &s) const
{
  return s.m_J;
}

int Lorentz_FS_Split::Construct(Splitting &s,const int mode) const
{
  if (mode&1) return Update(s,mode);
  int nk(0);
  Vec4D qa(s.p_c->Mom()), Kt;
  const Amplitude &a(*s.p_c->Ampl());
  for (size_t i(0);i<a.size();++i)
    if ((s.m_rcl[i]&2) && s.p_c!=a[i]) {
      Kt+=a[i]->Mom();
      ++nk;
    }
  double Q2((qa+Kt).Abs2()), Kt2(Kt.Abs2());
  s.m_y=s.m_t/(Q2-s.m_mi2-s.m_mj2-Kt2);
  s.m_x=s.m_z;
  s.m_J=1.0-s.m_y;
  if (s.m_y<0.0 || s.m_y>1.0) return -1;
  Kin_Args ff(s.m_y,s.m_x,s.m_phi);
  if (ConstructFFDipole
      (s.m_mi2,s.m_mj2,s.m_mij2,Kt2,qa,Kt,ff)<0)
    return -1;
  s.m_pi=ff.m_pi;
  s.m_pj=ff.m_pj;
  s.m_pk=ff.m_pk;
  s.m_q2=(qa+Kt).Abs2();
  s.m_zi=s.m_x;
  if (p_sk->PS()->KernelScheme()&1)
    s.m_zi=(ff.m_pi*ff.m_nb)/((ff.m_pi+ff.m_pj)*ff.m_nb);
  s.m_mk2=p_ms->Mass2(s.p_s->Flav());
  if (nk>1) {
    Poincare oldcms(Kt), newcms(ff.m_pk);
    newcms.Invert();
    s.m_pk=newcms*(oldcms*s.p_s->Mom());
    s.m_p.clear();
    s.m_p.reserve(a.size());
    for (size_t i(0);i<a.size();++i) {
      if ((s.m_rcl[i]&2)==0) s.m_p.push_back(a[i]->Mom());
      else s.m_p.push_back(newcms*(oldcms*a[i]->Mom()));
    }
    s.m_mk2=ff.m_pk.Abs2();
  }
  double q2(s.m_q2-s.m_mi2-s.m_mj2-s.m_mk2);
  s.m_J*=q2/sqrt(Lam(s.m_q2,s.m_mij2,s.m_mk2));
  s.m_J/=(1.0+(s.m_mi2+s.m_mj2-s.m_mij2)/(s.m_y*q2));
  s.m_kt2=s.m_t;
  if (s.m_kt2<s.m_t0) return -1;
  return 1;
}

bool Lorentz_FS_Split::Cluster
(Splitting &s,const int mode,
 const int i, const int j,const int k,
 const ATOOLS::Cluster_Amplitude *a) const
{
  int nk(0);
  Vec4D pij(s.p_c->Mom()+s.p_n->Mom()), K;
  for (size_t l(0);l<a->Legs().size();++l)
    if ((s.m_rcl[l]&2) && l!=i && l!=j) {
      K+=a->Leg(l)->Mom();
      ++nk;
    }
  double Q2((pij+K).Abs2()), K2(K.Abs2());
  Kin_Args ff=ClusterFFDipole
    (s.m_mi2,s.m_mj2,s.m_mij2,K2,
     s.p_c->Mom(),s.p_n->Mom(),K,mode);
  if (ff.m_stat<0) return false;
  SetParams(s);
  s.m_z=s.m_x=ff.m_z;
  s.m_y=ff.m_y;
  s.m_J=1.0;
  s.m_t=s.m_y*(Q2-s.m_mi2-s.m_mj2-K2);
  s.m_phi=ff.m_phi;
  s.m_p=a->Momenta();
  if (nk>1) {
    Poincare oldcms(K), newcms(ff.m_pk);
    newcms.Invert();
    for (size_t l(0);l<s.m_p.size();++l)
      if ((s.m_rcl[l]&2) && l!=i && l!=j)
	s.m_p[l]=newcms*(oldcms*s.m_p[l]);
  }
  s.m_p[i]=ff.m_pi;
  s.m_pk=s.m_p[k];
  s.m_kt2=s.m_t;
  s.m_zi=s.m_x;
  if (p_sk->PS()->KernelScheme()&1)
    s.m_zi=(s.m_pi*ff.m_nb)/((s.m_pi+s.m_pj)*ff.m_nb);
  s.m_p.erase(s.m_p.begin()+j);
  return true;
}

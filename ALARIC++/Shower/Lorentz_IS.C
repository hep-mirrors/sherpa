#include "ALARIC++/Shower/Lorentz_IS.H"

#include "ALARIC++/Shower/Shower.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace ALARIC;
using namespace PHASIC;
using namespace ATOOLS;

Lorentz_IS::Lorentz_IS(const Kernel_Key &k):
  Lorentz(k,1)
{
}

double Lorentz_IS::Jacobian(const Splitting &s) const
{
  if (s.m_clu&1) return 1.0;
  double t(p_sk->PS()->MuF2(s));
  double fo=p_sk->PS()->GetXPDF(s.m_eta,t,m_fl[0],s.p_c->Beam()-1);
  double fn=p_sk->PS()->GetXPDF(s.m_eta/s.m_z,t,m_fl[1],s.p_c->Beam()-1);
  if (dabs(fo)<p_sk->PS()->PDFMin(0)*
      log(1.0-s.m_eta)/log(1.0-p_sk->PS()->PDFMin(1))) return 0.0; 
  return fn/fo*s.m_J;
}

double Lorentz_IS::PDFEstimate(const Splitting &s) const
{
  double fo=p_sk->PS()->GetXPDF
    (s.m_eta,Min(s.m_t1,s.m_Q2),m_fl[0],s.p_c->Beam()-1);
  double fn=p_sk->PS()->GetXPDF
    (s.m_eta,Min(s.m_t1,s.m_Q2),m_fl[1],s.p_c->Beam()-1);
  if (m_fl[1].Mass(true)<1.0 && m_fl[0].Mass(true)>=1.0) {
    double tcut(Max(s.m_t0,sqr(2.0*m_fl[0].Mass(true))));
    double fo0=p_sk->PS()->GetXPDF(s.m_eta,tcut,m_fl[0],s.p_c->Beam()-1);
    double fn0=p_sk->PS()->GetXPDF(0.2,tcut,m_fl[1],s.p_c->Beam()-1);
    if (fo0 && dabs(fo0)<dabs(fo)) fo=fo0;
    if (dabs(fn0)>dabs(fn)) fn=fn0;
  }
  double min=p_sk->PS()->PDFMin(0)*
    log(1.0-s.m_eta)/log(1.0-p_sk->PS()->PDFMin(1));
  if (dabs(fo)<min) return 0.0;
  if (dabs(fn)<min) fn=fo;
  return dabs(fn/fo);
}

int Lorentz_IS::Construct(Splitting &s,const int mode) const
{
  if (mode&1) return Update(s,mode);
  size_t ij, k;
  Vec4D qa(s.p_c->Mom()), Kt;
  const Amplitude &a(*s.p_c->Ampl());
  for (size_t i(0);i<a.size();++i) {
    if (s.m_rcl[i]&2) Kt+=a[i]->Mom();
    if (a[i]==s.p_c) ij=i;
    if (a[i]==s.p_s) k=i;
  }
  double gam(2.*qa*Kt), t(s.m_t/gam);
  s.m_y=t/(1.0-1.0/s.m_z+t);
  s.m_J=(1.-1./s.m_z)/(1.-1./s.m_z+t);
  s.m_x=s.m_z;
  Ant_Args ff(s.m_y,s.m_x,s.m_phi,1);
  ff.m_b=s.m_rcl;
  ff.m_p.reserve(a.size());
  for (size_t i(0);i<a.size();++i) {
    if (a[i]==s.p_c) ff.m_ij=ff.m_p.size();
    else if (a[i]==s.p_s) ff.m_k=ff.m_p.size();
    ff.m_p.push_back(a[i]->Mom());
  }
  if (ConstructAntenna(ff,ij,k,s.m_mi2,s.m_mj2)<0) return -1;
  s.m_pi=ff.m_pi;
  s.m_pj=ff.m_pj;
  s.m_pk=ff.m_pk;
  s.m_kap=Kt.Abs2()/gam;
  s.m_kt2=-gam*s.m_y*((1.0-s.m_x)*(1.0-s.m_y)-s.m_y*s.m_kap);
  if (s.m_kt2<s.m_t0) return -1;
  s.m_p=ff.m_p;
  return 1;
}

bool Lorentz_IS::Cluster
(Splitting &s,const int mode,
 const int i, const int j,const int k,
 const ATOOLS::Cluster_Amplitude *a) const
{
  Ant_Args ff;
  ff.m_b=s.m_rcl;
  Vec4D K;
  for (size_t l(0);l<a->Legs().size();++l) {
    if ((s.m_rcl[l]&2) && l!=j) K+=a->Leg(l)->Mom();
    ff.m_p.push_back(a->Leg(l)->Mom());
  }
  ff.m_mode=1;
  ClusterAntenna(ff,i,j,k,s.m_mij2);
  if (ff.m_stat<0) return false;
  SetParams(s);
  double gam(2.*ff.m_pijt*ff.m_Kt);
  s.m_t=2.*(s.m_pi*s.m_pj)*(s.m_pj*K)/(s.m_pi*K);
  s.m_x=s.m_z=ff.m_z;
  s.m_phi=ff.m_phi;
  s.m_kap=K.Abs2()/gam;
  s.m_kt2=-gam*ff.m_y*((1.0-ff.m_z)*(1.0-ff.m_y)-ff.m_y*s.m_kap);
  s.m_p=ff.m_p;
  s.m_Kt=ff.m_Kt;
  return true;
}

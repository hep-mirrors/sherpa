#include "DIRE/Shower/Lorentz_FF.H"

#include "DIRE/Shower/Shower.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace DIRE;
using namespace PHASIC;
using namespace ATOOLS;

Lorentz_FF::Lorentz_FF(const Kernel_Key &k):
  Lorentz(k,0)
{
}

double Lorentz_FF::Jacobian(const Splitting &s) const
{
  double Q2(s.m_Q2+s.m_mi2+s.m_mj2+s.m_mk2);
  return s.m_Q2/sqrt(Lam(Q2,s.m_mij2,s.m_mk2));
}

int Lorentz_FF::Construct(Splitting &s,const int mode) const
{
  Kin_Args ff(s.m_y,s.m_x,s.m_phi);
  if (ConstructFFDipole
      (s.m_mi2,s.m_mj2,s.m_mij2,
       s.m_mk2,s.p_c->Mom(),s.p_s->Mom(),ff)<0)
    return -1;
  return Update(s,ff,mode);
}

bool Lorentz_FF::Cluster(Splitting &s,const int mode) const
{
  Kin_Args ff=ClusterFFDipole
    (s.m_mi2,s.m_mj2,s.m_mij2,s.m_mk2,
       s.p_c->Mom(),s.p_n->Mom(),s.p_s->Mom(),mode);
  if (ff.m_stat<0) return false;
  SetParams(s,ff);
  s.m_t=s.m_Q2*s.m_y*(1.0-s.m_y)*(1.0-s.m_x);
  s.m_z=1.0-(1.0-s.m_x)*(1.0-s.m_y);
  return true;
}

bool Lorentz_FF::Compute(Splitting &s,const int mode) const
{
  if (mode) {
    s.m_y=s.m_t/s.m_Q2/(1.0-s.m_z);
    s.m_x=(s.m_z-s.m_y)/(1.0-s.m_y);
  }
  if (s.m_mi2==0.0 && s.m_mj2==0.0 && s.m_mk2==0.0)
    return s.m_x>0.0 && s.m_x<1.0
      && s.m_y>0.0 && s.m_y<1.0;
  double nui2(s.m_mi2/s.m_Q2), nuj2(s.m_mj2/s.m_Q2);
  double nuk2(s.m_mk2/s.m_Q2);
  double viji=sqr(s.m_y)-4.0*nui2*nuj2;
  double vijk=sqr(1.0-s.m_y)-4.0*(s.m_y+nui2+nuj2)*nuk2;
  if (viji<0.0 || vijk<0.0) return false;
  viji=sqrt(viji)/(s.m_y+2.0*nui2);
  vijk=sqrt(vijk)/(1.0-s.m_y);
  double frac=(2.0*nui2+s.m_y)/(2.0*(nui2+nuj2+s.m_y));
  double zm=frac*(1.0-viji*vijk), zp=frac*(1.0+viji*vijk);
  return s.m_x>zm && s.m_x<zp
    && s.m_y>0.0 && s.m_y<1.0;
}

Lorentz_FF_123::Lorentz_FF_123(const Kernel_Key &k):
  Lorentz_FF(k)
{
}

double Lorentz_FF_123::Jacobian(const Splitting &s) const
{
  double Q2(s.m_Q2+s.m_mij2+s.m_mk2), J(1.0);
  if (s.m_s) J=sqrt(Lam(s.m_s,s.m_mj2,s.m_ml2))/s.m_s;
  double q2(Q2-s.m_s-s.m_mj2-s.m_mk2);
  return J*(1.0-s.m_y)*sqr(q2)/sqrt(Lam(Q2,s.m_mij2,s.m_mk2))*
    s.m_y/(q2*s.m_y+s.m_s+s.m_mj2-s.m_mij2);
}

int Lorentz_FF_123::Construct(Splitting &s,const int mode) const
{
  if (s.m_s<rpa->gen.SqrtAccu()) s.m_s=0.0;
  Kin_Args ff(s.m_y,s.m_x,s.m_phi);
  if (ConstructFFDipole
      (s.m_s,s.m_mj2,s.m_mij2,
       s.m_mk2,s.p_c->Mom(),s.p_s->Mom(),ff)<0) return -1;
  double y2(2.0*(ff.m_pi*ff.m_pk)/(s.m_s-s.m_mi2-s.m_ml2));
  Kin_Args ff2(s.m_s?1.0/(1.0+y2):0.0,s.m_z2,s.m_phi2);
  if (ConstructFFDipole
      (s.m_mi2,s.m_ml2,s.m_s,
       s.m_mk2,ff.m_pi,ff.m_pk,ff2)<0) return -1;
  ff.m_pi=ff2.m_pi;
  return Update(s,ff,mode,ff2.m_pj);
}

bool Lorentz_FF_123::Cluster(Splitting &s,const int mode) const
{
  return false;
}

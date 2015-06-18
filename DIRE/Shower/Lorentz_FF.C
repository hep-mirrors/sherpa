#include "DIRE/Shower/Lorentz_FF.H"

#include "DIRE/Shower/Shower.H"
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

bool Lorentz_FF::Compute(Splitting &s) const
{
  s.m_y=s.m_t/s.m_Q2/(1.0-s.m_z);
  s.m_x=(s.m_z-s.m_y)/(1.0-s.m_y);
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

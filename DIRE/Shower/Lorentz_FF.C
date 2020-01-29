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
  double Q2(s.m_q2-s.m_mi2-s.m_mj2-s.m_mk2);
  double J(Q2/sqrt(Lam(s.m_q2,s.m_mij2,s.m_mk2)));
  return J/(1.0+(s.m_mi2+s.m_mj2-s.m_mij2)/(s.m_y*Q2));
}

int Lorentz_FF::Construct(Splitting &s,const int mode) const
{
  if (mode&1) return Update(s,mode);
  if (s.m_q2<sqr(sqrt(s.m_mi2)+sqrt(s.m_mj2)+sqrt(s.m_mk2))) return -1;
  s.m_y=s.m_t/(s.m_q2-s.m_mi2-s.m_mj2-s.m_mk2)/(1.0-s.m_z);
  s.m_x=(s.m_z-s.m_y)/(1.0-s.m_y);
  Kin_Args ff(s.m_y,s.m_x,s.m_phi);
  if (ConstructFFDipole
      (s.m_mi2,s.m_mj2,s.m_mij2,
       s.m_mk2,s.p_c->Mom(),s.p_s->Mom(),ff)<0)
    return -1;
  s.m_pi=ff.m_pi;
  s.m_pj=ff.m_pj;
  s.m_pk=ff.m_pk;
  return 1;
}

bool Lorentz_FF::Cluster(Splitting &s,const int mode) const
{
  Kin_Args ff=ClusterFFDipole
    (s.m_mi2,s.m_mj2,s.m_mij2,s.m_mk2,
       s.p_c->Mom(),s.p_n->Mom(),s.p_s->Mom(),mode);
  if (ff.m_stat<0) return false;
  SetParams(s,ff);
  double Q2(s.m_q2-s.m_mi2-s.m_mj2-s.m_mk2);
  s.m_t=Q2*s.m_y*(1.0-s.m_y)*(1.0-s.m_x);
  s.m_z=1.0-(1.0-s.m_x)*(1.0-s.m_y);
  return true;
}

Lorentz_FF_123::Lorentz_FF_123(const Kernel_Key &k):
  Lorentz_FF(k)
{
}

void Lorentz_FF_123::SetMS(ATOOLS::Mass_Selector *const ms)
{
  p_ms=ms;
  if (p_ms->Mass(m_fl[0]) || p_ms->Mass(m_fl[1])) p_sk->SetOn(0);
  else p_sk->SetOn(1);
}

bool Lorentz_FF_123::Allowed(const Splitting &s) const
{
  if (p_ms->Mass(s.p_s->Flav())) return false;
  return Lorentz::Allowed(s);
}

double Lorentz_FF_123::Jacobian(const Splitting &s) const
{
  double q2(s.m_q2-s.m_mij2-s.m_mk2);
  double J1(q2/sqrt(Lam(s.m_q2,s.m_mij2,s.m_mk2)));
  double saik(s.m_z/s.m_z2*q2+s.m_s+s.m_mk2);
  double J2((saik-s.m_s-s.m_mk2)/sqrt(Lam(saik,s.m_s,s.m_mk2)));
  return J1*J2/(1.0 + (s.m_s+s.m_mj2-s.m_mij2) / (s.m_t*s.m_z2/s.m_z));
}

int Lorentz_FF_123::Construct(Splitting &s,const int mode) const
{
  //msg_Out()<<METHOD<<" for mode = "<<mode<<".\n";
  if ((mode&1) && !(s.m_mode&1)) return Update(s,mode);
  if (s.m_s<rpa->gen.SqrtAccu()) s.m_s=0.0;
  if ((mode&1) && (s.m_mode&1)) s.m_s=0.0;
  double Q2(s.m_q2-s.m_s-s.m_mj2-s.m_mk2);
  s.m_y=s.m_t*s.m_z2/s.m_z/Q2;
  s.m_x=s.m_z/s.m_z2/(1.0-s.m_y)*(s.m_q2-s.m_mij2-s.m_mk2)/Q2;
  Kin_Args ff(s.m_y,s.m_x,s.m_phi);
  // msg_Out()<<"    *** before first FFDipole: y = "<<s.m_y<<", x = "<<s.m_x<<", "
  // 	   <<"phi = "<<s.m_phi<<" for t = "<<s.m_t<<".\n";
  if (ConstructFFDipole
      (s.m_s,s.m_mj2,s.m_mij2,
       s.m_mk2,s.p_c->Mom(),s.p_s->Mom(),ff)<0) return -1;
  // msg_Out()<<"  --> pai = "<<ff.m_pi<<"\n"
  // 	   <<"      pk  = "<<ff.m_pk<<"\n";
  double y2(2.0*(ff.m_pi*ff.m_pk)/(s.m_s-s.m_mi2-s.m_ml2));
  Kin_Args ff2(s.m_s?1.0/(1.0+y2):0.0,s.m_z2,s.m_phi2);
  // msg_Out()<<"    *** before second FFDipole for "
  // 	   <<"sai = "<<s.m_s<<" y2 = "<<y2<<", z2 = "<<s.m_z2<<", phi2 = "<<s.m_phi2<<".\n";
  if (ConstructFFDipole
      (s.m_mi2,s.m_ml2,s.m_s,
       s.m_mk2,ff.m_pi,ff.m_pk,ff2)<0) return -1;
  s.m_pi=ff2.m_pi;
  s.m_pl=ff2.m_pj;
  s.m_pj=ff.m_pj;
  s.m_pk=ff.m_pk;
  // msg_Out()<<"    *** before final update:\n"
  // 	   <<"  pa = "<<s.m_pi<<"\n"
  // 	   <<"  pi = "<<s.m_pl<<"\n"
  // 	   <<"  pj = "<<s.m_pj<<"\n"
  // 	   <<"  pk = "<<s.m_pk<<"\n";
  return (mode&1)?Update(s,mode):1;
}

bool Lorentz_FF_123::Cluster(Splitting &s,const int mode) const
{
  return false;
}

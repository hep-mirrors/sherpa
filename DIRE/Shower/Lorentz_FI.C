#include "DIRE/Shower/Lorentz_FI.H"

#include "DIRE/Shower/Shower.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace DIRE;
using namespace PHASIC;
using namespace ATOOLS;

Lorentz_FI::Lorentz_FI(const Kernel_Key &k):
  Lorentz(k,2)
{
}

double Lorentz_FI::Jacobian(const Splitting &s) const
{
  if (s.m_clu&1) return 1.0;
  double eta(s.p_s->GetXB());
  double y(s.m_y*(1.0+(s.m_mij2-s.m_mi2-s.m_mj2)/s.m_Q2));
  double fo=p_sk->PS()->GetXPDF(eta,s.m_t,s.p_s->Flav(),s.p_s->Beam()-1);
  double fn=p_sk->PS()->GetXPDF(eta/y,s.m_t,s.p_s->Flav(),s.p_s->Beam()-1);
  if (dabs(fo)<p_sk->PS()->PDFMin(0)*
      log(1.0-eta)/log(1.0-p_sk->PS()->PDFMin(1))) return 0.0; 
  return (1.0-s.m_y)/(1.0-y)*fn/fo;
}

double Lorentz_FI::PDFEstimate(const Splitting &s) const
{
  return 1.0;
}

int Lorentz_FI::Construct(Splitting &s,const int mode) const
{
  if (mode&1) return Update(s,mode);
  s.m_y=1.0/(1.0+s.m_t/s.m_Q2/(1.0-s.m_z));
  s.m_x=s.m_z;
  Kin_Args ff(1.0-s.m_y,s.m_x,s.m_phi,1|8);
  if (ConstructFIDipole
      (s.m_mi2,s.m_mj2,s.m_mij2,
       s.m_mk2,s.p_c->Mom(),-s.p_s->Mom(),ff)<0)
    return -1;
  s.m_pi=ff.m_pi;
  s.m_pj=ff.m_pj;
  s.m_pk=-ff.m_pk;
  return 1;
}

bool Lorentz_FI::Cluster(Splitting &s,const int mode) const
{
  Kin_Args ff=ClusterFIDipole
    (s.m_mi2,s.m_mj2,s.m_mij2,s.m_mk2,
     s.p_c->Mom(),s.p_n->Mom(),-s.p_s->Mom(),mode|8);
  if (ff.m_stat<0) return false;
  ff.m_y=1.0-ff.m_y;
  SetParams(s,ff);
  s.m_t=s.m_Q2*(1.0-s.m_y)/s.m_y*(1.0-s.m_x);
  s.m_z=s.m_x;
  return true;
}

Lorentz_FI_123::Lorentz_FI_123(const Kernel_Key &k):
  Lorentz_FI(k)
{
}

void Lorentz_FI_123::SetMS(ATOOLS::Mass_Selector *const ms)
{
  p_ms=ms;
  if (p_ms->Mass(m_fl[0]) || p_ms->Mass(m_fl[1])) p_sk->SetOn(0);
  else p_sk->SetOn(1);
}

bool Lorentz_FI_123::Allowed(const Splitting &s) const
{
  if (p_ms->Mass(s.p_s->Flav())) return false;
  return Lorentz::Allowed(s);
}

double Lorentz_FI_123::Jacobian(const Splitting &s) const
{
  double eta(s.p_s->GetXB());
  double y(s.m_y*(1.0-(s.m_mij2-s.m_s-s.m_mj2)/
		  (s.m_q2-s.m_s-s.m_mj2-s.m_mk2)));
  double fo=p_sk->PS()->GetXPDF(eta,s.m_t,s.p_s->Flav(),s.p_s->Beam()-1);
  double fn=p_sk->PS()->GetXPDF(eta/y,s.m_t,s.p_s->Flav(),s.p_s->Beam()-1);
  if (dabs(fo)<p_sk->PS()->PDFMin(0)*
      log(1.0-eta)/log(1.0-p_sk->PS()->PDFMin(1))) return 0.0; 
  double saij(s.m_t*s.m_z2/s.m_z+s.m_s+s.m_mj2);
  double x((s.m_q2-s.m_mij2-s.m_mk2)/(s.m_q2-saij-s.m_mk2)), rho(x);
  double J1(rho/x*(saij+s.m_mk2-s.m_q2)/sqrt(Lam(saij,s.m_mk2,s.m_q2)));
  double sbai(s.m_z/s.m_z2*(s.m_q2-saij-s.m_mk2)+s.m_s+s.m_mk2);
  double J2((s.m_s+s.m_mk2-sbai)/sqrt(Lam(s.m_s,s.m_mk2,sbai)));
  return J1*J2*fn/fo/(1.0+(s.m_s+s.m_mj2-s.m_mij2)/(s.m_t*s.m_z2/s.m_z));
}

int Lorentz_FI_123::Construct(Splitting &s,const int mode) const
{
  if ((mode&1) && !(s.m_mode&1)) return Update(s,mode);
  if (s.m_s<rpa->gen.SqrtAccu()) s.m_s=0.0;
  if ((mode&1) && (s.m_mode&1)) s.m_s=0.0;
  s.m_x=s.m_z/s.m_z2;
  s.m_y=1.0/(1.0-s.m_t/s.m_x/(s.m_q2-s.m_s-s.m_mj2-s.m_mk2));
  double y(s.m_y*(1.0-(s.m_mij2-s.m_s-s.m_mj2)/
		  (s.m_q2-s.m_s-s.m_mj2-s.m_mk2)));
  if (y<s.p_s->GetXB()) return -1;
  Kin_Args ff(1.0-s.m_y,s.m_x,s.m_phi,1|8);
  if (ConstructFIDipole
      (s.m_s,s.m_mj2,s.m_mij2,
       s.m_mk2,s.p_c->Mom(),-s.p_s->Mom(),ff)<0) return -1;
  double y2(2.0*(ff.m_pi*ff.m_pk)/(s.m_s-s.m_mi2-s.m_ml2));
  Kin_Args ff2(s.m_s?1.0/(1.0-y2):0.0,s.m_z2,s.m_phi2);
  if (ConstructFFDipole
      (s.m_mi2,s.m_ml2,s.m_s,
       s.m_mk2,ff.m_pi,-ff.m_pk,ff2)<0) return -1;
  s.m_pk=-ff.m_pk;
  s.m_pi=ff2.m_pi;
  s.m_pl=ff2.m_pj;
  s.m_pj=ff.m_pj;
  return (mode&1)?Update(s,mode):1;
}

bool Lorentz_FI_123::Cluster(Splitting &s,const int mode) const
{
  return false;
}

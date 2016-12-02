#include "DIRE/Shower/Lorentz_II.H"

#include "DIRE/Shower/Shower.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace DIRE;
using namespace PHASIC;
using namespace ATOOLS;

Lorentz_II::Lorentz_II(const Kernel_Key &k):
  Lorentz(k,3)
{
}

double Lorentz_II::Jacobian(const Splitting &s) const
{
  if (s.m_clu&1) return 1.0;
  double fo=p_sk->PS()->GetXPDF(s.m_eta,s.m_t,m_fl[0],s.p_c->Beam()-1);
  double fn=p_sk->PS()->GetXPDF(s.m_eta/s.m_x,s.m_t,m_fl[1],s.p_c->Beam()-1);
  if (dabs(fo)<p_sk->PS()->PDFMin()) return 0.0; 
  return fn/fo;
}

double Lorentz_II::PDFEstimate(const Splitting &s) const
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
  if (fo<p_sk->PS()->PDFMin()) return 0.0;
  return dabs(fn/fo);
}

int Lorentz_II::Construct(Splitting &s,const int mode) const
{
  Kin_Args ff(s.m_y,s.m_x,s.m_phi,s.m_kin);
  if (ConstructIIDipole
      (s.m_mi2,s.m_mj2,s.m_mij2,
       s.m_mk2,-s.p_c->Mom(),-s.p_s->Mom(),ff)<0)
    return -1;
  ff.m_pi=-ff.m_pi;
  ff.m_pk=-ff.m_pk;
  return Update(s,ff,mode);
}

bool Lorentz_II::Cluster(Splitting &s,const int mode) const
{
  Kin_Args ff=ClusterIIDipole
    (s.m_mi2,s.m_mj2,s.m_mij2,s.m_mk2,
     -s.p_c->Mom(),s.p_n->Mom(),-s.p_s->Mom(),mode);
  if (ff.m_stat<0) return false;
  SetParams(s,ff);
  s.m_t=s.m_Q2*s.m_y*(1.0-s.m_x-s.m_y);
  s.m_z=s.m_x+s.m_y;
  return true;
}

bool Lorentz_II::Compute(Splitting &s,const int mode) const
{
  if (mode) {
    s.m_y=s.m_t/s.m_Q2/(1.0-s.m_z);
    s.m_x=s.m_z-s.m_t/s.m_Q2/(1.0-s.m_z);
  }
  return s.m_x>s.m_eta && s.m_y>=0.0 && s.m_x+s.m_y<1.0;
}

Lorentz_II_123::Lorentz_II_123(const Kernel_Key &k):
  Lorentz_II(k)
{
}

double Lorentz_II_123::Jacobian(const Splitting &s) const
{
  Splitting c(s);
  c.m_x=s.m_z;
  return Lorentz_II::Jacobian(c);
}

int Lorentz_II_123::Construct(Splitting &s,const int mode) const
{
  if (s.m_s<rpa->gen.SqrtAccu()) s.m_s=0.0;
  double Q2((s.p_c->Mom()+s.p_s->Mom()).Abs2()-s.m_mij2-s.m_mk2);
  double q2(Q2+s.m_mij2-s.m_mj2+s.m_s), xa(s.m_z2), za(s.m_z);
  Kin_Args ff(s.m_t/Q2*za/xa,q2/(xa/za*Q2),s.m_phi,s.m_kin);
  if (ConstructIIDipole
      (-s.m_s,s.m_mj2,s.m_mij2,s.m_mk2,
       -s.p_c->Mom(),-s.p_s->Mom(),ff)<0)
    return -1;
  double v((s.m_s+s.m_mi2+s.m_ml2)/(Q2/za));
  Kin_Args ff2(v,xa-v,s.m_phi2,s.m_kin);
  if (ConstructIIDipole
      (s.m_mi2,s.m_ml2,-s.m_s,s.m_mk2,
       ff.m_pi,ff.m_pk,ff2)<0) return -1;
  ff.m_pj=ff2.m_lam*ff.m_pj;
  ff.m_lam.insert(ff.m_lam.end(),ff2.m_lam.begin(),ff2.m_lam.end());
  ff.m_pk=-ff2.m_pk;
  ff.m_pi=-ff2.m_pi;
  return Update(s,ff,mode,ff2.m_pj);
}

bool Lorentz_II_123::Cluster(Splitting &s,const int mode) const
{
  return false;
}

#include "DIRE/Shower/Lorentz_II.H"

#include "DIRE/Shower/Shower.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
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
  if (s.m_t<sqr(2.0*m_fl[0].Mass(true)) && fn>fo) fn=fo;
  return fn/fo;
}

double Lorentz_II::PDFEstimate(const Splitting &s) const
{
  double fo=p_sk->PS()->GetXPDF
    (s.m_eta,Min(s.m_t1,s.m_Q2),m_fl[0],s.p_c->Beam()-1);
  double fn=p_sk->PS()->GetXPDF
    (s.m_eta,Min(s.m_t1,s.m_Q2),m_fl[1],s.p_c->Beam()-1);
  if (m_fl[1]==Flavour(kf_u).Bar()||m_fl[1]==Flavour(kf_d).Bar()) {
    double fo0=p_sk->PS()->GetXPDF(s.m_eta,s.m_t0,m_fl[0],s.p_c->Beam()-1);
    double fn0=p_sk->PS()->GetXPDF(0.2,s.m_t0,m_fl[1],s.p_c->Beam()-1);
    if (fo0 && dabs(fo0)<dabs(fo)) fo=fo0;
    if (dabs(fn0)>dabs(fn)) fn=fn0;
  }
  if (fo==0.0) return 0.0;
  return dabs(fn/fo)*pow(1.0e6,s.m_eta*s.m_t0/Min(s.m_t1,s.m_Q2));
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

bool Lorentz_II::Compute(Splitting &s) const
{
  s.m_y=s.m_t/s.m_Q2/(1.0-s.m_z);
  s.m_x=s.m_z-s.m_t/s.m_Q2/(1.0-s.m_z);
  return s.m_x>s.m_eta && s.m_y>0.0 && s.m_x+s.m_y<1.0;
}

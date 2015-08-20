#include "DIRE/Shower/Lorentz_FI.H"

#include "DIRE/Shower/Shower.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
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
  if (dabs(fo)<p_sk->PS()->PDFMin()) return 0.0; 
  return y*fn/fo;
}

double Lorentz_FI::PDFEstimate(const Splitting &s) const
{
  return pow(1.0e6,s.p_s->GetXB()*s.m_t0/Min(s.m_t1,s.m_Q2));
}

int Lorentz_FI::Construct(Splitting &s,const int mode) const
{
  Kin_Args ff(1.0-s.m_y,s.m_x,s.m_phi,1|8);
  if (ConstructFIDipole
      (s.m_mi2,s.m_mj2,s.m_mij2,
       s.m_mk2,s.p_c->Mom(),-s.p_s->Mom(),ff)<0)
    return -1;
  ff.m_pk=-ff.m_pk;
  return Update(s,ff,mode);
}

bool Lorentz_FI::Cluster(Splitting &s,const int mode) const
{
  Kin_Args ff=ClusterFIDipole
    (s.m_mi2,s.m_mj2,s.m_mij2,s.m_mk2,
     s.p_c->Mom(),s.p_n->Mom(),-s.p_s->Mom(),mode|8);
  if (ff.m_stat<0) return false;
  ff.m_y=1.0-ff.m_y;
  SetParams(s,ff);
  s.m_t=s.m_Q2*(1.0-s.m_y)*(1.0-s.m_x);
  s.m_z=s.m_x;
  return true;
}

bool Lorentz_FI::Compute(Splitting &s) const
{
  s.m_y=1.0-s.m_t/s.m_Q2/(1.0-s.m_z);
  s.m_x=s.m_z;
  if (s.m_mi2==0.0 && s.m_mj2==0.0)
    return s.m_y>s.p_s->GetXB();
  double nui2(s.m_mi2/s.m_Q2*s.m_y), nuj2(s.m_mj2/s.m_Q2*s.m_y);
  double viji=sqr(1.0-s.m_y)-4.0*nui2*nuj2;
  if (viji<0.0 || s.m_y>1.0) return false;
  viji=sqrt(viji)/(1.0-s.m_y+2.0*nui2);
  double frac=(1.0-s.m_y+2.0*nui2)/(2.0*(1.0-s.m_y+nui2+nuj2));
  double zm=frac*(1.0-viji), zp=frac*(1.0+viji);
  double y=s.m_y*(1.0+(s.m_mij2-s.m_mi2-s.m_mj2)/s.m_Q2);
  return s.m_x>zm && s.m_x<zp && y>s.p_s->GetXB();
}

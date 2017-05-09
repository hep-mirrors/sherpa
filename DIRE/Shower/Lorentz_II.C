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
  if (dabs(fo)<p_sk->PS()->PDFMin()) return 0.0;
  if (dabs(fn)<p_sk->PS()->PDFMin()) fn=fo;
  return dabs(fn/fo);
}

int Lorentz_II::Construct(Splitting &s,const int mode) const
{
  if (mode&1) return Update(s,mode);
  s.m_y=s.m_t/s.m_Q2/(1.0-s.m_z);
  s.m_x=s.m_z-s.m_t/s.m_Q2/(1.0-s.m_z);
  Kin_Args ff(s.m_y,s.m_x,s.m_phi,s.m_kin);
  if (ConstructIIDipole
      (s.m_mi2,s.m_mj2,s.m_mij2,
       s.m_mk2,-s.p_c->Mom(),-s.p_s->Mom(),ff)<0)
    return -1;
  s.m_pi=-ff.m_pi;
  s.m_pj=ff.m_pj;
  s.m_pk=-ff.m_pk;
  s.m_lam=ff.m_lam;
  return 1;
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

Lorentz_II_123::Lorentz_II_123(const Kernel_Key &k):
  Lorentz_II(k)
{
}

void Lorentz_II_123::SetMS(ATOOLS::Mass_Selector *const ms)
{
  p_ms=ms;
  if (p_ms->Mass(m_fl[0]) || p_ms->Mass(m_fl[1])) p_sk->SetOn(0);
  else p_sk->SetOn(1);
}

bool Lorentz_II_123::Allowed(const Splitting &s) const
{
  if (p_ms->Mass(s.p_s->Flav())) return false;
  return Lorentz::Allowed(s);
}

double Lorentz_II_123::Jacobian(const Splitting &s) const
{
  double fo=p_sk->PS()->GetXPDF(s.m_eta,s.m_t,m_fl[0],s.p_c->Beam()-1);
  double fn=p_sk->PS()->GetXPDF(s.m_eta/s.m_x,s.m_t,m_fl[1],s.p_c->Beam()-1);
  if (dabs(fo)<p_sk->PS()->PDFMin()) return 0.0;
  double sab(s.m_q2/s.m_z+s.m_mi2+s.m_mk2);
  double J1((sab-s.m_mi2-s.m_mk2)/sqrt(Lam(sab,s.m_mi2,s.m_mk2)));
  double sjq(s.m_q2*s.m_z2/s.m_z-s.m_s+s.m_mk2);
  double J2((sjq+s.m_s-s.m_mk2)/sqrt(Lam(sjq,-s.m_s,s.m_mk2)));
  return J1*J2*fn/fo/(1.0+(-s.m_s+s.m_mj2-s.m_mij2)/(-s.m_t/s.m_z2));
}

int Lorentz_II_123::Construct(Splitting &s,const int mode) const
{
  if ((mode&1) && !(s.m_mode&1)) return Update(s,mode);
  if (s.m_s<rpa->gen.SqrtAccu()) s.m_s=0.0;
  if ((mode&1) && (s.m_mode&1)) s.m_s=0.0;
  s.m_x=s.m_z*(s.m_q2-s.m_mi2-s.m_ml2-s.m_mk2)/s.m_q2;
  s.m_y=s.m_z*(s.m_s+s.m_mi2+s.m_ml2)/s.m_q2;
  double sai(-s.m_s-s.m_mi2-s.m_ml2), xa(s.m_z2), za(s.m_z);
  double x(za*(s.m_q2-s.m_mi2-s.m_ml2-s.m_mk2)/s.m_q2);
  Kin_Args ff(-za*sai/s.m_q2,x,s.m_phi,s.m_kin);
  ff.m_x=xa;
  if (ConstructIIDipole
      (s.m_mi2,s.m_ml2,s.m_mij2,s.m_mk2,
       -s.p_c->Mom(),-s.p_s->Mom(),ff)<0) return -1;
  Vec4D pjq(ff.m_pi-ff.m_pj+ff.m_pk);
  double sjq(s.m_q2*xa/za-s.m_s+s.m_mk2);
  double y2((s.m_q2*xa/za-2.0*s.m_s)/(sjq-s.m_mj2-s.m_q2));
  double z2(s.m_t/xa/(s.m_q2*xa/za-2.0*s.m_s));
  if (y2<0.0 || z2<0.0) return -1;
  Kin_Args ff2(1.0/(1.0+y2),z2,s.m_phi2);
  if (ConstructFFDipole
      (s.m_mj2,s.m_q2,sjq,-s.m_s,
       pjq,ff.m_pi-ff.m_pj,ff2)<0) return -1;
  if (mode<0) return 1;
  Vec4D pj(ff2.m_pi), q(ff2.m_pj), pai(ff2.m_pk);
  s.m_lam.clear();
  if (s.m_kin==0) {
    Vec4D pbh(ff.m_pk);
    s.m_lam.push_back(ff.m_pi-ff.m_pj-ff2.m_pi+ff.m_pk);
    s.m_lam.back().Boost(pbh);
    s.m_lam.push_back(Poincare(pbh,ff.m_pk));
    s.m_lam.push_back(Poincare(-s.p_c->Mom()+ff.m_pk));
    s.m_lam.back().Invert();
    s.m_lam.Invert();
  }
  else {
    s.m_lam.push_back
      (Poincare(-s.p_c->Mom()-s.p_s->Mom(),
		ff.m_pi-ff.m_pj-ff2.m_pi+ff.m_pk,1));
  }
  s.m_pi=-ff.m_pi;
  s.m_pk=-ff.m_pk;
  s.m_pl=ff.m_pj;
  s.m_pj=ff2.m_pi;
  return (mode&1)?Update(s,mode):1;
}

bool Lorentz_II_123::Cluster(Splitting &s,const int mode) const
{
  return false;
}

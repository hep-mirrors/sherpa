#include "DIRE/Shower/Lorentz_IF.H"

#include "DIRE/Shower/Shower.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace DIRE;
using namespace PHASIC;
using namespace ATOOLS;

Lorentz_IF::Lorentz_IF(const Kernel_Key &k):
  Lorentz(k,1)
{
}

double Lorentz_IF::Jacobian(const Splitting &s) const
{
  if (s.m_clu&1) return 1.0;
  double fo=p_sk->PS()->GetXPDF(s.m_eta,s.m_t,m_fl[0],s.p_c->Beam()-1);
  double fn=p_sk->PS()->GetXPDF(s.m_eta/s.m_x,s.m_t,m_fl[1],s.p_c->Beam()-1);
  if (dabs(fo)<p_sk->PS()->PDFMin(0)*
      log(1.0-s.m_eta)/log(1.0-p_sk->PS()->PDFMin(1))) return 0.0; 
  return fn/fo;
}

double Lorentz_IF::PDFEstimate(const Splitting &s) const
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

int Lorentz_IF::Construct(Splitting &s,const int mode) const
{
  if (mode&1) return Update(s,mode);
  s.m_y=s.m_t/s.m_Q2/(1.0-s.m_z);
  s.m_x=s.m_z;
  Parton *b(NULL);
  if (s.m_kin==0)
    for (size_t i(0);i<s.p_c->Ampl()->size();++i)
      if ((*s.p_c->Ampl())[i]->Beam()==3-s.p_c->Beam()) {
	b=(*s.p_c->Ampl())[i];
	break;
      }
  Kin_Args ff(s.m_y,s.m_x,s.m_phi,s.m_kin);
  if (ConstructIFDipole
      (s.m_mi2,s.m_mj2,s.m_mij2,s.m_mk2,b?p_ms->Mass2(b->Flav()):0.0,
       -s.p_c->Mom(),s.p_s->Mom(),b?-b->Mom():Vec4D(),ff)<0)
    return -1;
  s.m_pi=-ff.m_pi;
  s.m_pj=ff.m_pj;
  s.m_pk=ff.m_pk;
  if (b && p_sk->PS()->RemnantTest(s.p_c,s.m_pi)<0) return -1; 
  s.m_lam=ff.m_lam;
  return 1;
}

bool Lorentz_IF::Cluster(Splitting &s,const int mode) const
{
  Kin_Args ff=ClusterIFDipole
    (s.m_mi2,s.m_mj2,s.m_mij2,s.m_mk2,0.0,
     -s.p_c->Mom(),s.p_n->Mom(),s.p_s->Mom(),Vec4D(),
     mode|(s.m_kin?4:0));
  if (ff.m_stat<0) return false;
  SetParams(s,ff);
  s.m_t=s.m_Q2*s.m_y*(1.0-s.m_x);
  s.m_z=s.m_x;
  return true;
}

Lorentz_IF_123::Lorentz_IF_123(const Kernel_Key &k):
  Lorentz_IF(k)
{
}

void Lorentz_IF_123::SetMS(ATOOLS::Mass_Selector *const ms)
{
  p_ms=ms;
  if (p_ms->Mass(m_fl[0]) || p_ms->Mass(m_fl[1])) p_sk->SetOn(0);
  else p_sk->SetOn(1);
}

bool Lorentz_IF_123::Allowed(const Splitting &s) const
{
  if (p_ms->Mass(s.p_s->Flav())) return false;
  return Lorentz::Allowed(s);
}

double Lorentz_IF_123::Jacobian(const Splitting &s) const
{
  double fo=p_sk->PS()->GetXPDF(s.m_eta,s.m_t,m_fl[0],s.p_c->Beam()-1);
  double fn=p_sk->PS()->GetXPDF(s.m_eta/s.m_z,s.m_t,m_fl[1],s.p_c->Beam()-1);
  if (dabs(fo)<p_sk->PS()->PDFMin(0)*
      log(1.0-s.m_eta)/log(1.0-p_sk->PS()->PDFMin(1))) return 0.0; 
  double sijk(s.m_q2*(1.0-1.0/s.m_z)-s.m_mi2), rho(s.m_z);
  double J1(rho/s.m_z*(sijk+s.m_mi2-s.m_q2)/sqrt(Lam(sijk,s.m_mi2,s.m_q2)));
  double sjk(-s.m_q2*(s.m_z2/s.m_z-1.0)+s.m_t/s.m_z2+s.m_s);
  double J2(-s.m_q2*s.m_z2/s.m_z/sqrt(Lam(sjk,-s.m_s,s.m_q2)));
  return J1*J2*fn/fo/(1.0+(-s.m_s+s.m_mj2-s.m_mij2)/(-s.m_t/s.m_z2));
}

int Lorentz_IF_123::Construct(Splitting &s,const int mode) const
{
  if ((mode&1) && !(s.m_mode&1)) return Update(s,mode);
  Parton *b(NULL);
  if (s.m_kin==0)
    for (size_t i(0);i<s.p_c->Ampl()->size();++i)
      if ((*s.p_c->Ampl())[i]->Beam()==3-s.p_c->Beam()) {
	b=(*s.p_c->Ampl())[i];
	break;
      }
  if (s.m_s<rpa->gen.SqrtAccu()) s.m_s=0.0;
  if ((mode&1) && (s.m_mode&1)) s.m_s=0.0;
  s.m_y=(-s.m_s-s.m_mi2-s.m_ml2)*s.m_z/s.m_q2;
  s.m_x=s.m_y+s.m_z2-s.m_t*s.m_z/(s.m_q2*s.m_z2);
  double sai(s.m_s+s.m_mi2+s.m_ml2), xa(s.m_z2), za(s.m_z);
  double Q2(-(s.p_c->Mom()+s.p_s->Mom()).Abs2());
  Kin_Args ff(sai*za/Q2,xa+sai*za/Q2+s.m_t*za/(Q2*xa),s.m_phi,s.m_kin);
  if (ff.m_z>1.0) return -1;
  ff.m_mk2=Q2*(xa/za-1.0)+s.m_s+s.m_t/xa;
  if (ff.m_mk2<sqr(sqrt(s.m_mj2)+sqrt(s.m_mk2))) return -1;
  if (ConstructIFDipole
      (s.m_mi2,s.m_ml2,s.m_mij2,s.m_mk2,b?p_ms->Mass2(b->Flav()):0.0,
       -s.p_c->Mom(),s.p_s->Mom(),b?-b->Mom():Vec4D(),ff)<0)
    return -1;
  double y2((s.m_t/xa+Q2*xa/za)/(ff.m_mk2-s.m_mj2-s.m_mk2));
  Kin_Args ff2(1.0/(1.0-y2),1.0/(1.0+Q2/s.m_t*xa*xa/za),s.m_phi2);
  if (ConstructFFDipole
      (s.m_mj2,s.m_mk2,ff.m_mk2,
       -s.m_s,ff.m_pk,-ff.m_pi+ff.m_pj,ff2)<0) return -1;
  s.m_pk=ff2.m_pj;
  s.m_pi=-ff.m_pi;
  s.m_pj=ff2.m_pi;
  s.m_pl=ff.m_pj;
  s.m_lam=ff.m_lam;
  return (mode&1)?Update(s,mode):1;
}

bool Lorentz_IF_123::Cluster(Splitting &s,const int mode) const
{
  return false;
}

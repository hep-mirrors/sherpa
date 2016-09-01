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
  if (dabs(fo)<p_sk->PS()->PDFMin()) return 0.0; 
  return fn/fo;
}

double Lorentz_IF::PDFEstimate(const Splitting &s) const
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

int Lorentz_IF::Construct(Splitting &s,const int mode) const
{
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
  ff.m_pi=-ff.m_pi;
  if (b && p_sk->PS()->RemnantTest(s.p_c,ff.m_pi)<0) return -1; 
  return Update(s,ff,mode);
}

bool Lorentz_IF::Cluster(Splitting &s,const int mode) const
{
  Kin_Args ff=ClusterIFDipole
    (s.m_mi2,s.m_mj2,s.m_mij2,s.m_mk2,0.0,
     -s.p_c->Mom(),s.p_n->Mom(),s.p_s->Mom(),Vec4D(),mode|4);
  if (ff.m_stat<0) return false;
  SetParams(s,ff);
  s.m_t=s.m_Q2*s.m_y*(1.0-s.m_x);
  s.m_z=s.m_x;
  return true;
}

bool Lorentz_IF::Compute(Splitting &s,const int mode) const
{
  s.m_y=s.m_t/s.m_Q2/(1.0-s.m_z);
  s.m_x=s.m_z;
  if (s.m_mk2==0.0)
    return s.m_y>0.0 && s.m_y<1.0;
  double muk2(s.m_mk2*s.m_z/s.m_Q2);
  double zp((1.0-s.m_z)/(1.0-s.m_z+muk2));
  return s.m_y>0.0 && s.m_y<zp;
}

Lorentz_IF_123::Lorentz_IF_123(const Kernel_Key &k):
  Lorentz_IF(k)
{
}

double Lorentz_IF_123::Jacobian(const Splitting &s) const
{
  double za(s.m_x), xa(s.m_z2);
  double u(s.m_s*za/s.m_Q2), x(xa+s.m_s/s.m_Q2*za/xa);
  double ct(1.0-2.0/(1.0+xa*xa/za/(s.m_s/s.m_Q2-u*(x/za-1.0))));
  double J((1.0-ct*ct)/xa+(x-za)/(2.0*xa*xa)*sqr(1.0+ct)*((x-u)/xa-2.0));
  return J/2.0*Lorentz_IF::Jacobian(s);
}

int Lorentz_IF_123::Construct(Splitting &s,const int mode) const
{
  Parton *b(NULL);
  if (s.m_kin==0)
    for (size_t i(0);i<s.p_c->Ampl()->size();++i)
      if ((*s.p_c->Ampl())[i]->Beam()==3-s.p_c->Beam()) {
	b=(*s.p_c->Ampl())[i];
	break;
      }
  if (s.m_s<rpa->gen.SqrtAccu()) s.m_s=0.0;
  double xa(s.m_z2), za(s.m_x);
  double Q2(s.m_mij2+s.m_mk2-(s.p_c->Mom()+s.p_s->Mom()).Abs2());
  Kin_Args ff(s.m_s*za/Q2,xa+s.m_s*za/Q2+s.m_t*za/(Q2*xa),s.m_phi,s.m_kin);
  if (ff.m_z>1.0) return -1.0;
  ff.m_mk2=Q2*(xa/za-1.0)+s.m_s+s.m_t/xa;
  if (ConstructIFDipole
      (s.m_mi2,s.m_ml2,s.m_mij2,s.m_mk2,b?p_ms->Mass2(b->Flav()):0.0,
       -s.p_c->Mom(),s.p_s->Mom(),b?-b->Mom():Vec4D(),ff)<0)
    return -1;
  double y2(2.0*(ff.m_pk*(ff.m_pi-ff.m_pj))/(ff.m_mk2-s.m_mk2-s.m_mj2));
  Kin_Args ff2(ff.m_mk2?1.0/(1.0+y2):0.0,1.0-2.0/
	       (1.0+xa*xa/za/(s.m_t/Q2-ff.m_y*(ff.m_z/za-1.0))),s.m_phi2);
  if (ConstructFFDipole
      (s.m_mk2,s.m_mj2,ff.m_mk2,-s.m_s,
       ff.m_pk,ff.m_pi-ff.m_pj,ff2)<0) return -1;
  ff.m_pk=ff2.m_pi;
  ff.m_pi=-ff.m_pi;
  return Update(s,ff,mode,ff2.m_pj);
}

bool Lorentz_IF_123::Cluster(Splitting &s,const int mode) const
{
  return false;
}

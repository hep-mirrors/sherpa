#include "DIM/Shower/Lorentz_FF.H"

#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "DIM/Shower/Shower.H"

using namespace DIM;
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
  Kin_Args ff;
  switch(m_dipole_case){
  case EXTAMP::IDa:
    {
    const double mw2 = sqr(Flavour(24).Mass());
    ATOOLS::Vec4D paitilde = s.p_c->Mom();
    ATOOLS::Vec4D pb       = s.p_s->Mom();
    ATOOLS::Vec4D pwtilde  = s.p_kinspec->Mom();
    /* boost into pai+p_ rest-frame in order to be able to use phi_ib and construct additional
       momenta */
          ATOOLS::Vec4D pminus = pwtilde - mw2/s.m_Qprime2*paitilde;
    const ATOOLS::Vec4D pboost = paitilde+pminus;
    Poincare bst(pboost);
    bst.Boost(paitilde);
    bst.Boost(pwtilde);
    bst.Boost(pminus);
    bst.Boost(pb);
    ff = Kin_Args(1.-s.m_z,1.-s.m_vi,s.m_phi);
    ff.m_res    = true;
    ff.m_pb     = pb;
    ff.m_pminus = pminus;
    if (ConstructFFDipole(mw2,0.,mw2,0.,pwtilde,paitilde,ff)<0)
      return -1;
    bst.BoostBack(ff.m_pj); // pi
    bst.BoostBack(ff.m_pk); // pa
    bst.BoostBack(ff.m_pi); // pw
    break;
    }
  default:
    ff = Kin_Args(s.m_y,s.m_x,s.m_phi);
    if (ConstructFFDipole(s.m_mi2,s.m_mj2,s.m_mij2,
         s.m_mk2,s.p_c->Mom(),s.p_s->Mom(),ff)<0)
      return -1;
  }
  return Update(s,ff,mode);
}

double Lorentz_FF::KT2(const Splitting &s, const double &vi) const {
  const double mw2     = sqr(Flavour(24).Mass());
  const double Qprime2 = s.m_Qprime2;
  const double paipb   = s.m_paipb;
  const double alpha   = s.m_alpha;
  const double cosphi  = cos(s.m_phi);
  const double A       = alpha*vi*Qprime2/2.;
  const double B       = paipb*(vi*(s.m_z-1.-mw2/Qprime2)+1-s.m_z);
  const double kt2     = vi*Qprime2*(A+B-2.*sqrt(A*B)*cosphi)/(2.*paipb);
  return kt2;
}

double Lorentz_FF::GetVimax(const Splitting &s) const
{
  const ATOOLS::Vec4D paitilde = s.p_c->Mom();
  const ATOOLS::Vec4D pwtilde  = s.p_kinspec->Mom();
  const ATOOLS::Vec4D pa       = s.m_z*paitilde;
  const ATOOLS::Vec4D n        = pwtilde - pa + paitilde;
  return 2.*pa*n/(n*n)*(1.-s.m_z)/s.m_z;
}

double Lorentz_FF::GetVi(const Splitting &s, const double &vimax) const
{
  std::vector<double> vi;
  const double mw2     = sqr(Flavour(24).Mass());
  const double Q2      = s.m_Qprime2+mw2;
  const double paipb   = s.m_paipb;
  const double alpha   = s.m_alpha;
  const double kt2     = s.m_t;
  const double z       = s.m_z;
  const double cos2phi = sqr(cos(s.m_phi));

  /* Do not set up Kinematics_FF each time the method is called, but do the following
     to avoid calling time-consuming default-constructor */
  MCATNLO::Kinematics_FF ff = m_ff;
  ff.m_calcV.calculate_helpers(Q2,kt2,z,paipb,alpha,cos2phi);
  const double vi1 = ff.m_calcV.GetV1();
  const double vi2 = ff.m_calcV.GetV2();
  const double vi3 = ff.m_calcV.GetV3();
  const double vi4 = ff.m_calcV.GetV4();

//DEBUG_VAR(m_vimax);
//DEBUG_VAR(vi1);
//DEBUG_VAR(vi2);
//DEBUG_VAR(vi3);
//DEBUG_VAR(vi4);
//DEBUG_VAR(m_kperp2);
//DEBUG_VAR(KT2(vi1));
//DEBUG_VAR(KT2(vi2));
//DEBUG_VAR(KT2(vi3));
//DEBUG_VAR(KT2(vi4));
  if(!IsBad(vi1) && (vi1>0. && vi1<vimax) && IsEqual(KT2(s,vi1),kt2,1e-2))
    vi.push_back(vi1);
  if(!IsBad(vi2) && (vi2>0. && vi2<vimax) && IsEqual(KT2(s,vi2),kt2,1e-2))
    vi.push_back(vi2);
  if(!IsBad(vi3) && (vi3>0. && vi3<vimax) && IsEqual(KT2(s,vi3),kt2,1e-2))
    vi.push_back(vi3);
  if(!IsBad(vi4) && (vi4>0. && vi4<vimax) && IsEqual(KT2(s,vi4),kt2,1e-2))
    vi.push_back(vi4);
//DEBUG_VAR(m_vi);

  /* select one element randomly */
  const double n = vi.size();
  if(n==0) return -1;
  const double r = ran->Get();

  if(r<=1./n)              return vi[0];
  else if(n>=2 && r<=2./n) return vi[1];
  else if(n>=3 && r<=3./n) return vi[2];
  else if(n>=4 && r<=4./n) return vi[3];

  THROW(fatal_error, "Upps.\n");
}

double Lorentz_FF::GetViab(const Splitting &s) const
{
        ATOOLS::Vec4D paitilde = s.p_c->Mom();
        ATOOLS::Vec4D pb       = s.p_s->Mom();
        ATOOLS::Vec4D pwtilde  = s.p_kinspec->Mom();
  const double mw2             = sqr(Flavour(24).Mass());
  const double y_wia           = 1.-s.m_z;
  const double ztilde_w        = 1.-s.m_vi;
  Kin_Args ff(y_wia,ztilde_w,s.m_phi);

  /* boost into pai+p_ rest-frame in order to be able to use phi_ib and construct additional
     momentum */
        ATOOLS::Vec4D pminus   = pwtilde - mw2/s.m_Qprime2*paitilde;
  const ATOOLS::Vec4D pboost   = paitilde + pminus;
  Poincare bst(pboost);
  bst.Boost(paitilde);
  bst.Boost(pwtilde);
  bst.Boost(pminus);
  bst.Boost(pb);
  ff.m_res    = true;
  ff.m_pb     = pb;
  ff.m_pminus = pminus;
  if (ConstructFFDipole(mw2,0.,mw2,0.,pwtilde,paitilde,ff)<0) return -1.;
  const ATOOLS::Vec4D pi = ff.m_pj;
  const ATOOLS::Vec4D pa = ff.m_pk;

  const double kt2 = (pa*pi)*(pb*pi)/(pa*pb);
  if(IsEqual(s.m_t, kt2, 1e-3)) return pa*pb / (pi*(pa+pb));
  return -1.;
}

bool Lorentz_FF::Compute(Splitting &s) const
{
  switch(m_dipole_case){
  case EXTAMP::IDa:
    {
    const double vimax = GetVimax(s);
    s.m_vi             = GetVi(s, vimax);    if(s.m_vi<0.)   return false;
    s.m_viab           = GetViab(s);         if(s.m_viab<0.) return false;
    return s.m_z>0. && s.m_z<1. && s.m_viab>0.;
    }

  default:
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

double Lorentz_FF::MEPSWeight(const Splitting &s) const
{
  switch(m_dipole_case){
  case EXTAMP::IDa:
    return (8.0*M_PI)/(2.*s.p_n->Mom()*s.p_c->Mom())/Jacobian(s);
  default:
    return (8.0*M_PI)/(s.m_Q2*s.m_y)/Jacobian(s);
  }
}

#include "DIRE/Shower/Lorentz_FF.H"

#include "MODEL/Main/Single_Vertex.H"
#include "DIRE/Shower/Kernel.H"
#include "ATOOLS/Math/Random.H"

using namespace ATOOLS;

namespace DIRE {
  
  class FFV_FF: public Lorentz_FF {
  public:

    inline FFV_FF(const Kernel_Key &key):
      Lorentz_FF(key) {}

    double Value(const Splitting &s) const
    {
      double z(s.m_z), y(s.m_y);
      double A1=2.0*(1.0-z)/(sqr(1.0-z)+s.m_t/s.m_Q2);
      double B1=-(1.0+z);
      if (s.m_mij2==0.0 && s.m_mi2==0.0 && s.m_mk2==0.0)
	return A1*(1.0+p_sk->GF()->K(s))+B1;
      double Q2(s.m_Q2+s.m_mi2+s.m_mj2+s.m_mk2);
      double muij2(s.m_mij2/Q2), mui2(s.m_mi2/Q2);
      double muk2(s.m_mk2/Q2), vtijk=Lam(1.0,muij2,muk2);
      double vijk=sqr(2.0*muk2+(1.0-mui2-muk2)*(1.0-y))-4.0*muk2;
      if (vtijk<0.0 || vijk<0.0) return 0.0;
      vtijk=sqrt(vtijk)/(1.0-muij2-muk2);
      vijk=sqrt(vijk)/((1.0-mui2-muk2)*(1.0-y));
      double pipj=s.m_Q2*s.m_y/2.0;
      B1=vtijk/vijk*(B1-s.m_mi2/pipj);
      return A1*(1.0+p_sk->GF()->K(s))+B1;
    }

    double Integral(const Splitting &s) const
    {
      double I=log(1.0+s.m_Q2/s.m_t0);
      return I*(1.0+p_sk->GF()->KMax(s));
    }

    double Estimate(const Splitting &s) const
    {
      double E=2.0*(1.0-s.m_z)/(sqr(1.0-s.m_z)+s.m_t0/s.m_Q2);
      return E*(1.0+p_sk->GF()->KMax(s));
    }

    bool GeneratePoint(Splitting &s) const
    {
      s.m_z=1.0-sqrt(s.m_t0/s.m_Q2*(pow(1.0+s.m_Q2/s.m_t0,ran->Get())-1.0));
      s.m_phi=2.0*M_PI*ran->Get();
      return true;
    }

  };// end of class FFV_FF

  class VFF_FF: public Lorentz_FF {
  public:

    inline VFF_FF(const Kernel_Key &key):
      Lorentz_FF(key) {}

    double Value(const Splitting &s) const
    {
      double z(s.m_z), y(s.m_y);
      if (s.m_mi2==0.0 && s.m_mj2==0.0 && s.m_mk2==0.0) {
	double V=1.0-2.0*z*(1.0-z);
	return V;
      }
      double nui2(s.m_mi2/s.m_Q2), nuk2(s.m_mk2/s.m_Q2);
      double vijj=sqr(s.m_y)-4.0*nui2*nui2;
      double vijk=sqr(1.0-s.m_y)-4.0*(s.m_y+2.0*nui2)*nuk2;
      if (vijj<0.0 || vijk<0.0) return 0.0;
      vijj=sqrt(vijj)/(s.m_y+2.0*nui2);
      vijk=sqrt(vijk)/(1.0-s.m_y);
      double zm=1.0-0.5*(1.0-y)*(1.0+vijj*vijk);
      double zp=1.0-0.5*(1.0-y)*(1.0-vijj*vijk);
      double V=vijk*(1.0-2.0*(z*(1.0-z)-2.0*zp*zm));
      V/=1.0+2.0*nui2/s.m_y;
      return V;
    }

    double Integral(const Splitting &s) const
    {
      return 1.0;
    }

    double Estimate(const Splitting &s) const
    {
      return 1.0;
    }

    bool GeneratePoint(Splitting &s) const
    {
      s.m_z=ran->Get();
      s.m_phi=2.0*M_PI*ran->Get();
      return true;
    }

  };// end of class VFF_FF

}// end of namespace DIRE

using namespace DIRE;

DECLARE_GETTER(FFV_FF,"FF_FFV",Lorentz,Kernel_Key);

Lorentz *ATOOLS::Getter<Lorentz,Kernel_Key,FFV_FF>::
operator()(const Parameter_Type &args) const
{
  if (args.m_type!=0) return NULL;
  if (args.p_v->in[0].IntSpin()==1 &&
      args.p_v->in[1+args.m_mode].IntSpin()==1 &&
      args.p_v->in[2-args.m_mode].IntSpin()==2) {
    return new FFV_FF(args);
  }
  if (args.p_v->in[0].IntSpin()==2 &&
      args.p_v->in[1].IntSpin()==1 &&
      args.p_v->in[2].IntSpin()==1) {
    return new VFF_FF(args);
  }
  return NULL;
}

void ATOOLS::Getter<Lorentz,Kernel_Key,FFV_FF>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"FFV Lorentz Function";
}

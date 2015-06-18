#include "DIRE/Shower/Lorentz_FI.H"

#include "MODEL/Main/Single_Vertex.H"
#include "DIRE/Shower/Shower.H"
#include "ATOOLS/Math/Random.H"

using namespace ATOOLS;

namespace DIRE {
  
  class FFV_FI: public Lorentz_FI {
  private:

    double m_jmax;

  public:

    inline FFV_FI(const Kernel_Key &key):
      Lorentz_FI(key), m_jmax(5.0) {}

    double Value(const Splitting &s) const
    {
      double A1=2.0*(1.0-s.m_z)/(sqr(1.0-s.m_z)+s.m_t/s.m_Q2);
      double B1=-(1.0+s.m_z);
      if (s.m_mij2==0.0 && s.m_mi2==0.0)
	return A1*(1.0+p_sk->GF()->K(s))+B1;
      double pipj=s.m_Q2*(1.0-s.m_y)/s.m_y/2.0;
      B1=B1-s.m_mi2/pipj;
      return A1*(1.0+p_sk->GF()->K(s))+B1;
    }

    double Integral(const Splitting &s) const
    {
      double I=log(1.0+s.m_Q2/s.m_t0);
      return I*(1.0+p_sk->GF()->KMax(s))*m_jmax;
    }

    double Estimate(const Splitting &s) const
    {
      double E=2.0*(1.0-s.m_z)/(sqr(1.0-s.m_z)+s.m_t0/s.m_Q2);
      return E*(1.0+p_sk->GF()->KMax(s))*m_jmax;
    }

    bool GeneratePoint(Splitting &s) const
    {
      s.m_z=1.0-sqrt(s.m_t0/s.m_Q2*(pow(1.0+s.m_Q2/s.m_t0,ran->Get())-1.0));
      s.m_phi=2.0*M_PI*ran->Get();
      return true;
    }

  };// end of class FFV_FI

  class VFF_FI: public Lorentz_FI {
  private:

    double m_jmax;

  public:

    inline VFF_FI(const Kernel_Key &key):
      Lorentz_FI(key), m_jmax(5.0) {}

    double Value(const Splitting &s) const
    {
      double V=1.0-2.0*s.m_z*(1.0-s.m_z);
      if (s.m_mi2==0.0 && s.m_mj2==0.0) return V;
      double nui2(s.m_mi2/s.m_Q2*s.m_y);
      double viji=sqr(1.0-s.m_y)-4.0*nui2*nui2;
      if (viji<0.0 || s.m_y>1.0) return false;
      viji=sqrt(viji)/(1.0-s.m_y+2.0*nui2);
      double zm=0.5*(1.0-viji), zp=0.5*(1.0+viji);
      V=1.0-2.0*(zp-s.m_z)*(s.m_z-zm);
      return V;
    }

    double Integral(const Splitting &s) const
    {
      return m_jmax;
    }

    double Estimate(const Splitting &s) const
    {
      return m_jmax;
    }

    bool GeneratePoint(Splitting &s) const
    {
      s.m_z=ran->Get();
      s.m_phi=2.0*M_PI*ran->Get();
      return true;
    }

  };// end of class VFF_FI

}// end of namespace DIRE

using namespace DIRE;

DECLARE_GETTER(FFV_FI,"FI_FFV",Lorentz,Kernel_Key);

Lorentz *ATOOLS::Getter<Lorentz,Kernel_Key,FFV_FI>::
operator()(const Parameter_Type &args) const
{
  if (args.m_type!=2) return NULL;
  if (args.p_v->in[0].IntSpin()==1 &&
      args.p_v->in[1+args.m_mode].IntSpin()==1 &&
      args.p_v->in[2-args.m_mode].IntSpin()==2) {
    return new FFV_FI(args);
  }
  if (args.p_v->in[0].IntSpin()==2 &&
      args.p_v->in[1].IntSpin()==1 &&
      args.p_v->in[2].IntSpin()==1) {
    return new VFF_FI(args);
  }
  return NULL;
}

void ATOOLS::Getter<Lorentz,Kernel_Key,FFV_FI>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"FFV Lorentz Function";
}

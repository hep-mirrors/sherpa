#include "DIM/Shower/Lorentz_FF.H"

#include "MODEL/Main/Single_Vertex.H"
#include "DIM/Shower/Shower.H"
#include "ATOOLS/Math/Random.H"

using namespace ATOOLS;

namespace DIM {
  
  class VVV_FF: public Lorentz_FF {
  private:

    int m_mode;

  public:

    inline VVV_FF(const Kernel_Key &key,const int mode):
      Lorentz_FF(key), m_mode(mode) {}

    double Value(const Splitting &s) const
    {
      double z(s.m_z), y(s.m_y);
      double A1=2.0*(1.0-z)/(sqr(1.0-z)+s.m_t/s.m_Q2);
      double B1=-2.0+z*(1.0-z);
      if (s.m_mk2==0.0) return A1*(1.0+p_sk->GF()->K(s))+B1;
      double nuk2(s.m_mk2/s.m_Q2), vijk=sqr(1.0-y)-4.0*y*nuk2;
      if (vijk<0.0) return 0.0;
      vijk=sqrt(vijk)/(1.0-y);
      B1=(-2.0+z*(1.0-z))/vijk;
      return A1*(1.0+p_sk->GF()->K(s))+B1;
    }

    double AsymmetryFactor(const Splitting &s) const
    {
      double y(s.m_y), zi(s.m_z), zj(1.0-s.m_x*(1.0-y));
      double A11=2.0*(1.0-zi)/(sqr(1.0-zi)+s.m_t/s.m_Q2), A12=2.0/zi;
      double B11=-2.0+zi*(1.0-zi), B12=-2.0+zj*(1.0-zj);
      if (s.m_mk2==0.0)
	return (A11*(1.0+p_sk->GF()->K(s))+B11)/
	  ((A11+A12)*(1.0+p_sk->GF()->K(s))+(B11+B12));
      double nuk2(s.m_mk2/s.m_Q2), vijk=sqr(1.0-y)-4.0*y*nuk2;
      if (vijk<0.0) return 0.0;
      vijk=sqrt(vijk)/(1.0-y);
      double zm=0.5*(1.0-vijk), zp=0.5*(1.0+vijk);
      B11=(-2.0+zi*(1.0-zi))/vijk;
      B12=(-2.0+zj*(1.0-zj))/vijk;
      return (A11*(1.0+p_sk->GF()->K(s))+B11)/
	((A11+A12)*(1.0+p_sk->GF()->K(s))+(B11+B12));
    }

    double Integral(const Splitting &s) const
    {
      double I=log(1.0+s.m_Q2/s.m_t0);
      return I*(1.0+p_sk->GF()->KMax(s));
    }

    double Estimate(const Splitting &s) const
    {
      double z(s.m_z);
      double E=2.0*(1.0-z)/(sqr(1.0-z)+s.m_t0/s.m_Q2);
      return E*(1.0+p_sk->GF()->KMax(s));
    }

    bool GeneratePoint(Splitting &s) const
    {
      s.m_z=1.0-sqrt(s.m_t0/s.m_Q2*(pow(1.0+s.m_Q2/s.m_t0,ran->Get())-1.0));
      s.m_phi=2.0*M_PI*ran->Get();
      return true;
    }

  };// end of class VVV_FF

}// end of namespace DIM

using namespace DIM;

DECLARE_GETTER(VVV_FF,"FF_VVV",Lorentz,Kernel_Key);

Lorentz *ATOOLS::Getter<Lorentz,Kernel_Key,VVV_FF>::
operator()(const Parameter_Type &args) const
{
  if (args.m_type!=0) return NULL;
  if (args.p_v->in[0].IntSpin()==2 &&
      args.p_v->in[1].IntSpin()==2 &&
      args.p_v->in[2].IntSpin()==2) {
    return new VVV_FF(args,args.m_mode);
  }
  return NULL;
}

void ATOOLS::Getter<Lorentz,Kernel_Key,VVV_FF>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"VVV Lorentz Function";
}

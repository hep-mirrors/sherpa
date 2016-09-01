#include "DIRE/Shower/Lorentz_FF.H"

#include "MODEL/Main/Single_Vertex.H"
#include "DIRE/Shower/Shower.H"
#include "ATOOLS/Math/Random.H"

using namespace ATOOLS;

namespace DIRE {
  
  class VVV_FF: public Lorentz_FF {
  public:

    inline VVV_FF(const Kernel_Key &key):
      Lorentz_FF(key) {}

    double Value(const Splitting &s) const
    {
      double z(s.m_z), y(s.m_y);
      double A=2.0*(1.0-z)/(sqr(1.0-z)+s.m_t/s.m_Q2);
      double B=-2.0+z*(1.0-z);
      if (s.m_kfac&2) {
	double CF=4./3., CA=3., TF=.5*p_sk->GF()->Nf(s), x=p_sk->Mode()?1.0-z:z;
	double B2=TF*(4*(-1+x)*(-23+x*(6+x*(10+x*(4+23*x))))+24*(1+x)*(2+(-1+x)*x*(3+x*(-3+2*x)))*log(x))+
	  (CF*TF*(-12*(1+x)*(8+x*(7-x*(2+x)*(-3+8*x)))*log(x)-8*(1+x)*(23+x*(14+41*x))*sqr(-1+x)+
		  36*(-1+x)*x*sqr(1+x)*sqr(log(x))))/CA+72*CA*(-1+x)*DiLog(1/(1+x))*sqr(1+x+sqr(x))+
	  CA*(-6*(1+x)*(-22+x*(11+x*(30+x*(-19+22*x))))*log(x)+
	      (1-x)*(x*(1+x)*(25+109*x)+6*(2+x*(1+2*x*(1+x)))*sqr(M_PI))-72*(1+x)*log(1-x)*log(x)*sqr(1+(-1+x)*x)+
	      36*(2+x*(1+(-4+x)*(-1+x)*x*(1+x)))*sqr(log(x))+36*(-1+x)*sqr(log(1+x))*sqr(1+x+sqr(x)));
	B2-=2.*(x*x-1.0)*40*TF/(1.0+x*x/(s.m_t/s.m_Q2));
	B+=p_sk->GF()->Coupling(s)/(2.0*M_PI)*B2/(18.0*x*(x*x-1.0))/2.0;
      }
      if (s.m_mk2==0.0)
	return (p_sk->Mode()?1.0-z:z)*(A*(1.0+p_sk->GF()->K(s))+B);
      double nuk2(s.m_mk2/s.m_Q2), vijk=sqr(1.0-y)-4.0*y*nuk2;
      if (vijk<0.0) return 0.0;
      vijk=sqrt(vijk)/(1.0-y);
      B=(-2.0+z*(1.0-z))/vijk;
      return (p_sk->Mode()?1.0-z:z)*(A*(1.0+p_sk->GF()->K(s))+B);
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

}// end of namespace DIRE

using namespace DIRE;

DECLARE_GETTER(VVV_FF,"FF_VVV",Lorentz,Kernel_Key);

Lorentz *ATOOLS::Getter<Lorentz,Kernel_Key,VVV_FF>::
operator()(const Parameter_Type &args) const
{
  if (args.m_type!=0) return NULL;
  if (args.m_swap) return NULL;
  if (args.p_v->in[0].IntSpin()==2 &&
      args.p_v->in[1].IntSpin()==2 &&
      args.p_v->in[2].IntSpin()==2) {
    return new VVV_FF(args);
  }
  return NULL;
}

void ATOOLS::Getter<Lorentz,Kernel_Key,VVV_FF>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"VVV Lorentz Function";
}

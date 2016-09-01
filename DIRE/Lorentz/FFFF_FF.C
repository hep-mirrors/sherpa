#include "DIRE/Shower/Lorentz_FF.H"

#include "MODEL/Main/Single_Vertex.H"
#include "DIRE/Shower/Kernel.H"
#include "ATOOLS/Math/Random.H"

using namespace ATOOLS;

namespace DIRE {
  
  class FFFF_FF: public Lorentz_FF_123 {
  public:

    inline FFFF_FF(const Kernel_Key &key):
      Lorentz_FF_123(key) {}

    double Value(const Splitting &s) const
    {
      if (m_fl[1].Kfcode()>p_sk->GF()->Nf(s)) return 0.0; 
      double x(s.m_z), CF(4.0/3.0), CA(3.0), TF(0.5), B2(0.0);
      if (m_fl[1]==m_fl[0].Bar()) {
	B2=TF*((4*(-1+x)*(5+x*(23+14*x)))/(9.*x)+log(x)*(-5-9*x+(1+x)*log(x)-(8*sqr(x))/3.))+
	  (-CA/2.+CF)*(4-4*x+2*(1+x)*log(x)+(1+sqr(x))*(-12*DiLog(1/(1+x))+sqr(M_PI)+
							3*sqr(log(x))-6*sqr(log(1+x)))/(3.*(1+x)));
      }
      else {
	B2=TF*((4*(-1+x)*(5+x*(23+14*x)))/(9.*x)
	       +log(x)*(-5-9*x+(1+x)*log(x)-(8*sqr(x))/3.));
      }
      B2+=20/(9*x)*TF/(1.0+x*x/(s.m_t/s.m_Q2));
      B2*=p_sk->GF()->Coupling(s)/(2.0*M_PI);
      return x*B2;
    }

    double Integral(const Splitting &s) const
    {
      double I=20.0/9.0*0.5*0.5*log((s.m_Q2+s.m_t0)/s.m_t0);
      return I*p_sk->GF()->CplMax(s)/(2.0*M_PI);
    }

    double Estimate(const Splitting &s) const
    {
      double E=20.0/9.0*0.5*s.m_z/(sqr(s.m_z)+s.m_t0/s.m_Q2);
      return E*p_sk->GF()->CplMax(s)/(2.0*M_PI);
    }

    bool GeneratePoint(Splitting &s) const
    {
      double k2(s.m_t0/s.m_Q2);
      s.m_z=sqrt(pow((1.0+k2)/k2,-ran->Get())*(1.0+k2)-k2);
      s.m_phi=2.0*M_PI*ran->Get();
      do s.m_z2=pow(s.m_z,ran->Get());
      while (sqr(s.m_z2)+sqr(1.0-s.m_z2)<ran->Get());
      s.m_phi2=2.*M_PI*ran->Get();
      return true;
    }

    bool Compute(Splitting &s,const int mode) const
    {
      s.m_x=s.m_z/s.m_z2;
      s.m_y=s.m_t/(s.m_Q2*s.m_x*(1.0-s.m_x));
      s.m_s=sqr(p_ms->Mass(m_fl[1])+p_ms->Mass(m_fl[3]));
      if (!Lorentz_FF::Compute(s,0)) return false;
      Splitting c(s);
      double amax(p_sk->GF()->CplMax(s));
      double I(amax/(2.0*M_PI)*0.5*(log(1.0/s.m_z)-sqr(1.0-s.m_z)));
      for (s.m_s=s.m_Q2*pow(ran->Get(),Max(1.0/I,1.0e-3));
      	   s.m_s>s.m_t0;s.m_s*=pow(ran->Get(),Max(1.0/I,1.0e-3))) {
      	s.m_mi2=c.m_t=s.m_s;
      	if (p_sk->GF()->Coupling(c)/amax<ran->Get()) continue;
      	if (Lorentz_FF::Compute(s,0)) {
      	  s.m_mi2=p_ms->Mass2(m_fl[1]);
      	  return true;
      	}
      }
      s.m_mi2=p_ms->Mass2(m_fl[1]);
      s.m_s=sqr(p_ms->Mass(m_fl[1])+p_ms->Mass(m_fl[3]));
      return Lorentz_FF::Compute(s,0);
    }

  };// end of class FFFF_FF

}// end of namespace DIRE

using namespace DIRE;

DECLARE_GETTER(FFFF_FF,"FF_FFFF",Lorentz,Kernel_Key);

Lorentz *ATOOLS::Getter<Lorentz,Kernel_Key,FFFF_FF>::
operator()(const Parameter_Type &args) const
{
  if (args.m_type!=0 || args.p_v) return NULL;
  if (args.m_lfid=="FFFF") {
    return new FFFF_FF(args);
  }
  return NULL;
}

void ATOOLS::Getter<Lorentz,Kernel_Key,FFFF_FF>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"FFFF Lorentz Function";
}

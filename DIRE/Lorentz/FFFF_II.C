#include "DIRE/Shower/Lorentz_II.H"

#include "MODEL/Main/Single_Vertex.H"
#include "DIRE/Shower/Kernel.H"
#include "ATOOLS/Math/Random.H"

using namespace ATOOLS;

namespace DIRE {
  
  class FFFF_II: public Lorentz_II_123 {
  private:

    double m_jmax;

  public:

    inline FFFF_II(const Kernel_Key &key):
      Lorentz_II_123(key), m_jmax(5.0) {}

    double Value(const Splitting &s) const
    {
      if (m_fl[1].Kfcode()>p_sk->GF()->Nf(s)) return 0.0; 
      double x(s.m_z), CF(4.0/3.0), CA(3.0), TF(0.5), B2(0.0);
      if (m_fl[1]==m_fl[0].Bar()) {
	B2=TF*(3*x*log(x)*(3+x*(15+8*x)-3*(1+x)*log(x))-2*(-1+x)*(10+x+28*sqr(x)))/(9.*x)+
	  (-CA/2.+CF)*(4-4*x+2*(1+x)*log(x)+(1+sqr(x))*(-12*DiLog(1/(1+x))+sqr(M_PI)+
							3*sqr(log(x))-6*sqr(log(1+x)))/(3.*(1+x)));
      }
      else {
	B2=TF*(3*x*log(x)*(3+x*(15+8*x)-3*(1+x)*log(x))-2*(-1+x)*(10+x+28*sqr(x)))/(9.*x);
      }
      B2+=20/(9*x)*TF/(1.0+x*x/(s.m_t/s.m_Q2));
      B2*=p_sk->GF()->Coupling(s)/(2.0*M_PI);
      return x*B2;
    }

    double Integral(const Splitting &s) const
    {
      double I=20.0/9.0*0.5*0.5*log((s.m_Q2+s.m_t0)/(s.m_Q2*sqr(s.m_eta)+s.m_t0));
      return I*p_sk->GF()->CplMax(s)/(2.0*M_PI)*m_jmax;
    }

    double Estimate(const Splitting &s) const
    {
      double E=20.0/9.0*0.5*s.m_z/(sqr(s.m_z)+s.m_t0/s.m_Q2);
      return E*p_sk->GF()->CplMax(s)/(2.0*M_PI)*m_jmax;
    }

    bool GeneratePoint(Splitting &s) const
    {
      double k2(s.m_t0/s.m_Q2);
      s.m_z=sqrt(pow((1.0+k2)/(sqr(s.m_eta)+k2),-ran->Get())*(1.0+k2)-k2);
      s.m_phi=2.0*M_PI*ran->Get();
      do s.m_z2=pow(s.m_z,ran->Get());
      while (1.0+sqr(1.0-s.m_z2)<2.0*ran->Get());
      s.m_phi2=2.*M_PI*ran->Get();
      return true;
    }

    bool Compute(Splitting &s,const int mode) const
    {
      s.m_x=s.m_z;
      s.m_y=s.m_t/s.m_Q2*sqr(s.m_z2/s.m_x);
      s.m_s=s.m_z2*(p_ms->Mass(m_fl[1])-p_ms->Mass(m_fl[3])/(1.0-s.m_z2));
      if (!Lorentz_II::Compute(s,0)) return false;
      Splitting c(s);
      double amax(p_sk->GF()->CplMax(s));
      double I(amax/(2.0*M_PI)*4.0/3.0*
	       (2.0*log(1.0/s.m_z)-(3.0-s.m_z)*(1.0-s.m_z)/2.0));
      for (s.m_s=s.m_Q2*pow(ran->Get(),Max(1.0/I,1.0e-3));
      	   s.m_s>s.m_t0;s.m_s*=pow(ran->Get(),Max(1.0/I,1.0e-3))) {
      	if (p_sk->GF()->Coupling(c)/amax<ran->Get()) continue;
      	if (Lorentz_II::Compute(s,0)) return true;
      }
      s.m_s=s.m_z2*(p_ms->Mass(m_fl[1])-p_ms->Mass(m_fl[3])/(1.0-s.m_z2));
      return Lorentz_II::Compute(s,0);
    }

  };// end of class FFFF_II

}// end of namespace DIRE

using namespace DIRE;

DECLARE_GETTER(FFFF_II,"II_FFFF",Lorentz,Kernel_Key);

Lorentz *ATOOLS::Getter<Lorentz,Kernel_Key,FFFF_II>::
operator()(const Parameter_Type &args) const
{
  if (args.m_type!=3 || args.p_v) return NULL;
  if (args.m_lfid=="FFFF") {
    return new FFFF_II(args);
  }
  return NULL;
}

void ATOOLS::Getter<Lorentz,Kernel_Key,FFFF_II>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"FFFF Lorentz Function";
}

#include "DIRE/Shower/Lorentz_II.H"

#include "MODEL/Main/Single_Vertex.H"
#include "DIRE/Shower/Shower.H"
#include "ATOOLS/Math/Random.H"

using namespace ATOOLS;

namespace DIRE {
  
  class FFV_II: public Lorentz_II {
  private:

    double m_jmax;

  public:

    inline FFV_II(const Kernel_Key &key):
      Lorentz_II(key), m_jmax(m_fl[0].Kfcode()<3?5.0:2.0) {}

    double Value(const Splitting &s) const
    {
      double A=2.0*(1.0-s.m_z)/(sqr(1.0-s.m_z)+s.m_t/s.m_Q2);
      double B=-(1.0+s.m_z);
      if (s.m_kfac&2) {
	double CF=4./3., CA=3., TF=.5*p_sk->GF()->Nf(s), x=s.m_z;
	double B2=(-1+x)*(-8*TF*(-5+(-1+x)*x*(-5+14*x))+x*(90*CF*(-1+x)+CA*(53-187*x+3*(1+x)*sqr(M_PI))))+
	  3*x*log(x)*(-2*(TF+CF*(-9+6*(-1+x)*x)+TF*x*(12-x*(9+8*x)))+12*CF*log(1-x)*(1+sqr(x))-CA*(17+5*sqr(x)))-
	  9*x*(CA-CF-2*TF+(CA+CF+2*TF)*sqr(x))*sqr(log(x));
	B2-=(x-1.)*40*TF/(1.0+x*x/(s.m_t/s.m_Q2));
	B+=p_sk->GF()->Coupling(s)/(2.0*M_PI)*B2/(18.*x*(x-1.));
      }
      return A*(1.0+p_sk->GF()->K(s))+B;
    }

    double Integral(const Splitting &s) const
    {
      double I=log(1.0+sqr(1.0-s.m_eta)*s.m_Q2/s.m_t0);
      return I*(1.0+p_sk->GF()->KMax(s))*m_jmax;
    }

    double Estimate(const Splitting &s) const
    {
      double E=2.0*(1.0-s.m_z)/(sqr(1.0-s.m_z)+s.m_t0/s.m_Q2);
      return E*(1.0+p_sk->GF()->KMax(s))*m_jmax;
    }

    bool GeneratePoint(Splitting &s) const
    {
      double k2(s.m_t0/s.m_Q2);
      s.m_z=1.0-sqrt(k2*(pow(1.0+sqr(1.0-s.m_eta)/k2,ran->Get())-1.0));
      s.m_phi=2.0*M_PI*ran->Get();
      return true;
    }

  };// end of class FFV_II

  class FVF_II: public Lorentz_II {
  private:

    double m_jmax;

  public:

    inline FVF_II(const Kernel_Key &key):
      Lorentz_II(key), m_jmax(5.0) {}

    double Value(const Splitting &s) const
    {
      double B=2.0*s.m_z/(sqr(s.m_z)+s.m_t/s.m_Q2)-(2.0-s.m_z);
      if (s.m_kfac&2) {
	double CF=4./3., CA=3., TF=.5*p_sk->GF()->Nf(s), x=s.m_z;
	double B2=-9*CF*x*(5+7*x)-16*TF*(5+x*(-5+4*x))+36*CA*(2+x*(2+x))*DiLog(1/(1+x))+
	  2*CA*(9+x*(19+x*(37+44*x))-3*sqr(M_PI)*(2+sqr(x)))+
	  3*(-2*log(1-x)*(CA*(-22+(22-17*x)*x)+4*TF*(2+(-2+x)*x)+3*CF*(6+x*(-6+5*x))+6*CA*(2+(-2+x)*x)*log(x))+
	     x*log(x)*(3*CF*(4+7*x)-2*CA*(36+x*(15+8*x))+3*(CF*(-2+x)+2*CA*(2+x))*log(x))+
	     6*(CA-CF)*(2+(-2+x)*x)*sqr(log(1-x))+6*CA*(2+x*(2+x))*sqr(log(1+x)));
	B2+=2.*40.*TF/(1.0+x*x/(s.m_t/s.m_Q2));
	B+=p_sk->GF()->Coupling(s)/(2.0*M_PI)*B2/(18.*x);
      }
      return B;
    }

    double Integral(const Splitting &s) const
    {
      double I=log((s.m_Q2+s.m_t0)/(s.m_Q2*sqr(s.m_eta)+s.m_t0));
      return I*m_jmax*PDFEstimate(s);
    }

    double Estimate(const Splitting &s) const
    {
      double E=2.0*s.m_z/(sqr(s.m_z)+s.m_t0/s.m_Q2);
      return E*m_jmax*PDFEstimate(s);
    }

    bool GeneratePoint(Splitting &s) const
    {
      double k2(s.m_t0/s.m_Q2);
      s.m_z=sqrt(pow((1.0+k2)/(sqr(s.m_eta)+k2),-ran->Get())*(1.0+k2)-k2);
      s.m_phi=2.0*M_PI*ran->Get();
      return true;
    }

  };// end of class FVF_II

  class VFF_II: public Lorentz_II {
  private:

    double m_jmax;

  public:

    inline VFF_II(const Kernel_Key &key):
      Lorentz_II(key), m_jmax(5.0) {}

    double Value(const Splitting &s) const
    {
      double B=1.0-2.0*s.m_z*(1.0-s.m_z);
      if (s.m_kfac&2) {
	double CF=4./3., CA=3., x=s.m_z;
	double B2=CF*(4-9*x+4*log(1-x)+(-1+4*x)*log(x)-(2*(1+2*(-1+x)*x)*(-15-3*(-2+log(-1+1/x))*log(-1+1/x)+sqr(M_PI)))/3.+(-1+2*x)*sqr(log(x)))
	  +(2*CA*(20-18*x*(1+2*x*(1+x))*DiLog(1/(1+x))+x*(-18+(225-218*x)*x+sqr(M_PI)*(3+6*sqr(x)))+
		  3*x*(12*(-1+x)*x*log(1-x)+log(x)*(3+4*x*(6+11*x)-3*(1+2*x)*log(x))+(-3-6*(-1+x)*x)*sqr(log(1-x))-
		       3*(1+2*x*(1+x))*sqr(log(1+x)))))/(9.*x);
	B2-=40*CA/(9.*x)/(1.0+x*x/(s.m_t/s.m_Q2));
	B+=p_sk->GF()->Coupling(s)/(2.0*M_PI)*B2/2.0;
      }
      return B;
    }

    double Integral(const Splitting &s) const
    {
      return (1.0-s.m_eta)*m_jmax*PDFEstimate(s);
    }

    double Estimate(const Splitting &s) const
    {
      return m_jmax*PDFEstimate(s);
    }

    bool GeneratePoint(Splitting &s) const
    {
      s.m_z=s.m_eta+(1.0-s.m_eta)*ran->Get();
      s.m_phi=2.0*M_PI*ran->Get();
      return true;
    }

  };// end of class VFF_II

}// end of namespace DIRE

using namespace DIRE;

DECLARE_GETTER(FFV_II,"II_FFV",Lorentz,Kernel_Key);

Lorentz *ATOOLS::Getter<Lorentz,Kernel_Key,FFV_II>::
operator()(const Parameter_Type &args) const
{
  if (args.m_type!=3) return NULL;
  if ((args.m_mode==0 &&
       args.p_v->in[0].IntSpin()==1 &&
       args.p_v->in[1].IntSpin()==1 &&
       args.p_v->in[2].IntSpin()==2) ||
      (args.m_mode==1 &&
       args.p_v->in[0].IntSpin()==1 &&
       args.p_v->in[2].IntSpin()==1 &&
       args.p_v->in[1].IntSpin()==2)) {
    return new FFV_II(args);
  }
  if ((args.m_mode==0 &&
       args.p_v->in[0].IntSpin()==1 &&
       args.p_v->in[1].IntSpin()==2 &&
       args.p_v->in[2].IntSpin()==1) ||
      (args.m_mode==1 &&
       args.p_v->in[0].IntSpin()==1 &&
       args.p_v->in[2].IntSpin()==2 &&
       args.p_v->in[1].IntSpin()==1)) {
    return new VFF_II(args);
  }
  if (args.p_v->in[0].IntSpin()==2 &&
      args.p_v->in[1].IntSpin()==1 &&
      args.p_v->in[2].IntSpin()==1) {
    return new FVF_II(args);
  }
  return NULL;
}

void ATOOLS::Getter<Lorentz,Kernel_Key,FFV_II>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"FFV Lorentz Function";
}

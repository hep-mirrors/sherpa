#include "DIRE/Shower/Lorentz_FI.H"

#include "MODEL/Main/Single_Vertex.H"
#include "DIRE/Shower/Shower.H"
#include "ATOOLS/Math/Random.H"

using namespace ATOOLS;

namespace DIRE {
  
  class FFV_FI: public Lorentz_FI {
  private:

    int    m_swap;
    double m_jmax;

  public:

    inline FFV_FI(const Kernel_Key &key):
      Lorentz_FI(key), m_swap(key.m_swap), m_jmax(5.0) {}

    double Value(const Splitting &s) const
    {
      double z(s.m_z);
      double A=2.0*(1.0-z)/(sqr(1.0-z)+s.m_t/s.m_Q2);
      double B=-(1.0+z);
      if (s.m_kfac&2) {
	if (m_swap) {
	  B+=A;
	  A=0.0;
	  double CF=4./3., CA=3., TF=.5*p_sk->GF()->Nf(s), x=1.0-s.m_z;
	  double B2=9*CF*x*(-1+9*x)+144*(CA-CF)*(2+(-2+x)*x)*DiLog(x)+36*CA*(2+x*(2+x))*DiLog(1/(1+x))-
	    2*CA*(-17+9*(-5+x)*x+44*pow(x,3)+3*sqr(M_PI)*(2+sqr(x)))+
	    3*(12*log(1-x)*((3*CA-2*CF)*(2+(-2+x)*x)*log(x)+(-CA+CF)*sqr(x))+
	       log(x)*(3*CF*(-16+x)*x+2*CA*(-18+x*(24+x*(27+8*x)))-3*log(x)*(CF*(-2+x)*x+CA*(8+4*x+6*sqr(x))))-
	       6*(CA-CF)*(2+(-2+x)*x)*sqr(log(1-x))+6*CA*(2+x*(2+x))*sqr(log(1+x)));
	  B2-=40*TF/(1.0+x*x/(s.m_t/s.m_Q2));
	  B+=p_sk->GF()->Coupling(s)/(2.0*M_PI)*B2/(18.*x);
	}
	else {
	  double CF=4./3., CA=3., TF=.5*p_sk->GF()->Nf(s), x=s.m_z;
	  double B2=(-1+x)*(4*TF*(-10+x*(-37+x*(29+28*x)))+x*(90*CF*(-1+x)+CA*(53-187*x+3*(1+x)*sqr(M_PI))))+
	    3*x*log(x)*(34*TF+12*(CF-CF*x+2*TF*x)-2*(9*CF+TF*(17+8*x))*sqr(x)-12*CF*log(1-x)*(1+sqr(x))-
			CA*(17+5*sqr(x))-3*log(x)*(CA-3*CF+2*TF+(CA-5*CF-2*TF)*sqr(x)));
	  B2+=(x-1.)*40*TF/(1.0+x*x/(s.m_t/s.m_Q2));
	  B+=p_sk->GF()->Coupling(s)/(2.0*M_PI)*B2/(18.*x*(x-1.0));
	}
      }
      if (s.m_mij2==0.0 && s.m_mi2==0.0)
	return (m_swap?1.0-z:z)*(A*(1.0+p_sk->GF()->K(s))+B);
      double pipj=s.m_Q2*(1.0-s.m_y)/s.m_y/2.0;
      B=B-s.m_mi2/pipj;
      return (m_swap?1.0-z:z)*(A*(1.0+p_sk->GF()->K(s))+B);
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

    int    m_swap;
    double m_jmax;

  public:

    inline VFF_FI(const Kernel_Key &key):
      Lorentz_FI(key), m_swap(key.m_swap), m_jmax(5.0) {}

    double Value(const Splitting &s) const
    {
      double z(s.m_z);
      double B=1.0-2.0*z*(1.0-z);
      if (s.m_mi2==0.0 && s.m_mj2==0.0) {
	if (s.m_kfac&2) {
	  double CF=4./3., CA=3., TF=.5*p_sk->GF()->Nf(s), x=m_swap?1.0-z:z;
	  double B2=TF*(-8/3.-(8*(1+2*(-1+x)*x)*(2+3*log(1-x)+3*log(x)))/9.)+
	    CF*(-2+3*x-4*log(1-x)+(-7+8*x)*log(x)+(1-2*x)*sqr(log(x))-
		(2*(1+2*(-1+x)*x)*(15-24*DiLog(x)+3*log(-1+1/x)-24*log(1-x)*log(x)+sqr(M_PI)+3*sqr(log(-(-1+x)*x))))/3.)+
	    (CA*(-152-40/x+166*x+36*log(1-x)-12*(1+19*x)*log(x)+
		 (1+2*(-1+x)*x)*(178-144*DiLog(x)+log(1-x)*(30-72*log(x))-3*log(x)*(4+3*log(x))+3*sqr(M_PI)+
				 18*sqr(log(1-x)))+9*(2+8*x)*sqr(log(x))+
		 3*(1+2*x*(1+x))*(-12*DiLog(1/(1+x))+sqr(M_PI)+3*sqr(log(x))-6*sqr(log(1+x)))))/9.;
	  B2+=40*CA/(9.*x)/(1.0+x*x/(s.m_t/s.m_Q2));
	  B+=p_sk->GF()->Coupling(s)/(2.0*M_PI)*B2/2.0;
	}
	return (m_swap?1.0-z:z)*B;
      }
      double nui2(s.m_mi2/s.m_Q2*s.m_y);
      B=1.0-2.0*s.m_z*(1.0-s.m_z)+nui2/((1.0-s.m_y)/2.0+nui2);
      return (m_swap?1.0-z:z)*B;
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
  if (args.m_swap) return NULL;
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

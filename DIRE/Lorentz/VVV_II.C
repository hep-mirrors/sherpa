#include "DIRE/Shower/Lorentz_II.H"

#include "MODEL/Main/Single_Vertex.H"
#include "DIRE/Shower/Shower.H"
#include "ATOOLS/Math/Random.H"

using namespace ATOOLS;

namespace DIRE {
  
  class VVV_II: public Lorentz_II {
  private:

    double m_jmax;

    int m_mode;

  public:

    inline VVV_II(const Kernel_Key &key,const int mode):
      Lorentz_II(key), m_jmax(1.0), m_mode(mode) {}

    double Value(const Splitting &s) const
    {
      double A=0.0, B=0.0;
      if (m_mode) B=2.0*s.m_z*(1.0-s.m_z)+s.m_z/(sqr(s.m_z)+s.m_t/s.m_Q2)-1.0;
      else {
	A=2.0*(1.0-s.m_z)/(sqr(1.0-s.m_z)+s.m_t/s.m_Q2);
	B=-2.0+s.m_z/(sqr(s.m_z)+s.m_t/s.m_Q2)-1.0;
	A*=1.0+p_sk->GF()->K(s);
      }
      if (s.m_kfac&2) {
	double CF=4./3., CA=3., TF=.5*p_sk->GF()->Nf(s), x=s.m_z;
	double B2=TF*(4*(-1+x)*(-23+x*(6+x*(10+x*(4+23*x))))+24*(1-x)*x*log(x)*sqr(1+x))+
	  (CF*TF*(36*(1-x)*x*(1+x)*(3+5*x)*log(x)+24*(1+x)*(-1+x*(11+5*x))*sqr(-1+x)-36*(-1+x)*x*sqr(1+x)*sqr(log(x))))/CA+
	  CA*(6*(1-x)*x*(1+x)*(25+11*x*(-1+4*x))*log(x)+(1-x)*(x*(1+x)*(25+109*x)+6*(2+x*(1+2*x*(1+x)))*sqr(M_PI))+
	      72*(1+x)*log(1-x)*log(x)*sqr(1+(-1+x)*x)-36*x*sqr(log(x))*sqr(1+x-sqr(x))+36*(-1+x)*sqr(log(1+x))*sqr(1+x+sqr(x)));
	B2+=72*CA*(x-1)*DiLog(1/(1+x))*sqr(1+x+x*x);
	B2+=2.*(x*x-1.0)*40.*TF/(1.0+x*x/(s.m_t/s.m_Q2));
	B+=p_sk->GF()->Coupling(s)/(2.0*M_PI)*B2/(18.0*x*(x*x-1.0))/2.0;
      }
      return A+B;
    }

    double Integral(const Splitting &s) const
    {
      if (m_mode) 
	return 0.5*log((s.m_Q2+s.m_t0)/(s.m_Q2*sqr(s.m_eta)+s.m_t0))*m_jmax;
      double k2=s.m_t0/s.m_Q2;
      double I=log((k2+sqr(1.0-s.m_eta))/(s.m_eta*k2));
      return I*(1.0+p_sk->GF()->KMax(s))*m_jmax;
    }

    double Estimate(const Splitting &s) const
    {
      if (m_mode) return s.m_z/(sqr(s.m_z)+s.m_t0/s.m_Q2)*m_jmax;
      double E=2.0*(1.0-s.m_z)/(sqr(1.0-s.m_z)+s.m_t0/s.m_Q2)+1.0/s.m_z;
      return E*(1.0+p_sk->GF()->KMax(s))*m_jmax;
    }

    bool GeneratePoint(Splitting &s) const
    {
      double k2(s.m_t0/s.m_Q2);
      if (m_mode) {
	s.m_z=sqrt(pow((1.0+k2)/(sqr(s.m_eta)+k2),-ran->Get())*(1.0+k2)-k2);
      }
      else {
	double FF((k2+sqr(1.0-s.m_eta))/(s.m_eta*k2));
	FF=1.0+k2/2.0*pow(FF,ran->Get());
	s.m_z=FF-sqrt(FF*FF-(1.0+k2));
      }
      s.m_phi=2.0*M_PI*ran->Get();
      return true;
    }

  };// end of class VVV_II

}// end of namespace DIRE

using namespace DIRE;

DECLARE_GETTER(VVV_II,"II_VVV",Lorentz,Kernel_Key);

Lorentz *ATOOLS::Getter<Lorentz,Kernel_Key,VVV_II>::
operator()(const Parameter_Type &args) const
{
  if (args.m_type!=3) return NULL;
  if (args.m_swap) return NULL;
  if (args.p_v->in[0].IntSpin()==2 &&
      args.p_v->in[1].IntSpin()==2 &&
      args.p_v->in[2].IntSpin()==2) {
    return new VVV_II(args,args.m_mode);
  }
  return NULL;
}

void ATOOLS::Getter<Lorentz,Kernel_Key,VVV_II>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"VVV Lorentz Function";
}

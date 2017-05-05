#include "DIRE/Shower/Lorentz_FI.H"

#include "MODEL/Main/Single_Vertex.H"
#include "DIRE/Shower/Kernel.H"
#include "ATOOLS/Math/Random.H"

using namespace ATOOLS;

namespace DIRE {
  
  class FFFF_FI: public Lorentz_FI_123 {
  private:

    double m_jmax;

  public:

    inline FFFF_FI(const Kernel_Key &key):
      Lorentz_FI_123(key), m_jmax(5.0) {}

    double Value(const Splitting &s) const
    {
      if (m_fl[1].Kfcode()>p_sk->GF()->Nf(s)) return 0.0; 
      double TR(0.5), B2(0.0);
      double s123(s.m_t*s.m_z2/s.m_z+s.m_s+s.m_mj2);
      double z1(s.m_z), z2(s.m_z/s.m_z2-z1), z3(1.0-z1-z2);
      if (s.m_mode==0) {// differential kernels
	if (IsZero(s.m_s)) return 0.0;
	double s13((s.m_pi+s.m_pj).Abs2()), s23((s.m_pl+s.m_pj).Abs2());
	double s12(s.m_s), t123((2.0*(z1*s23-z2*s13)+(z1-z2)*s12)/(z1+z2));
	double cp13(CosPhi(s.m_pi,s.m_pj,s.m_pl,s.m_pk));
	B2=0.5*TR*s123/s12*(-sqr(t123)/s12/s123+(4.0*z3+sqr(z1-z2))/(z1+z2)+(z1+z2-s12/s123));
	B2-=TR*s123/s12*(1.0+z3*z3)/(1.0-z3)*(1-2.0*z1*z2/sqr(z1+z2));
	B2-=TR*s123/s12*4.0*z1*z2*z3/(1.0-z3)/sqr(z1+z2)*(1.0-2.0*sqr(cp13));
	if (m_fl[1]==m_fl[0].Bar()) {
	  double t132((2.0*(z1*s23-z3*s12)+(z1-z3)*s13)/(z1+z3));
	  double cp12(CosPhi(s.m_pi,s.m_pl,s.m_pj,s.m_pk)), CF(4.0/3.0), CA(3.0);
	  B2+=0.5*TR*s123/s13*(-sqr(t132)/s13/s123+(4.0*z2+sqr(z1-z3))/(z1+z3)+(z1+z3-s13/s123));
	  B2-=TR*s123/s13*(1.0+z2*z2)/(1.0-z2)*(1-2.0*z1*z3/sqr(z1+z3));
	  B2-=TR*s123/s13*4.0*z1*z2*z3/(1.0-z2)/sqr(z1+z3)*(1.0-2.0*sqr(cp12));
	  B2+=(CF-0.5*CA)*(2.*s23/s12+s123/s12*((1+z1*z1)/(1-z2)-2.0*z2/(1.0-z3))
			   +2.*s23/s13+s123/s13*((1+z1*z1)/(1-z3)-2.0*z3/(1.0-z2))
			   -s123*s123/(s12*s13)*z1*(1.0+z1*z1)/(1.0-z2)/(1.0-z3));
	}
      }
      else {// subtraction terms
	B2=TR*((1.0+z3*z3)/(1.0-z3)+(1.0-2.0*z1*z2/sqr(z1+z2))*(1.0-z3+(1.0+z3*z3)/(1.0-z3)*(log(z1*z2*z3)-1.0)));
	B2-=2.0*TR*((1.0+z3*z3)/(1.0-z3)*log(z3*(1.0-z3))+1.0-z3)*(1.0-2.0*z1*z2/sqr(z1+z2));
	if (m_fl[1]==m_fl[0].Bar()) {
	  B2+=TR*((1.0+z2*z2)/(1.0-z2)+(1.0-2.0*z1*z3/sqr(z1+z3))*(1.0-z2+(1.0+z2*z2)/(1.0-z2)*(log(z1*z2*z3)-1.0)));
	  B2-=2.0*TR*((1.0+z2*z2)/(1.0-z2)*log(z2*(1.0-z2))+1.0-z2)*(1.0-2.0*z1*z3/sqr(z1+z3));
	}
      }
      // summation and phase-space weight
      B2*=2.0*log(1.0/z1)/(1.0-s.m_s/s123);
      B2*=p_sk->GF()->Coupling(s)/(2.0*M_PI);
      // desymmetrization
      if (m_fl[1]==m_fl[0].Bar()) B2*=(s.m_z/s.m_z2-s.m_z)/(1.0-s.m_z);
      return s.m_z*B2;
    }

    double Integral(const Splitting &s) const
    {
      double I=20.0/9.0*0.5*log((s.m_Q2+s.m_t0)/s.m_t0);
      return I*p_sk->GF()->CplMax(s)/(2.0*M_PI)*m_jmax;
    }

    double Estimate(const Splitting &s) const
    {
      double E=20.0/9.0*0.5/(s.m_z+s.m_t0/s.m_Q2);
      return E*p_sk->GF()->CplMax(s)/(2.0*M_PI)*m_jmax;
    }

    bool GeneratePoint(Splitting &s) const
    {
      double k2(s.m_t0/s.m_Q2);
      s.m_z=pow((1.0+k2)/k2,-ran->Get())*(1.0+k2)-k2;
      s.m_phi=2.0*M_PI*ran->Get();
      s.m_z2=pow(s.m_z,ran->Get());
      double v(ran->Get());
      s.m_s=v/(1.0-v)*(s.m_t*s.m_z2/s.m_z+s.m_mj2);
      s.m_phi2=2.0*M_PI*ran->Get();
      s.m_mode=ran->Get()>0.5;
      return true;
    }

  };// end of class FFFF_FI

}// end of namespace DIRE

using namespace DIRE;

DECLARE_GETTER(FFFF_FI,"FI_FFFF",Lorentz,Kernel_Key);

Lorentz *ATOOLS::Getter<Lorentz,Kernel_Key,FFFF_FI>::
operator()(const Parameter_Type &args) const
{
  if (args.m_type!=2 || args.p_v) return NULL;
  if (args.m_lfid=="FFFF") {
    return new FFFF_FI(args);
  }
  return NULL;
}

void ATOOLS::Getter<Lorentz,Kernel_Key,FFFF_FI>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"FFFF Lorentz Function";
}

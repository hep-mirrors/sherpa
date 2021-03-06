#include "DIM/Shower/Lorentz_II.H"

#include "MODEL/Main/Single_Vertex.H"
#include "DIM/Shower/Shower.H"
#include "ATOOLS/Math/Random.H"

using namespace ATOOLS;

namespace DIM {
  
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
      if (m_mode) B=2.0*s.m_z*(1.0-s.m_z)+(1.0-s.m_z)/s.m_z;
      else {
	A=2.0*(1.0-s.m_z)/(sqr(1.0-s.m_z)+s.m_t/s.m_Q2);
	B=-2.0+(1.0-s.m_z)/s.m_z;
	A*=1.0+p_sk->GF()->K(s)+p_sk->GF()->RenCT(s);
      }
      return A+B;
    }

    double AsymmetryFactor(const Splitting &s) const
    {
      double V;
      if (m_mode) V=2.0*s.m_z*(1.0-s.m_z)+(1.0-s.m_z)/s.m_z;
      else {
	double A1=2.0*(1.0-s.m_z)/(sqr(1.0-s.m_z)+s.m_t/s.m_Q2);
	double B1=-2.0+(1.0-s.m_z)/s.m_z;
	V=A1*(1.0+p_sk->GF()->K(s)+p_sk->GF()->RenCT(s))+B1;
      }
      double A1=2.0*(1.0-s.m_z)/(sqr(1.0-s.m_z)+s.m_t/s.m_Q2);
      double B11=-2.0+(1.0-s.m_z)/s.m_z;
      double B12=2.0*s.m_z*(1.0-s.m_z)+(1.0-s.m_z)/s.m_z;
      return V/(A1*(1.0+p_sk->GF()->K(s)+p_sk->GF()->RenCT(s))+(B11+B12));
    }

    double Integral(const Splitting &s) const
    {
      if (m_mode) return log(1.0/s.m_eta)*m_jmax;
      double k2=s.m_t0/s.m_Q2;
      double I=log((k2+sqr(1.0-s.m_eta))/(s.m_eta*k2));
      return I*(1.0+p_sk->GF()->KMax(s))*m_jmax;
    }

    double Estimate(const Splitting &s) const
    {
      if (m_mode) return 1.0/s.m_z*m_jmax;
      double E=2.0*(1.0-s.m_z)/(sqr(1.0-s.m_z)+s.m_t0/s.m_Q2)+1.0/s.m_z;
      return E*(1.0+p_sk->GF()->KMax(s))*m_jmax;
    }

    bool GeneratePoint(Splitting &s) const
    {
      double k2(s.m_t0/s.m_Q2);
      if (m_mode) {
	s.m_z=pow(s.m_eta,ran->Get());
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

}// end of namespace DIM

using namespace DIM;

DECLARE_GETTER(VVV_II,"II_VVV",Lorentz,Kernel_Key);

Lorentz *ATOOLS::Getter<Lorentz,Kernel_Key,VVV_II>::
operator()(const Parameter_Type &args) const
{
  if (args.m_type!=3) return NULL;
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

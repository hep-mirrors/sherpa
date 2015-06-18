#include "DIRE/Shower/Lorentz_FI.H"

#include "MODEL/Main/Single_Vertex.H"
#include "DIRE/Shower/Shower.H"
#include "ATOOLS/Math/Random.H"

using namespace ATOOLS;

namespace DIRE {
  
  class VVV_FI: public Lorentz_FI {
  private:

    int m_mode;

  public:

    inline VVV_FI(const Kernel_Key &key,const int mode):
      Lorentz_FI(key), m_mode(mode) {}

    double Value(const Splitting &s) const
    {
      double z(s.m_z);
      double A1=2.0*(1.0-z)/(sqr(1.0-z)+s.m_t/s.m_Q2);
      double B1=-2.0+z*(1.0-z);
      return A1*(1.0+p_sk->GF()->K(s))+B1;
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

  };// end of class VVV_FI

}// end of namespace DIRE

using namespace DIRE;

DECLARE_GETTER(VVV_FI,"FI_VVV",Lorentz,Kernel_Key);

Lorentz *ATOOLS::Getter<Lorentz,Kernel_Key,VVV_FI>::
operator()(const Parameter_Type &args) const
{
  if (args.m_type!=2) return NULL;
  if (args.p_v->in[0].IntSpin()==2 &&
      args.p_v->in[1].IntSpin()==2 &&
      args.p_v->in[2].IntSpin()==2) {
    return new VVV_FI(args,args.m_mode);
  }
  return NULL;
}

void ATOOLS::Getter<Lorentz,Kernel_Key,VVV_FI>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"VVV Lorentz Function";
}

#include "ALARIC++/Shower/Lorentz_IS.H"

#include "MODEL/Main/Single_Vertex.H"
#include "ALARIC++/Shower/Shower.H"
#include "ATOOLS/Math/Random.H"

using namespace ATOOLS;

namespace ALARIC {
  
  class VVV_IS: public Lorentz_IS {
  private:

    double m_jmax;

    int m_mode;

  public:

    inline VVV_IS(const Kernel_Key &key,const int mode):
      Lorentz_IS(key), m_jmax(1.0), m_mode(mode) {}

    double Value(const Splitting &s) const
    {
      double xi=s.m_x+s.m_y-s.m_x*s.m_y*(1.0+s.m_kap);
      double B=(1.0-xi)/s.m_x;
      if (m_mode) B+=2.0*s.m_x*(1.0-s.m_x);
      return B;
    }

    double Integral(const Splitting &s) const
    {
      double I=log(1.0/s.m_eta);
      return I*m_jmax;
    }

    double Estimate(const Splitting &s) const
    {
      double E=1.0/s.m_z;
      return E*m_jmax;
    }

    bool GeneratePoint(Splitting &s) const
    {
      s.m_z=pow(s.m_eta,ran->Get());
      s.m_phi=2.0*M_PI*ran->Get();
      return true;
    }

  };// end of class VVV_IS

}// end of namespace ALARIC

using namespace ALARIC;

DECLARE_GETTER(VVV_IS,"IS_VVV",Lorentz,Kernel_Key);

Lorentz *ATOOLS::Getter<Lorentz,Kernel_Key,VVV_IS>::
operator()(const Parameter_Type &args) const
{
  if (args.m_type!=1) return NULL;
  if (args.m_swap) return NULL;
  if (args.p_v->in[0].IntSpin()==2 &&
      args.p_v->in[1].IntSpin()==2 &&
      args.p_v->in[2].IntSpin()==2) {
    return new VVV_IS(args,args.m_mode);
  }
  return NULL;
}

void ATOOLS::Getter<Lorentz,Kernel_Key,VVV_IS>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"VVV Lorentz Function";
}

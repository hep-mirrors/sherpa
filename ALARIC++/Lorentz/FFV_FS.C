#include "ALARIC++/Shower/Lorentz_FS.H"

#include "MODEL/Main/Single_Vertex.H"
#include "ALARIC++/Shower/Kernel.H"
#include "ATOOLS/Math/Random.H"

using namespace ATOOLS;

namespace ALARIC {
  
  class FFV_FS: public Lorentz_FS_Split {
  private:

    int m_swap;

  public:

    inline FFV_FS(const Kernel_Key &key):
      Lorentz_FS_Split(key), m_swap(key.m_swap) {}

    double Value(const Splitting &s) const
    {
      double B(1.0-s.m_z);
      B*=1.0+p_sk->GF()->K(s);
      return (s.m_clu?1.0:(m_swap?1.0-s.m_z:s.m_z))*B;
    }

    double Integral(const Splitting &s) const
    {
      return 1.0;
    }

    double Estimate(const Splitting &s) const
    {
      return 1.0;
    }

    bool GeneratePoint(Splitting &s) const
    {
      s.m_z=ran->Get();
      s.m_phi=2.0*M_PI*ran->Get();
      return true;
    }

  };// end of class FFV_FS

  class VFF_FS: public Lorentz_FS_Split {
  private:

    int m_swap;

  public:

    inline VFF_FS(const Kernel_Key &key):
      Lorentz_FS_Split(key), m_swap(key.m_swap) {}

    double Value(const Splitting &s) const
    {
      if (s.m_t<2.*sqr(m_fl[1].Mass(true))) return 0.;
      double B(1.0-2.0*s.m_z*(1.0-s.m_z));
      return (s.m_clu?1.0:(m_swap?1.0-s.m_z:s.m_z))*B;
    }

    double Integral(const Splitting &s) const
    {
      return 1.0;
    }

    double Estimate(const Splitting &s) const
    {
      return 1.0;
    }

    bool GeneratePoint(Splitting &s) const
    {
      s.m_z=ran->Get();
      s.m_phi=2.0*M_PI*ran->Get();
      return true;
    }

  };// end of class VFF_FS

}// end of namespace ALARIC

using namespace ALARIC;

DECLARE_GETTER(FFV_FS,"FS_FFV",Lorentz,Kernel_Key);

Lorentz *ATOOLS::Getter<Lorentz,Kernel_Key,FFV_FS>::
operator()(const Parameter_Type &args) const
{
  if (args.m_type!=0) return NULL;
  if (args.p_v->in[0].IntSpin()==1 &&
      args.p_v->in[1+args.m_mode].IntSpin()==1 &&
      args.p_v->in[2-args.m_mode].IntSpin()==2) {
    return new FFV_FS(args);
  }
  if (args.m_mode) return NULL;
  if (args.p_v->in[0].IntSpin()==2 &&
      args.p_v->in[1].IntSpin()==1 &&
      args.p_v->in[2].IntSpin()==1) {
    return new VFF_FS(args);
  }
  return NULL;
}

void ATOOLS::Getter<Lorentz,Kernel_Key,FFV_FS>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"FFV Lorentz Function";
}

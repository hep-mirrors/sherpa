#include "ALARIC++/Shower/Lorentz_FS.H"

#include "MODEL/Main/Single_Vertex.H"
#include "ALARIC++/Shower/Shower.H"
#include "ATOOLS/Math/Random.H"

using namespace ATOOLS;

namespace ALARIC {
  
  class VVV_FS: public Lorentz_FS_Split {
  public:

    inline VVV_FS(const Kernel_Key &key):
      Lorentz_FS_Split(key) {}

    double Value(const Splitting &s) const
    {
      double B(s.m_z*(1.0-s.m_z));
      double sf(1.0);
      if (!s.m_clu) sf=p_sk->Mode()?1.0-s.m_z:s.m_z;
      return sf*B*(1.0+p_sk->GF()->K(s));
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

  };// end of class VVV_FS

}// end of namespace ALARIC

using namespace ALARIC;

DECLARE_GETTER(VVV_FS,"FS_VVV",Lorentz,Kernel_Key);

Lorentz *ATOOLS::Getter<Lorentz,Kernel_Key,VVV_FS>::
operator()(const Parameter_Type &args) const
{
  if (args.m_type!=0) return NULL;
  if (args.m_swap) return NULL;
  if (args.p_v->in[0].IntSpin()==2 &&
      args.p_v->in[1].IntSpin()==2 &&
      args.p_v->in[2].IntSpin()==2) {
    return new VVV_FS(args);
  }
  return NULL;
}

void ATOOLS::Getter<Lorentz,Kernel_Key,VVV_FS>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"VVV Lorentz Function";
}

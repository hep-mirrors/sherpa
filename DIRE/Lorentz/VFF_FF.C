#include "DIRE/Shower/Lorentz_FF.H"

#include "MODEL/Interaction_Models/Single_Vertex.H"
#include "ATOOLS/Math/Random.H"

using namespace ATOOLS;

namespace DIRE {
  
  class FFV_FF: public Lorentz_FF {
  public:

    inline FFV_FF(const Kernel_Key &key):
      Lorentz_FF(key) {}

    double Value(const Splitting &s) const
    {
      double z(s.m_z), y(s.m_y);
      if (s.m_mij2==0.0 && s.m_mi2==0.0 && s.m_mk2==0.0) {
	double V=2.0/(1.0-z*(1.0-y))-(1.0+z);
	return V;
      }
      double muij2(s.m_mij2/s.m_Q2), mui2(s.m_mi2/s.m_Q2);
      double muk2(s.m_mk2/s.m_Q2), vtijk=Lam(1.0,muij2,muk2);
      double vijk=sqr(2.0*muk2+(1.0-mui2-muk2)*(1.0-y))-4.0*muk2;
      if (vtijk<0.0 || vijk<0.0) return 0.0;
      vtijk=sqrt(vtijk)/(1.0-muij2-muk2);
      vijk=sqrt(vijk)/((1.0-mui2-muk2)*(1.0-y));
      double pipj=s.m_Q2*(1.0-mui2-muk2)*s.m_y/2.0;
      double V=2.0/(1.0-z*(1.0-y))-vtijk/vijk*(1.0+z+s.m_mi2/pipj);
      V/=(1.0-mui2-muk2)+1.0/y*(mui2-muij2);
      return V;
    }

    double Integral(const Splitting &s) const
    {
      return 2.0*log((1.0-s.m_zmin)/(1.0-s.m_zmax));
    }

    double Estimate(const Splitting &s) const
    {
      return 2.0/(1.0-s.m_z);
    }

    bool GeneratePoint(Splitting &s) const
    {
      s.m_z=1.0-(1.0-s.m_zmin)*
	pow((1.0-s.m_zmax)/(1.0-s.m_zmin),ATOOLS::ran->Get());
      s.m_phi=2.0*M_PI*ran->Get();
      return true;
    }

  };// end of class FFV_FF

}// end of namespace DIRE

using namespace DIRE;

DECLARE_GETTER(FFV_FF,"Gamma",Lorentz,Kernel_Key);

Lorentz *ATOOLS::Getter<Lorentz,Kernel_Key,VFF_FF>::
operator()(const Parameter_Type &args) const
{
  if (args.p_v->in[0].IntSpin()==2 &&
      args.p_v->in[1].IntSpin()==1 &&
      args.p_v->in[2].IntSpin()==1) {
    switch (args.m_type) {
    case 0: return new VFF_FF(args);
    }
  }
  return NULL;
}

void ATOOLS::Getter<Lorentz,Kernel_Key,VFF_FF>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"VFF Lorentz Function";
}

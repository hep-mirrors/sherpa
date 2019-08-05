#include "CFPSHOWER++/Calculators/FF/SF_FF12.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class VFF_FF : public SF_FF12 {
  private:
    double B1(const double & z,const double & kappa2) const;
  public:
    VFF_FF(const Kernel_Info & info);
    double operator()(const Splitting & split) const;
    double Integral(const Splitting & split) const;
    double OverEstimate(const Splitting & split) const;
    void   GeneratePoint(Splitting & split) const;
  };
}


using namespace CFPSHOWER;
using namespace ATOOLS;

VFF_FF::VFF_FF(const Kernel_Info & info) : SF_FF12(info) {
  SetName("V->FF");
}

double VFF_FF::operator()(const Splitting & split) const {
  double mi2(split.m2(1)), mj2(split.m2(2)), mspect2(split.mspect2());
  double z(split.z(0)), y(split.y()), Q2(split.Q2()), kappa2(split.t()/Q2);
  double value = 0.;
  if (mi2==0. && mj2==0. && mspect2==0.) value += B1(z,kappa2);
  else {
    double mui2(mi2/Q2);
    double vijl  = sqr(1.-y)-4.*(y+2.*mui2)*(mspect2/Q2);
    if (vijl<0.) return 0.;
    value       += (1.-y)/sqrt(vijl) * (B1(z,kappa2)+2.*mui2/(y+2.*mui2));
  }
  if (split.Clustered()==0) value *= split.z(m_tagsequence[0]);
  return value;
}

double VFF_FF::Integral(const Splitting & split)     const { return 1.; }

double VFF_FF::OverEstimate(const Splitting & split) const { return 1.; }

void VFF_FF::GeneratePoint(Splitting & split) const {
  double z = ran->Get();
  split.Set_z(0,z);
  split.Set_z(1,1.-z);
  split.Set_phi();
}

double VFF_FF::B1(const double & z,const double & kappa2) const {
  return sqr(z)+sqr(1.-z);
}

DECLARE_GETTER(VFF_FF,"FF_VFF",SF_Base,Kernel_Info);

SF_Base *
ATOOLS::Getter<SF_Base,Kernel_Info,VFF_FF>::
operator()(const Parameter_Type & info) const
{
  if (info.Type()==kernel_type::FF &&
      (info.GetSplit().IsVector() &&
       info.GetFlavs()[0].IsFermion() &&
       info.GetFlavs()[0].IsAnti() && 
       info.GetFlavs()[1].IsFermion())) {
    return new VFF_FF(info);
  }
  return NULL;
}

void ATOOLS::Getter<SF_Base,Kernel_Info,VFF_FF>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"VFF Splitting Function (FF)";
}

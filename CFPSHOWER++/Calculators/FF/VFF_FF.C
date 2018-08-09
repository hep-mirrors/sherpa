#include "CFPSHOWER++/Calculators/FF/SF_FF.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

//#include "CFPSHOWER++/Calculators/GQQ.C"

namespace CFPSHOWER {
  class VFF_FF : public SF_FF {
  private:
    double B1(const double & z,const double & kappa2) const;
    double B2(const double & z,const double & kappa2) const;
  public:
    VFF_FF(const Kernel_Info & info);
    double operator()(const Splitting & split) const;
    double Integral(const Splitting & split) const;
    double OverEstimate(const Splitting & split) const;
  };
}


using namespace CFPSHOWER;
using namespace ATOOLS;

VFF_FF::VFF_FF(const Kernel_Info & info) :
  SF_FF(info) {
  SetName("V->FF");
}

double VFF_FF::operator()(const Splitting & split) const {
  double mij2(split.mij2()), mi2(split.mi2()), mj2(split.mj2()), mk2(split.mk2());
  double z(split.Z()), y(split.Y()), Q2(split.Q2()), kappa2(split.T()/Q2);
  double value = 0.;
  if (mi2==0. && mj2==0. && mk2==0.) {
    value += B1(z,kappa2);
  }
  else {
    double mui2(mi2/Q2), muk2(mk2/Q2);
    double vijk  = sqr(1.-y)-4.*(y+2.*mui2)*muk2;
    if (vijk<0.) return 0.;
    vijk         = sqrt(vijk)/(1.-y);
    value       += 1./vijk*(B1(z,kappa2)+2.*mui2/(y+2.*mui2));
  }
  if (split.Clustered()==0) value *= m_swap?(1.-z):z;
  return value;
}

double VFF_FF::Integral(const Splitting & split) const { return 1.; }

double VFF_FF::OverEstimate(const Splitting & split) const { return 1.; }

double VFF_FF::B1(const double & z,const double & kappa2) const {
  return sqr(z)+sqr(1.-z);
}

double VFF_FF::B2(const double & z,const double & kappa2) const {
  double b2 = 0.0;
  return b2;
}

DECLARE_GETTER(VFF_FF,"FF_VFF",SF_Base,Kernel_Info);

SF_Base *
ATOOLS::Getter<SF_Base,Kernel_Info,VFF_FF>::
operator()(const Parameter_Type & info) const
{
  return NULL;
  if (info.Type()==kernel_type::FF &&
      (info.GetFlavs()[0].IsVector() &&
       info.GetFlavs()[1].IsFermion() && info.GetFlavs()[1].IsAnti() && 
       info.GetFlavs()[2].IsFermion())) {
    return new VFF_FF(info);
  }
  return NULL;
}

void ATOOLS::Getter<SF_Base,Kernel_Info,VFF_FF>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"VFF Splitting Function (FF)";
}

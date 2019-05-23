#include "CFPSHOWER++/Calculators/FF/SF_FF.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

//#include "CFPSHOWER++/Calculators/GGG.C"

namespace CFPSHOWER {
  class VVV_FF : public SF_FF {
  private:
    double A1(const double & z,const double & kappa2) const;
    double B1(const double & z,const double & kappa2) const;
    double B2(const double & z,const double & kappa2) const;
  public:
    VVV_FF(const Kernel_Info & info);
    double operator()(const Splitting & split) const;
    double Integral(const Splitting & split) const;
    double OverEstimate(const Splitting & split) const;
    void   GeneratePoint(Splitting & split) const;
  };
}


using namespace CFPSHOWER;
using namespace ATOOLS;

VVV_FF::VVV_FF(const Kernel_Info & info) :
  SF_FF(info) {
  SetName("V->VV");
}

double VVV_FF::operator()(const Splitting & split) const {
  double mij2(split.mij2()), mi2(split.mi2()), mj2(split.mj2()), mk2(split.mk2());
  double z(split.Z()), y(split.Y()), Q2(split.Q2()), kappa2(split.T()/Q2);
  // Start with the soft term only, including possible K factors
  // (cusp anomalous dimensions), obtained from the gauge part of the kernel
  double hofactor = (1.+split.GetKernel()->GetGauge()->K(split));
  double value    = A1(z,kappa2) * hofactor;
  // All massless: just add the collinear parts.
  // TODO: Add the B2 parts as soon as the LL shower is validated.
  if (mi2==0. && mj2==0. && mk2==0.) {
    if (m_orderB>0) value += B1(z,kappa2);
  }
  else {
    double muk2(split.mk2()/Q2);
    double vijk  = sqr(1.-y)-4.*y*muk2;
    if (vijk<0.) return 0.;
    vijk         = sqrt(vijk)/(1.-y);
    value       += B1(z,kappa2)/vijk;
    value       -= 2.*(muk2*y)/((1.-z)*(1.-z+y));
  }
  if (split.Clustered()==0) value *= m_swap?(1.-z):z;
  return value;
}

double VVV_FF::Integral(const Splitting & split) const {
  double homax = (1.+split.GetKernel()->GetGauge()->KMax(split));
  return log(1.0+split.Q2()/split.T0()) * homax;
}

double VVV_FF::OverEstimate(const Splitting & split) const {
  double homax = (1.+split.GetKernel()->GetGauge()->KMax(split));
  return A1(split.Z(),split.T0()/split.Q2()) * homax;
}

void VVV_FF::GeneratePoint(Splitting & split) const {
  double kappa2 = split.T0()/split.Q2();
  double z      = 1.-sqrt(kappa2 * (pow((1.+1./kappa2),ran->Get())-1.)); 
  split.SetZ(z);
  split.Setphi();
}

double VVV_FF::A1(const double & z,const double & kappa2) const {
  return 2.*(1.-z)/(sqr(1.-z)+kappa2);
}

double VVV_FF::B1(const double & z,const double & kappa2) const {
  return -2.+z*(1.-z);
}

double VVV_FF::B2(const double & z,const double & kappa2) const {
  double b2 = 0.0;
  return b2;
}

DECLARE_GETTER(VVV_FF,"FF_VVV",SF_Base,Kernel_Info);

SF_Base * ATOOLS::Getter<SF_Base,Kernel_Info,VVV_FF>::
operator()(const Parameter_Type & info) const
{
  return NULL;
  if (info.Type()==kernel_type::FF &&
      (info.GetFlavs()[0].IsVector() &&
       info.GetFlavs()[1].IsVector() &&
       info.GetFlavs()[2].IsVector())) return new VVV_FF(info);
  return NULL;
}

void ATOOLS::Getter<SF_Base,Kernel_Info,VVV_FF>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"VVV Splitting Function (FF)";
}

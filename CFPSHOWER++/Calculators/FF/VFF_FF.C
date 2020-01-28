#include "CFPSHOWER++/Calculators/FF/SF_FF12.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class VFF_FF : public SF_FF12 {
  private:
    double B1(const double & z,const double & kappa2) const;
  public:
    VFF_FF(const Kernel_Info & info);
    double operator()(const Splitting & split);
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

double VFF_FF::operator()(const Splitting & split) {
  double z(split.z()), kappa2(split.t()/split.Q2red());
  double value = 0.;
  if (split.IsMassive()) {
    double mi2(split.m2(0)), mk2(split.mspect2());
    double Q2      = split.Q2(), sij = split.y()*(Q2-mk2);
    double v2_ij_k = Lambda2(Q2,sij,mk2);
    if (v2_ij_k<0.) return 0.;
    value += sqrt(v2_ij_k)/(Q2-sij-mk2) * (B1(z,kappa2) + 2.*mi2/sij);
  }
  else value += B1(z,kappa2);
  if (split.Clustered()==0) value *= (m_tags[0]==0) ? z : 1-z;
  return value;
}

double VFF_FF::Integral(const Splitting & split)     const { return 1.; }

double VFF_FF::OverEstimate(const Splitting & split) const { return 1.; }

void VFF_FF::GeneratePoint(Splitting & split) const {
  double z = ran->Get();
  split.Set_z(z);
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

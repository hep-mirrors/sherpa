#include "CFPSHOWER++/Calculators/SF_Base.H"
#include "CFPSHOWER++/Calculators/SF_Base.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class VVV_FF : public SF_Base {
  private:
    double A1(const double & z,const double & kappa2) const;
    double B1(const double & z,const double & kappa2) const;
  public:
    VVV_FF(const Kernel_Info & info);
    double operator()(const Splitting & split);
    double Integral(const Splitting & split) const;
    double OverEstimate(const Splitting & split) const;
    void   GeneratePoint(Splitting & split) const;
  };
}


using namespace CFPSHOWER;
using namespace ATOOLS;

VVV_FF::VVV_FF(const Kernel_Info & info) : SF_Base(info) {
  SetName("V->VV");
}

double VVV_FF::operator()(const Splitting & split) {
  double mspect2(split.mspect2());
  double z(split.z(0)), kappa2(split.t(0)/split.Q2red());
  // Start with the soft term only, including possible K factors (cusp anomalous
  // dimensions), obtained from the gauge part of the kernel
  double Kfactor = m_CMW==1 ? (1.+split.GetKernel()->GetGauge()->K(split)) : 1.;
  double value   = A1(z,kappa2) * Kfactor;
  if (split.IsMassive()) {
    // Massive spectator adjustments
    // directly return 0 if the splitting is kinematically not viable
    double mk2     = split.mspect2();
    double Q2      = split.Q2(), y = split.y(), sij = y*(Q2-mk2);
    double v2_ij_k = Lambda2(Q2,sij,mk2);
    if (v2_ij_k<0.) return 0.;
    double v_ij_k  = sqrt(v2_ij_k)/(Q2-sij-mk2);
    // Terms ~(1.-z)/(1-z+y) come from reshuffling the mass-dependent terms in the eikonal
    // to ensure correct soft behaviour
    value += B1(z,kappa2)/v_ij_k - 2.*(mk2*y)/((Q2-mk2)*(1.-z)*(1.-z+y));
  }
  else value += B1(z,kappa2);
  if (split.Clustered()==0) value *= split.ztilde(m_tags[0]);
  return value;
}

double VVV_FF::Integral(const Splitting & split) const {
  double Kmax = (m_CMW==1.) ? (1.+split.GetKernel()->GetGauge()->KMax(split)) : 1.;
  return log(1.0+split.Q2()/split.tcut()) * Kmax;
}

double VVV_FF::OverEstimate(const Splitting & split) const {
  double Kmax = (m_CMW==1.) ? (1.+split.GetKernel()->GetGauge()->KMax(split)) : 1.;
  return A1(split.z(0),split.tcut()/split.Q2()) * Kmax;
}

void VVV_FF::GeneratePoint(Splitting & split) const {
  double kappa2 = split.tcut()/split.Q2();
  double z     = 1.-sqrt(kappa2 * (pow((1.+1./kappa2),ran->Get())-1.)); 
  split.Set_z(0,z);
  split.Set_phi(0);
}

double VVV_FF::A1(const double & z,const double & kappa2) const {
  return 2.*(1.-z)/(sqr(1.-z)+kappa2);
}

double VVV_FF::B1(const double & z,const double & kappa2) const {
  return -2.+z*(1.-z);
}

DECLARE_GETTER(VVV_FF,"FF_VVV",SF_Base,Kernel_Info);

SF_Base * ATOOLS::Getter<SF_Base,Kernel_Info,VVV_FF>::
operator()(const Parameter_Type & info) const
{
  return NULL;
  if (info.Type()==kernel_type::FF &&
      (info.GetSplit().IsVector() &&
       info.GetFlavs()[0].IsVector() &&
       info.GetFlavs()[1].IsVector())) {
    return new VVV_FF(info);
  }
  return NULL;
}

void ATOOLS::Getter<SF_Base,Kernel_Info,VVV_FF>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"VVV Splitting Function (FF)";
}

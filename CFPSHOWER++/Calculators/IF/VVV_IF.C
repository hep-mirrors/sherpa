#include "CFPSHOWER++/Calculators/IF/SF_IF12.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class VVV_IF : public SF_IF12 {
  private:    
    double A1(const double & z,const double & kappa2) const;
    double B1(const double & z,const double & kappa2) const;
  public:
    VVV_IF(const Kernel_Info & info);
    double operator()(const Splitting & split);
    double Integral(const Splitting & split) const;
    double OverEstimate(const Splitting & split) const;
    void   GeneratePoint(Splitting & split) const;
  };
}


using namespace CFPSHOWER;
using namespace ATOOLS;

VVV_IF::VVV_IF(const Kernel_Info & info) : SF_IF12(info) {
  SetName("V->VV");
}

double VVV_IF::operator()(const Splitting & split) {
  double z(split.z()), kappa2(split.t()/split.Q2red());
  // Start with the soft term only, including possible K factors
  // (cusp anomalous dimensions), obtained from the gauge part of the kernel
  double Kfactor = m_CMW==1 ? (1.+split.GetKernel()->GetGauge()->K(split)) : 1.;
  double value   = A1(z,kappa2) * Kfactor;
  // All massless: just add the collinear parts.
  value += B1(z,kappa2);
  if (split.IsMassive()) {
    double mk2 = split.mspect2(), y = split.y(), pipj = split.Q2red()*(1.-y)/(2.*y);
    value -= mk2/pipj;
  }
  return value;
}

double VVV_IF::Integral(const Splitting & split) const {
  if (m_tags[0]==1) return log(1./split.eta());
  double Kmax = (m_CMW==1.) ? (1.+split.GetKernel()->GetGauge()->KMax(split)) : 1.;
  double eta  = split.eta(), kappa2 = split.tcut()/split.Q2red();
  double I    = log((kappa2+sqr(1.-eta))/(kappa2*eta));
  return I * Kmax;
}

double VVV_IF::OverEstimate(const Splitting & split) const {
  if (m_tags[0]==1) return 1./split.z();
  double Kmax = (m_CMW==1.) ? (1.+split.GetKernel()->GetGauge()->KMax(split)) : 1.;
  return (A1(split.z(),split.tcut()/split.Q2red()) * Kmax + 1./split.z());
}

void VVV_IF::GeneratePoint(Splitting & split) const {
  if (m_tags[0]==1) split.Set_z(pow(split.eta(),ran->Get()));
  else {
    double eta  = split.eta(), kappa2 = split.tcut()/split.Q2();
    double arg  = (kappa2+sqr(1.-eta))/(eta*kappa2);
    double help = 1.+kappa2/2.*pow(arg,ran->Get());
    split.Set_z(help-sqrt(sqr(help)-(1.+kappa2)));
  }
  split.Set_phi();
}

double VVV_IF::A1(const double & z,const double & kappa2) const {
  return m_tags[0]==1 ? 0.: 2.*(1.-z)/(sqr(1.-z)+kappa2);
}

double VVV_IF::B1(const double & z,const double & kappa2) const {
  return m_tags[0]==1 ? 2.*z*(1.-z)+(1.-z)/z : -2.+(1.-z)/z;
}

DECLARE_GETTER(VVV_IF,"IF_VVV",SF_Base,Kernel_Info);

SF_Base * ATOOLS::Getter<SF_Base,Kernel_Info,VVV_IF>::
operator()(const Parameter_Type & info) const
{
  return NULL;
  if (info.Type()==kernel_type::IF &&
      (info.GetSplit().IsVector() &&
       info.GetFlavs()[0].IsVector() &&
       info.GetFlavs()[1].IsVector())) return new VVV_IF(info);
  return NULL;
}

void ATOOLS::Getter<SF_Base,Kernel_Info,VVV_IF>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"VVV Splitting Function (IF)";
}

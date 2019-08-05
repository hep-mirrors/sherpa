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
    double operator()(const Splitting & split) const;
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

double VVV_IF::operator()(const Splitting & split) const {
  double mspect2(split.mspect2());
  double z(split.z(0)), y(split.y()), Q2(split.Q2()), kappa2(split.t()/Q2);
  // Start with the soft term only, including possible K factors
  // (cusp anomalous dimensions), obtained from the gauge part of the kernel
  double Kfactor = m_CMW==1 ? (1.+split.GetKernel()->GetGauge()->K(split)) : 1.;
  double value   = A1(z,kappa2) * Kfactor;
  // All massless: just add the collinear parts.
  value += B1(z,kappa2);
  if (mspect2>0.) value -= mspect2*y/(Q2*(1.-y)); 
  return value;
}

double VVV_IF::Integral(const Splitting & split) const {
  if (m_tagsequence[0]==2) return log(1./split.eta());
  double Kmax = (m_CMW==1.) ? (1.+split.GetKernel()->GetGauge()->KMax(split)) : 1.;
  double kappa2 = split.t0()/split.Q2();
  double I      = log((kappa2+sqr(1.-split.eta()))/(kappa2*split.eta()));
  return I * Kmax;
}

double VVV_IF::OverEstimate(const Splitting & split) const {
  if (m_tagsequence[0]==2) return 1./split.z(0);
  double Kmax = (m_CMW==1.) ? (1.+split.GetKernel()->GetGauge()->KMax(split)) : 1.;
  return (A1(split.z(0),split.t0()/split.Q2()) * Kmax + 1./split.z(0));
}

void VVV_IF::GeneratePoint(Splitting & split) const {
  if (m_tagsequence[0]==1) split.Set_z(0,pow(split.eta(),ran->Get()));
  else {
    double kappa2 = split.t0()/split.Q2();
    double arg    = (kappa2+sqr(1.-split.eta()))/(split.eta()*kappa2);
    double help   = 1.+kappa2/2.*pow(arg,ran->Get());
    split.Set_z(0,help-sqrt(sqr(help)-(1.+kappa2)));
  }
  split.Set_z(1,1.-split.z(0));
  split.Set_phi();
}

double VVV_IF::A1(const double & z,const double & kappa2) const {
  return m_tagsequence[0]==2 ? 0.: 2.*(1.-z)/(sqr(1.-z)+kappa2);
}

double VVV_IF::B1(const double & z,const double & kappa2) const {
  return m_tagsequence[0]==2 ? 2.*z*(1.-z)+(1.-z)/z : -2.+(1.-z)/z;
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

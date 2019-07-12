#include "CFPSHOWER++/Calculators/IF/SF_IF.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class VVV_IF : public SF_IF {
  private:
    double m_jmax;
    
    double A1(const double & z,const double & kappa2) const;
    double B1(const double & z,const double & kappa2) const;
    double B2(const double & z,const double & kappa2) const;
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

VVV_IF::VVV_IF(const Kernel_Info & info) :
  SF_IF(info), m_jmax(1.) {
  SetName("V->VV");
}

double VVV_IF::operator()(const Splitting & split) const {
  double mij2(split.mij2()), mi2(split.mi2()), mj2(split.mj2()), mk2(split.mk2());
  double z(split.Z()), y(split.Y()), Q2(split.Q2()), kappa2(split.T()/Q2);
  // Start with the soft term only, including possible K factors
  // (cusp anomalous dimensions), obtained from the gauge part of the kernel
  double hofactor = (1.+split.GetKernel()->GetGauge()->K(split));
  double value    = A1(z,kappa2) * hofactor;
  // All massless: just add the collinear parts.
  // TODO: Add the B2 parts as soon as the LL shower is validated.
  if (m_orderB>0) {
    value += B1(z,kappa2);
    if (mk2>0.) value -= mk2*y/(Q2*(1.-y));
  }
  //msg_Out()<<"   * value = "<<(A1(z,kappa2)*hofactor)<<"+"<<B1(z,kappa2)<<" = "<<value<<"\n";
  return value;
}

double VVV_IF::Integral(const Splitting & split) const {
  if (m_swap) return log(1./split.Eta()) * m_jmax;
  double homax  = (1.+split.GetKernel()->GetGauge()->KMax(split));
  double kappa2 = split.T0()/split.Q2();
  double I      = log((kappa2+sqr(1.-split.Eta()))/(kappa2*split.Eta()));
  return I * homax * m_jmax;
}

double VVV_IF::OverEstimate(const Splitting & split) const {
  if (m_swap) return 1./split.Z() * m_jmax;
  double homax = (1.+split.GetKernel()->GetGauge()->KMax(split));
  return (A1(split.Z(),split.T0()/split.Q2()) * homax + 1./split.Z()) * m_jmax;
}

void VVV_IF::GeneratePoint(Splitting & split) const {
  if (m_swap) split.SetZ(pow(split.Eta(),ran->Get()));
  else {
    double kappa2 = split.T0()/split.Q2();
    double arg    = (kappa2+sqr(1.-split.Eta()))/(split.Eta()*kappa2);
    double help   = 1.+kappa2/2.*pow(arg,ran->Get());
    //std::cout<<"   * generate z with eta = "<<split.Eta()<<" -> FF = "<<help<<"\n";
    split.SetZ(help-sqrt(sqr(help)-(1.+kappa2)));
  }
  split.Setphi();
}

double VVV_IF::A1(const double & z,const double & kappa2) const {
  return m_swap ? 0.: 2.*(1.-z)/(sqr(1.-z)+kappa2);
}

double VVV_IF::B1(const double & z,const double & kappa2) const {
  return m_swap ? 2.*z*(1.-z)+(1.-z)/z : -2.+(1.-z)/z;
}

double VVV_IF::B2(const double & z,const double & kappa2) const {
  double b2 = 0.0;
  return b2;
}

DECLARE_GETTER(VVV_IF,"IF_VVV",SF_Base,Kernel_Info);

SF_Base * ATOOLS::Getter<SF_Base,Kernel_Info,VVV_IF>::
operator()(const Parameter_Type & info) const
{
  return NULL;
  //if (info.Swapped()) return NULL;
  if (info.Type()==kernel_type::IF &&
      (info.GetFlavs()[0].IsVector() &&
       info.GetFlavs()[1].IsVector() &&
       info.GetFlavs()[2].IsVector())) return new VVV_IF(info);
  return NULL;
}

void ATOOLS::Getter<SF_Base,Kernel_Info,VVV_IF>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"VVV Splitting Function (IF)";
}

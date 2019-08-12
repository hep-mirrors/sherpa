#include "CFPSHOWER++/Calculators/FI/SF_FI12.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class VVV_FI : public SF_FI12 {
  private:
    double A1(const double & z,const double & kappa2) const;
    double B1(const double & z,const double & kappa2) const;
  public:
    VVV_FI(const Kernel_Info & info);
    double operator()(const Splitting & split) const;
    double Integral(const Splitting & split) const;
    double OverEstimate(const Splitting & split) const;
    void   GeneratePoint(Splitting & split) const;
  };
}


using namespace CFPSHOWER;
using namespace ATOOLS;

VVV_FI::VVV_FI(const Kernel_Info & info) : SF_FI12(info) {
  SetName("V->VV");
}

double VVV_FI::operator()(const Splitting & split) const {
  double z(split.z()), kappa2(split.tcut()/split.Q2red());
  // Start with the soft term only, including possible K factors
  // (cusp anomalous dimensions), obtained from the gauge part of the kernel
  double Kfactor = m_CMW==1 ? (1.+split.GetKernel()->GetGauge()->K(split)) : 1.;
  double value    = A1(z,kappa2) * Kfactor;
  // All massless: just add the collinear parts.
  if (split.IsMassive()) {
    msg_Error()<<"Error in "<<METHOD<<": did not expect massive spectator in IS.\n"
	       <<"   Will exit the run.\n";
    exit(1);
  }
  else {
    value += B1(z,kappa2);
  }
  if (split.Clustered()==0) value *= (m_tags[0]==0) ? z : 1-z;
  return value;
}

double VVV_FI::Integral(const Splitting & split) const {
  double Kmax = (m_CMW==1.) ? (1.+split.GetKernel()->GetGauge()->KMax(split)) : 1.;
  return log(1.0+split.Q2red()/split.tcut()) * Kmax;
}

double VVV_FI::OverEstimate(const Splitting & split) const {
  double Kmax = (m_CMW==1.) ? (1.+split.GetKernel()->GetGauge()->KMax(split)) : 1.;
  return A1(split.z(),split.tcut()/split.Q2red()) * Kmax;
}

void VVV_FI::GeneratePoint(Splitting & split) const {
  double kappa2 = split.tcut()/split.Q2red();
  split.Set_z(1.-sqrt(kappa2 * (pow((1.+1./kappa2),ran->Get())-1.)));
  split.Set_phi();
}

double VVV_FI::A1(const double & z,const double & kappa2) const {
  return 2.*(1.-z)/(sqr(1.-z)+kappa2);
}

double VVV_FI::B1(const double & z,const double & kappa2) const {
  return -2.+z*(1.-z);
}

DECLARE_GETTER(VVV_FI,"FI_VVV",SF_Base,Kernel_Info);

SF_Base * ATOOLS::Getter<SF_Base,Kernel_Info,VVV_FI>::
operator()(const Parameter_Type & info) const
{
  return NULL;
  if (info.Type()==kernel_type::FI &&
      (info.GetSplit().IsVector() &&
       info.GetFlavs()[0].IsVector() &&
       info.GetFlavs()[1].IsVector())) return new VVV_FI(info);
  return NULL;
}

void ATOOLS::Getter<SF_Base,Kernel_Info,VVV_FI>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"VVV Splitting Function (FI)";
}

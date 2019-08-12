#include "CFPSHOWER++/Calculators/FI/SF_FI12.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class FFV_FI : public SF_FI12 {
  private:
    double m_jmax;
    
    double A1(const double & z,const double & kappa2) const;
    double B1(const double & z,const double & kappa2) const;
  public:
    FFV_FI(const Kernel_Info & info);
    double operator()(const Splitting & split) const;
    double Integral(const Splitting & split) const;
    double OverEstimate(const Splitting & split) const;
    void   GeneratePoint(Splitting & split) const;
  };
}

using namespace CFPSHOWER;
using namespace ATOOLS;

FFV_FI::FFV_FI(const Kernel_Info & info) : SF_FI12(info), m_jmax(5.) {
  SetName(m_tags[0]==1?"F->VF":"F->FV");
}

double FFV_FI::operator()(const Splitting & split) const {
  double z(split.z()), kappa2(split.t()/split.Q2red());
  // Start with the soft term only, including possible K factors
  // (cusp anomalous dimensions), obtained from the gauge part of the kernel
  double Kfactor = m_CMW==1 ? (1.+split.GetKernel()->GetGauge()->K(split)) : 1.;
  double value   = A1(z,kappa2) * Kfactor;
  // All massless: just add the collinear parts.
  // TODO: Add the B2 parts as soon as the LL shower is validated.
  if (split.mspect2()>1.e-12) {
    msg_Error()<<"Error in "<<METHOD<<": did not expect massive spectator in IS.\n"
	       <<"   Will exit the run.\n";
    exit(1);
  }
  value += B1(z,kappa2);
  if (split.IsMassive()) {
    double mi2 = split.m2(0), y = split.y(), pipj = split.Q2red()*(1.-y)/(2.*y);
    value -= mi2/pipj;
  }
  if (split.Clustered()==0) value *= (m_tags[0]==0) ? z : 1-z;
  return value;
}

double FFV_FI::Integral(const Splitting & split) const {
  double Kmax = (m_CMW==1.) ? (1.+split.GetKernel()->GetGauge()->KMax(split)) : 1.;
  return log(1.0+split.Q2red()/split.tcut()) * Kmax * m_jmax;
}

double FFV_FI::OverEstimate(const Splitting & split) const {
  double Kmax = (m_CMW==1.) ? (1.+split.GetKernel()->GetGauge()->KMax(split)) : 1.;
  return A1(split.z(),split.tcut()/split.Q2red()) * Kmax * m_jmax;
}

void FFV_FI::GeneratePoint(Splitting & split) const {
  double kappa2 = split.tcut()/split.Q2red();
  split.Set_z(1.-sqrt(kappa2 * (pow((1.+1./kappa2),ran->Get())-1.)));
  split.Set_phi();
}

double FFV_FI::A1(const double & z,const double & kappa2) const {
  return 2.*(1.-z)/(sqr(1.-z)+kappa2);
}

double FFV_FI::B1(const double & z,const double & kappa2) const {
  return -(1.+z);
}

DECLARE_GETTER(FFV_FI,"FI_FFV",SF_Base,Kernel_Info);

SF_Base * ATOOLS::Getter<SF_Base,Kernel_Info,FFV_FI>::
operator()(const Parameter_Type & info) const
{
  return NULL;
  if (info.Type()==kernel_type::FI &&
      info.GetSplit().IsFermion() && 
      info.GetFlavs()[0].IsFermion() &&
      info.GetFlavs()[1].IsVector()) {
    return new FFV_FI(info);
  }
  return NULL;
}

void ATOOLS::Getter<SF_Base,Kernel_Info,FFV_FI>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"FFV Splitting Function (FI)";
}

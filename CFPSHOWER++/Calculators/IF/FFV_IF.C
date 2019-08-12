#include "CFPSHOWER++/Calculators/IF/SF_IF12.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class FFV_IF : public SF_IF12 {
    // somewhat tricky nomenclature to keep track off.
    // incoming fermion F (fl[0], mass mij), splitting into an incoming fermion 
    // F (fl[1], mass mi) and an outgoing vector boson V (fl[2], mass mj);
    // spectator has mass mk.
  private:
    double m_jmax;
    
    double A1(const double & z,const double & kappa2) const;
    double B1(const double & z,const double & kappa2) const;
  public:
    FFV_IF(const Kernel_Info & info);
    double operator()(const Splitting & split) const;
    double Integral(const Splitting & split) const;
    double OverEstimate(const Splitting & split) const;
    void   GeneratePoint(Splitting & split) const;
  };
}


using namespace CFPSHOWER;
using namespace ATOOLS;


// over-estimator: 5 for potential valence quarks, 2 for all others.
// TODO: may need to monitor this.
FFV_IF::FFV_IF(const Kernel_Info & info) :
  SF_IF12(info), m_jmax(info.GetSplit().Kfcode() < 3 ? 5.:2.)
{
  SetName("F->FV");
}

double FFV_IF::operator()(const Splitting & split) const {
  double z(split.z()), kappa2(split.t()/split.Q2red());
  // Start with the soft term only, including possible K factors
  // (cusp anomalous dimensions), obtained from the gauge part of the kernel
  double Kfactor = m_CMW==1 ? 1.+split.GetKernel()->GetGauge()->K(split) : 1.;
  double value   = A1(z,kappa2) * Kfactor;
  // All massless: just add the collinear parts.
  // TODO: Add the DIS ME correction
  value += B1(z,kappa2);
  if (split.IsMassive()) {
    double mi2(split.m2(0)), mj2(split.m2(1)), mspect2(split.mspect2());
    if (mi2!=0. || mj2!=0.) {
      msg_Debugging()<<METHOD<<" without mass corrections yet.\n";
    }
    double y = split.y(), pipj = split.Q2red()*(1.-y)/(2.*y); 
    if (mspect2>0.) value -= mspect2/pipj;
  }
  return value;
}

double FFV_IF::Integral(const Splitting & split) const {
  double Kmax = (m_CMW==1.) ? (1.+split.GetKernel()->GetGauge()->KMax(split)) : 1.;
  return log(1.+sqr(1.-split.eta())*split.Q2red()/split.tcut()) * Kmax * m_jmax;
}

double FFV_IF::OverEstimate(const Splitting & split) const {
  double Kmax = (m_CMW==1.) ? (1.+split.GetKernel()->GetGauge()->KMax(split)) : 1.;
  return A1(split.z(),split.tcut()/split.Q2red()) * Kmax * m_jmax;
}

void FFV_IF::GeneratePoint(Splitting & split) const {
  double kappa2 = split.tcut()/split.Q2red();
  split.Set_z(1.-sqrt(kappa2 * (pow(1.+sqr(1.-split.eta())/kappa2,ran->Get())-1.)));
  split.Set_phi();
}

double FFV_IF::A1(const double & z,const double & kappa2) const {
  return 2.*(1.-z)/(sqr(1.-z)+kappa2);
}

double FFV_IF::B1(const double & z,const double & kappa2) const {
  return -(1.+z);
}

DECLARE_GETTER(FFV_IF,"IF_FFV",SF_Base,Kernel_Info);

SF_Base * ATOOLS::Getter<SF_Base,Kernel_Info,FFV_IF>::
operator()(const Parameter_Type & info) const
{
  return NULL;
  if (info.Type()==kernel_type::IF &&
      info.GetSplit().IsFermion() && 
      info.GetFlavs()[0].IsFermion() &&
      info.GetFlavs()[1].IsVector()) {
    return new FFV_IF(info);
  }
  return NULL;
}

void ATOOLS::Getter<SF_Base,Kernel_Info,FFV_IF>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"FFV splitting Function (IF)";
}

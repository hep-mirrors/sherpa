#include "CFPSHOWER++/Calculators/IF/SF_IF12.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class VFF_IF : public SF_IF12 {
  private:
    double m_jmax;
    
    double B1(const double & z,const double & kappa2) const;
  public:
    VFF_IF(const Kernel_Info & info);
    double operator()(const Splitting & split);
    double Integral(const Splitting & split) const;
    double OverEstimate(const Splitting & split) const;
    void   GeneratePoint(Splitting & split) const;
  };
}


using namespace CFPSHOWER;
using namespace ATOOLS;

// TODO: may need to monitor the overestimator
VFF_IF::VFF_IF(const Kernel_Info & info) : SF_IF12(info), m_jmax(5.)
{
  SetName("F->VF");
}

double VFF_IF::operator()(const Splitting & split) {
  double z(split.z(0)), kappa2(split.t(0)/split.Q2red());
  // No LL term, so no A1 term and no HO factor
  // TODO: Add the DIS ME correction
  double value = B1(z,kappa2);
  if (split.IsMassive()) {
    double mi2(split.m2(0)), mj2(split.m2(1)), mk2(split.mspect2());
    if (mi2!=0. || mj2!=0.) {
      msg_Debugging()<<METHOD<<" without mass corrections yet.\n";
    }
    double y = split.y(), pipj = split.Q2red()*(1.-y)/(2.*y);
    value -= mk2/pipj;
  }
  return value;
}

double VFF_IF::Integral(const Splitting & split) const {
  return 2.*log(1.0/split.eta()) * m_jmax * PDFEstimate(split);
}

double VFF_IF::OverEstimate(const Splitting & split) const {
  return 2./split.z(0) * m_jmax * PDFEstimate(split);
}

void VFF_IF::GeneratePoint(Splitting & split) const {
  double z = pow(split.eta(),ran->Get());
  split.Set_z(0, z);
  split.Set_phi(0);
}

double VFF_IF::B1(const double & z,const double & kappa2) const {
  return 2./z-(2.-z);
}

DECLARE_GETTER(VFF_IF,"IF_VFF",SF_Base,Kernel_Info);

SF_Base * ATOOLS::Getter<SF_Base,Kernel_Info,VFF_IF>::
operator()(const Parameter_Type & info) const
{
  return NULL;
  if (info.Type()==kernel_type::IF &&
      info.GetSplit().IsVector() && 
      info.GetFlavs()[0].IsFermion() &&
      info.GetFlavs()[1].IsFermion()) {
    return new VFF_IF(info);
  }
  return NULL;
}

void ATOOLS::Getter<SF_Base,Kernel_Info,VFF_IF>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"VFF splitting Function (IF)";
}

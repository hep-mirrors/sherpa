#include "CFPSHOWER++/Calculators/IF/SF_IF.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class FVF_IF : public SF_IF {
  private:
    double m_jmax;
    
    double A1(const double & z,const double & kappa2) const;
    double B1(const double & z,const double & kappa2) const;
    double B2(const double & z,const double & kappa2) const;
  public:
    FVF_IF(const Kernel_Info & info);
    double operator()(const Splitting & split) const;
    double Integral(const Splitting & split) const;
    double OverEstimate(const Splitting & split) const;
    void   GeneratePoint(Splitting & split) const;
  };
}


using namespace CFPSHOWER;
using namespace ATOOLS;

FVF_IF::FVF_IF(const Kernel_Info & info) :
  SF_IF(info),
  // TODO: may need to monitor the overestimator
  m_jmax(5.)
{
  SetName("F->FV");
}

double FVF_IF::operator()(const Splitting & split) const {
  double mij2(split.mij2()), mi2(split.mi2()), mj2(split.mj2()), mk2(split.mk2());
  double z(split.Z()), y(split.Y()), Q2(split.Q2()), kappa2(split.T()/Q2);
  // No LL term, so no A1 term and no HO factor
  double value    = 0.;
  // All massless: just add the collinear parts.
  // TODO: Add the B2 parts as soon as the LL shower is validated.
  // TODO: Add the DIS ME correction
  if (m_orderB>0) {
    value += B1(z,kappa2);
    if (mi2!=0. || mj2!=0.) {
      msg_Debugging()<<METHOD<<" without mass corrections yet.\n";
    }
    if (mk2>0.) value -= 2.*mk2/Q2*y/(1.-y);
  }
  return value;
}

double FVF_IF::Integral(const Splitting & split) const {
  return log(1.0/split.Eta()) * m_jmax * PDFEstimate(split);
}

double FVF_IF::OverEstimate(const Splitting & split) const {
  return 2./split.Z() * m_jmax * PDFEstimate(split);
}

void FVF_IF::GeneratePoint(Splitting & split) const {
  double z      = pow(split.Eta(),ran->Get());
  split.SetZ(z);
  split.Setphi();
}

double FVF_IF::A1(const double & z,const double & kappa2) const {
  return 0.;
}

double FVF_IF::B1(const double & z,const double & kappa2) const {
  return 2./z-(2.-z);
}

double FVF_IF::B2(const double & z,const double & kappa2) const {
  double b2 = 0.0;
  return b2;
}

DECLARE_GETTER(FVF_IF,"IF_FVF",SF_Base,Kernel_Info);

SF_Base * ATOOLS::Getter<SF_Base,Kernel_Info,FVF_IF>::
operator()(const Parameter_Type & info) const
{
  return NULL;
  if (info.Type()==kernel_type::IF &&
      info.GetFlavs()[0].IsFermion() && 
      info.GetFlavs()[1].IsVector() &&
      info.GetFlavs()[2].IsFermion()) {
    return new FVF_IF(info);
  }
  return NULL;
}

void ATOOLS::Getter<SF_Base,Kernel_Info,FVF_IF>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"FVF splitting Function (IF)";
}

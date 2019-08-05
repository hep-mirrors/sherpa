#include "CFPSHOWER++/Calculators/IF/SF_IF12.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class FVF_IF : public SF_IF12 {
  private:
    double m_jmax;
    double B1(const double & z,const double & kappa2) const;
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

FVF_IF::FVF_IF(const Kernel_Info & info) : SF_IF12(info), m_jmax(1.) {
  SetName("V->FF");
}

double FVF_IF::operator()(const Splitting & split) const {
  // still have to fill this in.
  return 0.;
}

double FVF_IF::Integral(const Splitting & split) const {
  return (1.0 - split.eta()) * m_jmax * PDFEstimate(split);
}

double FVF_IF::OverEstimate(const Splitting & split) const {
  return m_jmax * PDFEstimate(split);
}

void FVF_IF::GeneratePoint(Splitting & split) const {
  double eta = split.eta(), z = eta+(1.0-eta)*ran->Get();
  split.Set_z(0,z);
  split.Set_z(1,1.-z);
  split.Set_phi();
}

double FVF_IF::B1(const double & z,const double & kappa2) const {
  return 1.-2.*z*(1.-z);
}

DECLARE_GETTER(FVF_IF,"IF_FVF",SF_Base,Kernel_Info);

SF_Base * ATOOLS::Getter<SF_Base,Kernel_Info,FVF_IF>::
operator()(const Parameter_Type & info) const
{
  return NULL;
  if (info.Type()!=kernel_type::IF) return NULL;
  if ((info.GetSplit().IsFermion() &&
       info.GetFlavs()[0].IsVector() &&
       info.GetFlavs()[1].IsFermion()) ||
      (info.GetSplit().IsFermion() &&
       info.GetFlavs()[0].IsFermion() &&
       info.GetFlavs()[1].IsVector())) return new FVF_IF(info);
  return NULL;
}

void ATOOLS::Getter<SF_Base,Kernel_Info,FVF_IF>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"FVF Splitting Function (IF)";
}

#include "CFPSHOWER++/SplittingFunctions/SF_Base.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class FI_Coll_VFF : public SF_Base {
    bool m_swapped;
    double B1(const Splitting & split) const;
  public:
    FI_Coll_VFF(const Kernel_Info & info);
    double operator()(const Splitting & split);
    double Integral(const Splitting & split) const;
    double OverEstimate(const Splitting & split) const;
    void   GeneratePoint(Splitting & split) const;
  };
}

using namespace CFPSHOWER;
using namespace ATOOLS;

FI_Coll_VFF::FI_Coll_VFF(const Kernel_Info & info)  :
  SF_Base(info),
  m_swapped(info.TagSequence()[0]!=0)
{
  m_name = std::string("FI: collinear V->FF");
}

double FI_Coll_VFF::operator()(const Splitting & split) {
  return B1(split) * (m_swapped ? 1.-split.Z() : split.Z());;
}

double FI_Coll_VFF::Integral(const Splitting & split) const {
  return 1.;
}

double FI_Coll_VFF::OverEstimate(const Splitting & split) const {
  return 1.;
}

void FI_Coll_VFF::GeneratePoint(Splitting & split) const {
  split.SetZ(ran->Get());
  split.SetPhi(2.*M_PI*ran->Get());
}

double FI_Coll_VFF::B1(const Splitting & split) const {
  return sqr(1.-split.Z())+sqr(split.Z());
  double pin = split.Mom(0)*split.GetKinSpect();
  double pjn = split.Mom(1)*split.GetKinSpect();
  return (sqr(pin)+sqr(pjn))/sqr(pin+pjn);
}

DECLARE_GETTER(FI_Coll_VFF,"FI_Coll_VFF",SF_Base,Kernel_Info);

SF_Base * ATOOLS::Getter<SF_Base,Kernel_Info,FI_Coll_VFF>::
operator()(const Parameter_Type & info) const
{
  if (info.Type()==kernel_type::FI &&
      info.LogType()==log_type::coll &&
      int(info.SFType() & 4)>0 &&
      info.GetSplit().IsVector() && 
      info.GetFlavs().size()==2 &&
      info.GetFlavs()[0].IsFermion() && info.GetFlavs()[0].IsAnti() &&
      info.GetFlavs()[1].IsFermion() && !info.GetFlavs()[1].IsAnti()) {
    return new FI_Coll_VFF(info);
  }
  return NULL;
}

void ATOOLS::Getter<SF_Base,Kernel_Info,FI_Coll_VFF>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"VFF splitting function: coll only (FI)";
}



#include "CFPSHOWER++/SplittingFunctions/SF_Base.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class II_Coll_FFV : public SF_Base {
    bool m_swapped;
    double B1(const Splitting & split) const;
  public:
    II_Coll_FFV(const Kernel_Info & info);
    double operator()(const Splitting & split);
    double Integral(const Splitting & split) const;
    double OverEstimate(const Splitting & split) const;
    void   GeneratePoint(Splitting & split) const;
  };
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

using namespace CFPSHOWER;
using namespace ATOOLS;

II_Coll_FFV::II_Coll_FFV(const Kernel_Info & info)  :
  SF_Base(info),
  m_swapped(info.TagSequence()[0]!=0)
{
  m_name = std::string("II: collinear F->FV");
}

double II_Coll_FFV::operator()(const Splitting & split) {
  return B1(split) * (m_swapped ? 1.-split.Z() : split.Z());
}

double II_Coll_FFV::Integral(const Splitting & split) const {
  return 1.;
}

double II_Coll_FFV::OverEstimate(const Splitting & split) const {
  return 1.;
}

void II_Coll_FFV::GeneratePoint(Splitting & split) const {
  split.SetZ(ran->Get());
  split.SetPhi(2.*M_PI*ran->Get());
}

double II_Coll_FFV::B1(const Splitting & split) const {
  return 1.-split.Z();
  double pin = split.Mom(0)*split.GetKinSpect();
  double pjn = split.Mom(1)*split.GetKinSpect();
  return (pjn)/(pin+pjn);
}


DECLARE_GETTER(II_Coll_FFV,"II_Coll_FFV",SF_Base,Kernel_Info);

SF_Base * ATOOLS::Getter<SF_Base,Kernel_Info,II_Coll_FFV>::
operator()(const Parameter_Type & info) const
{
  if (info.Type()==kernel_type::II &&
      info.LogType()==log_type::coll &&
      int(info.SFType() & 1)>0 &&
      info.GetSplit().IsFermion() &&
      info.GetFlavs().size()==2 &&
      info.GetFlavs()[0].IsFermion() && info.GetFlavs()[1].IsVector()) {
    return new II_Coll_FFV(info);
  }
  return NULL;
}

void ATOOLS::Getter<SF_Base,Kernel_Info,II_Coll_FFV>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Collinear FFV splitting function (II)";
}



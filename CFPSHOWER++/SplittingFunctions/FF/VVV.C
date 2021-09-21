#include "CFPSHOWER++/SplittingFunctions/SF_Base.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class FF_Coll_VVV : public SF_Base {
    bool   m_swapped;
    double B1(const Splitting & split) const;
  public:
    FF_Coll_VVV(const Kernel_Info & info);
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

FF_Coll_VVV::FF_Coll_VVV(const Kernel_Info & info)  :
  SF_Base(info),
  m_swapped(info.TagSequence()[0]!=0)
{
  m_name = std::string("FF: collinear V->VV");
}

double FF_Coll_VVV::operator()(const Splitting & split) {
  return B1(split) * (m_swapped ? 1.-split.Z() : split.Z());
}

double FF_Coll_VVV::Integral(const Splitting & split) const {
  return 1./4.;
}

double FF_Coll_VVV::OverEstimate(const Splitting & split) const {
  return 1./4.;
}

void FF_Coll_VVV::GeneratePoint(Splitting & split) const {
  split.SetZ(ran->Get());
  split.SetPhi(2.*M_PI*ran->Get());
}

double FF_Coll_VVV::B1(const Splitting & split) const {
  double pin = split.Mom(0)*split.GetKinSpect();
  double pjn = split.Mom(1)*split.GetKinSpect();
  return pin*pjn/sqr(pin+pjn);
}


DECLARE_GETTER(FF_Coll_VVV,"FF_Coll_VVV",SF_Base,Kernel_Info);

SF_Base * ATOOLS::Getter<SF_Base,Kernel_Info,FF_Coll_VVV>::
operator()(const Parameter_Type & info) const
{
  if (info.Type()==kernel_type::FF &&
      info.LogType()==log_type::coll  &&
      int(info.SFType() & 2)>0 &&
      info.GetSplit().IsVector() &&
      info.GetFlavs().size()==2 &&
      info.GetFlavs()[0].IsVector() && info.GetFlavs()[1].IsVector()) {
    return new FF_Coll_VVV(info);
  }
  return NULL;
}

void ATOOLS::Getter<SF_Base,Kernel_Info,FF_Coll_VVV>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Collinear VVV splitting function (FF)";
}






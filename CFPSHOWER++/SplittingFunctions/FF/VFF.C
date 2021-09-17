#include "CFPSHOWER++/SplittingFunctions/SF_Base.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class FF_Coll_VFF : public SF_Base {
    bool m_swapped;
    double B1(const Splitting & split) const;
  public:
    FF_Coll_VFF(const Kernel_Info & info);
    double operator()(const Splitting & split);
    double Integral(const Splitting & split) const;
    double OverEstimate(const Splitting & split) const;
    void   GeneratePoint(Splitting & split) const;
  };
}

using namespace CFPSHOWER;
using namespace ATOOLS;

FF_Coll_VFF::FF_Coll_VFF(const Kernel_Info & info)  :
  SF_Base(info),
  m_swapped(info.TagSequence()[0]!=0)
{
  m_name = std::string("FF: collinear V->FF");
}

double FF_Coll_VFF::operator()(const Splitting & split) {
  return B1(split) * (m_swapped ? 1.-split.Z() : split.Z());;
}

double FF_Coll_VFF::Integral(const Splitting & split) const {
  return 1.;
}

double FF_Coll_VFF::OverEstimate(const Splitting & split) const {
  return 1.;
}

void FF_Coll_VFF::GeneratePoint(Splitting & split) const {
  split.SetZ(ran->Get());
  split.SetPhi(2.*M_PI*ran->Get());
}

double FF_Coll_VFF::B1(const Splitting & split) const {
  double pin = split.Mom(0)*split.GetKinSpect();
  double pjn = split.Mom(1)*split.GetKinSpect();
  return (sqr(pin)+sqr(pjn))/sqr(pin+pjn);
}

DECLARE_GETTER(FF_Coll_VFF,"FF_Coll_VFF",SF_Base,Kernel_Info);

SF_Base * ATOOLS::Getter<SF_Base,Kernel_Info,FF_Coll_VFF>::
operator()(const Parameter_Type & info) const
{
  return NULL;
  //msg_Out()<<"*** VFF Getter for "
  //	   <<info.GetSplit()<<" --> "
  //	   <<info.GetFlavs()[0]<<" ("<<info.TagSequence()[0]<<") "
  //	   <<info.GetFlavs()[1]<<" ("<<info.TagSequence()[1]<<")\n";
  if (info.Type()==kernel_type::FF &&
      info.LogType()==log_type::coll &&
      info.KinType()==kin_type::CataniSeymour &&
      info.GetSplit().IsVector() && 
      info.GetFlavs().size()==2 &&
      info.GetFlavs()[0].IsFermion() && info.GetFlavs()[0].IsAnti() &&
      info.GetFlavs()[1].IsFermion() && !info.GetFlavs()[1].IsAnti()) {
    return new FF_Coll_VFF(info);
  }
  return NULL;
}

void ATOOLS::Getter<SF_Base,Kernel_Info,FF_Coll_VFF>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"VFF splitting function: coll only (FF)";
}



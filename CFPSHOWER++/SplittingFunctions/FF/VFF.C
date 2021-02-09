#include "CFPSHOWER++/SplittingFunctions/SF_Base.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class FF_VFF_Coll : public SF_Base {
    double B1(const Splitting & split) const;
  public:
    FF_VFF_Coll(const Kernel_Info & info);
    double operator()(const Splitting & split);
    double Integral(const Splitting & split) const;
    double OverEstimate(const Splitting & split) const;
    void   GeneratePoint(Splitting & split) const;
  };
}

using namespace CFPSHOWER;
using namespace ATOOLS;

FF_VFF_Coll::FF_VFF_Coll(const Kernel_Info & info)  : SF_Base(info) {
  SetName("FF: V->FF (coll)");
}

double FF_VFF_Coll::operator()(const Splitting & split) {
  return B1(split);
}

double FF_VFF_Coll::Integral(const Splitting & split) const {
  double I = 1.;
  switch (split.GetKernel()->GetKinematics()->Scheme()) {
  case kin_type::CS:
    break;
  case kin_type::PanGlobal:
    I = 0.;
    break;
  default:
    break;
  }
  return I;
}

double FF_VFF_Coll::OverEstimate(const Splitting & split) const { return 1.; }

void FF_VFF_Coll::GeneratePoint(Splitting & split) const {
  switch (split.GetKernel()->GetKinematics()->Scheme()) {
  case kin_type::CS:
    split.SetZ(ran->Get());
    break;
  case kin_type::PanGlobal:
  default:
    msg_Error()<<"Error in "<<METHOD<<": wrong kinematics scheme for collinear part: "
	       <<split.GetKernel()->GetKinematics()->Scheme()<<"\n";
    exit(1);
    break;
  }
  split.SetPhi(2.*M_PI*ran->Get());
}

double FF_VFF_Coll::B1(const Splitting & split) const {
  return sqr(split.Z()) + sqr(1.-split.Z());
}

DECLARE_GETTER(FF_VFF_Coll,"FF_VFF_Coll",SF_Base,Kernel_Info);

SF_Base * ATOOLS::Getter<SF_Base,Kernel_Info,FF_VFF_Coll>::
operator()(const Parameter_Type & info) const
{
  if (info.Type()==kernel_type::FF &&
      info.LogType()==log_type::coll &&
      info.GetFlavs().size()==2 &&
      info.GetSplit().IsVector() && 
      info.GetFlavs()[0].IsFermion() &&
      info.GetFlavs()[0].IsAnti() &&
      info.TagSequence()[0]==0 &&
      info.GetFlavs()[1].IsFermion()) {
    return NULL; //new FF_VFF_Coll(info);
  }
  return NULL;
}

void ATOOLS::Getter<SF_Base,Kernel_Info,FF_VFF_Coll>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"VFF splitting function: coll only (FF)";
}



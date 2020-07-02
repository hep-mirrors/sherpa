#include "CFPSHOWER++/Calculators/FF/SF_FF2_Coll.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class FF_VVV_Coll : public SF_FF2_Coll {
    double B1(const Splitting & split) const;
  public:
    FF_VVV_Coll(const Kernel_Info & info);
    double operator()(const Splitting & split);
    double Integral(const Splitting & split) const;
    double OverEstimate(const Splitting & split) const;
    void   GeneratePoint(Splitting & split) const;
  };
}

using namespace CFPSHOWER;
using namespace ATOOLS;

FF_VVV_Coll::FF_VVV_Coll(const Kernel_Info & info)  : SF_FF2_Coll(info) {
  SetName("FF: V->VV(coll)");
}

double FF_VVV_Coll::operator()(const Splitting & split) {
  // Collinear part of the kernel, no K factor
  double value   = B1(split);
  if (!split.IsClustered()) value *= m_z[m_tags[0]];
  return value;
}

double FF_VVV_Coll::Integral(const Splitting & split) const { return 1./4.; }

double FF_VVV_Coll::OverEstimate(const Splitting & split) const { return 1./4.; }

void FF_VVV_Coll::GeneratePoint(Splitting & split) const {
  split.SetZ();
  split.SetPhi();
}

double FF_VVV_Coll::B1(const Splitting & split) const {
  return m_z[0]*(1.-m_z[0]);
}


DECLARE_GETTER(FF_VVV_Coll,"FF_VVV_Coll",SF_Base,Kernel_Info);

SF_Base * ATOOLS::Getter<SF_Base,Kernel_Info,FF_VVV_Coll>::
operator()(const Parameter_Type & info) const
{
  if (info.Type()==kernel_type::FF &&
      info.LogType()==log_type::coll &&
      info.GetFlavs().size()==2 &&
      info.GetSplit().IsVector() && 
      info.GetFlavs()[0].IsVector() &&
      info.GetFlavs()[1].IsVector()) {
    return new FF_VVV_Coll(info);
  }
  return NULL;
}

void ATOOLS::Getter<SF_Base,Kernel_Info,FF_VVV_Coll>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"VVV coll splitting Function (FF)";
}



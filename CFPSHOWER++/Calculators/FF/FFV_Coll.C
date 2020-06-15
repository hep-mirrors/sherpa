#include "CFPSHOWER++/Calculators/SF_Base.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class FF_FFV_Coll : public SF_Base {
    double B1(const ATOOLS::Vec4D & psplit,
	      const ATOOLS::Vec4D & pnew,
	      const ATOOLS::Vec4D & pspect) const;
  public:
    FF_FFV_Coll(const Kernel_Info & info);
    double operator()(const Splitting & split);
    double Integral(const Splitting & split) const;
    double OverEstimate(const Splitting & split) const;
    void   GeneratePoint(Splitting & split) const;
  };
}

using namespace CFPSHOWER;
using namespace ATOOLS;

FF_FFV_Coll::FF_FFV_Coll(const Kernel_Info & info)  : SF_Base(info) {
  SetName("FF: F->FV(coll)");
}

double FF_FFV_Coll::operator()(const Splitting & split) {
  // Collinear part - no K factor
  double value   = B1(split.Momentum(0),split.Momentum(1),split.SpectatorMomentum());
  if (split.Clustered()==0) value *= split.ztilde(m_tags[0]);
  return value;
}

double FF_FFV_Coll::Integral(const Splitting & split) const { return 1.; }

double FF_FFV_Coll::OverEstimate(const Splitting & split) const { return 1.; }

void FF_FFV_Coll::GeneratePoint(Splitting & split) const {
  split.Set_z(0,ran->Get());
  split.Set_phi(0);
}

double FF_FFV_Coll::B1(const ATOOLS::Vec4D & psplit,
		       const ATOOLS::Vec4D & pnew,
		       const ATOOLS::Vec4D & pspect) const {
  double pipk = psplit*pspect;
  double pjpk = pnew*pspect;
  return pjpk/(pipk+pjpk);
}


DECLARE_GETTER(FF_FFV_Coll,"FF_FFV_Coll",SF_Base,Kernel_Info);

SF_Base * ATOOLS::Getter<SF_Base,Kernel_Info,FF_FFV_Coll>::
operator()(const Parameter_Type & info) const
{
  if (info.Type()==kernel_type::FF &&
      info.LogType()==log_type::coll &&
      info.GetSplit().IsFermion() && 
      info.GetFlavs()[0].IsFermion() &&
      info.GetFlavs()[1].IsVector()) {
    return new FF_FFV_Coll(info);
  }
  return NULL;
}

void ATOOLS::Getter<SF_Base,Kernel_Info,FF_FFV_Coll>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"FFV coll splitting Function (FF)";
}



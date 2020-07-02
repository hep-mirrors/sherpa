#include "CFPSHOWER++/Calculators/FF/SF_FF2_Coll.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class FF_VFF_Coll : public SF_FF2_Coll {
    double B1(const Splitting & split) const;
    double Massive(const Splitting & split) const;
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

FF_VFF_Coll::FF_VFF_Coll(const Kernel_Info & info)  : SF_FF2_Coll(info) {
  SetName("FF: V->FF(coll)");
}

double FF_VFF_Coll::operator()(const Splitting & split) {
  // collinear part of the kernel
  double value = m_ismassive ? Massive(split) : B1(split);
  if (!split.IsClustered()) value *= m_z[m_tags[0]];
  return value;
}

double FF_VFF_Coll::Integral(const Splitting & split) const {
  return 1.;
}

double FF_VFF_Coll::OverEstimate(const Splitting & split) const {
  return 1.;
}

void FF_VFF_Coll::GeneratePoint(Splitting & split) const {
  split.SetZ(); 
  split.SetPhi();
}

double FF_VFF_Coll::B1(const Splitting & split) const {
  return sqr(m_z[0]) + sqr(1.-m_z[0]);
}

double FF_VFF_Coll::Massive(const Splitting & split) const {
  double value = 0.;
  double vijk2 = sqr(1.-m_y) - 4.*(m_y+2.*m_m2[0]/m_shat)*m_mspect2/m_shat;
  if (vijk2 > 0.0) {
    double s01 = (m_moms[0]+m_moms[1]).Abs2();
    value += (1.-m_y)/sqrt(vijk2) * (B1(split) + 2.*m_m2[0]/s01);
  }
  return value;
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
      info.GetFlavs()[1].IsFermion()) {
    return new FF_VFF_Coll(info);
  }
  return NULL;
}

void ATOOLS::Getter<SF_Base,Kernel_Info,FF_VFF_Coll>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"VFF coll splitting Function (FF)";
}



#include "CFPSHOWER++/Calculators/FF/SF_FF12.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class FFV_FF : public SF_FF12 {
    double A1(const double & z,const double & kappa2) const;
    double B1(const double & z,const double & kappa2) const;
  public:
    FFV_FF(const Kernel_Info & info);
    double operator()(const Splitting & split) const;
    double Integral(const Splitting & split) const;
    double OverEstimate(const Splitting & split) const;
    void   GeneratePoint(Splitting & split) const;
  };
}

using namespace CFPSHOWER;
using namespace ATOOLS;

FFV_FF::FFV_FF(const Kernel_Info & info)  : SF_FF12(info) {
  SetName(m_tagsequence[0]==1?"F->VF":"F->FV");
}

double FFV_FF::operator()(const Splitting & split) const {
  double mij2(split.m2(0)), mi2(split.m2(1)), mj2(split.m2(2)), mspect2(split.mspect2());
  double z(split.z(0)), y(split.y()), Q2(split.Q2()), kappa2(split.t()/Q2);
  // Start with the soft term only, including possible K factors (cusp anomalous
  // dimensions), obtained from the gauge part of the kernel
  double Kfactor = (m_CMW==1) ? (1.+split.GetKernel()->GetGauge()->K(split)) : 1.;
  double value   = A1(z,kappa2) * Kfactor;
  // All massless: just add the collinear parts.
  if (mi2==0 && mj2==0 && mspect2==0.) value += B1(z,kappa2);
  else {
    // Massive splitting adjustments
    // directly return 0 if the splitting is kinematically not viable 
    double s(split.s());
    double muij2(mij2/s), mui2(mi2/s), muspect2(mspect2/s);
    double vijl  = sqr(2.*muspect2+(1.-mui2-muspect2)*(1.-y))-4.*muspect2;
    double vtijl = Lambda2(1.,muij2,muspect2);
    if (vtijl<0. || vijl<0.) return 0.;
    double vtlji, vlji;
    if (muspect2>0.) {
      vtlji = 1.-4*muspect2*mui2/sqr(1.-muspect2-mui2);
      vlji  = 1.-4*(Q2*(1.-z)+mspect2)*mi2/sqr(Q2*z);
      if (vtlji<0. || vlji<0.) return 0.;
    }
    vtijl        = sqrt(vtijl)/(1.-muij2-muspect2);
    vijl         = sqrt(vijl)/((1.-muij2-muspect2) * (1.-y));
    value       += (vtijl/vijl) * ( B1(z,kappa2) -
				    (2.*(1.-z)*mi2)/((1.-z+y)*Q2*y));
    if (muspect2>0.) {
      value    -= 2.*sqrt(vtlji/vlji)*mspect2/((1.-z)*Q2) * y/((1.-z)+y);
    }
  }
  if (split.Clustered()==0) value *= split.z(m_tagsequence[0]);
  return value;
}

double FFV_FF::Integral(const Splitting & split) const {
  double Kmax = (m_CMW==1.) ? (1.+split.GetKernel()->GetGauge()->KMax(split)) : 1.;
  return log(1.0+split.Q2()/split.t0()) * Kmax;
}

double FFV_FF::OverEstimate(const Splitting & split) const {
  double Kmax = (m_CMW==1.) ? (1.+split.GetKernel()->GetGauge()->KMax(split)) : 1.;
  return A1(split.z(0),split.t0()/split.Q2()) * Kmax;
}

void FFV_FF::GeneratePoint(Splitting & split) const {
  double kappa2 = split.t0()/split.Q2();
  double z      = 1.-sqrt(kappa2 * (pow(1.+1./kappa2,ran->Get())-1.) );
  split.Set_z(0,z);
  split.Set_z(1,1.-z);
  split.Set_phi();
}

double FFV_FF::A1(const double & z,const double & kappa2) const {
  return 2.*(1.-z)/(ATOOLS::sqr(1.-z)+kappa2);
}

double FFV_FF::B1(const double & z,const double & kappa2) const {
  return -(1.+z);
}

DECLARE_GETTER(FFV_FF,"FF_FFV",SF_Base,Kernel_Info);

SF_Base * ATOOLS::Getter<SF_Base,Kernel_Info,FFV_FF>::
operator()(const Parameter_Type & info) const
{
  if (info.Type()==kernel_type::FF &&
      info.GetSplit().IsFermion() && 
      info.GetFlavs()[0].IsFermion() &&
      info.GetFlavs()[1].IsVector()) {
    return new FFV_FF(info);
  }
  return NULL;
}

void ATOOLS::Getter<SF_Base,Kernel_Info,FFV_FF>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"FFV splitting Function (FF)";
}

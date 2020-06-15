#include "CFPSHOWER++/Calculators/SF_Base.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class FF_FFV_Soft : public SF_Base {
    double A1(const ATOOLS::Vec4D & psplit,
	      const ATOOLS::Vec4D & pnew,
	      const ATOOLS::Vec4D & pspect) const;
  public:
    FF_FFV_Soft(const Kernel_Info & info);
    double operator()(const Splitting & split);
    double Integral(const Splitting & split) const;
    double OverEstimate(const Splitting & split) const;
    void   GeneratePoint(Splitting & split) const;
  };
}

using namespace CFPSHOWER;
using namespace ATOOLS;

FF_FFV_Soft::FF_FFV_Soft(const Kernel_Info & info)  : SF_Base(info) {
  SetName("FF: F->FV(soft)");
}

double FF_FFV_Soft::operator()(const Splitting & split) {
  // Soft term, including possible K factors (cusp anomalous
  // dimensions), obtained from the gauge part of the kernel
  double Kfactor = (m_CMW==1) ? (1.+split.GetKernel()->GetGauge()->K(split)) : 1.;
  double value   = A1(split.Momentum(0),split.Momentum(1),split.SpectatorMomentum()) * Kfactor;
  if (split.Clustered()==0) value *= split.ztilde(m_tags[0]);
  return value;
}

double FF_FFV_Soft::Integral(const Splitting & split) const {
  double Kmax = (m_CMW==1.) ? (1.+split.GetKernel()->GetGauge()->KMax(split)) : 1.;
  double I    = 2.*log((1.-split.zmin())/(1.-split.zmax()));
  return I * Kmax;
}

double FF_FFV_Soft::OverEstimate(const Splitting & split) const {
  double Kmax = (m_CMW==1.) ? (1.+split.GetKernel()->GetGauge()->KMax(split)) : 1.;
  double SF   = 2./(1.-split.z(0));
  return SF * Kmax;  
}

void FF_FFV_Soft::GeneratePoint(Splitting & split) const {
  double z = 1. - (1.-split.zmin()) * pow((1.-split.zmax())/(1.-split.zmin()), ran->Get());
  split.Set_z(0,z); 
  split.Set_phi(0);
}

double FF_FFV_Soft::A1(const ATOOLS::Vec4D & psplit,
		       const ATOOLS::Vec4D & pnew,
		       const ATOOLS::Vec4D & pspect) const {
  double pipk = psplit*pspect;
  double pjpk = pnew*pspect;
  double pipj = psplit*pnew;
  return 2.*pipk/(pipj+pjpk);
}



DECLARE_GETTER(FF_FFV_Soft,"FF_FFV_Soft",SF_Base,Kernel_Info);

SF_Base * ATOOLS::Getter<SF_Base,Kernel_Info,FF_FFV_Soft>::
operator()(const Parameter_Type & info) const
{
  if (info.Type()==kernel_type::FF &&
      info.LogType()==log_type::soft &&
      info.GetSplit().IsFermion() && 
      info.GetFlavs()[0].IsFermion() &&
      info.GetFlavs()[1].IsVector()) {
    return new FF_FFV_Soft(info);
  }
  return NULL;
}

void ATOOLS::Getter<SF_Base,Kernel_Info,FF_FFV_Soft>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"FFV soft splitting Function (FF)";
}



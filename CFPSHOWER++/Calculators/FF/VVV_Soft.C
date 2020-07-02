#include "CFPSHOWER++/Calculators/FF/SF_FF2_Soft.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class FF_VVV_Soft : public SF_FF2_Soft {
    double A1(const Splitting & split) const;
    bool   MassCheck(const Splitting & split) const;
    double MassCorrection(const Splitting & split) const;
  public:
    FF_VVV_Soft(const Kernel_Info & info);
    double operator()(const Splitting & split);
    double Integral(const Splitting & split) const;
    double OverEstimate(const Splitting & split) const;
    void   GeneratePoint(Splitting & split) const;
  };
}

using namespace CFPSHOWER;
using namespace ATOOLS;

FF_VVV_Soft::FF_VVV_Soft(const Kernel_Info & info)  : SF_FF2_Soft(info) {
  SetName("FF: V->VV(soft)");
  msg_Out()<<METHOD<<".\n";
}

double FF_VVV_Soft::operator()(const Splitting & split) {
  // Soft term, including possible K factors (cusp anomalous
  // dimensions), obtained from the gauge part of the kernel
  double value = 0.;
  if (MassCheck(split)) {
    double Kfactor = (m_CMW==1) ? (1.+split.GetKernel()->GetGauge()->K(split)) : 1.;
    value          = (A1(split) + MassCorrection(split)) * Kfactor;
    if (!split.IsClustered()) value *= m_z[m_tags[0]];
  }
  return value;
}

double FF_VVV_Soft::Integral(const Splitting & split) const {
  double Kmax = (m_CMW==1.) ? (1.+split.GetKernel()->GetGauge()->KMax(split)) : 1.;
  double I    = 2.*log((1.-split.Zmin())/(1.-split.Zmax()));
  return I * Kmax;
}

double FF_VVV_Soft::OverEstimate(const Splitting & split) const {
  double Kmax = (m_CMW==1.) ? (1.+split.GetKernel()->GetGauge()->KMax(split)) : 1.;
  double SF   = 2./(1.-split.Z());
  return SF * Kmax;  
}

void FF_VVV_Soft::GeneratePoint(Splitting & split) const {
  double z = 1. - (1.-split.Zmin()) * pow((1.-split.Zmax())/(1.-split.Zmin()), ran->Get());
  split.SetZ(0,z); 
  split.SetPhi();
}

bool FF_VVV_Soft::MassCheck(const Splitting & split) const {
  if (m_ismassive) {
    double s01      = (m_moms[0]+m_moms[1]).Abs2();
    double vt2_ij_k = Lambda2(m_shat,m_msplit2,m_mspect2);
    double v2_ij_k  = Lambda2(m_shat,s01,m_mspect2);
    return ( sqrt(s01)    > m_m[0]+m_m[1] &&
	     sqrt(m_shat) > sqrt(s01)+m_spect &&
	     vt2_ij_k     > 0. &&
	     v2_ij_k      > 0.);
  }
  return true;
}

double FF_VVV_Soft::A1(const Splitting & split) const {
  return 2.*m_z[0]*(1.-m_y)/(1.-m_z[0]*(1.-m_y));
}

double FF_VVV_Soft::MassCorrection(const Splitting & split) const {
  double mass_corr = 0.;
  if (m_ismassive && m_mspect2>0.) {
    mass_corr -= 2.*m_mspect2/m_pp[1][2] * (1.-m_y)*m_y/(1.-m_z[0]+m_y);
    //msg_Out()<<METHOD<<" yields mass correction = "<<mass_corr<<"\n";
  }
  return mass_corr;
}

DECLARE_GETTER(FF_VVV_Soft,"FF_VVV_Soft",SF_Base,Kernel_Info);

SF_Base * ATOOLS::Getter<SF_Base,Kernel_Info,FF_VVV_Soft>::
operator()(const Parameter_Type & info) const
{
  if (info.Type()==kernel_type::FF &&
      info.LogType()==log_type::soft &&
      info.GetFlavs().size()==2 &&
      info.GetSplit().IsVector() && 
      info.GetFlavs()[0].IsVector() &&
      info.GetFlavs()[1].IsVector()) {
    return new FF_VVV_Soft(info);
  }
  return NULL;
}

void ATOOLS::Getter<SF_Base,Kernel_Info,FF_VVV_Soft>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"VVV soft splitting Function (FF)";
}



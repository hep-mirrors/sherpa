#include "CFPSHOWER++/Calculators/FF/SF_FF12.H"
#include "CFPSHOWER++/Calculators/SF_Base.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class VVV_FF : public SF_FF12 {
  private:
    double A1(const double & z,const double & kappa2) const;
    double B1(const double & z,const double & kappa2) const;
  public:
    VVV_FF(const Kernel_Info & info);
    double operator()(const Splitting & split) const;
    double Integral(const Splitting & split) const;
    double OverEstimate(const Splitting & split) const;
    void   GeneratePoint(Splitting & split) const;
  };
}


using namespace CFPSHOWER;
using namespace ATOOLS;

VVV_FF::VVV_FF(const Kernel_Info & info) : SF_FF12(info) {
  SetName("V->VV");
}

double VVV_FF::operator()(const Splitting & split) const {
  double mspect2(split.mspect2());
  double z(split.z(0)), y(split.y()), Q2(split.Q2()), kappa2(split.t()/Q2);
  // Start with the soft term only, including possible K factors (cusp anomalous
  // dimensions), obtained from the gauge part of the kernel
  double Kfactor = m_CMW==1 ? (1.+split.GetKernel()->GetGauge()->K(split)) : 1.;
  double value   = A1(z,kappa2) * Kfactor;
  if (mspect2==0.) value += B1(z,kappa2);
  else {
    // Massive spectator adjustments
    // directly return 0 if the splitting is kinematically not viable 
    double muspect2(mspect2/Q2);
    double vijl  = sqr(1.-y)-4.*y*muspect2;
    if (vijl<0.) return 0.;
    value += (1.-y)*B1(z,kappa2)/sqrt(vijl) - 2.*(muspect2*y)/((1.-z)*(1.-z+y));
  }
  if (split.Clustered()==0) value *= split.z(m_tagsequence[0]);
  return value;
}

double VVV_FF::Integral(const Splitting & split) const {
  double Kmax = (m_CMW==1.) ? (1.+split.GetKernel()->GetGauge()->KMax(split)) : 1.;
  return log(1.0+split.Q2()/split.t0()) * Kmax;
}

double VVV_FF::OverEstimate(const Splitting & split) const {
  double Kmax = (m_CMW==1.) ? (1.+split.GetKernel()->GetGauge()->KMax(split)) : 1.;
  return A1(split.z(0),split.t0()/split.Q2()) * Kmax;
}

void VVV_FF::GeneratePoint(Splitting & split) const {
  double kappa2 = split.t0()/split.Q2();
  double z     = 1.-sqrt(kappa2 * (pow((1.+1./kappa2),ran->Get())-1.)); 
  split.Set_z(0,z);
  split.Set_z(1,1.-z);
  split.Set_phi();
}

double VVV_FF::A1(const double & z,const double & kappa2) const {
  return 2.*(1.-z)/(sqr(1.-z)+kappa2);
}

double VVV_FF::B1(const double & z,const double & kappa2) const {
  return -2.+z*(1.-z);
}

DECLARE_GETTER(VVV_FF,"FF_VVV",SF_Base,Kernel_Info);

SF_Base * ATOOLS::Getter<SF_Base,Kernel_Info,VVV_FF>::
operator()(const Parameter_Type & info) const
{
  msg_Out()<<METHOD<<" for "<<info<<"   * check: "
	   <<"[type = "<<info.Type()<<"] and "
	   <<info.GetSplit().IsVector()
	   <<info.GetFlavs()[0].IsVector()
	   <<info.GetFlavs()[1].IsVector()<<"\n"; 
  if (info.Type()==kernel_type::FF &&
      (info.GetSplit().IsVector() &&
       info.GetFlavs()[0].IsVector() &&
       info.GetFlavs()[1].IsVector())) {
    return new VVV_FF(info);
  }
  return NULL;
}

void ATOOLS::Getter<SF_Base,Kernel_Info,VVV_FF>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"VVV Splitting Function (FF)";
}

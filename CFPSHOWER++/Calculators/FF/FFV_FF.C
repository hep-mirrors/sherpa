#include "CFPSHOWER++/Calculators/FF/SF_FF.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

//#include "CFPSHOWER++/Calculators/QQG.C"

namespace CFPSHOWER {
  class FFV_FF : public SF_FF {
    double A1(const double & z,const double & kappa2) const;
    double B1(const double & z,const double & kappa2) const;
    double B2(const double & z,const double & kappa2) const;
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

FFV_FF::FFV_FF(const Kernel_Info & info)  : SF_FF(info) {
  SetName(m_swap?"F->VF":"F->FV");
}

double FFV_FF::operator()(const Splitting & split) const {
  double mij2(split.mij2()), mi2(split.mi2()), mj2(split.mj2()), mk2(split.mk2());
  double z(split.Z()), y(split.Y()), Q2(split.Q2()), kappa2(split.T()/Q2);
  // Start with the soft term only, including possible K factors (cusp anomalous
  // dimensions), obtained from the gauge part of the kernel
  double Kfactor = (1.+split.GetKernel()->GetGauge()->K(split));
  double value   = A1(z,kappa2) * Kfactor, Beff = 0.;
  // All massless: just add the collinear parts.
  if (mi2==0 && mj2==0 && mk2==0.) {
    if (m_orderB>0) value += Beff = B1(z,kappa2);
  }
  // Massive splitting - directly return 0 if the splitting is kinematically
  // not viable 
  else {
    double sijk(split.sijk());
    double muij2(mij2/sijk), mui2(mi2/sijk), muk2(mk2/sijk);
    double vijk  = sqr(2.*muk2+(1.-mui2-muk2)*(1.-y))-4.*muk2;
    double vtijk = Lambda2(1.,muij2,muk2);
    if (vtijk<0. || vijk<0.) return 0.;
    double vtkji, vkji;
    if (muk2>0.) {
      vtkji = 1.-4*muk2*mui2/sqr(1.-muk2-mui2);
      vkji  = 1.-4*(Q2*(1.-z)+mk2)*mi2/sqr(Q2*z);
      if (vtkji<0. || vkji<0.) return 0.;
    }
    vtijk        = sqrt(vtijk)/(1.-muij2-muk2);
    vijk         = sqrt(vijk)/((1.-muij2-muk2) * (1.-y));
    double pipj  = Q2*y/2.;
    value       += Beff = (vtijk/vijk) * (B1(z,kappa2) -
					  ((1.-z)*mi2)/((1.-z+y)*pipj));
    if (muk2>0.) {
      value    -= 2.*sqrt(vtkji/vkji)*mk2/((1.-z)*Q2)*y/((1.-z)+y);
      Beff     -= 2.*sqrt(vtkji/vkji)*mk2/((1.-z)*Q2)*y/((1.-z)+y);
    }
  }
  if (split.Clustered()==0) value *= (m_swap)?(1.-z):z;
  return value;
}

double FFV_FF::Integral(const Splitting & split) const {
  double Kmax = (1.+split.GetKernel()->GetGauge()->KMax(split));
  return log(1.0+split.Q2()/split.T0()) * Kmax;
}

double FFV_FF::OverEstimate(const Splitting & split) const {
  double Kmax = (1.+split.GetKernel()->GetGauge()->KMax(split));
  return A1(split.Z(),split.T0()/split.Q2()) * Kmax;
}

void FFV_FF::GeneratePoint(Splitting & split) const {
  double kappa2 = split.T0()/split.Q2();
  split.SetZ(1.-sqrt(kappa2 * (pow(1.+1./kappa2,ran->Get())-1.) ) );
  split.Setphi();
}

double FFV_FF::A1(const double & z,const double & kappa2) const {
  return 2.*(1.-z)/(ATOOLS::sqr(1.-z)+kappa2);
}

double FFV_FF::B1(const double & z,const double & kappa2) const {
  return -(1.+z);
}

double FFV_FF::B2(const double & z,const double & kappa2) const {
  double b2 = 0.0;
  return b2;
}


DECLARE_GETTER(FFV_FF,"FF_FFV",SF_Base,Kernel_Info);

SF_Base * ATOOLS::Getter<SF_Base,Kernel_Info,FFV_FF>::
operator()(const Parameter_Type & info) const
{
  return NULL;
  if (info.Type()==kernel_type::FF &&
      info.GetFlavs()[0].IsFermion() && 
      info.GetFlavs()[1].IsFermion() &&
      info.GetFlavs()[2].IsVector()) {
    return new FFV_FF(info);
  }
  return NULL;
}

void ATOOLS::Getter<SF_Base,Kernel_Info,FFV_FF>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"FFV splitting Function (FF)";
}

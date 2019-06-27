#include "CFPSHOWER++/Calculators/FF/SF_FF.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

//#include "CFPSHOWER++/Calculators/QGQ.C"

namespace CFPSHOWER {
  class FVF_FF : public SF_FF {
  private:
    double A1inv(const double & z,const double & kappa2) const;
    double B1(const double & z,const double & kappa2) const;
    double B2(const double & z,const double & kappa2) const;
  public:
    FVF_FF(const Kernel_Info & info);
    double operator()(const Splitting & split) const;
    double Integral(const Splitting & split) const;
    double OverEstimate(const Splitting & split) const;
    void   GeneratePoint(Splitting & split)  const;
  };
}


using namespace CFPSHOWER;
using namespace ATOOLS;

FVF_FF::FVF_FF(const Kernel_Info & info) : SF_FF(info) {
  SetName("F->VF");
}

double FVF_FF::operator()(const Splitting & split) const {
  double mij2(split.mij2()), mi2(split.mi2()), mj2(split.mj2()), mk2(split.mk2());
  double z(split.Z()), y(split.Y()), Q2(split.Q2()), kappa2(split.T()/split.Q2());
  double Kfactor  = (1.+split.GetKernel()->GetGauge()->K(split));
  double value    = A1inv(z,kappa2) * Kfactor;
  if (mi2==0. && mj2==0. && mk2==0.) {
    value += B1(z,kappa2);
  }
  else {
    double sijk(split.sijk());
    double muij2(mij2/sijk), muj2(mj2/sijk), muk2(mk2/sijk);
    double vijk  = ATOOLS::sqr(2.*muk2+(1.-muj2-muk2)*(1.-y))-4.*muk2;
    double vtijk = Lambda2(1.,muij2,muk2);
    if (vtijk<0. || vijk<0.) return 0.;
    double vtkji, vkji;
    if (muk2>0.) {
      vtkji = 1.-4*muk2*muj2/sqr(1.-muk2-muj2);
      vkji  = 1.-4*(Q2*z+mk2)*mj2/sqr(Q2*(1.-z)); // 1-z  <--> z below
      if (vtkji<0. || vkji<0.) return 0.;
    }
    vtijk       = sqrt(vtijk)/(1.-muij2-muk2);
    vijk        = sqrt(vtijk)/((1.-muj2-muk2) * (1.-y));
    double pipj = Q2*y/2.;
    value      += (vtijk/vijk) * (B1(z,kappa2) - (z*mj2)/((z+y)*pipj));
    if (muk2>0.) {
      value    -= 2.*sqrt(vtkji/vkji)*mk2/(z*Q2)*y/(z+y);
    }
  }
  if (split.Clustered()==0) value *= z;
  return value;
}

double FVF_FF::Integral(const Splitting & split) const {
  double Kmax = (1.+split.GetKernel()->GetGauge()->KMax(split));
  return log(1.0+split.Q2()/split.T0()) * Kmax;
}

double FVF_FF::OverEstimate(const Splitting & split) const {
  double Kmax = (1.+split.GetKernel()->GetGauge()->KMax(split));
  return A1inv(split.Z(),split.T0()/split.Q2()) * Kmax;
}

void FVF_FF::GeneratePoint(Splitting & split) const {
  double kappa2 = split.T0()/split.Q2();
  split.SetZ(sqrt(kappa2 * (pow((1.+1./kappa2),ran->Get())-1.) ));
  split.Setphi();
}

double FVF_FF::A1inv(const double & z,const double & kappa2) const {
  return 2.*z/(ATOOLS::sqr(z)+kappa2);
}
double FVF_FF::B1(const double & z,const double & kappa2) const {
  return -(2.-z);
}

double FVF_FF::B2(const double & z,const double & kappa2) const {
  double b2 = 0.0;
  return b2;
}

using namespace CFPSHOWER;

DECLARE_GETTER(FVF_FF,"FF_FVF",SF_Base,Kernel_Info);

SF_Base *
ATOOLS::Getter<SF_Base,Kernel_Info,FVF_FF>::
operator()(const Parameter_Type & info) const
{
  return NULL;
  if (info.Type()==kernel_type::FF &&
      info.GetFlavs()[0].IsFermion() && 
      info.GetFlavs()[1].IsVector() &&
      info.GetFlavs()[2].IsFermion()) {
    //return new FVF_FF(info);
  }
  return NULL;
}

void ATOOLS::Getter<SF_Base,Kernel_Info,FVF_FF>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"FVF Splitting Function (FF)";
}

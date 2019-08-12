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
  SetName(m_tags[0]==1?"F->VF":"F->FV");
}

double FFV_FF::operator()(const Splitting & split) const {
  double z(split.z()), Q2red(split.Q2red()), kappa2(split.t()/Q2red);
  // Start with the soft term only, including possible K factors (cusp anomalous
  // dimensions), obtained from the gauge part of the kernel
  double Kfactor = (m_CMW==1) ? (1.+split.GetKernel()->GetGauge()->K(split)) : 1.;
  double value   = A1(z,kappa2) * Kfactor;
  // All massless: just add the collinear parts.
  if (split.IsMassive()) {
    // Massive splitting adjustments
    // directly return 0 if the splitting is kinematically not viable
    double mi2(split.m2(m_invtags[0])), mj2(split.m2(m_invtags[1]));
    double mk2(split.mspect2()), mij2 = split.msplit2(), y = split.y();
    double Q2 = split.Q2(), sij = y*(Q2-mk2)+(1.-y)*(mi2+mj2);
    double vt2_ij_k = Lambda2(Q2,mij2,mk2), v2_ij_k  = Lambda2(Q2,sij,mk2);
    if ( sij<mi2+mj2 || Q2<sij+mk2 || vt2_ij_k<0. || v2_ij_k < 0.) return 0.;
    double vt_ij_k = sqrt(vt2_ij_k)/(Q2-mij2-mk2), v_ij_k = sqrt(v2_ij_k)/(Q2-sij-mk2);
    // Terms ~y/(1-z+y) come from reshuffling the mass-dependent terms in the eikonal
    // to ensure correct soft behaviour
    value += vt_ij_k/v_ij_k * (B1(z,kappa2) +
			       2.*mi2/(sij-mi2-mj2) * (1.- y/(1.-z+y)) );
    if (mk2>0.) {
      double vt_kj_i = Lambda(Q2,mk2,mi2)/(Q2-mk2-mi2);
      double v_kj_i  = sqrt(1.- 4.*mi2*(Q2*(1.-z)+mk2)/sqr(Q2*z) );
      value  -= (vt_kj_i/v_kj_i  *
		 2.*mk2/(Q2*(1.-z)) * y/(1.-z+y));
    }
  }
  else value += B1(z,kappa2);
  if (split.Clustered()==0) value *= (m_tags[0]==0) ? z : 1-z;
  return value;
}

double FFV_FF::Integral(const Splitting & split) const {
  double Kmax = (m_CMW==1.) ? (1.+split.GetKernel()->GetGauge()->KMax(split)) : 1.;
  return log(1.0+split.Q2red()/split.tcut()) * Kmax;
}

double FFV_FF::OverEstimate(const Splitting & split) const {
  double Kmax = (m_CMW==1.) ? (1.+split.GetKernel()->GetGauge()->KMax(split)) : 1.;
  return A1(split.z(),split.tcut()/split.Q2red()) * Kmax;
}

void FFV_FF::GeneratePoint(Splitting & split) const {
  double kappa2 = split.tcut()/split.Q2red();
  split.Set_z(1.-sqrt(kappa2 * (pow(1.+1./kappa2,ran->Get())-1.) ));
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

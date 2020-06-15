#include "CFPSHOWER++/Calculators/FF/SF_FF12.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class FFV_FF : public SF_FF12 {
    double A1(const double & z,const double & kappa2) const;
    double B1(const double & z,const double & kappa2) const;
    double SoftCorrelation(const Splitting & split) const;
  public:
    FFV_FF(const Kernel_Info & info);
    double operator()(const Splitting & split);
    double Integral(const Splitting & split) const;
    double OverEstimate(const Splitting & split) const;
    void   GeneratePoint(Splitting & split) const;
  };
}

using namespace CFPSHOWER;
using namespace ATOOLS;

FFV_FF::FFV_FF(const Kernel_Info & info)  : SF_FF12(info) {
  SetName("F->FV");
}

double FFV_FF::operator()(const Splitting & split) {
  double z(split.z(0)), kappa2(split.t(0)/split.Q2red());
  // Start with the soft term only, including possible K factors (cusp anomalous
  // dimensions), obtained from the gauge part of the kernel
  double Kfactor = (m_CMW==1) ? (1.+split.GetKernel()->GetGauge()->K(split)) : 1.;
  double value   = A1(z,kappa2) * Kfactor;
  if (split.IsMassive()) {
    // Massive splitting adjustments
    // directly return 0 if the splitting is kinematically not viable
    double mi2(split.m2(0)), mj2(split.m2(1));
    double mk2(split.mspect2()), mij2(split.msplit2()), y(split.y());
    double Q2(split.Q2()), sij(y*(Q2-mk2)+(1.-y)*(mi2+mj2)); // = Q2red + (mi2+mj2)
    double vt2_ij_k(Lambda2(Q2,mij2,mk2)), v2_ij_k(Lambda2(Q2,sij,mk2));
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
  // All massless: just add the collinear parts plus, possibly, soft reweighting.
  else {
    // soft corrections
    if (m_tags[0]==0) {
      // explicit correlation for "tagged" quark only if soft partners are set correctly
      if (m_softcorr>0 &&
	  split.GetSplitter()->SoftPartner(0)==split.GetSpectator() &&
	  split.GetSplitter()->SoftPartner(2)!=NULL)
	value += SoftCorrelation(split);
      // adding soft endpoint
      if (m_endpoint>0)
	value += A1(z,kappa2) * split.GetKernel()->GetGauge()->SoftEndpoint(split);
    }
    value += B1(z,kappa2);
  }
  if (split.Clustered()==0) value *= split.ztilde(m_tags[0]); // m_tags[0]==0?z:1.-z;
  /*msg_Out()<<METHOD<<" for t = "<<split.t(0)<<", z = "<<split.z(0)
  	   <<", y = "<<split.y()<<", ztilde = "<<split.ztilde(0)<<"\n"
  	   <<"*** "<<split
  	   <<"*** value = "<<value<<" from A = "<<A1(z,kappa2)
  	   <<", B = "<<B1(z,kappa2)
  	   <<", endpoint = "<<split.GetKernel()->GetGauge()->SoftEndpoint(split)<<".\n";
  */
  return value;
}

double FFV_FF::Integral(const Splitting & split) const {
  double Kmax = (m_CMW==1.) ? (1.+split.GetKernel()->GetGauge()->KMax(split)) : 1.;
  return log(1.0+split.Q2red()/split.tcut()) * Kmax;
}

double FFV_FF::OverEstimate(const Splitting & split) const {
  double Kmax = (m_CMW==1.) ? (1.+split.GetKernel()->GetGauge()->KMax(split)) : 1.;
  return A1(split.z(0),split.tcut()/split.Q2red()) * Kmax;
}

void FFV_FF::GeneratePoint(Splitting & split) const {
  double kappa2 = split.tcut()/split.Q2red();
  double z      = 1.-sqrt(kappa2 * (pow(1.+1./kappa2,ran->Get())-1.)); 
  split.Set_z(0,z);
  split.Set_phi(0);
}

double FFV_FF::A1(const double & z,const double & kappa2) const {
  return 2.*(1.-z)/(ATOOLS::sqr(1.-z)+kappa2);
}

double FFV_FF::B1(const double & z,const double & kappa2) const {
  return -(1.+z);
}

double FFV_FF::SoftCorrelation(const Splitting & split) const {
  Vec4D pi = m_moms[0], p1 = m_moms[1], p2 = m_specmom;
  Vec4D pj = split.GetSplitter()->SoftPartner(2)->Mom();
  double sij = 2.*pi*pj, si1 = 2.*pi*p1, si2 = 2.*pi*p2;
  double sj1 = 2.*pj*p1, sj2 = 2.*pj*p2, s12 = 2.*p1*p2;
  // Equations (49) and (54) in 1805.03752
  double w_ij_12    = 1.-sij*s12/((si1+si2)*(sj1+sj2));
  double wbar_ij_12 = ((si1+si2)*(sj1+sj2)-sij*s12)/(si1*sj1+si2*sj2);
  // sub-leading-colour terms (~si2) of Equations (57) and (58) in 1805.03752
  // the extra term subtracts the soft limit of the splitting function,
  // (essentially a subtraction of the LO A-term in the soft limit).
  double value = 2.*si2/(si1+si2) *
    ( (w_ij_12+wbar_ij_12)/2. * (s_CA/(2.*s_CF)-1. + 1.) - 1.);
  // sub-leading colour terms (~sij) of Equation (57) in 1805.03752
  // they are only invoked if the relative transverse momentum is large enough
  // and therefore the eikonal and alphaS are safe
  double kt2_1ij = (si1*sj1)/(sij+si1+sj1); 
  if (kt2_1ij>split.tcut()) {
    Splitting dummy(NULL,NULL,kt2_1ij,split.tcut());
    double as_ratio = ((*split.GetKernel()->GetGauge())(dummy)/
		       (*split.GetKernel()->GetGauge())(split));
    value += -2.*sij/(si1+sj1) * as_ratio * (w_ij_12+wbar_ij_12)/2. * (s_CA/(2.*s_CF)-1.);
  }
  return value;
}

/*
DECLARE_GETTER(FFV_FF,"FF_FFV_Soft",SF_Base,Kernel_Info);

SF_Base * ATOOLS::Getter<SF_Base,Kernel_Info,FFV_FF>::
operator()(const Parameter_Type & info) const
{
  if (info.Type()==kernel_type::FF &&
      info.LogType()==log_type::soft &&
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
*/



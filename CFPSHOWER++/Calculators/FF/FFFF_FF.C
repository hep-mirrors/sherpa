#include "CFPSHOWER++/Calculators/SF_Base.H"
#include "CFPSHOWER++/Calculators/FF/SF_FF13.H"
#include "CFPSHOWER++/Calculators/Functions.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

#include <vector>

namespace CFPSHOWER {
  class FFFF_FF : public SF_FF13 {
    ffff_mode::code m_mode;
    double m_Nc, m_CF, m_CA, m_TR;
    double m_zi, m_zj, m_zk;
    double m_sij, m_sijk, m_tijk, m_tikj;

    void SortFlavours();
    const double & B2(const Splitting & split) const;
    inline const double SymmetryFactor(const Splitting & split) const {
      double z = split.z(), xi = z/split.z2();
      switch (m_mode) {
      case ffff_mode::same: return (xi-z)/(1.-z);
      case ffff_mode::diff:
      default:
	break;
      }
      return 1.;
    }
  public:
    FFFF_FF(const Kernel_Info & info);
    double operator()(const Splitting & split) const;
    double Integral(const Splitting & split) const;
    double OverEstimate(const Splitting & split) const;
    void   GeneratePoint(Splitting & split) const;
  };
}

/*
  The logic in the 1->3 splittings is that they come in two modes:
  either fully differential (m_mode=0), or as subtraction of the end-point 
  only (m_mode=1).  The latter effectively being described by a 1->2 
  kinematics, where one of the two outgoing objects splits into two 
  entirely collinear partons.  Which of the two is picked is defined at 
  random, with equal probability, in "GeneratePoint".

  The splitting function consists of three parts (similar to MC@NLO,
  but without the B and V terms): a collinear subtraction I, a real
  part R, and a subtraction term S. 
 */


using namespace CFPSHOWER;
using namespace ATOOLS;

FFFF_FF::FFFF_FF(const Kernel_Info & info)  : SF_FF13(info) {
  SortFlavours();
  if (m_flavs[1].Kfcode()==m_split.Kfcode()) {
    m_mode     = ffff_mode::same;
    SetName("q->qqq");
  }
  else {
    m_mode = ffff_mode::diff;
    SetName("q->qQQ");
  }
  m_subtract = (m_split.Bar()==m_flavs[m_tags[0]] ?
		(m_mode==ffff_mode::same ?
		 subtract::both : subtract::soft) : 
		subtract::coll );
  
  msg_Out()<<"--> "<<METHOD<<"["<<m_mode<<", sub = "<<m_subtract<<"] "
	   <<"for "<<m_split.Bar()<<" -> ";
  for (size_t i=0;i<m_flavs.size();i++) msg_Out()<<m_flavs[i]<<" ";
  msg_Out()<<" { ";
  for (size_t i=0;i<m_tags.size();i++)
    msg_Out()<<m_flavs[m_tags[i]]<<" ("<<m_tags[i]<<") ";
  msg_Out()<<"}\n"
	   <<"   compared "<<m_split.Bar()<<" with "<<m_flavs[m_tags[0]]
	   <<" ("<<m_tags[0]<<")\n"
	   <<"----------------------------------------------------------\n"
	   <<"----------------------------------------------------------\n";
}

void FFFF_FF::SortFlavours() {
  for (size_t i=0;i<3;i++) {
    if (m_split==m_flavs[i]) {
      std::swap(m_flavs[i], m_flavs[0]);
      for (size_t j=0; j<m_tags.size();j++) {
	if      (m_tags[j]==i) m_tags[j]=1;
	else if (m_tags[j]==1) m_tags[j]=i;
      }
      break;
    }
  }
  if (m_flavs[1].IsAnti() && !m_flavs[2].IsAnti()) {
    std::swap(m_flavs[1], m_flavs[2]);
    for (size_t j=0; j<m_tags.size();j++) {
      if      (m_tags[j]==1) m_tags[j]=2;
      else if (m_tags[j]==2) m_tags[j]=1;
    }
  }
}

double FFFF_FF::operator()(const Splitting & split) const {
  double value = B2(split);
  //value *= -2.*log(z1)/denom;
  value *= (*split.GetKernel()->GetGauge())(split)/(2.*M_PI);
  value *= SymmetryFactor(split);
  return split.z() * value;
}

const double & FFFF_FF::B2(const Splitting & split) const {
  double term = 0;
  if (split.IsEndPoint()) {
    if (IsSubtract(subtract::coll)) {
      //      double xi   = split.z()/split.z2();
      //double sijk = split.t()/xi+split.s01()+split.m2(2), denom = 1.-sijk/split.t1();
      //double z0   = split.z(0)/denom, z1 = (xi-split.z(0))/denom, z2 = 1.-z0-z1;
      //term = (m_mode==ffff_mode::same ?
      //      Functions::DeltaI_qqp_F(z0,z1,z2) :
      //      Functions::DeltaI_qqbar_F(z0,z1,z2)); 
    }
  }
  else {
  }
}

double FFFF_FF::Integral(const Splitting & split) const {
  double invkappa = split.tcut()/split.Q2();
  double SFInt    = 20./9. * 1./2. * log(1.+invkappa);
  return SFInt * split.GetKernel()->GetGauge()->OverEstimate(split)/(2.*M_PI); 
}

double FFFF_FF::OverEstimate(const Splitting & split) const {
  double kappa = split.tcut()/split.Q2();
  double SFEst = 20./9. * 1./(2.*(split.z()+kappa));
  return SFEst * split.GetKernel()->GetGauge()->OverEstimate(split)/(2.*M_PI); 
}

void FFFF_FF::GeneratePoint(Splitting & split) const {
  double kappa = split.tcut()/split.Q2(), invkappa = 1./kappa;
  double z     = pow(1.+invkappa, -ran->Get()) * (1.+kappa) - kappa;
  double z2    = pow(z, ran->Get());
  bool   endp  = ran->Get()>0.5;
  if (endp) {
    double help  = ran->Get();
    double sai   = help/(1.-help) * (z2/z * split.t() + split.m2(1));
    split.Set_t2(sai);
  }
  else {
    split.Set_t2(0.);
  }
  split.Set_z(z);
  split.Set_z2(z2);
  split.Set_phi();
  split.Set_phi2();
}

DECLARE_GETTER(FFFF_FF,"FF_FFFF",SF_Base,Kernel_Info);

SF_Base * ATOOLS::Getter<SF_Base,Kernel_Info,FFFF_FF>::
operator()(const Parameter_Type & info) const
{
  if (info.Type()==kernel_type::FF &&
      info.GetSplit().IsFermion() && 
      info.GetFlavs()[0].IsFermion() &&
      info.GetFlavs()[1].IsFermion() &&
      info.GetFlavs()[2].IsFermion()) {
     return new FFFF_FF(info);
  }
  return NULL;
}

void ATOOLS::Getter<SF_Base,Kernel_Info,FFFF_FF>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"FFFF splitting Function (FF)";
}

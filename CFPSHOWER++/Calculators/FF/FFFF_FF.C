#include "CFPSHOWER++/Calculators/SF_Base.H"
#include "CFPSHOWER++/Calculators/FF/SF_FF13.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

#include <vector>

namespace CFPSHOWER {
  class FFFF_FF : public SF_FF13 {
    ffff_mode::code m_mode;
    subtract::code  m_subtract;
    double m_Nc, m_CF, m_CA, m_TR;
    double m_zi, m_zj, m_zk;
    double m_sij, m_sijk, m_tijk, m_tikj;

    void SortFlavours();
    double B2(const Splitting & split) const;
    double R_qqprime(const Splitting & split) const;
    double S_qqprime(const Splitting & split) const;
    double I(const double & zi, const double & zj,
	     const double & zk, const double & x) const;
    double i(const double & zi, const double & zj,
	     const double & zk, const double & x) const;
    double R_qqbar(const Splitting & split) const;
    double S_qqbar(const Splitting & split) const;

    inline double Pqq(const double & z) const {
      return (1.+ATOOLS::sqr(z))/(1.-z);
    } 
    inline double Pgq(const double & z) const {
      return ATOOLS::sqr(z)+ATOOLS::sqr(1.-z);
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
  m_subtract = (m_split.Bar()==m_flavs[m_tagsequence[0]] ?
		(m_mode==ffff_mode::same ?
		 subtract::both : subtract::soft) : 
		subtract::coll );
  
  msg_Out()<<"--> "<<METHOD<<"["<<m_mode<<", sub = "<<m_subtract<<"] "
	   <<"for "<<m_split.Bar()<<" -> ";
  for (size_t i=0;i<m_flavs.size();i++) msg_Out()<<m_flavs[i]<<" ";
  msg_Out()<<" { ";
  for (size_t i=0;i<m_tagsequence.size();i++)
    msg_Out()<<m_flavs[m_tagsequence[i]]<<" ("<<m_tagsequence[i]<<") ";
  msg_Out()<<"}\n"
	   <<"   compared "<<m_split.Bar()<<" with "<<m_flavs[m_tagsequence[0]]
	   <<" ("<<m_tagsequence[0]<<")\n"
	   <<"----------------------------------------------------------\n"
	   <<"----------------------------------------------------------\n";
}

void FFFF_FF::SortFlavours() {
  for (size_t i=0;i<3;i++) {
    if (m_split==m_flavs[i]) {
      std::swap(m_flavs[i], m_flavs[0]);
      for (size_t j=0; j<m_tagsequence.size();j++) {
	if      (m_tagsequence[j]==i) m_tagsequence[j]=1;
	else if (m_tagsequence[j]==1) m_tagsequence[j]=i;
      }
      break;
    }
  }
  if (m_flavs[1].IsAnti() && !m_flavs[2].IsAnti()) {
    std::swap(m_flavs[1], m_flavs[2]);
    for (size_t j=0; j<m_tagsequence.size();j++) {
      if      (m_tagsequence[j]==1) m_tagsequence[j]=2;
      else if (m_tagsequence[j]==2) m_tagsequence[j]=1;
    }
  }
}

double FFFF_FF::operator()(const Splitting & split) const {
  double value = B2(split);
}

double FFFF_FF::B2(const Splitting & split) const {
}

double FFFF_FF::R_qqprime(const Splitting & split) const {
  //double R = 1./2. * m_saij/m_sai * ( -sqr(m_taij)/(m_sai*m_saij) +
  //				    (4.*m_zj + sqr(m_za-m_zi))/(m_za+m_zi) +
  //				    (m_za + m_zi - m_sai/m_saij) );
  //return R;
}

double FFFF_FF::S_qqprime(const Splitting & split) const {
  // cosphi may be wrong - check notation.
  //double cosphi_ik_aj = CosPhi(split.Mom(0),split.Mom(1),split.Mom(2),split.Mom(3));
  //double S = m_saij/m_sai * ( (1.+m_zj)/(1.-m_zj) * (1.-(2.*m_za*m_zi)/sqr(m_za+m_zi)) +
  //			      4.*m_za*m_zi*m_zj/pow(1.-m_zj,3.) * (1.-2.*sqr(cosphi_ik_aj)));
  //return S;
}

double FFFF_FF::R_qqbar(const Splitting & split) const {
  //double R = ( 1./2. * m_saij/m_sai * (-sqr(m_taij)/(m_sai*m_saij) +
  //				       4.*m_zj+sqr(m_za+m_zi)/(m_za+m_zi) +
  //				       (m_za+m_zi - m_sai/m_saij) ) -
  //	       1./m_Nc * m_saij/m_sai * (2.*m_sij/m_saij +
  //					 (1.+sqr(m_za))/(1.-m_zi) -
  //					 2.*m_zi/(1.-m_zj) -
  //					 ( m_saij/m_saj * m_za/2. *
  //					   (1.+sqr(m_za))/((1.-m_zi)*(1.-m_zj)) ) ) +
  //	       1./2. * m_saij/m_saj * (-sqr(m_taji)/(m_saj*m_saij) +
  //				       4.*m_zi+sqr(m_za+m_zj)/(m_za+m_zj) +
  //				       (m_za+m_zj - m_saj/m_saij) ) -
  //	       1./m_Nc * m_saij/m_saj * (2.*m_sij/m_saij +
  //					 (1.+sqr(m_za))/(1.-m_zj) -
  //					 2.*m_zj/(1.-m_zi) -
  //					 ( m_saij/m_sai * m_za/2. *
  //					   (1.+sqr(m_za))/((1.-m_zj)*(1.-m_zi)) ) ) );
}

double FFFF_FF::S_qqbar(const Splitting & split) const {
  //double cosphi_ik_aj = CosPhi(split.Mom(0),split.Mom(1),split.Mom(2),split.Mom(3));
  //double cosphi_jk_ai = CosPhi(split.Mom(0),split.Mom(1),split.Mom(2),split.Mom(3));
  //double S = (m_saij/m_sai * ( (1.+m_zj)/(1.-m_zj) * (1.-(2.*m_za*m_zi)/sqr(m_za+m_zi)) +
  //			       4.*m_za*m_zi*m_zj/pow(1.-m_zj,3.) * (1.-2.*sqr(cosphi_ik_aj))) +
  //	      m_saij/m_saj * ( (1.+m_zi)/(1.-m_zi) * (1.-(2.*m_za*m_zj)/sqr(m_za+m_zj)) +
  //			       4.*m_za*m_zj*m_zi/pow(1.-m_zi,3.) * (1.-2.*sqr(cosphi_jk_ai))) );
  //return S;
}

double FFFF_FF::I(const double & zi, const double & zj,
		  const double & zk, const double & x) const {
  //double P = Pqq(zk);
  //double I = ( P +
  //	       (1. - (2.*zi*zj)/sqr(zi+zj)) *
  //	       (1. - zk + P) * (log(x*zj*zk) - 1.) );
  //return I;
}

double FFFF_FF::i(const double & zi, const double & zj,
		  const double & zk, const double & x) const {
  //double i   = 2.*( Pqq(zk) * log(x*zk) + (1.-zk) ) * Pgq(zi/(zi+zj)));
  //return i;
}

double FFFF_FF::Integral(const Splitting & split) const {
  
}

double FFFF_FF::OverEstimate(const Splitting & split) const {

}

void FFFF_FF::GeneratePoint(Splitting & split) const {
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

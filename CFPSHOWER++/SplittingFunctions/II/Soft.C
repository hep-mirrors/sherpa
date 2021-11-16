#include "CFPSHOWER++/SplittingFunctions/SF_Base.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class II_Soft : public SF_Base {
    bool           m_is_g2gg, m_swapped;
    kin_type::code m_kin;
    double A1(const Splitting & split) const;
  public:
    II_Soft(const Kernel_Info & info);
    double operator()(const Splitting & split);
    double Integral(const Splitting & split) const;
    double OverEstimate(const Splitting & split) const;
    void   GeneratePoint(Splitting & split) const;
  };
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

using namespace CFPSHOWER;
using namespace ATOOLS;

II_Soft::II_Soft(const Kernel_Info & info)  :
  SF_Base(info),
  m_is_g2gg(info.GetSplit().IsVector()), 
  m_swapped(info.TagSequence()[0]!=0),
  m_kin(info.KinType())
{
  m_name = std::string("II: Soft");
}

double II_Soft::operator()(const Splitting & split) {
  double Kfactor  = ( (m_CMW==1) ?
		      (1.+split.GetKernel()->GetGauge()->K(split)) : 1. );
  double idfactor = ( m_is_g2gg ?
		      ( m_swapped ? 1.-split.Z() : split.Z() ) : 1. );
  double A = A1(split);
  return Kfactor * A * idfactor;
}

double II_Soft::Integral(const Splitting & split) const {
  double Kmax = ((m_CMW==1.) ?
		 (1.+split.GetKernel()->GetGauge()->KMax(split)) : 1.);
  double I    =  2. * split.IntPhi() * log((1.-split.Zmin())/(1.-split.Zmax()));
  return I * Kmax;
}

double II_Soft::OverEstimate(const Splitting & split) const {
  double Kmax = ((m_CMW==1.) ?
		 (1.+split.GetKernel()->GetGauge()->KMax(split)) : 1.);
  double E    = 2. * split.IntPhi() / (1.-split.Z()); 
  return E * Kmax;
}

void II_Soft::GeneratePoint(Splitting & split) const {
  split.SetZ(1. - (1.-split.Zmin()) *
	     pow((1.-split.Zmax())/(1.-split.Zmin()), ran->Get()));
  split.SetPhi(2.*M_PI*ran->Get());
}

double II_Soft::A1(const Splitting & split) const {
  double pipj = split.Mom(0)*split.Mom(1);
  double pipk = split.Mom(0)*split.Mom(2);
  double pjpk = split.Mom(1)*split.Mom(2);
  if (m_kin==kin_type::Alaric) {
    double pin  = split.Mom(0)*split.GetKinSpect();
    double pjn  = split.Mom(1)*split.GetKinSpect();
    double pkn  = split.Mom(2)*split.GetKinSpect();
    return 2.*pipk*pin/(pipj*pkn+pjpk*pin);
  }
  else if (m_kin==kin_type::CataniSeymour)
    return 2.*pipk/(pipj+pjpk);
  return 0.;
}

DECLARE_GETTER(II_Soft,"II_Soft",SF_Base,Kernel_Info);

SF_Base * ATOOLS::Getter<SF_Base,Kernel_Info,II_Soft>::
operator()(const Parameter_Type & info) const
{
  if (info.Type()==kernel_type::II &&
      info.LogType()==log_type::soft &&
      info.GetFlavs().size()==2) {
    if ((info.GetSplit().IsFermion() &&
	 info.GetFlavs()[0].IsFermion() && info.TagSequence()[0]==0 &&
	 info.GetFlavs()[1].IsVector() && info.TagSequence()[1]==1 &&
	 int(info.SFType() & 1)>0 ) ||
	(info.GetSplit().IsVector() &&
	 info.GetFlavs()[0].IsVector() && info.GetFlavs()[1].IsVector() &&
	 int(info.SFType() & 2)>0 ) )
    return new II_Soft(info);
  }
  return NULL;
}

void ATOOLS::Getter<SF_Base,Kernel_Info,II_Soft>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Soft splitting function (II)";
}

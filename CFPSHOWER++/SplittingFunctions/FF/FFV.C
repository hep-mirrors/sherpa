#include "CFPSHOWER++/SplittingFunctions/SF_Base.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class FF_FFV_Soft : public SF_Base {
    double A1(const Splitting & split) const;
    double B1(const Splitting & split) const;
  public:
    FF_FFV_Soft(const Kernel_Info & info);
    double operator()(const Splitting & split);
    double Integral(const Splitting & split) const;
    double OverEstimate(const Splitting & split) const;
    void   GeneratePoint(Splitting & split) const;
  };

  class FF_FFV_Coll : public SF_Base {
    double A1(const Splitting & split) const;
    double B1(const Splitting & split) const;
  public:
    FF_FFV_Coll(const Kernel_Info & info);
    double operator()(const Splitting & split);
    double Integral(const Splitting & split) const;
    double OverEstimate(const Splitting & split) const;
    void   GeneratePoint(Splitting & split) const;
  };
}

using namespace CFPSHOWER;
using namespace ATOOLS;

FF_FFV_Soft::FF_FFV_Soft(const Kernel_Info & info)  : SF_Base(info) {
  SetName("FF: F->FV (soft)");
}

double FF_FFV_Soft::operator()(const Splitting & split) {
  double Kfactor = (m_CMW==1) ? (1.+split.GetKernel()->GetGauge()->K(split)) : 1.;
  //msg_Out()<<METHOD<<" A(z = "<<split.Z()<<", y = "<<split.Y()<<") = "<<A1(split)<<" / "
  //	   <<"over = "<<OverEstimate(split)<<" for "
  //	   <<"[kin = "<<split.GetKernel()->GetKinematics()->Scheme()<<"]\n";
  return Kfactor * A1(split) + B1(split);
}

double FF_FFV_Soft::Integral(const Splitting & split) const {
  double Kmax = (m_CMW==1.) ? (1.+split.GetKernel()->GetGauge()->KMax(split)) : 1., I = 0;
  switch (split.GetKernel()->GetKinematics()->Scheme()) {
  case kin_type::CS:
    I = 2.*log((1.-split.Zmin())/(1.-split.Zmax()));
    break;
  case kin_type::PanGlobal:
    I = 2.*log(split.Q2()/split.Tcut());
    break;
  default:
    break;
  }
  return I * Kmax;
}

double FF_FFV_Soft::OverEstimate(const Splitting & split) const {
  double Kmax = (m_CMW==1.) ? (1.+split.GetKernel()->GetGauge()->KMax(split)) : 1., E = 1.e12;
  switch (split.GetKernel()->GetKinematics()->Scheme()) {
  case kin_type::CS:
    E = 2./(1.-split.Z());
    break;
  case kin_type::PanGlobal:
    E = 2.;
    break;
  default:
    break;
  }
  return E * Kmax;
}

void FF_FFV_Soft::GeneratePoint(Splitting & split) const {
  switch (split.GetKernel()->GetKinematics()->Scheme()) {
  case kin_type::CS:
    split.SetZ(1. - (1.-split.Zmin()) * pow((1.-split.Zmax())/(1.-split.Zmin()), ran->Get()));
    break;
  case kin_type::PanGlobal:
    split.SetEta((2.*ran->Get()-1.)/2. * log(split.Q2()/split.Tcut()));
    break;
  default:
    break;
  }
  split.SetPhi(2.*M_PI*ran->Get());
}

double FF_FFV_Soft::A1(const Splitting & split) const {
  double pipj = split.Mom(0)*split.Mom(1), pjpk = split.Mom(1)*split.Mom(2);
  double A1   = 2.*pjpk/(pipj+pjpk);
  switch (split.GetKernel()->GetKinematics()->Scheme()) {
  case kin_type::CS:
    A1 = 2./(1.-split.Z()*(1.-split.Y())) - 2.; 
    break;
  case kin_type::PanGlobal:
  default:
    break;
  }
  return A1;
}

double FF_FFV_Soft::B1(const Splitting & split) const { return 0.; }


DECLARE_GETTER(FF_FFV_Soft,"FF_FFV_Soft",SF_Base,Kernel_Info);

SF_Base * ATOOLS::Getter<SF_Base,Kernel_Info,FF_FFV_Soft>::
operator()(const Parameter_Type & info) const
{
  if (info.Type()==kernel_type::FF &&
      info.LogType()==log_type::soft &&
      info.GetFlavs().size()==2 &&
      info.GetSplit().IsFermion() && 
      info.GetFlavs()[0].IsFermion() &&
      info.TagSequence()[0]==0 &&
      info.GetFlavs()[1].IsVector() &&
      info.TagSequence()[1]==1) {
    return new FF_FFV_Soft(info);
  }
  return NULL;
}

void ATOOLS::Getter<SF_Base,Kernel_Info,FF_FFV_Soft>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"FFV splitting function: soft part (FF)";
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FF_FFV_Coll::FF_FFV_Coll(const Kernel_Info & info)  : SF_Base(info) {
  SetName("FF: F->FV (coll)");
}

double FF_FFV_Coll::operator()(const Splitting & split) { return B1(split); }

double FF_FFV_Coll::Integral(const Splitting & split) const {
  double I = 1.;
  switch (split.GetKernel()->GetKinematics()->Scheme()) {
  case kin_type::CS:
    break;
  case kin_type::PanGlobal:
    I = 0.;
    break;
  default:
    break;
  }
  return I;
}

double FF_FFV_Coll::OverEstimate(const Splitting & split) const { return 1.; }

void FF_FFV_Coll::GeneratePoint(Splitting & split) const {
  switch (split.GetKernel()->GetKinematics()->Scheme()) {
  case kin_type::CS:
    split.SetZ(ran->Get());
    break;
  case kin_type::PanGlobal:
  default:
    msg_Error()<<"Error in "<<METHOD<<": wrong kinematics scheme for collinear part: "
	       <<split.GetKernel()->GetKinematics()->Scheme()<<"\n";
    exit(1);
    break;
  }
  split.SetPhi(2.*M_PI*ran->Get());
}

double FF_FFV_Coll::A1(const Splitting & split) const { return 0.; }

double FF_FFV_Coll::B1(const Splitting & split) const { return 1.-split.Z(); }


DECLARE_GETTER(FF_FFV_Coll,"FF_FFV_Coll",SF_Base,Kernel_Info);

SF_Base * ATOOLS::Getter<SF_Base,Kernel_Info,FF_FFV_Coll>::
operator()(const Parameter_Type & info) const
{
  if (info.Type()==kernel_type::FF &&
      info.LogType()==log_type::coll &&
      info.GetFlavs().size()==2 &&
      info.GetSplit().IsFermion() && 
      info.GetFlavs()[0].IsFermion() &&
      info.TagSequence()[0]==0 &&
      info.GetFlavs()[1].IsVector() &&
      info.TagSequence()[1]==1) {
    return new FF_FFV_Coll(info);
  }
  return NULL;
}

void ATOOLS::Getter<SF_Base,Kernel_Info,FF_FFV_Coll>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"FFV splitting function: coll part (FF)";
}



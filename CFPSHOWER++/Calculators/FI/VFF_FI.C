#include "CFPSHOWER++/Calculators/FI/SF_FI12.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class VFF_FI : public SF_FI12 {
  private:
    double m_jmax;
    
    double B1(const double & z,const double & kappa2) const;
  public:
    VFF_FI(const Kernel_Info & info);
    double operator()(const Splitting & split);
    double Integral(const Splitting & split) const;
    double OverEstimate(const Splitting & split) const;
    void   GeneratePoint(Splitting & split) const;
  };
}


using namespace CFPSHOWER;
using namespace ATOOLS;

VFF_FI::VFF_FI(const Kernel_Info & info) : SF_FI12(info), m_jmax(5.) { 
  SetName("V->FF");
}

double VFF_FI::operator()(const Splitting & split) {
  double z(split.z(0)), kappa2(split.tcut()/split.Q2red());
  // TODO: have to add ME correction for DIS
  if (split.mspect2()>1.e-12) {
    msg_Error()<<"Error in "<<METHOD<<": did not expect massive spectator in IS.\n"
	       <<"   Will exit the run.\n";
    exit(1);
  }
  double value = B1(z,kappa2);
  if (!split.IsMassive()) {
    double mui2(split.m2(0)/split.Q2red()), y = split.y();
    value += 2.*y*mui2/((1.-y)+2.*y*mui2);
  }
  if (split.Clustered()==0) value *= split.ztilde(m_tags[0]);
  return value;
}

double VFF_FI::Integral(const Splitting & split)     const { return m_jmax; }

double VFF_FI::OverEstimate(const Splitting & split) const { return m_jmax; }

void VFF_FI::GeneratePoint(Splitting & split) const {
  split.Set_z(0,ran->Get());
  split.Set_phi(0);
}

double VFF_FI::B1(const double & z,const double & kappa2) const {
  return sqr(z)+sqr(1.-z);
}

DECLARE_GETTER(VFF_FI,"FI_VFF",SF_Base,Kernel_Info);

SF_Base * ATOOLS::Getter<SF_Base,Kernel_Info,VFF_FI>::
operator()(const Parameter_Type & info) const
{
  return NULL;
  if (info.Type()==kernel_type::FI &&
      (info.GetSplit().IsVector() &&
       info.GetFlavs()[0].IsFermion() &&
       info.GetFlavs()[0].IsAnti() && 
       info.GetFlavs()[1].IsFermion())) {
    return new VFF_FI(info);
  }
  return NULL;
}

void ATOOLS::Getter<SF_Base,Kernel_Info,VFF_FI>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"VFF Splitting Function (FI)";
}

#include "CFPSHOWER++/Calculators/FI/SF_FI.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

//#include "CFPSHOWER++/Calculators/GQQ.C"

namespace CFPSHOWER {
  class VFF_FI : public SF_FI {
  private:
    double m_jmax;
    
    double B1(const double & z,const double & kappa2) const;
    double B2(const double & z,const double & kappa2) const;
  public:
    VFF_FI(const Kernel_Info & info);
    double operator()(const Splitting & split) const;
    double Integral(const Splitting & split) const;
    double OverEstimate(const Splitting & split) const;
    //void GeneratePoint(Splitting & split) const; // new
  };
}


using namespace CFPSHOWER;
using namespace ATOOLS;

VFF_FI::VFF_FI(const Kernel_Info & info) :
  SF_FI(info), m_jmax(5.) { 
  SetName("V->FF");
}

double VFF_FI::operator()(const Splitting & split) const {
  double mij2(split.mij2()), mi2(split.mi2()), mj2(split.mj2()), mk2(split.mk2());
  double z(split.Z()), y(split.Y()), Q2(split.Q2()), kappa2(split.T()/Q2);
  // TODO: have to add ME correction for DIS
  double value = 0.;
  if (m_orderB>0) {
    if (mk2>1.e-12) {
      msg_Error()<<"Error in "<<METHOD<<": did not expect massive spectator in IS.\n"
		 <<"   Will exit the run.\n";
      exit(1);
    }
    value += B1(z,kappa2);
    if (!(mi2==0. && mj2==0.)) {
      double mui2(mi2/Q2);
      value += 2.*y*mui2/((1.-y)+2.*y*mui2);
    }
  }
  if (split.Clustered()==0) value *= m_swap?(1.-z):z;
  return value;
}

double VFF_FI::Integral(const Splitting & split)     const { return m_jmax; }
double VFF_FI::OverEstimate(const Splitting & split) const { return m_jmax; }

double VFF_FI::B1(const double & z,const double & kappa2) const {
  return sqr(z)+sqr(1.-z);
}

double VFF_FI::B2(const double & z,const double & kappa2) const {
  double b2 = 0.0;
  return b2;
}


DECLARE_GETTER(VFF_FI,"FI_VFF",SF_Base,Kernel_Info);

SF_Base * ATOOLS::Getter<SF_Base,Kernel_Info,VFF_FI>::
operator()(const Parameter_Type & info) const
{
  return NULL;
  if (info.Type()==kernel_type::FI &&
      (info.GetFlavs()[0].IsVector() &&
       info.GetFlavs()[1].IsFermion() &&
       info.GetFlavs()[1].IsAnti() && 
       info.GetFlavs()[2].IsFermion())) {
    return new VFF_FI(info);
  }
  return NULL;
}

void ATOOLS::Getter<SF_Base,Kernel_Info,VFF_FI>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"VFF Splitting Function (FI)";
}

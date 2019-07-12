#include "CFPSHOWER++/Calculators/IF/SF_IF.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class FFV_IF : public SF_IF {
    // somewhat tricky nomenclature to keep track off.
    // incoming fermion F (fl[0], mass mij), splitting into an incoming fermion 
    // F (fl[1], mass mi) and an outgoing vector boson V (fl[2], mass mj);
    // spectator has mass mk.
  private:
    double m_jmax;
    
    double A1(const double & z,const double & kappa2) const;
    double B1(const double & z,const double & kappa2) const;
    double B2(const double & z,const double & kappa2) const;
  public:
    FFV_IF(const Kernel_Info & info);
    double operator()(const Splitting & split) const;
    double Integral(const Splitting & split) const;
    double OverEstimate(const Splitting & split) const;
    void   GeneratePoint(Splitting & split) const;
  };
}


using namespace CFPSHOWER;
using namespace ATOOLS;


FFV_IF::FFV_IF(const Kernel_Info & info) :
  SF_IF(info),
  // over-estimator: 5 for potential valence quarks, 2 for all others.
  // TODO: may need to monitor this.
  m_jmax(info.GetFlavs()[0].Kfcode() < 3 ? 5.:2.)
{
  SetName("F->FV");
}

double FFV_IF::operator()(const Splitting & split) const {
  double mij2(split.mij2()), mi2(split.mi2()), mj2(split.mj2()), mk2(split.mk2());
  double z(split.Z()), y(split.Y()), Q2(split.Q2()), kappa2(split.T()/Q2);
  // Start with the soft term only, including possible K factors
  // (cusp anomalous dimensions), obtained from the gauge part of the kernel
  double hofactor = 1.+split.GetKernel()->GetGauge()->K(split);
  double value    = A1(z,kappa2) * hofactor;
  // All massless: just add the collinear parts.
  // TODO: Add the B2 parts as soon as the LL shower is validated.
  // TODO: Add the DIS ME correction
  if (m_orderB>0) {
    value += B1(z,kappa2);
    if (mi2!=0. || mj2!=0.) {
      msg_Debugging()<<METHOD<<" without mass corrections yet.\n";
    }
  }
  return value;
}

double FFV_IF::Integral(const Splitting & split) const {
  double homax = (1.+split.GetKernel()->GetGauge()->KMax(split));
  return log(1.+sqr(1.-split.Eta())*split.Q2()/split.T0()) * homax * m_jmax;
}

double FFV_IF::OverEstimate(const Splitting & split) const {
  double homax = (1.+split.GetKernel()->GetGauge()->KMax(split));
  return A1(split.Z(),split.T0()/split.Q2()) * homax * m_jmax;
}

void FFV_IF::GeneratePoint(Splitting & split) const {
  double kappa2 = split.T0()/split.Q2();
  double z      = 1.-sqrt(kappa2 * (pow(1.+sqr(1.-split.Eta())/kappa2,ran->Get())-1.)); 
  split.SetZ(z);
  split.Setphi();
}

double FFV_IF::A1(const double & z,const double & kappa2) const {
  return 2.*(1.-z)/(sqr(1.-z)+kappa2);
}

double FFV_IF::B1(const double & z,const double & kappa2) const {
  return -(1.+z);
}

double FFV_IF::B2(const double & z,const double & kappa2) const {
  double b2 = 0.0;
  return b2;
}

DECLARE_GETTER(FFV_IF,"IF_FFV",SF_Base,Kernel_Info);

SF_Base * ATOOLS::Getter<SF_Base,Kernel_Info,FFV_IF>::
operator()(const Parameter_Type & info) const
{
  if (info.Type()==kernel_type::IF &&
      info.GetFlavs()[0].IsFermion() && 
      info.GetFlavs()[1].IsFermion() &&
      info.GetFlavs()[2].IsVector()) {
    return new FFV_IF(info);
  }
  return NULL;
}

void ATOOLS::Getter<SF_Base,Kernel_Info,FFV_IF>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"FFV splitting Function (IF)";
}

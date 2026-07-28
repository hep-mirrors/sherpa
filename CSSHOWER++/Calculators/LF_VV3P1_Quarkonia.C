#include "CSSHOWER++/Showers/Splitting_Function_Base.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Phys/LDME.H"
#include "MODEL/Main/Single_Vertex.H"

namespace CSSHOWER {

class LF_VV3P1_Quarkonia_FF : public SF_Lorentz {
public:
  inline LF_VV3P1_Quarkonia_FF(const SF_Key &key) : SF_Lorentz(key) {}

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};
} // namespace CSSHOWER

#include "ATOOLS/Math/Random.H"
#include "MODEL/Main/Single_Vertex.H"

using namespace CSSHOWER;
using namespace ATOOLS;

double LF_VV3P1_Quarkonia_FF::operator()(
    const double z, const double y, const double eta, const double scale,
    const double Q2) { 
    const double m = m_flavs[2].IsC_Hadron() ? ATOOLS::Flavour(kf_c).Mass(true) : 0;
    const double mj  = m_flavs[2].Mass();
    const double mk  = ATOOLS::Flavour(m_flspec.Kfcode()).Mass(true);
    const double sij = scale/z/(1-z) + sqr(mj)/(1-z);
    const double M = 2*m;
    const double muj2 = sqr(mj)/Q2, muk2 = sqr(mk)/Q2;
    double value = 6*sqr(sij);
    value *= ( sqr(sij - sqr(M)) - 2*z*((1-z)*sij - sqr(M))*(sij - 2*sqr(M)) );
    value /= sqr(sqr(( sij - sqr(M))));
    // msg_Out() << METHOD << " s=" << sij << ", " << (scale/z/(1-z) + sqr(mj)/1-z) << std::endl;
    // msg_Out() << METHOD << " s=" << sij << ", M=" << M << ", m= " << M/2 << ", z=" << z << ", value=" << value << std::endl;
    value *= sqr(m)/sij;
    value *= 1. / ( (1 - muj2 - muk2) + 1./ y * ( muj2 ) );
    value *= 1. / (1 + sqr(z) * sqr(mj) / scale);
    const double prefactor = 32./243*GetLDME(m_flavs[2].Kfcode())/pow(m,5);
    // the correct normalisation according to Phys. Rev. D 50,5 (1994) p. 3176 is 64/243, but there is a floating factor 2 elsewhere in the code?
    return prefactor * sqr(p_cf->Coupling(scale, 0)) * value * JFF(y, 0, muj2, muk2, 0);
}

double LF_VV3P1_Quarkonia_FF::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  m_zmin = std::max(zmin, 0.15); //safety cut to catch singularity at z=0
  m_zmax = std::min(zmax, 1.);
  const double m = m_flavs[2].IsC_Hadron() ? ATOOLS::Flavour(kf_c).Mass(true) : 0;
  const double prefactor = 64./243*GetLDME(m_flavs[2].Kfcode())/pow(m,5);
  return prefactor * sqr(0.25 +0*p_cf->MaxCoupling(0)) * 123./4 * (zmax - zmin);
}

double LF_VV3P1_Quarkonia_FF::OverEstimated(const double z, const double y) {
  const double m = m_flavs[2].IsC_Hadron() ? ATOOLS::Flavour(kf_c).Mass(true) : 0;
  const double prefactor = 64./243*GetLDME(m_flavs[2].Kfcode())/pow(m,5);
  return prefactor * sqr(p_cf->MaxCoupling(0)) * 123/4.;
}

double LF_VV3P1_Quarkonia_FF::Z() {
  return m_zmin + (m_zmax - m_zmin) * ran->Get();
}


DECLARE_GETTER(LF_VV3P1_Quarkonia_FF, "VV3P1_Quarkonia", SF_Lorentz, SF_Key);

SF_Lorentz *ATOOLS::Getter<SF_Lorentz, SF_Key, LF_VV3P1_Quarkonia_FF>::operator()(
    const Parameter_Type &args) const {
  if (args.m_col < 0)
    return NULL;
  if ((args.m_mode == 0 && args.p_v->in[0].IntSpin() == 2 &&
       args.p_v->in[1].IntSpin() == 2 && args.p_v->in[2].IntSpin() == 2) ||
      (args.m_mode == 1 && args.p_v->in[0].IntSpin() == 2 &&
       args.p_v->in[1].IntSpin() == 2 && args.p_v->in[2].IntSpin() == 2)) {
    switch (args.m_type) {
    case cstp::FF:  return new LF_VV3P1_Quarkonia_FF(args);
    // case cstp::FI: return new LF_FFV_FI(args);
    // case cstp::IF: return new LF_FFV_IF(args);
    // case cstp::II: return new LF_FFV_II(args);
    case cstp::none:
      break;
    }
  }
  return NULL;
}

void ATOOLS::Getter<SF_Lorentz, SF_Key, LF_VV3P1_Quarkonia_FF>::PrintInfo(
    std::ostream &str, const size_t width) const {
  str << "VV3P1 lorentz functions";
}
#include "CSSHOWER++/Showers/Splitting_Function_Base.H"
#include "ATOOLS/Math/Random.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Phys/LDME.H"

namespace CSSHOWER {

class LF_VV1S0_Quarkonia_FF : public SF_Lorentz {
public:
  inline LF_VV1S0_Quarkonia_FF(const SF_Key &key) : SF_Lorentz(key) {}

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

double LF_VV1S0_Quarkonia_FF::operator()(
    const double z, const double y, const double eta, const double scale,
    const double Q2) { 
    const double mj  =  m_flavs[2].Mass();
    const double mi  =  m_flavs[0].Mass();
    const double mk  =  m_flspec.Mass();
    const double muj2 = sqr(m_flavs[2].Mass(true))/Q2;
    const double muk2 = sqr(m_flspec.Mass(true))/Q2;
    const double sij = y * (Q2 - sqr(mk)) + (1-y)*(sqr(mi) + sqr(mj));
    const double M = 2*Flavour(kf_c).Mass(true);
    double value = 0;
    value  = sqr(sij) + sqr(sqr(M)) - 2 * (1 - z) * (sij + sqr(M))*sij + 2*sqr(sij*(1-z));
    value /= sij*sqr(sij-sqr(M));
    value *= 1. / ( (1 - muj2 - muk2) + 1./ y * ( muj2 ) );
    value *= 1. / (1 + sqr(z) * sqr(mj) / scale);
    return sqr(p_cf->Coupling(scale, 0)) * value * JFF(y, 0, muj2, muk2, 0);
}

double LF_VV1S0_Quarkonia_FF::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  m_zmin = std::max(zmin,0.15);
  m_zmax = std::min(zmax,1.);
  const double m = Flavour(kf_c).Mass(true);
  const double prefactor = 4./27*GetLDME(m_flavs[2].Kfcode())/m;
  return prefactor * sqr(p_cf->MaxCoupling(0)) * 1./sqr(2*m) * (zmax - zmin);
}

double LF_VV1S0_Quarkonia_FF::OverEstimated(const double z, const double y) {
  const double m = Flavour(kf_c).Mass(true);
  return sqr(p_cf->MaxCoupling(0)) /sqr(2*m);
}

double LF_VV1S0_Quarkonia_FF::Z() {
  return m_zmin + (m_zmax - m_zmin) * ran->Get();
}


DECLARE_GETTER(LF_VV1S0_Quarkonia_FF, "VV1S0_Quarkonia", SF_Lorentz, SF_Key);

SF_Lorentz *ATOOLS::Getter<SF_Lorentz, SF_Key, LF_VV1S0_Quarkonia_FF>::operator()(
    const Parameter_Type &args) const {
  if (args.m_col < 0)
    return NULL;
  if ((args.m_mode == 0 && args.p_v->in[0].IntSpin() == 2 &&
       args.p_v->in[1].IntSpin() == 2 && args.p_v->in[2].IntSpin() == 0) ||
      (args.m_mode == 1 && args.p_v->in[0].IntSpin() == 2 &&
       args.p_v->in[1].IntSpin() == 0 && args.p_v->in[2].IntSpin() == 2)) {
    switch (args.m_type) {
    case cstp::FF:  return new LF_VV1S0_Quarkonia_FF(args);
    // case cstp::FI: return new LF_FFV_FI(args);
    // case cstp::IF: return new LF_FFV_IF(args);
    // case cstp::II: return new LF_FFV_II(args);
    case cstp::none:
      break;
    }
  }
  return NULL;
}

void ATOOLS::Getter<SF_Lorentz, SF_Key, LF_VV1S0_Quarkonia_FF>::PrintInfo(
    std::ostream &str, const size_t width) const {
  str << "VV1S0 lorentz functions";
}
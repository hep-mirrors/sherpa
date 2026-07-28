#include "CSSHOWER++/Showers/Splitting_Function_Base.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Phys/LDME.H"
#include "MODEL/Main/Single_Vertex.H"

namespace CSSHOWER {

class LF_VV3P0_Quarkonia_FF : public SF_Lorentz {
public:
  inline LF_VV3P0_Quarkonia_FF(const SF_Key &key) : SF_Lorentz(key) {}

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

double LF_VV3P0_Quarkonia_FF::operator()(
    const double z, const double y, const double eta, const double scale,
    const double Q2) { 
    const double m = m_flavs[2].IsC_Hadron() ? ATOOLS::Flavour(kf_c).Mass(true) : 0;
    const double mi = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true); 
    const double mj  = m_flavs[2].Mass();
    const double mk  = ATOOLS::Flavour(m_flspec.Kfcode()).Mass(true);
    const double sij = y * (Q2 - sqr(mk)) + (1-y)*(sqr(mi) + sqr(mj)); //scale/z/(1-z) + sqr(mj)/(1-z);
    const double M = 2*m;
    const double muj2 = sqr(mj)/Q2, muk2 = sqr(mk)/Q2;
    double value = sqr(sij - 3*sqr(M));
    value *= ( sqr(sij - sqr(M)) - 2*z*sij*((1-z)*sij - sqr(M)) );
    value /= sqr(sqr(( sij - sqr(M))));
      
    // msg_Out() << METHOD << " s=" << sij << ", " << scale/z/(1-z)+sqr(mj)/(1-z) << std::endl;
    // msg_Out() << METHOD << " t=" << scale << ", M=" << M << ", m= " << M/2 << ", z=" << z << ", mj= " << mj << ", value=" << value << std::endl;
    value *= sqr(M/2)/sij;
    value *= 1. / ( (1 - muj2 - muk2) + 1./ y * ( muj2 ) );
    value *= 1. / (1 + sqr(z) * sqr(mj) / scale);
    if (sqr(p_cf->Coupling(scale, 0)) * value * JFF(y, 0, muj2, muk2, 0) > sqr(p_cf->MaxCoupling(0)) * 0.2) {
          msg_Out() << METHOD << " too high " << std::endl;
    }
    const double prefactor = 64./243*GetLDME(m_flavs[2].Kfcode())/pow(m,5);
    return prefactor * sqr(p_cf->Coupling(scale, 0)) * value * JFF(y, 0, muj2, muk2, 0);
}

double LF_VV3P0_Quarkonia_FF::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  m_zmin = std::max(zmin, 0.); //safety cut to catch singularity at z=0
  m_zmax = std::min(zmax, 1.);
  const double m = m_flavs[2].IsC_Hadron() ? ATOOLS::Flavour(kf_c).Mass(true) : 0;
  const double prefactor = 64./243*GetLDME(m_flavs[2].Kfcode())/pow(m,5);
  return prefactor * sqr(p_cf->MaxCoupling(0)) * 0.2 * (m_zmax - m_zmin);  ;
}

double LF_VV3P0_Quarkonia_FF::OverEstimated(const double z, const double y) {
  const double mi  = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double mj  = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  const double mk  = ATOOLS::Flavour(m_flspec.Kfcode());
  const double m   = m_flavs[2].IsC_Hadron() ? ATOOLS::Flavour(kf_c).Mass(true) : 0;
  const double prefactor = 64./243*GetLDME(m_flavs[2].Kfcode())/pow(m,5);
  return prefactor * sqr(p_cf->MaxCoupling(0)) * 0.2;
}

double LF_VV3P0_Quarkonia_FF::Z() {
  // return 1. - (1. - m_zmin) * pow((1. - m_zmax) / (1. - m_zmin), ATOOLS::ran->Get());
  return m_zmin + (m_zmax - m_zmin) * ATOOLS::ran->Get();
}


DECLARE_GETTER(LF_VV3P0_Quarkonia_FF, "VV3P0_Quarkonia", SF_Lorentz, SF_Key);

SF_Lorentz *ATOOLS::Getter<SF_Lorentz, SF_Key, LF_VV3P0_Quarkonia_FF>::operator()(
    const Parameter_Type &args) const {
  if (args.m_col < 0)
    return NULL;
  if ((args.m_mode == 0 && args.p_v->in[0].IntSpin() == 2 &&
       args.p_v->in[1].IntSpin() == 2 && args.p_v->in[2].IntSpin() == 0) ||
      (args.m_mode == 1 && args.p_v->in[0].IntSpin() == 2 &&
       args.p_v->in[1].IntSpin() == 0 && args.p_v->in[2].IntSpin() == 2)) {
    switch (args.m_type) {
    case cstp::FF:  return new LF_VV3P0_Quarkonia_FF(args);
    // case cstp::FI: return new LF_FFV_FI(args);
    // case cstp::IF: return new LF_FFV_IF(args);
    // case cstp::II: return new LF_FFV_II(args);
    case cstp::none:
      break;
    }
  }
  return NULL;
}

void ATOOLS::Getter<SF_Lorentz, SF_Key, LF_VV3P0_Quarkonia_FF>::PrintInfo(
    std::ostream &str, const size_t width) const {
  str << "VV3P0 lorentz functions";
}
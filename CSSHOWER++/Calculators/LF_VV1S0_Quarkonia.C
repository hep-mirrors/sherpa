#include "CSSHOWER++/Showers/Splitting_Function_Base.H"
#include "ATOOLS/Math/Random.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Phys/LDME.H"

#define ZMIN 1e-9
#define ZMAX 1.

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

class LF_VV1S0_Quarkonia_FI : public SF_Lorentz {

protected:
  double m_Jmax;
public:
  inline LF_VV1S0_Quarkonia_FI(const SF_Key &key) : SF_Lorentz(key) {}

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_VV1S0_Quarkonia_IF : public SF_Lorentz {

protected:
  double m_Jmax;
public:
  inline LF_VV1S0_Quarkonia_IF(const SF_Key &key) : SF_Lorentz(key) {}

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_VV1S0_Quarkonia_II : public SF_Lorentz {

protected:
  double m_Jmax;

public:
  inline LF_VV1S0_Quarkonia_II(const SF_Key &key) : SF_Lorentz(key) {}

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_V1S0V_Quarkonia_FF : public SF_Lorentz {
public:
  inline LF_V1S0V_Quarkonia_FF(const SF_Key &key) : SF_Lorentz(key) {}

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_V1S0V_Quarkonia_FI : public SF_Lorentz {

protected:
  double m_Jmax;
public:
  inline LF_V1S0V_Quarkonia_FI(const SF_Key &key) : SF_Lorentz(key) {}

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
    const double mk  =  m_flspec.Mass();
    const double muj2 = sqr(m_flavs[2].Mass(true))/Q2;
    const double muk2 = sqr(m_flspec.Mass(true))/Q2;
    const double sij = y * (Q2 - sqr(mk)) + (1-y)*(sqr(mj));
    const double M = mj; //2*Flavour((m_flavs[2].Kfcode() / 10) % 10).Mass(true);
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
  m_zmin = std::max(zmin,ZMIN);
  m_zmax = std::min(zmax,1.);
  if (m_zmax <= m_zmin) return 0.;
  const double M = m_flavs[2].Mass(true);
  // Braaten-Yuan, PRL 71 (1993) 1673, Eq. (8), with
  // LDME = 9 |R(0)|^2/(2 pi). Includes compensation for CF_QCD's
  // colour factor, polarization averaging and the dt/(2 pi t) measure;
  // the inclusive parent-gluon rate sums both colour neighbours.
  const double prefactor = 64.*M_PI/27*GetLDME(m_flavs[2].Kfcode())/M;
  return prefactor * OverEstimated(0.5, 0.) * (m_zmax - m_zmin);
}

double LF_VV1S0_Quarkonia_FF::OverEstimated(const double z, const double y) {
  return sqr(p_cf->MaxCoupling(0) / m_flavs[2].Mass(true));
}

double LF_VV1S0_Quarkonia_FF::Z() {
  return m_zmin + (m_zmax - m_zmin) * ran->Get();
}

double LF_VV1S0_Quarkonia_FI::operator()(
  double z, const double y, const double eta, const double scale,
  const double Q2) {
  const double M = m_flavs[2].Mass(true);
  const double ma2 = sqr(p_ms->Mass(m_flspec.Kfcode()));
  const double yt = ((Q2+ma2)/(Q2+ma2+sqr(M))-(1-y))/(1-y);
  const double sij = sqr(M)*(1+yt)+yt*(Q2+ma2);
  if (sij < sqr(M)/(1-z)) return 0.;
  double value = sqr(sij) + sqr(sqr(M))
    - 2*(1-z)*(sij+sqr(M))*sij + 2*sqr(sij*(1-z));
  value /= sij*sqr(sij-sqr(M));
  // Convert the paper's ds measure to the shower's dt/t measure.
  value *= scale/(z*(1-z)*sij);
  return sqr(p_cf->Coupling(scale, 0))*value*JFI(y, eta, scale);
}

double LF_VV1S0_Quarkonia_FI::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  m_zmin = std::max(zmin,ZMIN);
  m_zmax = std::min(zmax,ZMAX);
  m_Jmax = 5.;
  if (m_zmax <= m_zmin) return 0.;
  const double M = m_flavs[2].Mass(true);
  // Same Braaten-Yuan normalization as FF: includes compensation for
  // Coupling=alpha_s/4, polarization averaging and dt/(2 pi t).
  // The inclusive parent-gluon rate sums both colour neighbours.
  const double prefactor = 64.*M_PI/27*GetLDME(m_flavs[2].Kfcode())/M;
  return prefactor * OverEstimated(0.5, 0.) * (m_zmax - m_zmin);
}

double LF_VV1S0_Quarkonia_FI::OverEstimated(const double z, const double y) {
  const double M = m_flavs[2].Mass(true);
  return sqr(p_cf->MaxCoupling(0)) /sqr(M) * m_Jmax;
}

double LF_VV1S0_Quarkonia_FI::Z() {
  return m_zmin + (m_zmax - m_zmin) * ran->Get();
}

double LF_VV1S0_Quarkonia_IF::operator()(
  const double z, const double y, const double eta, const double scale,
  const double Q2) {
  const double M2 = sqr(m_flavs[2].Mass(true));
  const double mk2 = p_ms->Mass2(m_flspec);
  // Spacelike invariant, reconstructed without selecting an evolution scheme.
  const double V = (Q2 + M2 + mk2)*y/z - M2;
  // Physical support of the collinear crossed expression.
  if (IsBad(V) || !(V >= z*M2/(1-z))) return 0.;
  const double r = (1-z) - z*M2/V;
  // Cross s*F_T with -z, then convert dV to dt/t; dV/dt=1/(1-z)
  // for the supported IF gluon evolution schemes. JIF supplies flux/PDF once.
  const double value = scale*(1+sqr(r))/(z*(1-z)*sqr(V+M2));
  return sqr(p_cf->Coupling(scale, 0))*value*JIF(z, y, eta, scale);
}

double LF_VV1S0_Quarkonia_IF::OverIntegrated(
  const double zmin, const double zmax, const double scale, const double xbj) {
  m_zmin = std::max(zmin,ZMIN);
  m_zmax = std::min(zmax,1.);
  m_Jmax = m_disablepdfrescaling ? 1. : 5.;
  const double M = m_flavs[2].Mass(true);
  // Same colour, polarization and dt/(2*pi*t) convention as FF/FI 1S0;
  const double prefactor = 64.*M_PI/27*GetLDME(m_flavs[2].Kfcode())/M;
  return prefactor * 0.25 * OverEstimated(0.5, 0.) * (m_zmax-m_zmin)/(m_zmin*m_zmax);
}

double LF_VV1S0_Quarkonia_IF::OverEstimated(const double z, const double y) {
  // Including PDF-free JIF=1/z, one envelope covers all supported gluon
  //  schemes: t <= (1-z)*(V+M^2)+M^2 and V >= z*M^2/(1-z).
  return 4.*sqr(p_cf->MaxCoupling(0)/m_flavs[2].Mass(true))*m_Jmax/sqr(z);
}

double LF_VV1S0_Quarkonia_IF::Z() {
  const double r = ran->Get();
  return 1./((1-r)/m_zmin+r/m_zmax);
}

double LF_VV1S0_Quarkonia_II::operator()(
  const double z, const double y, const double eta, const double scale, const double Q2) {
  if (!(z > 0. && z < 1. && y > 0. && y < 1-z && scale > 0.)) return 0.;
  const double M2 = sqr(m_flavs[2].Mass(true));
  const double mb2 = p_ms->Mass2(m_flspec);
  if (!(Q2 > M2+mb2)) return 0.;
  // ConstructIIDipole: V=-(pa-pH)^2=y*(sab-mb^2)-M^2.
  // Reconstruct the invariant without selecting an evolution scheme.
  const double V = (Q2-M2-mb2)*y/z-M2;
  if (IsBad(V) || !(V >= z*M2/(1-z))) return 0.;
  const double r = (1-z)-z*M2/V;
  // Same collinear -z crossing as IF. Convert dV to dt/t and call JII
  // once, following FFV_II's sampled-z flux/PDF convention.
  const double value = scale*(1+sqr(r))/(z*(1-z)*sqr(V+M2));
  return sqr(p_cf->Coupling(scale, 0))*value*JII(z, y, eta, scale);
}

double LF_VV1S0_Quarkonia_II::OverIntegrated(
  const double zmin, const double zmax,  const double scale, const double xbj) {
  m_zmin = std::max(zmin,ZMIN);
  m_zmax = std::min(zmax,1.);
  m_Jmax = m_disablepdfrescaling ? 1. : 5.;
  const double M = m_flavs[2].Mass(true);
  // Same per-colour-connection normalization as IF; do not count the
  // second incoming leg here. L, spin and evolution constants enter once.
  const double prefactor = 64.*M_PI/27*GetLDME(m_flavs[2].Kfcode())/M;
  return prefactor * 0.25 * OverEstimated(0.5, 0.)
    * (m_zmax-m_zmin)/(m_zmin*m_zmax);
}

double LF_VV1S0_Quarkonia_II::OverEstimated(const double z, const double y) {
  if (!(z > 0. && z < 1.)) return 0.;
  // All current gluon schemes obey t <= (1-z)*(V+M^2)+M^2.
  // Together with V >= z*M^2/(1-z), this bounds PDF-free JII=1/z.
  return 4.*sqr(p_cf->MaxCoupling(0)/m_flavs[2].Mass(true))*m_Jmax/sqr(z);
}

double LF_VV1S0_Quarkonia_II::Z() {
  const double r = ran->Get();
  return 1./((1-r)/m_zmin+r/m_zmax);
}

double LF_V1S0V_Quarkonia_FF::operator()(const double z, const double y, const double eta, const double scale, const double Q2) { 
    const double mi  =  m_flavs[1].Mass();
    const double mk  =  m_flspec.Mass();
    const double mui2 = mi/Q2;
    const double muk2 = sqr(m_flspec.Mass(true))/Q2;
    const double sij = y * (Q2 - sqr(mk)) + (1-y)*(sqr(mi));
    const double M = 2*Flavour((m_flavs[1].Kfcode() / 10) % 10).Mass(true);
    double value = 0;
    value  = sqr(sij) + sqr(sqr(M)) - 2 * (z) * (sij + sqr(M))*sij + 2*sqr(sij*(z));
    value /= sij*sqr(sij-sqr(M));
    value *= 1. / ( (1 - mui2 - muk2) + 1./ y * ( mui2 ) );
    value *= 1. / (1 + sqr(1-z) * sqr(mi) / scale);
    return sqr(p_cf->Coupling(scale, 0)) * value * JFF(y, 0, mui2, muk2, 0);
}

double LF_V1S0V_Quarkonia_FF::OverIntegrated(const double zmin, const double zmax, const double scale, const double xbj) {
  m_zmin = std::max(zmin,0.);
  m_zmax = std::min(zmax,1.-ZMIN);
  const double m = Flavour((m_flavs[2].Kfcode() / 10) % 10).Mass(true);
  const double prefactor = 4./27*GetLDME(m_flavs[1].Kfcode())/m;
  return prefactor * sqr(p_cf->MaxCoupling(0)) / sqr(2*m) * (m_zmax - m_zmin);
}

double LF_V1S0V_Quarkonia_FF::OverEstimated(const double z, const double y) {
  const double m = Flavour((m_flavs[1].Kfcode() / 10) % 10).Mass(true);
  return sqr(p_cf->MaxCoupling(0)) /sqr(2*m);
}

double LF_V1S0V_Quarkonia_FF::Z() {
  return m_zmin + (m_zmax - m_zmin) * ran->Get();
}

double LF_V1S0V_Quarkonia_FI::operator()(const double z, const double y, const double eta, const double scale, const double Q2) {
  const double mi  =  m_flavs[1].Mass();
  const double mk  =  m_flspec.Mass();
  const double mui2 = mi/Q2;
  const double muk2 = sqr(m_flspec.Mass(true))/Q2;
  const double sij = (Q2 + sqr(mi)) * y / (1 - y) + sqr(mi);
  const double M = 2*Flavour((m_flavs[1].Kfcode() / 10) % 10).Mass(true);
  double value = 0;
  value  = sqr(sij) + sqr(sqr(M)) - 2 * (z) * (sij + sqr(M))*sij + 2*sqr(sij*(z));
  value /= sij*sqr(sij-sqr(M));
  value *= 1. / ( (1 - mui2 - muk2) + 1./ y * ( mui2 ) );
  value *= 1. / (1 + sqr(1-z) * sqr(mi) / scale);
  return sqr(p_cf->Coupling(scale, 0)) * value * JFI(y, eta, scale);
}

double LF_V1S0V_Quarkonia_FI::OverIntegrated(const double zmin, const double zmax, const double scale, const double xbj) {
  m_zmin = std::max(zmin,ZMIN);
  m_zmax = std::min(zmax,1.);
  m_Jmax = 5.;//?
  const double m = Flavour((m_flavs[1].Kfcode() / 10) % 10).Mass(true);
  const double prefactor = 4./27*GetLDME(m_flavs[1].Kfcode())/m;
  return prefactor * sqr(p_cf->MaxCoupling(0)) / sqr(2*m) * (zmax - zmin) * m_Jmax;
}

double LF_V1S0V_Quarkonia_FI::OverEstimated(const double z, const double y) {
  const double m = Flavour((m_flavs[1].Kfcode() / 10) % 10).Mass(true);
  return sqr(p_cf->MaxCoupling(0)) /sqr(2*m) * m_Jmax;
}

double LF_V1S0V_Quarkonia_FI::Z() {
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
    case cstp::FF: return new LF_VV1S0_Quarkonia_FF(args);
    case cstp::FI: return new LF_VV1S0_Quarkonia_FI(args);
    case cstp::IF: return new LF_VV1S0_Quarkonia_IF(args);
    case cstp::II: return new LF_VV1S0_Quarkonia_II(args);
    case cstp::none:
      break;
    }
  }
if ((args.m_mode == 0 && args.p_v->in[0].IntSpin() == 2 &&
       args.p_v->in[1].IntSpin() == 0 && args.p_v->in[2].IntSpin() == 2) ||
      (args.m_mode == 1 && args.p_v->in[0].IntSpin() == 2 &&
       args.p_v->in[1].IntSpin() == 2 && args.p_v->in[2].IntSpin() == 0)) {
    switch (args.m_type) {
    case cstp::FF: return new LF_V1S0V_Quarkonia_FF(args);
    case cstp::FI: return new LF_V1S0V_Quarkonia_FI(args);
    case cstp::IF: NULL;
    case cstp::II: NULL;
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

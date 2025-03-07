#include "CSSHOWER++/Showers/Splitting_Function_Base.H"

namespace CSSHOWER {

class LF_FFV_Quarkonia_FF : public SF_Lorentz {
public:
  inline LF_FFV_Quarkonia_FF(const SF_Key &key) : SF_Lorentz(key) {}

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_FVF_Quarkonia_FF : public SF_Lorentz {
public:
  inline LF_FVF_Quarkonia_FF(const SF_Key &key) : SF_Lorentz(key) {}

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_VSV_Quarkonia_FF : public SF_Lorentz {
public:
  inline LF_VSV_Quarkonia_FF(const SF_Key &key) : SF_Lorentz(key) {}

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();

private:
  inline double f0(double r, double y);
  inline double f1(double r, double y);
  inline double f2(double r, double y);
  inline double g0(double r, double y);
  inline double g1(double r, double y);
  inline double g2(double r, double y);
};
} // namespace CSSHOWER

#include "ATOOLS/Math/Random.H"
#include "MODEL/Main/Single_Vertex.H"

using namespace CSSHOWER;
using namespace ATOOLS;

double LF_FFV_Quarkonia_FF::operator()(const double zz, const double y,
                                       const double eta, const double scale,
                                       const double Q2) {
  const double z = 1 - zz;
  double muij2 = sqr(p_ms->Mass(m_flavs[0])) / Q2;
  double mi2 = sqr(p_ms->Mass(m_flavs[1]));
  double mj2 = sqr(p_ms->Mass(m_flavs[2]));
  double mk2 = sqr(p_ms->Mass(m_flspec));
  double mui2 = mi2 / Q2, muj2 = mj2 / Q2, muk2 = mk2 / Q2;
  const double mjperp2 =
      (scale * z * (1. - z) - mi2 * z - mj2 * (1. - z)) + mj2;
  const double newscale = sqrt(scale * mjperp2);
  // the massless case
  if (muij2 == 0. || mui2 == 0. || muj2 == 0.) {
    msg_Error() << "Cannot make massless quarkonia emission" << std::endl;
    exit(EXIT_FAILURE);
  } else {
    // the massive case
    double value =
        1. / sqr(sqr(1 - muij2)) *
        ((1 - 2 * muij2 - 47 * sqr(muij2)) -
         z * (1 - muij2) * (1 - sqr(sqrt(mui2) + sqrt(muj2))) +
         4 * (z * (1 - z)) / (2 - z) * (1 - mui2) -
         4 * (8 - 7 * z - 5 * z * z) / (2 - z) * muij2 * (1 - muij2) +
         12 * (z * z * (1 - z)) / sqr(2 - z) * sqr(1 - muij2));
    msg_Debugging() << METHOD << "\tcpl max: " << p_cf->MaxCoupling(0)
                    << ", cpl: " << p_cf->Coupling(newscale, 0)
                    << ", newscale:  " << newscale << std::endl;
    msg_Debugging() << METHOD << "\treturn: "
                    << 16. / 27 * sqr(p_cf->Coupling(newscale, 0)) *
                           JFF(y, mui2, muj2, muk2, muij2)
                    << std::endl;
    return 16. / 27 * sqr(p_cf->Coupling(newscale, 0)) *
           JFF(y, mui2, muj2, muk2, muij2);
  }
}

double LF_FFV_Quarkonia_FF::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  m_zmin = zmin;
  m_zmax = zmax;
  return 16. / 27 * sqr(p_cf->MaxCoupling(0)) * (0.5) *
         log((1. - zmin) / (1. - zmax));
}

double LF_FFV_Quarkonia_FF::OverEstimated(const double z, const double y) {
  return 16. / 27 * sqr(p_cf->MaxCoupling(0)) * (0.5) / (1. - z);
}

double LF_FFV_Quarkonia_FF::Z() {
  return 1. -
         (1. - m_zmin) * pow((1. - m_zmax) / (1. - m_zmin), ATOOLS::ran->Get());
}

DECLARE_GETTER(LF_FFV_Quarkonia_FF, "FFV_Quarkonia", SF_Lorentz, SF_Key);

double LF_FVF_Quarkonia_FF::operator()(const double z, const double y,
                                       const double eta, const double scale,
                                       const double Q2) {
  double muij2 = sqr(p_ms->Mass(m_flavs[0])) / Q2;
  double mi2 = sqr(p_ms->Mass(m_flavs[1]));
  double mj2 = sqr(p_ms->Mass(m_flavs[2]));
  double mk2 = sqr(p_ms->Mass(m_flspec));
  double mui2 = mi2 / Q2, muj2 = mj2 / Q2, muk2 = mk2 / Q2;
  // the massless case
  if (muij2 == 0. || mui2 == 0. || muj2 == 0. || muk2 == 0.) {
    msg_Error() << "Cannot make massless quarkonia emission" << std::endl;
    exit(EXIT_FAILURE);
  } else {
    // the massive case
    double massive =
        1. / sqr(sqr(1 - muij2)) *
        ((1 - 2 * muij2 - 47 * sqr(muij2)) -
         z * (1 - muij2) * (1 - sqr(sqrt(mui2) + sqrt(muj2))) +
         4 * (z * (1 - z)) / (2 - z) * (1 - mui2) -
         4 * (8 - 7 * z - 5 * z * z) / (2 - z) * muij2 * (1 - muij2) +
         12 * (z * z * (1 - z)) / sqr(2 - z) * sqr(1 - muij2));
    // massive *= 1./((1.-mui2-muk2)+1./y*(mui2-muij2));
    double longpol = 0.5 * (1. - z);
    double value =
        8 / 27 / M_PI / (Q2 * mui2) * p_cf->Coupling(scale, 1) * massive +
        p_cf->Coupling(scale, 0) * longpol;
    return value * JFF(y, mui2, muj2, muk2, muij2);
  }
}

double LF_FVF_Quarkonia_FF::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  m_zmin = zmin;
  m_zmax = zmax;
  return (2 * p_cf->MaxCoupling(1)) * log((1. - zmin) / (1. - zmax));
}

double LF_FVF_Quarkonia_FF::OverEstimated(const double z, const double y) {
  const double zz = 1. - z;
  return (2 * p_cf->MaxCoupling(1)) / (1. - zz);
}

double LF_FVF_Quarkonia_FF::Z() {
  return 1. -
         (1. - m_zmin) * pow((1. - m_zmax) / (1. - m_zmin), ATOOLS::ran->Get());
}

SF_Lorentz *ATOOLS::Getter<SF_Lorentz, SF_Key, LF_FFV_Quarkonia_FF>::operator()(
    const Parameter_Type &args) const {
  if (args.m_col < 0)
    return NULL;
  if ((args.m_mode == 0 && args.p_v->in[0].IntSpin() == 1 &&
       args.p_v->in[1].IntSpin() == 1 && args.p_v->in[2].IntSpin() == 2) ||
      (args.m_mode == 1 && args.p_v->in[0].IntSpin() == 1 &&
       args.p_v->in[1].IntSpin() == 2 && args.p_v->in[2].IntSpin() == 1)) {
    switch (args.m_type) {
    case cstp::FF:
      return new LF_FFV_Quarkonia_FF(args);
    // case cstp::FI: return new LF_FFV_FI(args);
    // case cstp::IF: return new LF_FFV_IF(args);
    // case cstp::II: return new LF_FFV_II(args);
    case cstp::none:
      break;
    }
  }
  if ((args.m_mode == 0 && args.p_v->in[0].IntSpin() == 1 &&
       args.p_v->in[1].IntSpin() == 2 && args.p_v->in[2].IntSpin() == 1) ||
      (args.m_mode == 1 && args.p_v->in[0].IntSpin() == 1 &&
       args.p_v->in[2].IntSpin() == 2 && args.p_v->in[1].IntSpin() == 1)) {
    switch (args.m_type) {
    case cstp::FF:
      return new LF_FVF_Quarkonia_FF(args);
    // case cstp::FI: return new LF_FVF_FI(args);
    // case cstp::IF: return new LF_FVF_IF(args);
    // case cstp::II: return new LF_FVF_II(args);
    case cstp::none:
      break;
    }
  }
  if (args.p_v->in[0].IntSpin() == 2 && args.p_v->in[1].IntSpin() == 1 &&
      args.p_v->in[2].IntSpin() == 1) {
    switch (args.m_type) {
      return NULL;
    // case cstp::FF: return new LF_VFF_Quarkonia_FF(args);
    // case cstp::FI: return new LF_VFF_FI(args);
    // case cstp::IF: return new LF_VFF_IF(args);
    // case cstp::II: return new LF_VFF_II(args);
    case cstp::none:
      break;
    }
  }
  return NULL;
}

void ATOOLS::Getter<SF_Lorentz, SF_Key, LF_FFV_Quarkonia_FF>::PrintInfo(
    std::ostream &str, const size_t width) const {
  str << "ffv lorentz functions";
}

double LF_VSV_Quarkonia_FF::operator()(
    const double z, const double w, const double eta, const double scale,
    const double Q2) { // w replaces the usual y variable in CS dipoles
  double muij2 = sqr(p_ms->Mass(m_flavs[0])) / Q2;
  double mi2 = sqr(p_ms->Mass(m_flavs[1])); // this must be the J/psi
  double mj2 = sqr(p_ms->Mass(m_flavs[2])); // this must be gluon
  double mk2 =
      sqr(p_ms->Mass(m_flspec)); // this is whatever the recoil partner is
  double mui2 = mi2 / Q2, muj2 = mj2 / Q2, muk2 = mk2 / Q2;
  const double sij = (w * (Q2 - mk2) + (1. - w) * (mi2 + mj2));
  const double lam(sqr(Q2 + sij + mk2) - 4 * sij * mk2), beta(sqrt(lam)),
      gam(0.5 * (Q2 - sij - mk2) + beta);
  const double zb = (Q2 - sij - mk2) / beta *
                    (z - mk2 / gam * (sij - mk2) / (Q2 - sij - mk2));
  const double kt2 = zb * (1. - zb) * sij - zb * mj2;
  const double pkdotptij = 0.25 * (sqrt(lam) + Q2 + mk2 - sij);
  double const y = 1 / Q2 *
                   (Q2 * (1. - muk2) / 2 *
                    (gam * zb / beta - kt2 * mk2 / (zb * beta * gam) +
                     (pkdotptij) * (kt2 / zb / beta * (1 + mk2 / gam) -
                                    zb / beta * (sij + gam))));
  const double miperp2 =
      (scale * z * (1. - z) - mj2 * z - mi2 * (1. - z)) + mi2;
  //if ( y < (mui2 +z*z)/(2*z) || (1+mui2)/2.0 < y) return 0;
  const double newscale = sqrt(scale * miperp2);
  const double mq = ATOOLS::Flavour(kf_c).Mass();
  double value = (f0(mui2, y) + g0(mui2, y) * (1 + mui2 - 2 * y) /
                                    (2 * (y - mui2) * sqrt(y * y - mui2)) *
                                    log((y - mui2 + sqrt(y * y - mui2)) /
                                        (y - mui2 - sqrt(y * y - mui2))));
  msg_Debugging() << METHOD << ",  ->"
            << ((1 + mui2 - 2 * y) / (2 * (y - mui2) * sqrt(y * y - mui2)) *
                log((y - mui2 + sqrt(y * y - mui2)) /
                    (y - mui2 - sqrt(y * y - mui2))))
            << std::endl;
  value += z * (f1(mui2, y) + g1(mui2, y) * (1 + mui2 - 2 * y) /
                                  (2 * (y - mui2) * sqrt(y * y - mui2)) *
                                  log((y - mui2 + sqrt(y * y - mui2)) /
                                      (y - mui2 - sqrt(y * y - mui2))));
  value += sqr(z) * (f2(mui2, y) + g2(mui2, y) * (1 + mui2 - 2 * y) /
                                       (2 * (y - mui2) * sqrt(y * y - mui2)) *
                                       log((y - mui2 + sqrt(y * y - mui2)) /
                                           (y - mui2 - sqrt(y * y - mui2))));
  ;
  value *= 1. / (sqr(1 - y)) / sqr(y - mui2) / sqr(y * y - mui2);
  const double preF = (p_cf->Coupling(newscale, 0) / mq) *
                      (p_cf->Coupling(newscale, 0) / mq) *
                      (p_cf->Coupling(newscale, 0) / mq);
  // return 5. / 5184 / M_PI * preF * value * JFF(y, mui2, muj2, muk2, muij2);
  return 16. / 27 * sqr(p_cf->MaxCoupling(0)) * (0.5) / (1. - z) * 10;
}

double LF_VSV_Quarkonia_FF::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  m_zmin = zmin;
  m_zmax = zmax;
  return 100 * 16. / 27 * sqr(p_cf->MaxCoupling(0)) * (0.5) *
         log((1. - zmin) / (1. - zmax));
}

double LF_VSV_Quarkonia_FF::OverEstimated(const double z, const double y) {
  return 100 * 16. / 27 * sqr(p_cf->MaxCoupling(0)) * (0.5) / (1. - z);
}

double LF_VSV_Quarkonia_FF::Z() {
  return 1. -
         (1. - m_zmin) * pow((1. - m_zmax) / (1. - m_zmin), ATOOLS::ran->Get());
}

inline double LF_VSV_Quarkonia_FF::f0(double r, double y) {
  const double r2 = r * r, r3 = r2 * r;
  const double y2 = y * y, y3 = y2 * y, y4 = y3 * y;
  const double y5 = y4 * y, y6 = y5 * y, y7 = y6 * y;
  return r2 * (1 + r) * (3 + 12 * r + 13 * r2) -
         16 * r2 * (1 + r) * (1 + 3 * r) * y -
         2 * r * (3 - 9 * r - 21 * r2 + 7 * r3) * y2 +
         8 * r * (4 + 3 * r + 3 * r2) * y3 - 4 * r * (9 - 3 * r - 4 * r2) * y4 -
         16 * (1 + 3 * r + 3 * r2) * y5 + 8 * (6 + 7 * r) * y6 - 32 * y7;
}

// Equation (5)
inline double LF_VSV_Quarkonia_FF::f1(double r, double y) {
  const double r2 = r * r, r3 = r2 * r;
  const double y2 = y * y, y3 = y2 * y, y4 = y3 * y;
  const double y5 = y4 * y, y6 = y5 * y;
  return -2 * r * (1 + 5 * r + 19 * r2 + 7 * r3) * y + 96 * r2 * (1 + r) * y2 +
         8 * (1 - 5 * r - 22 * r2 - 2 * r3) * y3 + 16 * r * (7 + 3 * r) * y4 -
         8 * (5 + 7 * r) * y5 + 32 * y6;
}

// Equation (6)
inline double LF_VSV_Quarkonia_FF::f2(double r, double y) {
  const double r2 = r * r, r3 = r2 * r;
  const double y2 = y * y, y3 = y2 * y, y4 = y3 * y, y5 = y4 * y;
  return r * (1 + 5 * r + 19 * r2 + 7 * r3) - 48 * r2 * (1 + r) * y -
         4 * (1 - 5 * r - 22 * r2 - 2 * r3) * y2 - 8 * r * (7 + 3 * r) * y3 +
         4 * (5 + 7 * r) * y4 - 16 * y5;
}

// Equation (7)
inline double LF_VSV_Quarkonia_FF::g0(double r, double y) {
  const double r2 = r * r, r3 = r2 * r;
  const double y2 = y * y, y3 = y2 * y, y4 = y3 * y, y5 = y4 * y, y6 = y5 * y;
  return r3 * (1 - r) * (3 + 24 * r + 13 * r2) -
         4 * r3 * (7 - 3 * r - 12 * r2) * y -
         2 * r3 * (17 + 22 * r - 7 * r2) * y2 +
         4 * r2 * (13 + 5 * r - 6 * r2) * y3 -
         8 * r * (1 + 2 * r + 5 * r2 + 2 * r3) * y4 -
         8 * r * (3 - 11 * r - 6 * r2) * y5 + 8 * (1 - 2 * r - 5 * r2) * y6;
}

// Equation (8)
inline double LF_VSV_Quarkonia_FF::g1(double r, double y) {
  const double r2 = r * r, r3 = r2 * r;
  const double y2 = y * y, y3 = y2 * y, y4 = y3 * y, y5 = y4 * y;
  return -2 * r2 * (1 + r) * (1 - r) * (1 + 7 * r) * y +
         8 * r2 * (1 + 3 * r) * (1 - 4 * r) * y2 +
         4 * r * (1 + 10 * r + 57 * r2 + 4 * r3) * y3 -
         8 * r * (1 + 29 * r + 6 * r2) * y4 - 8 * (1 - 8 * r - 5 * r2) * y5;
}

// Equation (9)

inline double LF_VSV_Quarkonia_FF::g2(double r, double y) {
  const double r2 = r * r, r3 = r2 * r;
  const double y2 = y * y, y3 = y2 * y, y4 = y3 * y, y5 = y4 * y;
  return r2 * (1 + r) * (1 - r) * (1 + 7 * r) -
         4 * r2 * (1 + 3 * r) * (1 - 4 * r) * y -
         2 * r * (1 + 10 * r + 57 * r2 + 4 * r3) * y2 +
         4 * r * (1 + 29 * r + 6 * r2) * y2 + 4 * (1 - 8 * r - 5 * r2) * y4;
}

DECLARE_GETTER(LF_VSV_Quarkonia_FF, "VSV_Quarkonia", SF_Lorentz, SF_Key);

SF_Lorentz *ATOOLS::Getter<SF_Lorentz, SF_Key, LF_VSV_Quarkonia_FF>::operator()(
    const Parameter_Type &args) const {
  if (args.m_col < 0)
    return NULL;
  if ((args.m_mode == 0 && args.p_v->in[0].IntSpin() == 2 &&
       args.p_v->in[1].IntSpin() == 0 && args.p_v->in[2].IntSpin() == 2) ||
      (args.m_mode == 1 && args.p_v->in[0].IntSpin() == 2 &&
       args.p_v->in[1].IntSpin() == 2 && args.p_v->in[2].IntSpin() == 0)) {
    switch (args.m_type) {
    case cstp::FF:
      return new LF_VSV_Quarkonia_FF(args);
    // case cstp::FI: return new LF_FFV_FI(args);
    // case cstp::IF: return new LF_FFV_IF(args);
    // case cstp::II: return new LF_FFV_II(args);
    case cstp::none:
      break;
    }
  }
  if ((args.m_mode == 0 && args.p_v->in[0].IntSpin() == 1 &&
       args.p_v->in[1].IntSpin() == 2 && args.p_v->in[2].IntSpin() == 1) ||
      (args.m_mode == 1 && args.p_v->in[0].IntSpin() == 1 &&
       args.p_v->in[2].IntSpin() == 2 && args.p_v->in[1].IntSpin() == 1)) {
    switch (args.m_type) {
    case cstp::FF:
      return new LF_VSV_Quarkonia_FF(args);
    // case cstp::FI: return new LF_FVF_FI(args);
    // case cstp::IF: return new LF_FVF_IF(args);
    // case cstp::II: return new LF_FVF_II(args);
    case cstp::none:
      break;
    }
  }
  if (args.p_v->in[0].IntSpin() == 2 && args.p_v->in[1].IntSpin() == 1 &&
      args.p_v->in[2].IntSpin() == 1) {
    switch (args.m_type) {
      return NULL;
    // case cstp::FF: return new LF_VFF_Quarkonia_FF(args);
    // case cstp::FI: return new LF_VFF_FI(args);
    // case cstp::IF: return new LF_VFF_IF(args);
    // case cstp::II: return new LF_VFF_II(args);
    case cstp::none:
      break;
    }
  }
  return NULL;
}

void ATOOLS::Getter<SF_Lorentz, SF_Key, LF_VSV_Quarkonia_FF>::PrintInfo(
    std::ostream &str, const size_t width) const {
  str << "VSV lorentz functions";
}

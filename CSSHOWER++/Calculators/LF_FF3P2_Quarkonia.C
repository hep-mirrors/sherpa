#include "CSSHOWER++/Showers/Splitting_Function_Base.H"
#include "ATOOLS/Phys/LDME.H"

namespace CSSHOWER {

class LF_FF3P2_Quarkonia_FF : public SF_Lorentz {
public:
  inline LF_FF3P2_Quarkonia_FF(const SF_Key &key) : SF_Lorentz(key) {}

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_FF3P2_Quarkonia_FI : public SF_Lorentz {

protected:
  double m_Jmax;

public:
  inline LF_FF3P2_Quarkonia_FI(const SF_Key &key) : SF_Lorentz(key) {}

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_FF3P2_Quarkonia_IF : public SF_Lorentz {

protected:
  double m_Jmax;

public:
  inline LF_FF3P2_Quarkonia_IF(const SF_Key &key) : SF_Lorentz(key) {}

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_FF3P2_Quarkonia_II : public SF_Lorentz {

protected:
  double m_Jmax;

public:
  inline LF_FF3P2_Quarkonia_II(const SF_Key &key) : SF_Lorentz(key) {}

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

double LF_FF3P2_Quarkonia_FF::operator()(const double z, const double y,
                                        const double eta, const double scale,
                                        const double Q2) {
  // based on Eq. 23 of arXiv:hep-ph/9405348v1
  // again this is b -> c 3P2 but holds for any 3P2 singlet state
  const double mi  = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true); 
  const double mj  = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  const double mk  = ATOOLS::Flavour(m_flspec.Kfcode()).Mass(true);
  const double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double mui2 = sqr(mi) / Q2, muj2 = sqr(mj) / Q2, muk2 = sqr(mk) / Q2, muij2 = sqr(mij) / Q2;
  const double M = mi + mij;
  const double ri = mi / M, rij = mij / M;
  const double sij = y * (Q2 - sqr(mi) - sqr(mj) - sqr(mk)) + (sqr(mi) + sqr(mj));
  double value = sqr(sqr(sqr(M)))/sqr(sqr(sij - sqr(rij*M)))*( 320 * sqr(ri) * cube(rij) * sqr(sqr(1 - rij*(1-z))) );
  value += sqr(cube(M))/cube(sij - sqr(rij*M)) * 8 * ri * sqr(rij) * cube(1 - rij * (1 - z)) * (
    2 * (4 + 13 * ri) - (1 + 70 * ri - 26 * sqr(ri)) * (1 - z) - (7 + 8 * ri) * rij * sqr(1 - z)
  );
  value += sqr(sqr(M))/sqr(sij - sqr(rij*M)) * (-4) * sqr(rij) * sqr( 1 - rij * (1 - z)) * ( 
    4 * (1 + 4 * ri) - (7 + 12 * ri - 32 * sqr(ri)) * (1 - z) 
    + 2 * (1 + 13 * ri - 26 * sqr(ri) + 8 * cube(ri)) * sqr(1 - z) + (1 - 30 * ri - 5 * sqr(ri) + 4 * cube(ri)) * cube(1 - z)
  );
  value += sqr(M)/(sij - sqr(rij*M)) * 4 * sqr(rij) * z * (
    2 - 4 * (1 - 2 * ri) * (1 - z)
    + (5 - 8 * ri + 12 * sqr(ri)) * sqr(1 - z)
    - 2 * (1 - 2 * ri) * (3 + 2 * sqr(ri)) * cube(1 - z) 
    + (3 - 12 * ri + 12 * sqr(ri) + 2 * sqr(sqr(ri))) * sqr(sqr((1 - z)))
  );
  value *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
  value *= 1. / (1 + sqr( 1 - z) * sqr(mi) / scale + sqr(z) * sqr(mj) / scale);
  value *= sqr(p_cf->Coupling(scale, 0)) * ri * cube(rij) / sqr(sqr(1-rij*(1-z))) * JFF(y, mui2, muj2, muk2, muij2);
  // BEWARE: this is missing a global constant factor in front (so does the overestimate)
  return value;
}

double LF_FF3P2_Quarkonia_FF::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  const double mi = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double rij = mij / (mi + mij);
  const double ri = mi / (mi + mij);
  m_zmin = zmin; 
  m_zmax = zmax;
  const double prefactor = 4./27*GetLDME(m_flavs[2].Kfcode())/pow(mi,5);
  return prefactor * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.1561 * (m_zmax - m_zmin); ;
}

double LF_FF3P2_Quarkonia_FF::OverEstimated(const double z, const double y) {
  const double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double mi = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double rij = mij / (mi + mij);
  const double ri = mi / (mi + mij);
  return sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.1561;
}

double LF_FF3P2_Quarkonia_FF::Z() {
  return m_zmin + (m_zmax - m_zmin) * ATOOLS::ran->Get();
}

double LF_FF3P2_Quarkonia_FI::operator()(const double zz, const double y,
                                       const double eta, const double scale,
                                       const double Q2) {
  const double z = 1 - zz;
  double mi2 =
      sqr(p_ms->Mass(m_flavs[1])); // works with the mapping c -> J/Psi(1S) c
  double mj2 = sqr(p_ms->Mass(m_flavs[2]));
  double ma2 = sqr(p_ms->Mass(m_flspec));
  double mij2 = p_ms->Mass2(m_flavs[0]);
  const double sij = -y / (1 - y) * (Q2 - ma2) + (mi2 + mj2) / (1 - y);
  double muij2 = mij2 / sij;
  double mui2 = mi2 / sij, muj2 = mj2 / sij, mua2 = ma2 / sij;
  const double newscale = sqrt(scale);
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
    const double LDME = 1. / sqrt(mi2) / sqr(sij) *
                        (m_flavs[2].IsOctetMeson() ? (1.5E-02) : pow(0.82, 3));
    return 16. / 27 * LDME * sqr(p_cf->Coupling(newscale, 0)) * value /
           sqr(sij) * JFI(y, eta, scale);
  }
}

double LF_FF3P2_Quarkonia_FI::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  m_zmin = zmin;
  m_zmax = zmax;
  m_Jmax = 5.;
  return 16. / 27 * pow(0.82, 3) / Flavour(kf_c).Mass() *
         sqr(p_cf->MaxCoupling(0)) * (0.5) * log((1. - zmin) / (1. - zmax)) *
         m_Jmax;
}

double LF_FF3P2_Quarkonia_FI::OverEstimated(const double z, const double y) {
  return 16. / 27 * pow(0.82, 3) / Flavour(kf_c).Mass() *
         sqr(p_cf->MaxCoupling(0)) * (0.5) / (1. - z) * m_Jmax;
}

double LF_FF3P2_Quarkonia_FI::Z() {
  return 1. -
         (1. - m_zmin) * pow((1. - m_zmax) / (1. - m_zmin), ATOOLS::ran->Get());
}

double LF_FF3P2_Quarkonia_IF::operator()(const double zz, const double y,
                                       const double eta, const double scale,
                                       const double Q2) {
  // c --> c J/psi
  // (aj) - (k) --> a - j - k
  // zz = x_{jk,a};  y = u
  double ma2 = sqr(p_ms->Mass(m_flavs[1])); // quark
  double mj2 = sqr(p_ms->Mass(m_flavs[2])); // meson
  double mk2 = sqr(p_ms->Mass(m_flspec));
  double maj2 = p_ms->Mass2(m_flavs[0]); // quark
  const double pkpa = Q2 / 2 / zz * (1 - y);
  const double pkpj = 0.5 * (1 - zz) / zz * (mj2 + mk2 + ma2 - Q2);
  const double taj =
      y * Q2 + ma2 * (1 - y) + mj2 + y * (mj2 + mk2 + 2 * pkpj); // (pa - pj)^2
  const double z =
      pkpj / pkpa; // this is the momentum fraction of J/psi w.r.t. parent quark
  double muaj2 = maj2 / taj;
  double mua2 = ma2 / taj, muj2 = mj2 / taj, muk2 = mk2 / taj;
  const double newscale = sqrt(scale);
  // the massless case
  if (muaj2 == 0. || mua2 == 0. || muj2 == 0.) {
    msg_Error() << "Cannot make massless quarkonia emission" << std::endl;
    exit(EXIT_FAILURE);
  } else {
    // the massive case
    double value =
        1. / sqr(sqr(1 - muaj2)) *
        ((1 - 2 * muaj2 - 47 * sqr(muaj2)) -
         z * (1 - muaj2) * (1 - sqr(sqrt(fabs(mua2)) + sqrt(fabs(muj2)))) +
         4 * (z * (1 - z)) / (2 - z) * (1 - mua2) -
         4 * (8 - 7 * z - 5 * z * z) / (2 - z) * muaj2 * (1 - muaj2) +
         12 * (z * z * (1 - z)) / sqr(2 - z) * sqr(1 - muaj2));
    const double LDME = 1. / sqrt(ma2) / sqr(taj) *
                        (m_flavs[2].IsOctetMeson() ? (1.5E-02) : pow(0.82, 3));
    return 16. / 27 * LDME * sqr(p_cf->Coupling(newscale, 0)) * value /
           sqr(taj) * JFI(y, eta, scale);
  }
}

double LF_FF3P2_Quarkonia_IF::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  m_zmin = zmin;
  m_zmax = zmax;
  m_Jmax = m_flavs[0].Kfcode() < 3 ? 5. : 1.;
  return 16. / 27 * pow(0.82, 3) / Flavour(kf_c).Mass() *
         (2.0 * p_cf->MaxCoupling(0) * 2. + 0.5 * p_cf->MaxCoupling(1)) *
         log((1. - zmin) / (1. - zmax)) * m_Jmax;
}

double LF_FF3P2_Quarkonia_IF::OverEstimated(const double z, const double y) {
  return 16. / 27 * pow(0.82, 3) / Flavour(kf_c).Mass() *
         (2.0 * p_cf->MaxCoupling(0) * 2. + 0.5 * p_cf->MaxCoupling(1)) /
         (1. - z) * m_Jmax;
}

double LF_FF3P2_Quarkonia_IF::Z() {
  return 1. -
         (1. - m_zmin) * pow((1. - m_zmax) / (1. - m_zmin), ATOOLS::ran->Get());
}

double LF_FF3P2_Quarkonia_II::operator()(const double zz, const double y,
                                       const double eta, const double scale,
                                       const double Q2) {
  // works with the mapping [c aj] ->[c a] [J/psi j]
  // zz is effectively xjab and y vtildej
  double ma2 = sqr(p_ms->Mass(m_flavs[1]));
  double mj2 = sqr(p_ms->Mass(m_flavs[2]));
  double mb2 = sqr(p_ms->Mass(m_flspec));
  double maj2 = p_ms->Mass2(m_flavs[0]);
  const double sab = (Q2 - mj2) / zz - (1. - zz) / zz * (ma2 + mb2);
  const double z = 1 - (zz + y) + (ma2 + mb2 + mj2)/(sab - ma2 - mb2); // this is mom. fraction of J/psi w.r.t. to progenitor c
  const double taj =  ma2 + mj2 - 2 * y * (sab - ma2 - mb2); // this is (pa - pj)^2
  double muaj2 = maj2 / taj;
  double mua2 = ma2 / taj, muj2 = mj2 / taj, mub2 = mb2 / taj;
  const double newscale = scale;
  // the massless case
  if (muaj2 == 0. || mua2 == 0. || muj2 == 0.) {
    msg_Error() << "Cannot make massless quarkonia emission" << std::endl;
    exit(EXIT_FAILURE);
  } else {
    // the massive case
    double value =
        1. / sqr(sqr(1 - muaj2)) *
         ( (1 - 2 * muaj2 - 47 * sqr(muaj2)) -
         z * (1 - muaj2) * (1 - sqr(sqrt(fabs(mua2)) + sqrt(fabs(muj2)))) +
         4 * (z * (1 - z)) / (2 - z) * (1 - mua2) -
         4 * (8 - 7 * z - 5 * z * z) / (2 - z) * muaj2 * (1 - muaj2) +
         12 * (z * z * (1 - z)) / sqr(2 - z) * sqr(1 - muaj2));
    const double LDME =
        (m_flavs[2].IsOctetMeson() ? 1.5E-02 / M_PI_2
                                   : 9. / 2 * M_PI * pow(0.82, 3));
    return value * JII(z, y, eta, scale);
    return 16. / 27 / 9 * sqr(p_cf->Coupling(newscale, 0)) / sqrt(ma2) * LDME *
           value / sqr(taj) * JII(z, y, eta, scale);
  }
}

double LF_FF3P2_Quarkonia_II::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  m_zmin = zmin;
  m_zmax = zmax;
  m_Jmax = m_flavs[0].Kfcode() < 3 ? 5. : 1.;
  return (4.0 * p_cf->MaxCoupling(0) + 0.5 * p_cf->MaxCoupling(1)) *
         log((1. - zmin) / (1. - zmax)) * m_Jmax;
  return 16. / 27 * pow(0.82, 3) / Flavour(kf_c).Mass() *
         sqr(p_cf->MaxCoupling(0)) * (0.5) * log((1. - zmin) / (1. - zmax)) *
         m_Jmax;
}

double LF_FF3P2_Quarkonia_II::OverEstimated(const double z, const double y) {
  return (4.0 * p_cf->MaxCoupling(0) + 0.5 * p_cf->MaxCoupling(1)) / (1. - z) *
         m_Jmax;
  return 16. / 27 * pow(0.82, 3) / Flavour(kf_c).Mass() *
         sqr(p_cf->MaxCoupling(0)) * (0.5) / (1. - z) * m_Jmax;
}

double LF_FF3P2_Quarkonia_II::Z() {
  return 1. -
         (1. - m_zmin) * pow((1. - m_zmax) / (1. - m_zmin), ATOOLS::ran->Get());
}


DECLARE_GETTER(LF_FF3P2_Quarkonia_FF, "FF3P2_Quarkonia", SF_Lorentz, SF_Key);

SF_Lorentz *ATOOLS::Getter<SF_Lorentz, SF_Key, LF_FF3P2_Quarkonia_FF>::operator()(
    const Parameter_Type &args) const {
  if (args.m_col < 0)
    return NULL;
  if ((args.m_mode == 0 && args.p_v->in[0].IntSpin() == 1 &&
       args.p_v->in[1].IntSpin() == 1 && args.p_v->in[2].IntSpin() == 4) ||
      (args.m_mode == 1 && args.p_v->in[0].IntSpin() == 1 &&
       args.p_v->in[1].IntSpin() == 0 && args.p_v->in[2].IntSpin() == 4)) {
    switch (args.m_type) {
    case cstp::FF:
      return new LF_FF3P2_Quarkonia_FF(args);
    // case cstp::FI:
    //   return new LF_FF3P2_Quarkonia_FI(args);
    // case cstp::IF:
    //   return new LF_FF3P2_Quarkonia_IF(args);
    // case cstp::II: 
    //   return new LF_FF3P2_Quarkonia_II(args);
    // case cstp::none:
    //   break;
    }
  }
  if ((args.m_mode == 0 && args.p_v->in[0].IntSpin() == 1 &&
       args.p_v->in[1].IntSpin() == 2 && args.p_v->in[2].IntSpin() == 1) ||
      (args.m_mode == 1 && args.p_v->in[0].IntSpin() == 1 &&
       args.p_v->in[2].IntSpin() == 2 && args.p_v->in[1].IntSpin() == 1)) {
    switch (args.m_type) {
    case cstp::FF:
      // return new LF_FVF_Quarkonia_FF(args);
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

void ATOOLS::Getter<SF_Lorentz, SF_Key, LF_FF3P2_Quarkonia_FF>::PrintInfo(
    std::ostream &str, const size_t width) const {
  str << "FF3P2 lorentz functions";
}



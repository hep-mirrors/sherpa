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

class LF_FFV_Quarkonia_FI : public SF_Lorentz {

protected:
  double m_Jmax;

public:
  inline LF_FFV_Quarkonia_FI(const SF_Key &key) : SF_Lorentz(key) {}

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_FFV_Quarkonia_IF : public SF_Lorentz {

protected:
  double m_Jmax;

public:
  inline LF_FFV_Quarkonia_IF(const SF_Key &key) : SF_Lorentz(key) {}

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_FFV_Quarkonia_II : public SF_Lorentz {

protected:
  double m_Jmax;

public:
  inline LF_FFV_Quarkonia_II(const SF_Key &key) : SF_Lorentz(key) {}

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
} // namespace CSSHOWER

#include "ATOOLS/Math/Random.H"
#include "MODEL/Main/Single_Vertex.H"

using namespace CSSHOWER;
using namespace ATOOLS;

double LF_FFV_Quarkonia_FF::operator()(const double zz, const double y,
                                       const double eta, const double scale,
                                       const double Q2) {
  const double z = 1 - zz;
  double mi2 = sqr(p_ms->Mass(m_flavs[1])); // works with the mapping c -> c J/psi
  double mj2 = sqr(p_ms->Mass(m_flavs[2]));
  double mk2 = sqr(p_ms->Mass(m_flspec));
  double mij2 = p_ms->Mass2(m_flavs[0]);
  const double sij = y * (Q2 - mi2 - mj2 - mk2) + (mi2 + mj2);
  double muij2 = sqr(p_ms->Mass(m_flavs[0])) / Q2;
  double mui2 = mi2 / Q2, muj2 = mj2 / Q2, muk2 = mk2 / Q2;
  // const double mjperp2 =
  //     (scale * z * (1. - z) - mi2 * z - mj2 * (1. - z)) + mj2;
  const double newscale = 9*mij2;// sqrt(scale * mjperp2);
  // the massless case
  if (muij2 == 0. || mui2 == 0. || muj2 == 0.) {
    msg_Error() << "Cannot make massless quarkonia emission" << std::endl;
    exit(EXIT_FAILURE);
  } else {
    // the massive case
    double value = 1./(sqr(sij-mij2))*(
      (sqr(sij) - -2*mi2*sij -47*sqr(mi2)) - z*(sij-mij2)*(sij-sqr(sqrt(mi2) - sqrt(mj2))) +
      + 4*z*(1-z)/(2-z)*sij*(sij-mij2) - 4*(8-7*z-5*z*z)/(2-z)*mi2*(sij-mij2) + 
      12*z*z*(1-z)/sqr(2-z)*sqr(sij-mij2)
    );
    const double LDME = p_cf->Coupling(newscale,0)/(2*M_PI) * pow(0.82,2); // GeV^3 //(m_flavs[2].IsOctetMeson() ? (1.5E-02) : pow(0.82, 3));
    const double Jprop = 1./ (1 + (mui2 + muj2 - muij2)/y/(1-mui2-muj2-muk2) );
    return 16. / 27 / sqrt(mij2) * p_cf->Coupling(scale, 0) * LDME * value / (sij - mij2) * JFF(y, mui2, muj2, muk2, muij2);
  }
}

double LF_FFV_Quarkonia_FF::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  const double mij = p_ms->Mass(m_flavs[0]); // mass of heavy quark
  const double preF = ( 16. / 27 / mij * p_cf->MaxCoupling(0) )  * p_cf->Coupling(sqr(3*mij),0)/(2*M_PI) * pow(0.82,2);
  m_zmin = zmin; 
  m_zmax = zmax;
  return  preF * 236. / 100. / (8*sqr(mij));
}

double LF_FFV_Quarkonia_FF::OverEstimated(const double z, const double y) {
  const double mij = p_ms->Mass(m_flavs[0]);
  const double preF = ( 16. / 27 / mij * p_cf->MaxCoupling(0) )  * p_cf->Coupling(sqr(3*mij),0)/(2*M_PI) * pow(0.82,2);
  return preF * 236. / 100. / (8*sqr(mij));
}

double LF_FFV_Quarkonia_FF::Z() {
  return 1. - (1. - m_zmin) * pow((1. - m_zmax) / (1. - m_zmin), ATOOLS::ran->Get());
}

double LF_FFV_Quarkonia_FI::operator()(const double zz, const double y,
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

double LF_FFV_Quarkonia_FI::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  m_zmin = zmin;
  m_zmax = zmax;
  m_Jmax = 5.;
  return 16. / 27 * pow(0.82, 3) / Flavour(kf_c).Mass() *
         sqr(p_cf->MaxCoupling(0)) * (0.5) * log((1. - zmin) / (1. - zmax)) *
         m_Jmax;
}

double LF_FFV_Quarkonia_FI::OverEstimated(const double z, const double y) {
  return 16. / 27 * pow(0.82, 3) / Flavour(kf_c).Mass() *
         sqr(p_cf->MaxCoupling(0)) * (0.5) / (1. - z) * m_Jmax;
}

double LF_FFV_Quarkonia_FI::Z() {
  return 1. -
         (1. - m_zmin) * pow((1. - m_zmax) / (1. - m_zmin), ATOOLS::ran->Get());
}

double LF_FFV_Quarkonia_IF::operator()(const double zz, const double y,
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

double LF_FFV_Quarkonia_IF::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  m_zmin = zmin;
  m_zmax = zmax;
  m_Jmax = m_flavs[0].Kfcode() < 3 ? 5. : 1.;
  return 16. / 27 * pow(0.82, 3) / Flavour(kf_c).Mass() *
         (2.0 * p_cf->MaxCoupling(0) * 2. + 0.5 * p_cf->MaxCoupling(1)) *
         log((1. - zmin) / (1. - zmax)) * m_Jmax;
}

double LF_FFV_Quarkonia_IF::OverEstimated(const double z, const double y) {
  return 16. / 27 * pow(0.82, 3) / Flavour(kf_c).Mass() *
         (2.0 * p_cf->MaxCoupling(0) * 2. + 0.5 * p_cf->MaxCoupling(1)) /
         (1. - z) * m_Jmax;
}

double LF_FFV_Quarkonia_IF::Z() {
  return 1. -
         (1. - m_zmin) * pow((1. - m_zmax) / (1. - m_zmin), ATOOLS::ran->Get());
}

double LF_FFV_Quarkonia_II::operator()(const double zz, const double y,
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

double LF_FFV_Quarkonia_II::OverIntegrated(const double zmin, const double zmax,
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

double LF_FFV_Quarkonia_II::OverEstimated(const double z, const double y) {
  return (4.0 * p_cf->MaxCoupling(0) + 0.5 * p_cf->MaxCoupling(1)) / (1. - z) *
         m_Jmax;
  return 16. / 27 * pow(0.82, 3) / Flavour(kf_c).Mass() *
         sqr(p_cf->MaxCoupling(0)) * (0.5) / (1. - z) * m_Jmax;
}

double LF_FFV_Quarkonia_II::Z() {
  return 1. -
         (1. - m_zmin) * pow((1. - m_zmax) / (1. - m_zmin), ATOOLS::ran->Get());
}

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
    msg_Out() << "Called --> " << METHOD << std::endl;
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

DECLARE_GETTER(LF_FFV_Quarkonia_FF, "FFV_Quarkonia", SF_Lorentz, SF_Key);

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
    case cstp::FI:
      return new LF_FFV_Quarkonia_FI(args);
    case cstp::IF:
      return new LF_FFV_Quarkonia_IF(args);
    case cstp::II:
      return new LF_FFV_Quarkonia_II(args);
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



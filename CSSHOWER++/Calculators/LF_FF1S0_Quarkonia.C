#include "CSSHOWER++/Showers/Splitting_Function_Base.H"
#include "ATOOLS/Phys/LDME.H"

namespace CSSHOWER {

class LF_FF1S0_Quarkonia_FF : public SF_Lorentz {
public:
  inline LF_FF1S0_Quarkonia_FF(const SF_Key &key) : SF_Lorentz(key) {}

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_FF1S0_Quarkonia_FI : public SF_Lorentz {

protected:
  double m_Jmax;

public:
  inline LF_FF1S0_Quarkonia_FI(const SF_Key &key) : SF_Lorentz(key) {}

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_FF1S0_Quarkonia_IF : public SF_Lorentz {

protected:
  double m_Jmax;

public:
  inline LF_FF1S0_Quarkonia_IF(const SF_Key &key) : SF_Lorentz(key) {}

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_FF1S0_Quarkonia_II : public SF_Lorentz {

protected:
  double m_Jmax;

public:
  inline LF_FF1S0_Quarkonia_II(const SF_Key &key) : SF_Lorentz(key) {}

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_F1S0F_Quarkonia_FF : public SF_Lorentz {
public:
  inline LF_F1S0F_Quarkonia_FF(const SF_Key &key) : SF_Lorentz(key) {}

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_F1S0F_Quarkonia_FI : public SF_Lorentz {

protected:
  double m_Jmax;

public:
  inline LF_F1S0F_Quarkonia_FI(const SF_Key &key) : SF_Lorentz(key) {}

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

double LF_FF1S0_Quarkonia_FF::operator()(const double z, const double y,
                                       const double eta, const double scale,
                                       const double Q2) {
  // again this is c -> c eta_c
  double mi = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true); // works with the mapping c -> c J/psi
  double mj = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  double mk = ATOOLS::Flavour(m_flspec.Kfcode()).Mass(true);
  double mij = sqr(ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true));
  double mui2 = sqr(mi) / Q2, muj2 = sqr(mj) / Q2, muk2 = sqr(mk) / Q2, muij2 = sqr(mij) / Q2;
  const double sij = y * (Q2 - sqr(mi) - sqr(mj) - sqr(mk)) + (sqr(mi) + sqr(mj));
  const double M = mi + mij;
  const double ri = mi / M;
  const double rij = mij / M; 
  const double den = sij - sqr(rij*M);
  double value = 0;
  value += sqr(cube(M))/cube(den) * ( -4*ri*rij*sqr(1-rij*(1-z)) );
  value += sqr(sqr(M))/sqr(den) * ( -(1-rij*(1-z)))*( 2*(1-2*ri) - (3-4*ri+4*sqr(ri))*(1-z) + rij*(1-2*ri)*sqr(1-z) );
  value += sqr(M)/(den) * z * sqr( 1 + ri*z );
  value *= ri*cube(rij)/sqr(1-rij*(1-z));
  value *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
  value *= 1. / (1 + sqr( 1 - z) * sqr(mi) / scale + sqr(z) * sqr(mj) / scale);
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *= 4.0/27/cube(mi)  * (m_flavs[2].StrongCharge() == 0 ?  1. : 1.0/48);
  return prefactor * sqr(p_cf->Coupling(scale, 0)) * value * JFF(y, mui2, muj2, muk2, muij2); //

  // same mass case OLD
  //   double sameflavour(0.);
  //   sameflavour = 1. / sqr(sij - mij2) * (
  //      (sij+3*mij2) * (sij-5*mij2) 
  //     -(sij-mij2)*(sij-sqr(sqrt(mi2)+sqrt(mj2)))
  //     +4*(1-z)/(1+z)*(sij-mij2)*((sij-mij2) - (1-z)*(sij-3*mij2))
  //     +4*z*sqr(1-z)*sqr((sij-mij2)/(1+z)) 
  //   );
  //   sameflavour /= (sij - mij2);
  //   sameflavour *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
  //   sameflavour *= 1. / (1 + sqr( 1 - z) * mi2 / scale + sqr(z) * mj2 / scale);
  //   return sqr(p_cf->Coupling(scale, 0)) *  sameflavour * JFF(y, mui2, muj2, muk2, muij2);

}

double LF_FF1S0_Quarkonia_FF::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  double mi = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true); // works with the mapping c -> c J/psi
  double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double ri = mi / (mi + mij);
  const double rij = mij / (mi + mij);
  m_zmin = zmin; 
  m_zmax = zmax;
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *= m_flavs[2].StrongCharge() == 0 ?  4.0/27/cube(mi) : 4.0/27/cube(mi) * 1.0/48;
  return prefactor * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.1 * (zmax - zmin);
}

double LF_FF1S0_Quarkonia_FF::OverEstimated(const double z, const double y) {
  double mi = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double ri = mi / (mi + mij);
  const double rij = mij / (mi + mij);
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *= m_flavs[2].StrongCharge() == 0 ?  4.0/27/cube(mi) : 4.0/27/cube(mi) * 1.0/48;
  return prefactor * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.1 ;
}

double LF_FF1S0_Quarkonia_FF::Z() {
  return m_zmin + (m_zmax - m_zmin) * ATOOLS::ran->Get();
}

double LF_FF1S0_Quarkonia_FI::operator()(const double z, const double y,
                                       const double eta, const double scale,
                                       const double Q2) {
  // assuming the variables work in the same way as PHASIC::ClusterFIDipole
  double mi  = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true); // works with the mapping c -> c J/psi
  double mj  = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  double ma  = ATOOLS::Flavour(m_flspec.Kfcode()).Mass(true);
  double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  double mui2 = sqr(mi) / Q2, muj2 = sqr(mj) / Q2, muk2 = sqr(ma) / Q2, muij2 = sqr(mij) / Q2;
  const double yt = ((Q2 - sqr(ma) - sqr(mij)) / (Q2 - sqr(ma) - sqr(mi) - sqr(mj)) - (1-y))/(1-y);
  const double sij = (sqr(mi) + sqr(mj))*(1+yt) - yt * (Q2 - sqr(ma));
  const double M = mi + mij;
  const double ri = mi / M;
  const double rij = mij / M;
  const double den = sij - sqr(rij*M);
  double value = 0;
  value += sqr(cube(M))/cube(den) * ( -4*ri*rij*sqr(1-rij*(1-z)) );
  value += sqr(sqr(M))/sqr(den) * ( -(1-rij*(1-z)))*( 2*(1-2*ri) - (3-4*ri+4*sqr(ri))*(1-z) + rij*(1-2*ri)*sqr(1-z) );
  value += sqr(M)/(den) * z * sqr( 1 + ri*z );
  value *= ri*cube(rij)/sqr(1-rij*(1-z));
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *= 4.0 / 27 / cube(mi) * (m_flavs[2].StrongCharge() == 0 ? 1. :  1.0 / 48);
  return prefactor * sqr(p_cf->Coupling(scale, 0)) * value * JFI(y, eta, scale);

  // double mi2 =
  //     sqr(p_ms->Mass(m_flavs[1])); // works with the mapping c -> J/Psi(1S) c
  // double mj2 = sqr(p_ms->Mass(m_flavs[2]));
  // double ma2 = sqr(p_ms->Mass(m_flspec));
  // double mij2 = p_ms->Mass2(m_flavs[0]);
  // const double sij = -y / (1 - y) * (Q2 - ma2) + (mi2 + mj2) / (1 - y);
  // double muij2 = mij2 / sij;
  // double mui2 = mi2 / sij, muj2 = mj2 / sij, mua2 = ma2 / sij;
  // const double newscale = sqrt(scale);
  //   // the massive case
  //   double value =
  //       1. / sqr(sqr(1 - muij2)) *
  //       ((1 - 2 * muij2 - 47 * sqr(muij2)) -
  //        z * (1 - muij2) * (1 - sqr(sqrt(mui2) + sqrt(muj2))) +
  //        4 * (z * (1 - z)) / (2 - z) * (1 - mui2) -
  //        4 * (8 - 7 * z - 5 * z * z) / (2 - z) * muij2 * (1 - muij2) +
  //        12 * (z * z * (1 - z)) / sqr(2 - z) * sqr(1 - muij2));
  //   const double LDME = 1. / sqrt(mi2) / sqr(sij) *
  //                       (m_flavs[2].IsOctetMeson() ? (1.5E-02) : pow(0.82, 3));
  //   return 16. / 27 * LDME * sqr(p_cf->Coupling(newscale, 0)) * value /
  //          sqr(sij) * JFI(y, eta, scale);
}

double LF_FF1S0_Quarkonia_FI::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  m_zmin = zmin;
  m_zmax = zmax;
  m_Jmax = 5.;
  double mi = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true); // works with the mapping c -> c J/psi
  double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double ri = mi / (mi + mij);
  const double rij = mij / (mi + mij);
  double prefactor = GetLDME(m_flavs[2].Kfcode())*2;
  prefactor *= m_flavs[2].StrongCharge() == 0 ?  4.0/27/cube(mi) : 4.0/27/cube(mi) * 1.0/48;
  return prefactor * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.1 * (zmax - zmin) * m_Jmax;
}

double LF_FF1S0_Quarkonia_FI::OverEstimated(const double z, const double y) {
  double mi = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true); // works with the mapping c -> c J/psi
  double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double ri = mi / (mi + mij);
  const double rij = mij / (mi + mij);
  double prefactor = GetLDME(m_flavs[2].Kfcode())*2;
  prefactor *= m_flavs[2].StrongCharge() == 0 ?  4.0/27/cube(mi) : 4.0/27/cube(mi) * 1.0/48;
  return prefactor * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.1 * m_Jmax;
}

double LF_FF1S0_Quarkonia_FI::Z() {
  return m_zmin + (m_zmax - m_zmin) * ATOOLS::ran->Get();
}

double LF_FF1S0_Quarkonia_IF::operator()(const double z, const double y,
                                       const double eta, const double scale,
                                       const double Q2) {
  // c --> c 1S0
  double ma  = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true); // works with the mapping c -> c J/psi
  double mj  = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  double mk  = ATOOLS::Flavour(m_flspec.Kfcode()).Mass(true);
  double mi = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  double muj2 = sqr(mj) / Q2, mui2 = sqr(mi) / Q2, muk2 = sqr(mk) / Q2, mua2 = sqr(ma) / Q2;
  const double zt = (1-y)/(z-y);
  const double saj = (sqr(ma) + sqr(mj))*(1-y/z) + (Q2-sqr(mk))*y/z;
  const double M = mi + ma;
  const double ri = mi / M;
  const double rij = ma / M;
  const double den = saj - sqr(ri*M);
  double value = 0;
  value += sqr(cube(M))/cube(den) * ( -4*ri*rij*sqr(1-rij*(1-zt)) );
  value += sqr(sqr(M))/sqr(den) * ( -(1-rij*(1-zt)))*( 2*(1-2*ri) - (3-4*ri+4*sqr(ri))*(1-zt) + rij*(1-2*ri)*sqr(1-zt) );
  value += sqr(M)/(den) * zt * sqr( 1 + ri*zt );
  value *= ri*cube(rij)/sqr(1-rij*(1-zt));
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *= 4.0 / 27 / cube(rij*M) * (m_flavs[2].StrongCharge() == 0 ? 1. : 1.0 / 48);
  return prefactor * sqr(p_cf->Coupling(scale, 0)) * value * JIF(z, y, eta, scale);
  // (aj) - (k) --> a - j - k
  // zz = x_{jk,a};  y = u 
  // double ma2 = sqr(p_ms->Mass(m_flavs[1])); // quark
  // double mj2 = sqr(p_ms->Mass(m_flavs[2])); // meson
  // double mk2 = sqr(p_ms->Mass(m_flspec));
  // double maj2 = p_ms->Mass2(m_flavs[0]); // quark
  // const double pkpa = Q2 / 2 / zz * (1 - y);
  // const double pkpj = 0.5 * (1 - zz) / zz * (mj2 + mk2 + ma2 - Q2);
  // const double taj =
  //     y * Q2 + ma2 * (1 - y) + mj2 + y * (mj2 + mk2 + 2 * pkpj); // (pa - pj)^2
  // const double z =
  //     pkpj / pkpa; // this is the momentum fraction of J/psi w.r.t. parent quark
  // double muaj2 = maj2 / taj;
  // double mua2 = ma2 / taj, muj2 = mj2 / taj, muk2 = mk2 / taj;
  // const double newscale = sqrt(scale);
  // // the massless case
  // if (muaj2 == 0. || mua2 == 0. || muj2 == 0.) {
  //   msg_Error() << "Cannot make massless quarkonia emission" << std::endl;
  //   exit(EXIT_FAILURE);
  // } else {
  //   // the massive case
  //   double value =
  //       1. / sqr(sqr(1 - muaj2)) *
  //       ((1 - 2 * muaj2 - 47 * sqr(muaj2)) -
  //        z * (1 - muaj2) * (1 - sqr(sqrt(fabs(mua2)) + sqrt(fabs(muj2)))) +
  //        4 * (z * (1 - z)) / (2 - z) * (1 - mua2) -
  //        4 * (8 - 7 * z - 5 * z * z) / (2 - z) * muaj2 * (1 - muaj2) +
  //        12 * (z * z * (1 - z)) / sqr(2 - z) * sqr(1 - muaj2));
  //   const double LDME = 1. / sqrt(ma2) / sqr(taj) *
  //                       (m_flavs[2].IsOctetMeson() ? (1.5E-02) : pow(0.82, 3));
  //   return 16. / 27 * LDME * sqr(p_cf->Coupling(newscale, 0)) * value /
  //          sqr(taj) * JFI(y, eta, scale);
  // }
}

double LF_FF1S0_Quarkonia_IF::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  m_zmin = zmin;
  m_zmax = zmax;
  m_Jmax = m_flavs[0].Kfcode() < 3 ? 5. : 1.;
  double mi = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true); // works with the mapping c -> c J/psi
  double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double ri = mi / (mi + mij);
  const double rij = mij / (mi + mij);
  double prefactor = GetLDME(m_flavs[2].Kfcode())*2;
  prefactor *= m_flavs[2].StrongCharge() == 0 ?  4.0/27/cube(mi) : 4.0/27/cube(mi) * 1.0/48;
  return prefactor * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.1 * (zmax - zmin) * m_Jmax;
}

double LF_FF1S0_Quarkonia_IF::OverEstimated(const double z, const double y) {
  // works with the mapping c -> c J/psi
  double mai = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  double ma = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double ri = mai / (mai + ma);
  const double rij = ma / (mai + ma);
  double prefactor = GetLDME(m_flavs[2].Kfcode())*2;
  prefactor *= m_flavs[2].StrongCharge() == 0 ?  4.0/27/cube(mai) : 4.0/27/cube(mai) * 1.0/48;
  return prefactor * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.1 * m_Jmax;;
}

double LF_FF1S0_Quarkonia_IF::Z() {
  return m_zmin + (m_zmax - m_zmin) * ATOOLS::ran->Get();
}

double LF_FF1S0_Quarkonia_II::operator()(const double z, const double y,
                                       const double eta, const double scale,
                                       const double Q2) {
  // works with the mapping [c aj] ->[c a] [eta_c j]
  double ma  = p_ms->Mass(m_flavs[0]); 
  double mj  = p_ms->Mass(m_flavs[2]);
  double mb  = p_ms->Mass(m_flspec);
  double mai = p_ms->Mass(m_flavs[1]);
  double mua2 = sqr(ma) / Q2,  mub2 = sqr(mb) / Q2, muai2 = sqr(mai) / Q2;
  const double zt = 1./(z+y);
  const double sab = (Q2 - sqr(mj) - (1-z)*(sqr(ma) + sqr(mb))) / z;
  const double saj = sqr(ma) + sqr(mj) - y * (sab - sqr(ma) - sqr(mb));
  const double M = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true) + ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double ri = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true) / M;
  const double rij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true) / M;
  const double den = saj - sqr(rij*M);
  double value = 0;
  value += sqr(cube(M))/cube(den) * ( -4*ri*rij*sqr(1-rij*(1-zt)) );
  value += sqr(sqr(M))/sqr(den) * ( -(1-rij*(1-zt)))*( 2*(1-2*ri) - (3-4*ri+4*sqr(ri))*(1-zt) + rij*(1-2*ri)*sqr(1-zt) );
  value += sqr(M)/(den) * zt * sqr( 1 + ri*zt );
  value *= ri*cube(rij)/sqr(1-rij*(1-zt));
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *= (4.0 / 27) / cube(rij*M) * ( m_flavs[2].StrongCharge() == 0 ? 1. : 1.0 / 48);
  return prefactor * sqr(p_cf->Coupling(scale, 0)) * value * JII(z, y, eta, scale);
  
  // // zz is effectively xjab and y vtildej
  // double ma2 = sqr(p_ms->Mass(m_flavs[1]));
  // double mj2 = sqr(p_ms->Mass(m_flavs[2]));
  // double mb2 = sqr(p_ms->Mass(m_flspec));
  // double maj2 = p_ms->Mass2(m_flavs[0]);
  // const double sab = (Q2 - mj2) / zz - (1. - zz) / zz * (ma2 + mb2);
  // const double z = 1 - (zz + y) + (ma2 + mb2 + mj2)/(sab - ma2 - mb2); // this is mom. fraction of J/psi w.r.t. to progenitor c
  // const double taj =  ma2 + mj2 - 2 * y * (sab - ma2 - mb2); // this is (pa - pj)^2
  // double muaj2 = maj2 / taj;
  // double mua2 = ma2 / taj, muj2 = mj2 / taj, mub2 = mb2 / taj;
  // const double newscale = scale;
  // // the massless case
  // if (muaj2 == 0. || mua2 == 0. || muj2 == 0.) {
  //   msg_Error() << "Cannot make massless quarkonia emission" << std::endl;
  //   exit(EXIT_FAILURE);
  // } else {
  //   // the massive case
  //   double value =
  //       1. / sqr(sqr(1 - muaj2)) *
  //        ( (1 - 2 * muaj2 - 47 * sqr(muaj2)) -
  //        z * (1 - muaj2) * (1 - sqr(sqrt(fabs(mua2)) + sqrt(fabs(muj2)))) +
  //        4 * (z * (1 - z)) / (2 - z) * (1 - mua2) -
  //        4 * (8 - 7 * z - 5 * z * z) / (2 - z) * muaj2 * (1 - muaj2) +
  //        12 * (z * z * (1 - z)) / sqr(2 - z) * sqr(1 - muaj2));
  //   const double LDME =
  //       (m_flavs[2].IsOctetMeson() ? 1.5E-02 / M_PI_2
  //                                  : 9. / 2 * M_PI * pow(0.82, 3));
  //   return value * JII(z, y, eta, scale);
  //   return 16. / 27 / 9 * sqr(p_cf->Coupling(newscale, 0)) / sqrt(ma2) * LDME *
  //          value / sqr(taj) * JII(z, y, eta, scale);
  // }
}

double LF_FF1S0_Quarkonia_II::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  const double mai  = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double ma = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double ri = mai / (ma + mai);
  const double rij = ma / (ma + mai);
  m_zmin = zmin; m_zmax = zmax; 
  m_Jmax = m_flavs[0].Kfcode()<3?5.:1.;
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *= m_flavs[2].StrongCharge() == 0 ? 4.0 / 27 / cube(ma) : 4.0 / 27 / cube(ma) * 1.0 / 48;
  return prefactor * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.4 * m_Jmax * (m_zmax - m_zmin);}

double LF_FF1S0_Quarkonia_II::OverEstimated(const double z, const double y) {
  const double mai  = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double ma = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double ri = mai / (ma + mai);
  const double rij = ma / (ma + mai);
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *= m_flavs[2].StrongCharge() == 0 ? 4.0 / 27 / cube(ma) : 4.0 / 27 / cube(ma) * 1.0 / 48;
  return prefactor * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.4 * m_Jmax;
}

double LF_FF1S0_Quarkonia_II::Z() {
  return m_zmin + (m_zmax - m_zmin) * ATOOLS::ran->Get();
}


double LF_F1S0F_Quarkonia_FF::operator()(const double z, const double y,
                                       const double eta, const double scale,
                                       const double Q2) {
  // again this is c -> eta_c c
  double mi = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true); // works with the mapping c -> c J/psi
  double mj = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  double mk = ATOOLS::Flavour(m_flspec.Kfcode()).Mass(true);
  double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  double mui2 = sqr(mi) / Q2, muj2 = sqr(mj) / Q2, muk2 = sqr(mk) / Q2, muij2 = sqr(mij) / Q2;
  const double sij = (sqr(mi) + sqr(mj))*(1-y) + (Q2 - sqr(mk)) * y;
  const double M = mj + mij;
  const double rj = mj / M;
  const double rij = mij / M; 
  const double den = sij - sqr(rij*M);
  double value = 0;
  value += sqr(cube(M))/cube(den) * ( -4*rj*rij*sqr(1-rij*(z)) );
  value += sqr(sqr(M))/sqr(den) * ( -(1-rij*(z)))*( 2*(1-2*rj) - (3-4*rj+4*sqr(rj))*(z) + rij*(1-2*rj)*sqr(z) );
  value += sqr(M)/(den) * (1-z) * sqr( 1 + rj*(1-z) );
  value *= rj*cube(rij)/sqr(1-rij*(z));
  value *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
  value *= 1. / (1 + sqr(z) * sqr(mi) / scale + sqr(1-z) * sqr(mj) / scale);
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *= m_flavs[2].StrongCharge() == 0 ?  4.0/27 /cube(mj) : 4.0/27/cube(mj) * 1.0/48;
  return prefactor * sqr(p_cf->Coupling(scale, 0)) * value * JFF(y, mui2, muj2, muk2, muij2); //
}

double LF_F1S0F_Quarkonia_FF::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  double mj = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true); // works with the mapping c -> c J/psi
  double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double rj = mj / (mj + mij);
  const double rij = mij / (mj + mij);
  m_zmin = zmin; 
  m_zmax = zmax;
  double prefactor = GetLDME(m_flavs[2].Kfcode())*2;
  prefactor *= m_flavs[2].StrongCharge() == 0 ?  4.0/27/cube(mj) : 4.0/27/cube(mj) * 1.0/48;
  return prefactor * sqr(p_cf->MaxCoupling(0)) * rj * cube(rij) / sqr(sqr(1-rij)) * 0.1 * (zmax - zmin);
}

double LF_F1S0F_Quarkonia_FF::OverEstimated(const double z, const double y) {
  double mj = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double rj = mj / (mj + mij);
  const double rij = mij / (mj + mij);
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *= m_flavs[2].StrongCharge() == 0 ?  4.0/27/cube(mj) : 4.0/27/cube(mj) * 1.0/48;
  return prefactor * sqr(p_cf->MaxCoupling(0)) * rj * cube(rij) / sqr(sqr(1-rij)) * 0.1 ;
}

double LF_F1S0F_Quarkonia_FF::Z() {
  return m_zmin + (m_zmax - m_zmin) * ATOOLS::ran->Get();
}

double LF_F1S0F_Quarkonia_FI::operator()(const double z, const double y,
                                       const double eta, const double scale,
                                       const double Q2) {
  double mi  = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true); // works with the mapping c -> c J/psi
  double mj  = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  double ma  = ATOOLS::Flavour(m_flspec.Kfcode()).Mass(true);
  double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  double mui2 = sqr(mi) / Q2, muj2 = sqr(mj) / Q2, muk2 = sqr(ma) / Q2, muij2 = sqr(mij) / Q2;
  const double yt = ((Q2 - sqr(mij) - sqr(ma))/(Q2 - sqr(mi) - sqr(mj) - sqr(ma)) - (1-y))/(1-y);
  const double sij = (sqr(mi) + sqr(mj))*(1+yt) - yt * (Q2 - sqr(ma));
  const double M = mj + mij;
  const double rj = mj / M;
  const double rij = mij / M;
  const double den = sij - sqr(rij*M);
  double value = 0;
  value += sqr(cube(M))/cube(den) * ( -4*rj*rij*sqr(1-rij*(z)) );
  value += sqr(sqr(M))/sqr(den) * ( -(1-rij*(z)))*( 2*(1-2*rj) - (3-4*rj+4*sqr(rj))*(z) + rij*(1-2*rj)*sqr(z) );
  value += sqr(M)/(den) * (1-z) * sqr( 1 + rj*(1-z) );
  value *= rj*cube(rij)/sqr(1-rij*(z));
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *= m_flavs[2].StrongCharge() == 0 ? 4.0 / 27 / cube(mj) : 4.0 / 27 / cube(mj) * 1.0 / 48;
  return prefactor * sqr(p_cf->Coupling(scale, 0)) * value * JFI(y, eta, scale);
}

double LF_F1S0F_Quarkonia_FI::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  m_zmin = zmin;
  m_zmax = zmax;
  m_Jmax = 5.;
  double mj = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double rj = mj / (mj + mij);
  const double rij = mij / (mj + mij);
  double prefactor = GetLDME(m_flavs[2].Kfcode())*2;
  prefactor *= m_flavs[2].StrongCharge() == 0 ?  4.0/27/cube(mj) : 4.0/27/cube(mj) * 1.0/48;
  return prefactor * sqr(p_cf->MaxCoupling(0)) * rj * cube(rij) / sqr(sqr(1-rij)) * 0.1 * (zmax - zmin) * m_Jmax;
}

double LF_F1S0F_Quarkonia_FI::OverEstimated(const double z, const double y) {
  double mj = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double rj = mj / (mj + mij);
  const double rij = mij / (mj + mij);
  double prefactor = GetLDME(m_flavs[2].Kfcode())*2;
  prefactor *= m_flavs[2].StrongCharge() == 0 ?  4.0/27/cube(mj) : 4.0/27/cube(mj) * 1.0/48;
  return prefactor * sqr(p_cf->MaxCoupling(0)) * rj * cube(rij) / sqr(sqr(1-rij)) * 0.1 * m_Jmax;
}

double LF_F1S0F_Quarkonia_FI::Z() {
  return m_zmin + (m_zmax - m_zmin) * ATOOLS::ran->Get();
}

DECLARE_GETTER(LF_FF1S0_Quarkonia_FF, "FF1S0_Quarkonia", SF_Lorentz, SF_Key);

SF_Lorentz *ATOOLS::Getter<SF_Lorentz, SF_Key, LF_FF1S0_Quarkonia_FF>::operator()(
    const Parameter_Type &args) const {
  if (args.m_col < 0)
    return NULL;
  if ((args.m_mode == 0 && args.p_v->in[0].IntSpin() == 1 &&
       args.p_v->in[1].IntSpin() == 1 && args.p_v->in[2].IntSpin() == 0) ||
      (args.m_mode == 1 && args.p_v->in[0].IntSpin() == 1 &&
       args.p_v->in[1].IntSpin() == 0 && args.p_v->in[2].IntSpin() == 1)) {
    switch (args.m_type) {
    case cstp::FF:
      return new LF_FF1S0_Quarkonia_FF(args);
    case cstp::FI:
      return new LF_FF1S0_Quarkonia_FI(args);
    case cstp::IF:
      return new LF_FF1S0_Quarkonia_IF(args);
    case cstp::II:
      return new LF_FF1S0_Quarkonia_II(args);
    case cstp::none:
      break;
    }
  }
  if ((args.m_mode == 0 && args.p_v->in[0].IntSpin() == 1 &&
      args.p_v->in[1].IntSpin() == 0 && args.p_v->in[2].IntSpin() == 1) ||
      (args.m_mode == 1 && args.p_v->in[0].IntSpin() == 1 &&
      args.p_v->in[2].IntSpin() == 0 && args.p_v->in[1].IntSpin() == 1)) {
    switch (args.m_type) {
    case cstp::FF: return new LF_F1S0F_Quarkonia_FF(args);
    case cstp::FI: return new LF_F1S0F_Quarkonia_FI(args);
    case cstp::IF: return NULL;
    case cstp::II: return NULL;
    case cstp::none:
      break;
    }
  }
  return NULL;
}

void ATOOLS::Getter<SF_Lorentz, SF_Key, LF_FF1S0_Quarkonia_FF>::PrintInfo(
    std::ostream &str, const size_t width) const {
  str << "FF1S0 lorentz functions";
}



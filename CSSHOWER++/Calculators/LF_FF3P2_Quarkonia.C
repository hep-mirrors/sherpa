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

class LF_F3P2F_Quarkonia_FF : public SF_Lorentz {
public:
  inline LF_F3P2F_Quarkonia_FF(const SF_Key &key) : SF_Lorentz(key) {}

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_F3P2F_Quarkonia_FI : public SF_Lorentz {

protected:
  double m_Jmax;

public:
  inline LF_F3P2F_Quarkonia_FI(const SF_Key &key) : SF_Lorentz(key) {}

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
  const double den = sij - sqr(rij * M);
  double value = 0.0;
  value = sqr(sqr(sqr(M)))/sqr(sqr(den))*( 320 * sqr(ri) * cube(rij) * sqr(sqr(1 - rij*(1-z))) );
  value += sqr(cube(M))/cube(den) * 8 * ri * sqr(rij) * cube(1 - rij * (1 - z)) * (
    2 * (4 + 13 * ri) - (1 + 70 * ri - 26 * sqr(ri)) * (1 - z) - (7 + 8 * ri) * rij * sqr(1 - z)
  );
  value += sqr(sqr(M))/sqr(den) * (-4) * sqr(rij) * sqr( 1 - rij * (1 - z)) * ( 
    4 * (1 + 4 * ri) - (7 + 12 * ri - 32 * sqr(ri)) * (1 - z) 
    + 2 * (1 + 13 * ri - 26 * sqr(ri) + 8 * cube(ri)) * sqr(1 - z) + (1 - 30 * ri - 5 * sqr(ri) + 4 * cube(ri)) * cube(1 - z)
  );
  value += sqr(M)/(den) * 4 * sqr(rij) * z * (
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
  return prefactor * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.1561 * (m_zmax - m_zmin);
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

double LF_FF3P2_Quarkonia_FI::operator()(const double z, const double y,
                                       const double eta, const double scale,
                                       const double Q2) {
  double mi  = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true); // works with the mapping c -> c J/psi
  double mj  = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  double ma  = ATOOLS::Flavour(m_flspec.Kfcode()).Mass(true);
  double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  double mui2 = sqr(mi) / Q2, muj2 = sqr(mj) / Q2, muk2 = sqr(ma) / Q2, muij2 = sqr(mij) / Q2;
  // const double sij = (Q2*(1-y) + sqr(mij) )/y;
  const double sij =  y * (Q2 - sqr(mi) - sqr(mj) - sqr(ma)) + (sqr(mi) + sqr(mj));
  const double M = mi + mij;
  const double ri = mi / M;
  const double rij = mij / M;
  const double den = sij - sqr(rij * M);
  double value = 0;
  value = sqr(sqr(sqr(M)))/sqr(sqr(den))*( 320 * sqr(ri) * cube(rij) * sqr(sqr(1 - rij*(1-z))) );
  value += sqr(cube(M))/cube(den) * 8 * ri * sqr(rij) * cube(1 - rij * (1 - z)) * (
    2 * (4 + 13 * ri) - (1 + 70 * ri - 26 * sqr(ri)) * (1 - z) - (7 + 8 * ri) * rij * sqr(1 - z)
  );
  value += sqr(sqr(M))/sqr(den) * (-4) * sqr(rij) * sqr( 1 - rij * (1 - z)) * ( 
    4 * (1 + 4 * ri) - (7 + 12 * ri - 32 * sqr(ri)) * (1 - z) 
    + 2 * (1 + 13 * ri - 26 * sqr(ri) + 8 * cube(ri)) * sqr(1 - z) + (1 - 30 * ri - 5 * sqr(ri) + 4 * cube(ri)) * cube(1 - z)
  );
  value += sqr(M)/(den) * 4 * sqr(rij) * z * (
    2 - 4 * (1 - 2 * ri) * (1 - z)
    + (5 - 8 * ri + 12 * sqr(ri)) * sqr(1 - z)
    - 2 * (1 - 2 * ri) * (3 + 2 * sqr(ri)) * cube(1 - z) 
    + (3 - 12 * ri + 12 * sqr(ri) + 2 * sqr(sqr(ri))) * sqr(sqr((1 - z)))
  );
  value *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
  value *= 1. / (1 + sqr( 1 - z) * sqr(mi) / scale + sqr(z) * sqr(mj) / scale);
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *= 4.0 / 27 / pow(mi,5);
  return prefactor * sqr(p_cf->Coupling(scale, 0)) * value * JFI(y, eta, scale);
}

double LF_FF3P2_Quarkonia_FI::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  m_zmin = zmin;
  m_zmax = zmax;
  m_Jmax = 5.;
  const double mi = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double rij = mij / (mi + mij);
  const double ri = mi / (mi + mij);
  const double prefactor = 4./27*GetLDME(m_flavs[2].Kfcode())/pow(mi,5);
  return prefactor * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.1561 * (m_zmax - m_zmin);
}

double LF_FF3P2_Quarkonia_FI::OverEstimated(const double z, const double y) {
  const double mi = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double rij = mij / (mi + mij);
  const double ri = mi / (mi + mij);
  const double prefactor = 4./27*GetLDME(m_flavs[2].Kfcode())/pow(mi,5);
  return prefactor * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.1561 * m_Jmax;
}

double LF_FF3P2_Quarkonia_FI::Z() {
  return m_zmin + (m_zmax - m_zmin) * ATOOLS::ran->Get();
}

double LF_FF3P2_Quarkonia_IF::operator()(const double z, const double y,
                                       const double eta, const double scale,
                                       const double Q2) {
  // c --> c 3P0
  // morally y is xij,a
  double mi  = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true); // works with the mapping c -> c J/psi
  double mj  = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  double mk  = ATOOLS::Flavour(m_flspec.Kfcode()).Mass(true);
  double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  double mui2 = sqr(mi) / Q2, muj2 = sqr(mj) / Q2, muk2 = sqr(mk) / Q2, muij2 = sqr(mij) / Q2;
  // const double sij = (Q2*(1-y) + sqr(mij) )/y;
  const double tai = y * (Q2 - sqr(mi) - sqr(mj) - sqr(mk)) + (sqr(mi) + sqr(mj));
  const double M = mi + mij;
  const double ri = mi / M;
  const double rij = mij / M;
  const double den = tai - sqr(rij*M);
  double value = 0;
  value = sqr(sqr(sqr(M)))/sqr(sqr(den))*( 320 * sqr(ri) * cube(rij) * sqr(sqr(1 - rij*(1-z))) );
  value += sqr(cube(M))/cube(den) * 8 * ri * sqr(rij) * cube(1 - rij * (1 - z)) * (
    2 * (4 + 13 * ri) - (1 + 70 * ri - 26 * sqr(ri)) * (1 - z) - (7 + 8 * ri) * rij * sqr(1 - z)
  );
  value += sqr(sqr(M))/sqr(den) * (-4) * sqr(rij) * sqr( 1 - rij * (1 - z)) * ( 
    4 * (1 + 4 * ri) - (7 + 12 * ri - 32 * sqr(ri)) * (1 - z) 
    + 2 * (1 + 13 * ri - 26 * sqr(ri) + 8 * cube(ri)) * sqr(1 - z) + (1 - 30 * ri - 5 * sqr(ri) + 4 * cube(ri)) * cube(1 - z)
  );
  value += sqr(M)/(den) * 4 * sqr(rij) * z * (
    2 - 4 * (1 - 2 * ri) * (1 - z)
    + (5 - 8 * ri + 12 * sqr(ri)) * sqr(1 - z)
    - 2 * (1 - 2 * ri) * (3 + 2 * sqr(ri)) * cube(1 - z) 
    + (3 - 12 * ri + 12 * sqr(ri) + 2 * sqr(sqr(ri))) * sqr(sqr((1 - z)))
  );
  // value *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
  // value *= 1. / (1 + sqr( 1 - z) * sqr(mi) / scale + sqr(z) * sqr(mj) / scale);
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *=  4.0 / 27 / pow(mi,5);
  return prefactor * sqr(p_cf->Coupling(scale, 0)) * value * JFI(y, eta, scale);
}

double LF_FF3P2_Quarkonia_IF::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  m_zmin = zmin;
  m_zmax = zmax;
  m_Jmax = m_flavs[0].Kfcode() < 3 ? 5. : 1.;
  const double mi = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double rij = mij / (mi + mij);
  const double ri = mi / (mi + mij);
  const double prefactor = 4./27*GetLDME(m_flavs[2].Kfcode())/pow(mi,5);
  return prefactor * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.1561 * (m_zmax - m_zmin) * m_Jmax;
  // extra factor 2 is spurios but necessary
}

double LF_FF3P2_Quarkonia_IF::OverEstimated(const double z, const double y) {
  const double mi = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double rij = mij / (mi + mij);
  const double ri = mi / (mi + mij);
  const double prefactor = 4./27*GetLDME(m_flavs[2].Kfcode())/pow(mi,5);
  return prefactor * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.1561 * m_Jmax;
}

double LF_FF3P2_Quarkonia_IF::Z() {
  return m_zmin + (m_zmax - m_zmin) * ATOOLS::ran->Get();
}

double LF_FF3P2_Quarkonia_II::operator()(const double z, const double y,
                                       const double eta, const double scale,
                                       const double Q2) {
  double ma  = p_ms->Mass(m_flavs[1]); 
  double mi  = p_ms->Mass(m_flavs[2]);
  double mb  = p_ms->Mass(m_flspec);
  double mai = p_ms->Mass(m_flavs[0]);
  double mua2 = sqr(ma) / Q2,  mub2 = sqr(mb) / Q2, muai2 = sqr(mai) / Q2;
  // const double tai = sqr(ma) + sqr(mi) - y / z * (Q2 - sqr(mb) - sqr(ma) - sqr(mi));
  const double tai = y * (Q2 - sqr(ma) - sqr(mi) - sqr(mb)) + (sqr(ma) + sqr(mi));
  const double M = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true) + ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double ri = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true) / M;
  const double rij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true) / M;
  const double den = tai - sqr(rij*M);
  double value = 0;
  value = sqr(sqr(sqr(M)))/sqr(sqr(den))*( 320 * sqr(ri) * cube(rij) * sqr(sqr(1 - rij*(1-z))) );
  value += sqr(cube(M))/cube(den) * 8 * ri * sqr(rij) * cube(1 - rij * (1 - z)) * (
    2 * (4 + 13 * ri) - (1 + 70 * ri - 26 * sqr(ri)) * (1 - z) - (7 + 8 * ri) * rij * sqr(1 - z)
  );
  value += sqr(sqr(M))/sqr(den) * (-4) * sqr(rij) * sqr( 1 - rij * (1 - z)) * ( 
    4 * (1 + 4 * ri) - (7 + 12 * ri - 32 * sqr(ri)) * (1 - z) 
    + 2 * (1 + 13 * ri - 26 * sqr(ri) + 8 * cube(ri)) * sqr(1 - z) + (1 - 30 * ri - 5 * sqr(ri) + 4 * cube(ri)) * cube(1 - z)
  );
  value += sqr(M)/(den) * 4 * sqr(rij) * z * (
    2 - 4 * (1 - 2 * ri) * (1 - z)
    + (5 - 8 * ri + 12 * sqr(ri)) * sqr(1 - z)
    - 2 * (1 - 2 * ri) * (3 + 2 * sqr(ri)) * cube(1 - z) 
    + (3 - 12 * ri + 12 * sqr(ri) + 2 * sqr(sqr(ri))) * sqr(sqr((1 - z)))
  );
  // value *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
  // value *= 1. / (1 + sqr( 1 - z) * sqr(mi) / scale + sqr(z) * sqr(mj) / scale);
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *= 4.0 / 27 / pow(ma,5);
  return prefactor * sqr(p_cf->Coupling(scale, 0)) * value * JII(z, y, eta, scale);
}

double LF_FF3P2_Quarkonia_II::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  m_zmin = zmin;
  m_zmax = zmax;
  m_Jmax = m_flavs[0].Kfcode() < 3 ? 5. : 1.;
  const double mi = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double rij = mij / (mi + mij);
  const double ri = mi / (mi + mij);
  const double prefactor = 4./27*GetLDME(m_flavs[2].Kfcode())/pow(mi,5);
  return prefactor * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.1561 * (m_zmax - m_zmin) * m_Jmax;
}

double LF_FF3P2_Quarkonia_II::OverEstimated(const double z, const double y) {
  const double mi = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double rij = mij / (mi + mij);
  const double ri = mi / (mi + mij);
  const double prefactor = 4./27*GetLDME(m_flavs[2].Kfcode())/pow(mi,5);
  return prefactor * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.1561 * m_Jmax;
}

double LF_FF3P2_Quarkonia_II::Z() {
  return m_zmin + (m_zmax - m_zmin) * ATOOLS::ran->Get();
}


double LF_F3P2F_Quarkonia_FF::operator()(const double z, const double y,
                                       const double eta, const double scale,
                                       const double Q2) {
  // again this is c -> eta_c c
  double mi = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true); // works with the mapping c -> c J/psi
  double mj = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  double mk = ATOOLS::Flavour(m_flspec.Kfcode()).Mass(true);
  double mij = sqr(ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true));
  double mui2 = sqr(mi) / Q2, muj2 = sqr(mj) / Q2, muk2 = sqr(mk) / Q2, muij2 = sqr(mij) / Q2;
  const double sij = y * (Q2 - sqr(mi) - sqr(mj) - sqr(mk)) + (sqr(mi) + sqr(mj));
  const double M = mj + mij;
  const double ri = mj / M;
  const double rij = mij / M; 
  const double den = den;
  double value = 0;
  value = sqr(sqr(sqr(M)))/sqr(sqr(den))*( 320 * sqr(ri) * cube(rij) * sqr(sqr(1 - rij*(z))) );
  value += sqr(cube(M))/cube(den) * 8 * ri * sqr(rij) * cube(1 - rij * (z)) * (
    2 * (4 + 13 * ri) - (1 + 70 * ri - 26 * sqr(ri)) * (z) - (7 + 8 * ri) * rij * sqr(z)
  );
  value += sqr(sqr(M))/sqr(den) * (-4) * sqr(rij) * sqr( 1 - rij * (z)) * ( 
    4 * (1 + 4 * ri) - (7 + 12 * ri - 32 * sqr(ri)) * (z) 
    + 2 * (1 + 13 * ri - 26 * sqr(ri) + 8 * cube(ri)) * sqr(z) + (1 - 30 * ri - 5 * sqr(ri) + 4 * cube(ri)) * cube(z)
  );
  value += sqr(M)/(den) * 4 * sqr(rij) * (1-z) * (
    2 - 4 * (1 - 2 * ri) * (z)
    + (5 - 8 * ri + 12 * sqr(ri)) * sqr(z)
    - 2 * (1 - 2 * ri) * (3 + 2 * sqr(ri)) * cube(z) 
    + (3 - 12 * ri + 12 * sqr(ri) + 2 * sqr(sqr(ri))) * sqr(sqr((z)))
  );
  value *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
  value *= 1. / (1 + sqr(1-z) * sqr(mi) / scale + sqr(z) * sqr(mj) / scale);
  double prefactor = GetLDME(m_flavs[1].Kfcode());
  prefactor *= 4.0/27 / pow(mj,5);
  return prefactor * sqr(p_cf->Coupling(scale, 0)) * value * JFF(y, mui2, muj2, muk2, muij2); //
}

double LF_F3P2F_Quarkonia_FF::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  double mj = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true); // works with the mapping c -> c J/psi
  double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double ri = mj / (mj + mij);
  const double rij = mij / (mj + mij);
  m_zmin = zmin; 
  m_zmax = zmax;
  double prefactor = GetLDME(m_flavs[1].Kfcode())*2;
  prefactor *=  4.0/27/pow(mj,5);
  return prefactor * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.1561 * (m_zmax - m_zmin);
}

double LF_F3P2F_Quarkonia_FF::OverEstimated(const double z, const double y) {
  double mj = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double ri = mj / (mj + mij);
  const double rij = mij / (mj + mij);
  double prefactor = GetLDME(m_flavs[1].Kfcode());
  prefactor *= 4.0/27/pow(mj,5);
  return prefactor * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.1561;
}

double LF_F3P2F_Quarkonia_FF::Z() {
  return m_zmin + (m_zmax - m_zmin) * ATOOLS::ran->Get();
}

double LF_F3P2F_Quarkonia_FI::operator()(const double z, const double y,
                                       const double eta, const double scale,
                                       const double Q2) {
  double mi  = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true); // works with the mapping c -> c J/psi
  double mj  = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  double ma  = ATOOLS::Flavour(m_flspec.Kfcode()).Mass(true);
  double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  double mui2 = sqr(mi) / Q2, muj2 = sqr(mj) / Q2, muk2 = sqr(ma) / Q2, muij2 = sqr(mij) / Q2;
  // morally y is x ij,a
  // const double sij = (Q2*(1-y) + sqr(mij) )/y;
  const double sij =  y * (Q2 - sqr(mi) - sqr(mj) - sqr(ma)) + (sqr(mi) + sqr(mj));
  const double M = mj + mij;
  const double ri = mj / M;
  const double rij = mij / M;
  const double den = den;
  double value = 0;
  value = sqr(sqr(sqr(M)))/sqr(sqr(den))*( 320 * sqr(ri) * cube(rij) * sqr(sqr(1 - rij*(z))) );
  value += sqr(cube(M))/cube(den) * 8 * ri * sqr(rij) * cube(1 - rij * (z)) * (
    2 * (4 + 13 * ri) - (1 + 70 * ri - 26 * sqr(ri)) * (z) - (7 + 8 * ri) * rij * sqr(z)
  );
  value += sqr(sqr(M))/sqr(den) * (-4) * sqr(rij) * sqr( 1 - rij * (z)) * ( 
    4 * (1 + 4 * ri) - (7 + 12 * ri - 32 * sqr(ri)) * (z) 
    + 2 * (1 + 13 * ri - 26 * sqr(ri) + 8 * cube(ri)) * sqr(z) + (1 - 30 * ri - 5 * sqr(ri) + 4 * cube(ri)) * cube(z)
  );
  value += sqr(M)/(den) * 4 * sqr(rij) * (1-z) * (
    2 - 4 * (1 - 2 * ri) * (z)
    + (5 - 8 * ri + 12 * sqr(ri)) * sqr(z)
    - 2 * (1 - 2 * ri) * (3 + 2 * sqr(ri)) * cube(z) 
    + (3 - 12 * ri + 12 * sqr(ri) + 2 * sqr(sqr(ri))) * sqr(sqr((z)))
  );
  value *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
  value *= 1. / (1 + sqr(1-z) * sqr(mi) / scale + sqr(z) * sqr(mj) / scale);
  double prefactor = GetLDME(m_flavs[1].Kfcode());
  prefactor *= 4.0 / 27 / pow(mj,5);
  return prefactor * sqr(p_cf->Coupling(scale, 0)) * value * JFI(y, eta, scale);
}

double LF_F3P2F_Quarkonia_FI::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  m_zmin = zmin;
  m_zmax = zmax;
  m_Jmax = 5.;
  double mj = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double ri = mj / (mj + mij);
  const double rij = mij / (mj + mij);
  double prefactor = GetLDME(m_flavs[1].Kfcode())*2;
  prefactor *= 4.0/27/pow(mj,5);
  return prefactor * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.1561 * (m_zmax - m_zmin) * m_Jmax;
}

double LF_F3P2F_Quarkonia_FI::OverEstimated(const double z, const double y) {
  double mj = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double ri = mj / (mj + mij);
  const double rij = mij / (mj + mij);
  double prefactor = GetLDME(m_flavs[1].Kfcode())*2;
  prefactor *= 4.0/27/pow(mj,5);
  return prefactor * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.1561 * m_Jmax;
}

double LF_F3P2F_Quarkonia_FI::Z() {
  return m_zmin + (m_zmax - m_zmin) * ATOOLS::ran->Get();
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
    case cstp::FF: return new LF_FF3P2_Quarkonia_FF(args);
    case cstp::FI: return new LF_FF3P2_Quarkonia_FI(args);
    case cstp::IF: return new LF_FF3P2_Quarkonia_IF(args);
    case cstp::II: return new LF_FF3P2_Quarkonia_II(args);
    case cstp::none:
      break;
    }
  }
  if ((args.m_mode == 0 && args.p_v->in[0].IntSpin() == 1 &&
       args.p_v->in[1].IntSpin() == 2 && args.p_v->in[2].IntSpin() == 1) ||
      (args.m_mode == 1 && args.p_v->in[0].IntSpin() == 1 &&
       args.p_v->in[2].IntSpin() == 2 && args.p_v->in[1].IntSpin() == 1)) {
    switch (args.m_type) {
    case cstp::FF: return new LF_F3P2F_Quarkonia_FF(args);
    case cstp::FI: return new LF_F3P2F_Quarkonia_FI(args);
    case cstp::IF: NULL;
    case cstp::II: NULL;
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



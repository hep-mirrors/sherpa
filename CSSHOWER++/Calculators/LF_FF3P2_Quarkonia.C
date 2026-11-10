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

namespace {
// Equal-flavour constituent expression: M0=2*mc; the external state keeps
// its physical mass. The common r*rbar^3=1/16 lives in OverIntegrated.
// hep-ph/9405348v1 Eqs. (18)--(35), normalized to (102)--(105).
double FF3P2Reduced(const double z, const double mc2, const double den)
{
  const double x = 1.-z, u = 1.-0.5*x, rho = 4.*mc2/den;
  const double c0 = 10.*sqr(sqr(u));
  const double c1 = -0.5*cube(u)*(x+6.)*(11.*x-7.);
  const double c2 = 0.25*sqr(u)*(59.*cube(x)-16.*sqr(x)+20.*x-48.);
  const double c3 = z*(sqr(sqr(x))+32.*sqr(x)+16.)/8.;
  return rho*(c3+rho*(c2+rho*(c1+rho*c0)))/sqr(sqr(u));
}
}


double LF_FF3P2_Quarkonia_FF::operator()(const double z, const double y,
                                        const double eta, const double scale,
                                        const double Q2) {
  if (m_flavs[0] == m_flavs[1]) {
    const double mi = m_flavs[1].Mass(true);
    const double mj = m_flavs[2].Mass(true);
    const double mk = m_flspec.Mass(true);
    const double mij = m_flavs[0].Mass(true);
    const double mi2 = sqr(mi), mij2 = sqr(mij);
    const double mui2 = mi2 / Q2, muj2 = sqr(mj) / Q2,
                 muk2 = sqr(mk) / Q2, muij2 = mij2 / Q2;
    const double sij = y * (Q2 - mi2 - sqr(mj) - sqr(mk)) + mi2 + sqr(mj);
    const double den = sij - mij2;
    double value = FF3P2Reduced(z, mi2, den);
    value *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
    value *= 1. / (1 + sqr( 1 - z) * mi2 / scale + sqr(z) * sqr(mj) / scale);
    return sqr(p_cf->Coupling(scale, 0)) * value * JFF(y, mui2, muj2, muk2, muij2);
  }
  // Unequal-flavour channel retains its legacy prescription.

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
  if (m_flavs[0] == m_flavs[1]) {
    const double mi = m_flavs[1].Mass(true);
    const double mij = m_flavs[0].Mass(true);
    const double ri = mi / (mi + mij);
    const double rij = mij / (mi + mij);
    m_zmin = zmin;
    m_zmax = zmax;
    // Paper normalization, including CF_QCD's C_F^2, parent polarization
    // averaging and dt/(2*pi*t). A quark has one colour connection.
    // All fixed factors common to evaluation and the bound live here.
    const double prefactor = 16.*M_PI/27 * (GetLDME(m_flavs[2].Kfcode())/5.)
        / pow(mi,5) * ri*cube(rij)
        * (m_flavs[2].StrongCharge() == 0 ? 1. : 1./48);
    return prefactor * OverEstimated(0.5, 0.) * (m_zmax - m_zmin);
  }
  // Unequal-flavour channel retains its legacy prescription.

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
  if (m_flavs[0] == m_flavs[1]) {
    // Analytic constituent-mass bound; proof in notes/check_ff_quarkonia_states.py.
    return 5.*sqr(p_cf->MaxCoupling(0));
  }
  // Unequal-flavour channel retains its legacy prescription.

  const double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double mi = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double rij = mij / (mi + mij);
  const double ri = mi / (mi + mij);
  return sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.1561;
}

double LF_FF3P2_Quarkonia_FF::Z() {
  return m_zmin + (m_zmax-m_zmin) * ATOOLS::ran->Get();
}

double LF_FF3P2_Quarkonia_FI::operator()(const double z, const double y,
                                       const double eta, const double scale,
                                       const double Q2) {
  if (m_flavs[0] == m_flavs[1]) {
    if (!(z > 0. && z < 1. && y > 0. && y < 1. && scale > 0.)) return 0.;
    const double mi = m_flavs[1].Mass(true);
    const double mj = m_flavs[2].Mass(true);
    const double mij = m_flavs[0].Mass(true);
    const double mi2 = sqr(mi), mij2 = sqr(mij);
    const double ma2 = p_ms->Mass2(m_flspec);
    // Q2 is the positive spacelike magnitude. y is the rescaling variable
    // after Sudakov's transformation of Kinematics_FI::GetY.
    const double sij = (mij2 + y*(Q2+ma2))/(1.-y);
    const double den = sij - mij2;
    if (IsBad(sij) || !(den > 0.) || sij < mi2/z+sqr(mj)/(1.-z)) return 0.;
    double value = FF3P2Reduced(z, mi2, den);
    value *= scale/(z*(1.-z)*den);
    return sqr(p_cf->Coupling(scale, 0)) * value * JFI(y, eta, scale);
  }
  // Unequal-flavour channel retains its legacy prescription.

  double mi  = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true); // works with the mapping c -> c J/psi
  double mj  = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  double ma  = ATOOLS::Flavour(m_flspec.Kfcode()).Mass(true);
  double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  double mui2 = sqr(mi) / Q2, muj2 = sqr(mj) / Q2, muk2 = sqr(ma) / Q2, muij2 = sqr(mij) / Q2;
  const double yt = ((Q2 - sqr(ma) - sqr(mij)) / (Q2 - sqr(ma) - sqr(mi) - sqr(mj)) - (1.0-y)) / (1.0-y);
  const double sij = (sqr(mi) + sqr(mj)) * (1.0 + yt) - yt * (Q2 - sqr(ma));
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
  // value *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
  // value *= 1. / (1 + sqr( 1 - z) * sqr(mi) / scale + sqr(z) * sqr(mj) / scale);
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *= 4.0 / 27 / pow(mi,5);
  return prefactor * sqr(p_cf->Coupling(scale, 0)) * value * JFI(y, eta, scale);
}

double LF_FF3P2_Quarkonia_FI::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  if (m_flavs[0] == m_flavs[1]) {
    m_zmin = std::max(zmin,0.);
    m_zmax = std::min(zmax,1.);
    m_Jmax = m_disablepdfrescaling ? 1. : 5.;
    if (m_zmax <= m_zmin) return 0.;
    const double mi = m_flavs[1].Mass(true);
    const double mij = m_flavs[0].Mass(true);
    const double ri = mi/(mi+mij), rij = mij/(mi+mij);
    // Per-quark paper normalization; LDME spin multiplicity removed where needed.
    const double prefactor = 16.*M_PI/27 * (GetLDME(m_flavs[2].Kfcode())/5.)
        / pow(mi,5) * ri*cube(rij)
        * (m_flavs[2].StrongCharge() == 0 ? 1. : 1./48);
    return prefactor * OverEstimated(0.5,0.) * (m_zmax-m_zmin);
  }
  // Unequal-flavour channel retains its legacy prescription.

  m_zmin = zmin;
  m_zmax = zmax;
  m_Jmax = 1.;
  const double mi = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double rij = mij / (mi + mij);
  const double ri = mi / (mi + mij);
  const double prefactor = 4./27*GetLDME(m_flavs[2].Kfcode())/pow(mi,5);
  return prefactor * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.1561 * (m_zmax - m_zmin);
}

double LF_FF3P2_Quarkonia_FI::OverEstimated(const double z, const double y) {
  if (m_flavs[0] == m_flavs[1]) {
    // Analytic constituent-mass bound; proof in notes/check_ff_quarkonia_states.py.
    return 5.*sqr(p_cf->MaxCoupling(0)) * m_Jmax;
  }
  // Unequal-flavour channel retains its legacy prescription.

  const double mi = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double rij = mij / (mi + mij);
  const double ri = mi / (mi + mij);
  const double prefactor = 4./27*GetLDME(m_flavs[2].Kfcode())/pow(mi,5);
  return prefactor * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.1561 * m_Jmax;
}

double LF_FF3P2_Quarkonia_FI::Z() {
  return m_zmin + (m_zmax-m_zmin) * ATOOLS::ran->Get();
}

double LF_FF3P2_Quarkonia_IF::operator()(const double z, const double y,
                                       const double eta, const double scale,
                                       const double Q2) {
  if (m_flavs[0] == m_flavs[1]) {
    if (!(z > 0. && z < 1. && y > 0. && y < 1. && scale > 0. && Q2 > 0.)) return 0.;
    const double mc2 = sqr(m_flavs[1].Mass(true));
    const double ma2 = p_ms->Mass2(m_flavs[0]);
    const double mj2 = sqr(m_flavs[2].Mass(true));
    const double ms2 = p_ms->Mass2(m_flspec);
    const double V = (Q2+ma2+mj2+ms2)*y/z-ma2-mj2;
    if (IsBad(V) || !(V > 0.) ||
        V < z*mj2/(1.-z)-z*ma2) return 0.;
    // Continue the FF reduced kernel with s=-V and z_FF=1/z.
    // The -z crossing factor follows FFV; JIF supplies its flux once.
    const double den = -V-mc2, zt = 1./z;
    const double crossed = -z*FF3P2Reduced(zt, mc2, den);
    // All current incoming-quark -> quark + quarkonium schemes have dV/dt=1/(1-z).
    return sqr(p_cf->Coupling(scale,0))*crossed*scale/((1.-z)*(V+mc2))
        * JIF(z,y,eta,scale);
  }
  // Unequal-flavour channel retains its legacy prescription.

  // c --> c 3P0
  // morally y is xij,a
  double mi  = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true); // works with the mapping c -> c J/psi
  double mj  = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  double mk  = ATOOLS::Flavour(m_flspec.Kfcode()).Mass(true);
  double ma = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  double mui2 = sqr(mi) / Q2, muj2 = sqr(mj) / Q2, muk2 = sqr(mk) / Q2, muij2 = sqr(ma) / Q2;
  const double zt = (1.0 - y) / (z - y);
  const double saj = (sqr(mi) + sqr(mj)) * (1.0 - y / z) + (Q2 - sqr(mk)) * (y / z);
  const double M = mi + ma;
  const double ri = mi / M;
  const double rij = ma / M;
  const double den = saj - sqr(rij * M);
  double value = 0;
  value = sqr(sqr(sqr(M)))/sqr(sqr(den))*( 320 * sqr(ri) * cube(rij) * sqr(sqr(1 - rij*(1-zt))) );
  value += sqr(cube(M))/cube(den) * 8 * ri * sqr(rij) * cube(1 - rij * (1 - zt)) * (
    2 * (4 + 13 * ri) - (1 + 70 * ri - 26 * sqr(ri)) * (1 - zt) - (7 + 8 * ri) * rij * sqr(1 - zt)
  );
  value += sqr(sqr(M))/sqr(den) * (-4) * sqr(rij) * sqr( 1 - rij * (1 - zt)) * ( 
    4 * (1 + 4 * ri) - (7 + 12 * ri - 32 * sqr(ri)) * (1 - zt) 
    + 2 * (1 + 13 * ri - 26 * sqr(ri) + 8 * cube(ri)) * sqr(1 - zt) + (1 - 30 * ri - 5 * sqr(ri) + 4 * cube(ri)) * cube(1 - zt)
  );
  value += sqr(M)/(den) * 4 * sqr(rij) * zt * (
    2 - 4 * (1 - 2 * ri) * (1 - zt)
    + (5 - 8 * ri + 12 * sqr(ri)) * sqr(1 - zt)
    - 2 * (1 - 2 * ri) * (3 + 2 * sqr(ri)) * cube(1 - zt) 
    + (3 - 12 * ri + 12 * sqr(ri) + 2 * sqr(sqr(ri))) * sqr(sqr((1 - zt)))
  );
  // value *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
  // value *= 1. / (1 + sqr( 1 - z) * sqr(mi) / scale + sqr(z) * sqr(mj) / scale);
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *=  4.0 / 27 / pow(rij * M,5);
  return prefactor * sqr(p_cf->Coupling(scale, 0)) * value * JIF(z, y, eta, scale);
}

double LF_FF3P2_Quarkonia_IF::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  if (m_flavs[0] == m_flavs[1]) {
    m_zmin = Max(0.,zmin); m_zmax = Min(1.,zmax);
    m_Jmax = m_disablepdfrescaling ? 1. : 5.;
    if (!(m_zmin > 0. && m_zmax > m_zmin)) return 0.;
    const double mc = m_flavs[1].Mass(true);
    const double prefactor = M_PI/27*(GetLDME(m_flavs[2].Kfcode())/5.)/pow(mc,5)
        * (m_flavs[2].StrongCharge() == 0 ? 1. : 1./48);
    return prefactor*OverEstimated(0.5,0.)/4.
        * (m_zmax-m_zmin)/(m_zmin*m_zmax);
  }
  // Unequal-flavour channel retains its legacy prescription.

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
  if (m_flavs[0] == m_flavs[1]) {
    // Common bound for the current evolution schemes, M_H >= 2*m_c,
    // and an incoming kinematic mass no larger than the constituent mass.
    return 256.*sqr(p_cf->MaxCoupling(0))*m_Jmax/sqr(z);
  }
  // Unequal-flavour channel retains its legacy prescription.

  const double mi = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double rij = mij / (mi + mij);
  const double ri = mi / (mi + mij);
  const double prefactor = 4./27*GetLDME(m_flavs[2].Kfcode())/pow(mi,5);
  return prefactor * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.1561 * m_Jmax;
}

double LF_FF3P2_Quarkonia_IF::Z() {
  if (m_flavs[0] == m_flavs[1]) {
    const double r = ATOOLS::ran->Get();
    return 1./((1.-r)/m_zmin+r/m_zmax);
  }
  // Unequal-flavour channel retains its legacy prescription.

  return m_zmin + (m_zmax - m_zmin) * ATOOLS::ran->Get();
}

double LF_FF3P2_Quarkonia_II::operator()(const double z, const double y,
                                       const double eta, const double scale,
                                       const double Q2) {
  if (m_flavs[0] == m_flavs[1]) {
    if (!(z > 0. && z < 1. && y > 0. && y < 1.-z && scale > 0. && Q2 > 0.)) return 0.;
    const double mc2 = sqr(m_flavs[1].Mass(true));
    const double ma2 = p_ms->Mass2(m_flavs[0]);
    const double mj2 = sqr(m_flavs[2].Mass(true));
    const double ms2 = p_ms->Mass2(m_flspec);
    const double V = (Q2-ma2-mj2-ms2)*y/z-ma2-mj2;
    if (IsBad(V) || !(V > 0.) ||
        V < z*mj2/(1.-z)-z*ma2) return 0.;
    // Continue the FF reduced kernel with s=-V and z_FF=1/z.
    // The -z crossing factor follows FFV; JII supplies its flux once.
    const double den = -V-mc2, zt = 1./z;
    const double crossed = -z*FF3P2Reduced(zt, mc2, den);
    // All current incoming-quark -> quark + quarkonium schemes have dV/dt=1/(1-z).
    return sqr(p_cf->Coupling(scale,0))*crossed*scale/((1.-z)*(V+mc2))
        * JII(z,y,eta,scale);
  }
  // Unequal-flavour channel retains its legacy prescription.

  double ma  = p_ms->Mass(m_flavs[1]); 
  double mi  = p_ms->Mass(m_flavs[2]);
  double mb  = p_ms->Mass(m_flspec);
  double mai = p_ms->Mass(m_flavs[0]);
  double mua2 = sqr(ma) / Q2,  mub2 = sqr(mb) / Q2, muai2 = sqr(mai) / Q2;
  const double zt = 1.0 / (z + y);
  const double sab = (Q2 - sqr(mi) - (1.0 - z) * (sqr(ma) + sqr(mb))) / z;
  const double saj = sqr(ma) + sqr(mi) - y * (sab - sqr(ma) - sqr(mb));
  const double M = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true) + ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double ri = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true) / M;
  const double rij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true) / M;
  const double den = saj - sqr(rij * M);
  double value = 0;
  value = sqr(sqr(sqr(M)))/sqr(sqr(den))*( 320 * sqr(ri) * cube(rij) * sqr(sqr(1 - rij*(1-zt))) );
  value += sqr(cube(M))/cube(den) * 8 * ri * sqr(rij) * cube(1 - rij * (1 - zt)) * (
    2 * (4 + 13 * ri) - (1 + 70 * ri - 26 * sqr(ri)) * (1 - zt) - (7 + 8 * ri) * rij * sqr(1 - zt)
  );
  value += sqr(sqr(M))/sqr(den) * (-4) * sqr(rij) * sqr( 1 - rij * (1 - zt)) * ( 
    4 * (1 + 4 * ri) - (7 + 12 * ri - 32 * sqr(ri)) * (1 - zt) 
    + 2 * (1 + 13 * ri - 26 * sqr(ri) + 8 * cube(ri)) * sqr(1 - zt) + (1 - 30 * ri - 5 * sqr(ri) + 4 * cube(ri)) * cube(1 - zt)
  );
  value += sqr(M)/(den) * 4 * sqr(rij) * zt * (
    2 - 4 * (1 - 2 * ri) * (1 - zt)
    + (5 - 8 * ri + 12 * sqr(ri)) * sqr(1 - zt)
    - 2 * (1 - 2 * ri) * (3 + 2 * sqr(ri)) * cube(1 - zt) 
    + (3 - 12 * ri + 12 * sqr(ri) + 2 * sqr(sqr(ri))) * sqr(sqr((1 - zt)))
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
  if (m_flavs[0] == m_flavs[1]) {
    m_zmin = Max(0.,zmin); m_zmax = Min(1.,zmax);
    m_Jmax = m_disablepdfrescaling ? 1. : 5.;
    if (!(m_zmin > 0. && m_zmax > m_zmin)) return 0.;
    const double mc = m_flavs[1].Mass(true);
    const double prefactor = M_PI/27*(GetLDME(m_flavs[2].Kfcode())/5.)/pow(mc,5)
        * (m_flavs[2].StrongCharge() == 0 ? 1. : 1./48);
    return prefactor*OverEstimated(0.5,0.)/4.
        * (m_zmax-m_zmin)/(m_zmin*m_zmax);
  }
  // Unequal-flavour channel retains its legacy prescription.

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
  if (m_flavs[0] == m_flavs[1]) {
    // Common bound for the current evolution schemes, M_H >= 2*m_c,
    // and an incoming kinematic mass no larger than the constituent mass.
    return 256.*sqr(p_cf->MaxCoupling(0))*m_Jmax/sqr(z);
  }
  // Unequal-flavour channel retains its legacy prescription.

  const double mi = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double rij = mij / (mi + mij);
  const double ri = mi / (mi + mij);
  const double prefactor = 4./27*GetLDME(m_flavs[2].Kfcode())/pow(mi,5);
  return prefactor * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.1561 * m_Jmax;
}

double LF_FF3P2_Quarkonia_II::Z() {
  if (m_flavs[0] == m_flavs[1]) {
    const double r = ATOOLS::ran->Get();
    return 1./((1.-r)/m_zmin+r/m_zmax);
  }
  // Unequal-flavour channel retains its legacy prescription.

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
  const double sij = (Q2 + sqr(mi)) * y / (1 - y) + sqr(mi) + sqr(mj);
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



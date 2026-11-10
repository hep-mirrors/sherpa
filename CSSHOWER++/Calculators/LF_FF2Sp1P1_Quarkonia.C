#include "CSSHOWER++/Showers/Splitting_Function_Base.H"
#include "ATOOLS/Phys/LDME.H"

namespace CSSHOWER {

class LF_FF2Sp1P1_Quarkonia_FF : public SF_Lorentz {
  protected:
    double m_ctheta, m_stheta;
public:
  inline LF_FF2Sp1P1_Quarkonia_FF(const SF_Key &key) : SF_Lorentz(key) {}

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_FF2Sp1P1_Quarkonia_FI : public SF_Lorentz {

protected:
  double m_Jmax;
  double m_ctheta, m_stheta;
public:
  inline LF_FF2Sp1P1_Quarkonia_FI(const SF_Key &key) : SF_Lorentz(key) {}

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_FF2Sp1P1_Quarkonia_IF : public SF_Lorentz {

protected:
  double m_Jmax;
  double m_ctheta, m_stheta;
public:
  inline LF_FF2Sp1P1_Quarkonia_IF(const SF_Key &key) : SF_Lorentz(key) {}

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_FF2Sp1P1_Quarkonia_II : public SF_Lorentz {

protected:
  double m_Jmax;
  double m_ctheta, m_stheta;

public:
  inline LF_FF2Sp1P1_Quarkonia_II(const SF_Key &key) : SF_Lorentz(key) {}

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_F2Sp1P1F_Quarkonia_FF : public SF_Lorentz {
  protected:
    double m_ctheta, m_stheta;
public:
  inline LF_F2Sp1P1F_Quarkonia_FF(const SF_Key &key) : SF_Lorentz(key) {}

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_F2Sp1P1F_Quarkonia_FI : public SF_Lorentz {

protected:
  double m_Jmax;
  double m_ctheta, m_stheta;
public:
  inline LF_F2Sp1P1F_Quarkonia_FI(const SF_Key &key) : SF_Lorentz(key) {}

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
double FF2Sp1P1Reduced(const double z, const double mc2, const double den, const bool spin_singlet)
{
  const double x = 1.-z, u = 1.-0.5*x, rho = 4.*mc2/den;
  if (spin_singlet) {
    // The x^4 coefficient is 9/16: this reproduces Eq. (102),
    // also independently confirmed by arXiv:2602.05615 Eq. (B.1).
    const double c0 = 2.*sqr(sqr(u));
    const double c1 = 1.5*cube(u)*(2.-5.*x);
    const double c2 = 0.5*sqr(u)*(15.*sqr(x)-8.*x-4.);
    const double c3 = z*(1.+2.5*sqr(x)-2.*cube(x)+9./16*sqr(sqr(x)));
    return rho*(c3+rho*(c2+rho*(c1+rho*c0)))/sqr(sqr(u));
  }
  const double c0 = 6.*sqr(sqr(u));
  const double c1 = 1.5*cube(u)*(sqr(x)-15.*x+2.);
  const double c2 = 0.75*sqr(u)*(cube(x)+4.*sqr(x)+44.*x-32.);
  const double c3 = 3./8*z*(sqr(sqr(x))+16.);
  return rho*(c3+rho*(c2+rho*(c1+rho*c0)))/sqr(sqr(u));
}
}



double D_FF1P1(const double z, const double s, const double ri, const double rij, const double M) {
  const std::vector<double> fn = { 
    64*sqr(ri)*cube(rij)*sqr(sqr(1-rij*(1-z))),
    8*ri*rij*cube(1-rij*(1-z))*( 
      3 - 2*ri - 2*sqr(ri)
      - 2*rij*(1-z)*( 2 + 4*ri - sqr(ri))+sqr(rij*(1-z))*(1 - 2*ri)),
    -sqr(1-rij*(1-z))*( 2*(1 - 2*ri + 4*sqr(ri))-(1-z)*( 3 - 42*ri + 64*sqr(ri) - 16*cube(ri))
      - 2*ri*rij*sqr(1-z)*( 23 - 14*ri - 4*sqr(ri)) + sqr(rij)*cube(1-z)*(1 + 12*ri)*(1 - 2*ri) ),
    z*(1 - 2*(1-2*ri)*(1-z) + sqr(1-z)*(3-2*ri+2*sqr(ri))
      - 2*rij*cube(1-z)*(2+ri-2*sqr(ri)) + cube(rij)*sqr(sqr(1-z))*(2+sqr(ri)) )
    };
  const double denom = s - sqr(rij * M);
  double value = 0;
  for (size_t n = 0; n < fn.size(); ++n) {
    value += fn[n] * std::pow(M, 8.0 - 2.0 * n) / std::pow(denom, 4.0 - n);
  }
  return value;
}

double D_FF3P1(const double z, const double s, const double ri, const double rij, const double M) {
  const std::vector<double> fn = { 
    192*sqr(ri)*cube(rij)*sqr(sqr(1-rij*(1-z))),
    24*ri*rij*cube(1-rij*(1-z))*(
      2*(1-ri-sqr(ri))
      -rij*(1-z)*(3+10*ri-2*sqr(ri))
      +sqr(rij*(1-z))
      ),
    -6*sqr(1-rij*(1-z))*( 2*(1+2*ri)
        -(1-z)*(5-2*ri+6*sqr(ri))
        +2*rij*sqr(1-z)*(2-3*ri-4*sqr(ri))
        -sqr(rij)*cube(1-z)*(1-2*ri+2*sqr(ri))
      ),
      6*z*( 1-2*(1-2*ri)*(1-z) 
        + sqr(1-z)*(1-4*ri)*(1-2*ri)
        + 2*ri*rij*cube(1-z)*(1-2*ri)
        + sqr(ri*rij*sqr(1-z)))
  };
  double value = 0;
  for (size_t n = 0; n < fn.size(); ++n) {
    value += fn[n] * std::pow(M, 8.0 - 2.0 * n) / std::pow(s - sqr(rij * M), 4.0 - n);
  }
  return value;
}

double D_FFMixingP1(const double z, const double s, const double ri, const double rij, const double M) {
  const std::vector<double> fn = { 
    -8*ri*rij*cube(1-rij*(1-z))*(2-3*ri-rij*(1+ri)*(1-z)),
    (1-rij*(1-z))*(
      2*(1+2*sqr(ri)) - (1-z)*(5-24*ri+36*sqr(ri)-8*cube(ri))
      +2*rij*sqr(1-z)*(2-15*ri+10*sqr(ri)-2*cube(ri))
      -sqr(rij)*cube(1-z)*(1-8*ri-4*sqr(ri))
      ),
      -z*(1- (1-z)*(1-3*ri)+ri*(2+ri)*sqr(1-z)+rij*sqr(ri)*cube(1-z))
    };
    double value = 0;
    for (size_t n = 0; n < fn.size(); ++n) {
      value += fn[n] * std::pow(M, 6.0 - 2.0 * n) / std::pow(s - sqr(rij * M), 3.0 - n);
    };
  return value;
}

double LF_FF2Sp1P1_Quarkonia_FF::operator()(const double z, const double y,
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
    double value = FF2Sp1P1Reduced(z, mi2, den, m_flavs[2].Kfcode()/10000 == 1);
    value *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
    value *= 1. / (1 + sqr( 1 - z) * mi2 / scale + sqr(z) * sqr(mj) / scale);
    return sqr(p_cf->Coupling(scale, 0)) * value * JFF(y, mui2, muj2, muk2, muij2);
  }
  // Unequal-flavour channel retains its legacy prescription.

  // based on Eq. 45 of arXiv:hep-ph/9405348v1
  // again this is b -> c (1P1 + 3P1) 
  const double mi  = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true); 
  const double mj  = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  const double mk  = ATOOLS::Flavour(m_flspec.Kfcode()).Mass(true);
  const double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double mui2 = sqr(mi) / Q2, muj2 = sqr(mj) / Q2, muk2 = sqr(mk) / Q2, muij2 = sqr(mij) / Q2;
  const double M = mi + mij;
  const double ri = mi / M, rij = mij / M;
  const double sij = y * (Q2 - sqr(mi) - sqr(mj) - sqr(mk)) + (sqr(mi) + sqr(mj));
  double value = 0.0;
  value  = sqr(m_ctheta) * D_FF1P1(z, sij, ri, rij, M) + sqr(m_stheta) * D_FF3P1(z, sij, ri, rij, M); 
  //+ 2 * m_ctheta * m_stheta * (1-rij*(1-z)) * D_FFMixingP1(z, sij, ri, rij, M); // missing factor in front of mixing (not relevant for heavy quarks)
  value *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
  value *= 1. / (1 + sqr( 1 - z) * sqr(mi) / scale + sqr(z) * sqr(mj) / scale);
  value *= sqr(p_cf->Coupling(scale, 0)) * ri * cube(rij) / sqr(sqr(1-rij*(1-z))) * JFF(y, mui2, muj2, muk2, muij2);
  // BEWARE: this is missing a global constant factor in front (so does the overestimate)
  return value;
}

double LF_FF2Sp1P1_Quarkonia_FF::OverIntegrated(const double zmin, const double zmax,
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
    const double prefactor = 16.*M_PI/27 * (GetLDME(m_flavs[2].Kfcode())/(m_flavs[2].Kfcode()/10000 == 1 ? 1. : 3.))
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
  const kf_code kfc = m_flavs[2].Kfcode();
  double prefactor = GetLDME(m_flavs[2].Kfcode()) * 4.0/27/pow(mi,5);
  if ( kfc % 2 == 1 && kfc / 10000 == 1 ) { // 1P1 state
    m_ctheta = 1.0;
    m_stheta = 0.0;
    prefactor *= .5; // 1P1 state
  } else if ( kfc % 2 == 1 && kfc / 10000 == 2 ) { // 3P1 state
    m_ctheta = 0.0;
    m_stheta = 1.0;
  } else {
    throw std::runtime_error("LF_FF2Sp1P1_Quarkonia_FF::OverIntegrated: unrecognized kf_code for quarkonium state.");
  } // this is a bit of a hack, but it works as long as mixing is not relevant

  return prefactor * sqr(p_cf->MaxCoupling(0))  * ri * cube(rij) / sqr(sqr(1-rij)) * (m_zmax - m_zmin);
}

double LF_FF2Sp1P1_Quarkonia_FF::OverEstimated(const double z, const double y) {
  if (m_flavs[0] == m_flavs[1]) {
    // Analytic constituent-mass bound; proof in notes/check_ff_quarkonia_states.py.
    return (m_flavs[2].Kfcode()/10000 == 1 ? 2. : 10.)*sqr(p_cf->MaxCoupling(0));
  }
  // Unequal-flavour channel retains its legacy prescription.

  const double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double mi = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double rij = mij / (mi + mij);
  const double ri = mi / (mi + mij);
  return sqr(p_cf->MaxCoupling(0))  * ri * cube(rij) / sqr(sqr(1-rij));
}

double LF_FF2Sp1P1_Quarkonia_FF::Z() {
  return m_zmin + (m_zmax-m_zmin) * ATOOLS::ran->Get();
}

double LF_FF2Sp1P1_Quarkonia_FI::operator()(const double z, const double y,
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
    double value = FF2Sp1P1Reduced(z, mi2, den, m_flavs[2].Kfcode()/10000 == 1);
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
  value  = sqr(m_ctheta) * D_FF1P1(z, sij, ri, rij, M) + sqr(m_stheta) * D_FF3P1(z, sij, ri, rij, M);
  // value *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
  // value *= 1. / (1 + sqr( 1 - z) * sqr(mi) / scale + sqr(z) * sqr(mj) / scale);
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *= 4.0 / 27 / pow(mi,5);
  return prefactor * sqr(p_cf->Coupling(scale, 0)) * value * JFI(y, eta, scale);
}

double LF_FF2Sp1P1_Quarkonia_FI::OverIntegrated(const double zmin, const double zmax,
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
    const double prefactor = 16.*M_PI/27 * (GetLDME(m_flavs[2].Kfcode())/(m_flavs[2].Kfcode()/10000 == 1 ? 1. : 3.))
        / pow(mi,5) * ri*cube(rij)
        * (m_flavs[2].StrongCharge() == 0 ? 1. : 1./48);
    return prefactor * OverEstimated(0.5,0.) * (m_zmax-m_zmin);
  }
  // Unequal-flavour channel retains its legacy prescription.

  m_zmin = zmin;
  m_zmax = zmax;
  m_Jmax = 5.;
  const double mi = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double rij = mij / (mi + mij);
  const double ri = mi / (mi + mij);
  const kf_code kfc = m_flavs[2].Kfcode();
  if ( kfc % 2 == 1 && kfc / 10000 == 1 ) { // 1P1 state
    m_ctheta = 1.0;
    m_stheta = 0.0;
  } else if ( kfc % 2 == 1 && kfc / 10000 == 2 ) { // 3P1 state
    m_ctheta = 0.0;
    m_stheta = 1.0;
  } else {
    throw std::runtime_error("LF_FF2Sp1P1_Quarkonia_FF::OverIntegrated: unrecognized kf_code for quarkonium state.");
  } // this is a bit of a hack, but it works as long as mixing is not relevant
  const double prefactor = GetLDME(m_flavs[2].Kfcode()) * 4.0/27/pow(mi,5);
  return prefactor * sqr(p_cf->MaxCoupling(0))  * ri * cube(rij) / sqr(sqr(1-rij)) * (m_zmax - m_zmin) * m_Jmax;

}

double LF_FF2Sp1P1_Quarkonia_FI::OverEstimated(const double z, const double y) {
  if (m_flavs[0] == m_flavs[1]) {
    // Analytic constituent-mass bound; proof in notes/check_ff_quarkonia_states.py.
    return (m_flavs[2].Kfcode()/10000 == 1 ? 2. : 10.)*sqr(p_cf->MaxCoupling(0)) * m_Jmax;
  }
  // Unequal-flavour channel retains its legacy prescription.

  const double mi = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double rij = mij / (mi + mij);
  const double ri = mi / (mi + mij);
  const double prefactor = GetLDME(m_flavs[2].Kfcode()) * 4.0/27/pow(mi,5);
  return prefactor * sqr(p_cf->MaxCoupling(0))  * ri * cube(rij) / sqr(sqr(1-rij)) * m_Jmax;
}

double LF_FF2Sp1P1_Quarkonia_FI::Z() {
  return m_zmin + (m_zmax-m_zmin) * ATOOLS::ran->Get();
}

double LF_FF2Sp1P1_Quarkonia_IF::operator()(const double z, const double y,
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
    const double crossed = -z*FF2Sp1P1Reduced(zt, mc2, den, m_flavs[2].Kfcode()/10000 == 1);
    // All current incoming-quark -> quark + quarkonium schemes have dV/dt=1/(1-z).
    return sqr(p_cf->Coupling(scale,0))*crossed*scale/((1.-z)*(V+mc2))
        * JIF(z,y,eta,scale);
  }
  // Unequal-flavour channel retains its legacy prescription.

  // c --> c (2S+1)P1 
  double mi  = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true); // works with the mapping c -> c J/psi
  double mj  = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  double mk  = ATOOLS::Flavour(m_flspec.Kfcode()).Mass(true);
  double ma = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  double muj2 = sqr(mj) / Q2, mui2 = sqr(mi) / Q2, muk2 = sqr(mk) / Q2, mua2 = sqr(ma) / Q2;
  const double zt = (1.0 - y) / (z - y);
  const double saj = (sqr(ma) + sqr(mj)) * (1.0 - y / z) + (Q2 - sqr(mk)) * (y / z);
  const double M = mi + ma;
  const double ri = mi / M;
  const double rij = ma / M;
  const double den = saj - sqr(rij * M);
  double value = 0;
  value  = sqr(m_ctheta) * D_FF1P1(zt, saj, ri, rij, M) + sqr(m_stheta) * D_FF3P1(zt, saj, ri, rij, M);
  // value *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
  // value *= 1. / (1 + sqr( 1 - z) * sqr(mi) / scale + sqr(z) * sqr(mj) / scale);
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *=  4.0 / 27 / pow(ri * M,5);
  return prefactor * sqr(p_cf->Coupling(scale, 0)) * value * JIF(z, y, eta, scale);
}

double LF_FF2Sp1P1_Quarkonia_IF::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  if (m_flavs[0] == m_flavs[1]) {
    m_zmin = Max(0.,zmin); m_zmax = Min(1.,zmax);
    m_Jmax = m_disablepdfrescaling ? 1. : 5.;
    if (!(m_zmin > 0. && m_zmax > m_zmin)) return 0.;
    const double mc = m_flavs[1].Mass(true);
    const double prefactor = M_PI/27*(GetLDME(m_flavs[2].Kfcode())/(m_flavs[2].Kfcode()/10000 == 1 ? 1. : 3.))/pow(mc,5)
        * (m_flavs[2].StrongCharge() == 0 ? 1. : 1./48);
    return prefactor*OverEstimated(0.5,0.)/4.
        * (m_zmax-m_zmin)/(m_zmin*m_zmax);
  }
  // Unequal-flavour channel retains its legacy prescription.

  m_zmin = zmin;
  m_zmax = zmax;
  m_Jmax = m_flavs[0].Kfcode() < 3 ? 5. : 1.;
  const double mai = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double ma = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double rij = ma / (mai + ma);
  const double ri = mai / (mai + ma);
  const kf_code kfc = m_flavs[2].Kfcode();
  if ( kfc % 2 == 1 && kfc / 10000 == 1 ) { // 1P1 state
    m_ctheta = 1.0;
    m_stheta = 0.0;
  } else if ( kfc % 2 == 1 && kfc / 10000 == 2 ) { // 3P1 state
    m_ctheta = 0.0;
    m_stheta = 1.0;
  } else {
    throw std::runtime_error("LF_FF2Sp1P1_Quarkonia_FF::OverIntegrated: unrecognized kf_code for quarkonium state.");
  } // this is a bit of a hack, but it works as long as mixing is not relevant
  const double prefactor = GetLDME(m_flavs[2].Kfcode()) * 4.0/27/pow(mai,5);
  return prefactor * sqr(p_cf->MaxCoupling(0))  * ri * cube(rij) / sqr(sqr(1-rij)) * (m_zmax - m_zmin) * m_Jmax;
}

double LF_FF2Sp1P1_Quarkonia_IF::OverEstimated(const double z, const double y) {
  if (m_flavs[0] == m_flavs[1]) {
    // Common bound for the current evolution schemes, M_H >= 2*m_c,
    // and an incoming kinematic mass no larger than the constituent mass.
    return (m_flavs[2].Kfcode()/10000 == 1 ? 24. : 144.)*sqr(p_cf->MaxCoupling(0))*m_Jmax/sqr(z);
  }
  // Unequal-flavour channel retains its legacy prescription.

  const double mai = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double ma = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double rij = ma / (mai + ma);
  const double ri = mai / (mai + ma);
  const double prefactor = GetLDME(m_flavs[2].Kfcode()) * 4.0/27/pow(mai,5);
  return prefactor * sqr(p_cf->MaxCoupling(0))  * ri * cube(rij) / sqr(sqr(1-rij)) * m_Jmax;
}

double LF_FF2Sp1P1_Quarkonia_IF::Z() {
  if (m_flavs[0] == m_flavs[1]) {
    const double r = ATOOLS::ran->Get();
    return 1./((1.-r)/m_zmin+r/m_zmax);
  }
  // Unequal-flavour channel retains its legacy prescription.

  return m_zmin + (m_zmax - m_zmin) * ATOOLS::ran->Get();
}

double LF_FF2Sp1P1_Quarkonia_II::operator()(const double z, const double y,
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
    const double crossed = -z*FF2Sp1P1Reduced(zt, mc2, den, m_flavs[2].Kfcode()/10000 == 1);
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
  value  = sqr(m_ctheta) * D_FF1P1(zt, saj, ri, rij, M) + sqr(m_stheta) * D_FF3P1(zt, saj, ri, rij, M);
  // value *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
  // value *= 1. / (1 + sqr( 1 - z) * sqr(mi) / scale + sqr(z) * sqr(mj) / scale);
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *= 4.0 / 27 / pow(ma,5);
  return prefactor * sqr(p_cf->Coupling(scale, 0)) * value * JII(z, y, eta, scale);
}

double LF_FF2Sp1P1_Quarkonia_II::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  if (m_flavs[0] == m_flavs[1]) {
    m_zmin = Max(0.,zmin); m_zmax = Min(1.,zmax);
    m_Jmax = m_disablepdfrescaling ? 1. : 5.;
    if (!(m_zmin > 0. && m_zmax > m_zmin)) return 0.;
    const double mc = m_flavs[1].Mass(true);
    const double prefactor = M_PI/27*(GetLDME(m_flavs[2].Kfcode())/(m_flavs[2].Kfcode()/10000 == 1 ? 1. : 3.))/pow(mc,5)
        * (m_flavs[2].StrongCharge() == 0 ? 1. : 1./48);
    return prefactor*OverEstimated(0.5,0.)/4.
        * (m_zmax-m_zmin)/(m_zmin*m_zmax);
  }
  // Unequal-flavour channel retains its legacy prescription.

  m_zmin = zmin;
  m_zmax = zmax;
  m_Jmax = m_flavs[0].Kfcode() < 3 ? 5. : 1.;
  const double mai = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double ma = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double rij = ma / (mai + ma);
  const double ri = mai / (mai + ma);
  const kf_code kfc = m_flavs[2].Kfcode();
  if ( kfc % 2 == 1 && kfc / 10000 == 1 ) { // 1P1 state
    m_ctheta = 1.0;
    m_stheta = 0.0;
  } else if ( kfc % 2 == 1 && kfc / 10000 == 2 ) { // 3P1 state
    m_ctheta = 0.0;
    m_stheta = 1.0;
  } else {
    throw std::runtime_error("LF_FF2Sp1P1_Quarkonia_FF::OverIntegrated: unrecognized kf_code for quarkonium state.");
  } // this is a bit of a hack, but it works as long as mixing is not relevant
  const double prefactor = GetLDME(m_flavs[2].Kfcode()) * 4.0/27/pow(mai,5);
  return prefactor * sqr(p_cf->MaxCoupling(0))  * ri * cube(rij) / sqr(sqr(1-rij)) * (m_zmax - m_zmin) * m_Jmax;
}

double LF_FF2Sp1P1_Quarkonia_II::OverEstimated(const double z, const double y) {
  if (m_flavs[0] == m_flavs[1]) {
    // Common bound for the current evolution schemes, M_H >= 2*m_c,
    // and an incoming kinematic mass no larger than the constituent mass.
    return (m_flavs[2].Kfcode()/10000 == 1 ? 24. : 144.)*sqr(p_cf->MaxCoupling(0))*m_Jmax/sqr(z);
  }
  // Unequal-flavour channel retains its legacy prescription.

  const double mai = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double ma = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double rij = ma / (mai + ma);
  const double ri = mai / (mai + ma);
  const double prefactor = GetLDME(m_flavs[2].Kfcode()) * 4.0/27/pow(mai,5);
  return prefactor * sqr(p_cf->MaxCoupling(0))  * ri * cube(rij) / sqr(sqr(1-rij)) * m_Jmax;
}

double LF_FF2Sp1P1_Quarkonia_II::Z() {
  if (m_flavs[0] == m_flavs[1]) {
    const double r = ATOOLS::ran->Get();
    return 1./((1.-r)/m_zmin+r/m_zmax);
  }
  // Unequal-flavour channel retains its legacy prescription.

  return m_zmin + (m_zmax - m_zmin) * ATOOLS::ran->Get();
}


double LF_F2Sp1P1F_Quarkonia_FF::operator()(const double z, const double y,
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
  value  = sqr(m_ctheta) * D_FF1P1(1-z, sij, ri, rij, M) + sqr(m_stheta) * D_FF3P1(1-z, sij, ri, rij, M);
  value *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
  value *= 1. / (1 + sqr(1-z) * sqr(mi) / scale + sqr(z) * sqr(mj) / scale);
  double prefactor = GetLDME(m_flavs[1].Kfcode());
  prefactor *= 4.0/27 / pow(mj,5);
  return prefactor * sqr(p_cf->Coupling(scale, 0)) * value * JFF(y, mui2, muj2, muk2, muij2); //
}

double LF_F2Sp1P1F_Quarkonia_FF::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  double mj = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true); // works with the mapping c -> c J/psi
  double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double ri = mj / (mj + mij);
  const double rij = mij / (mj + mij);
  m_zmin = zmin; 
  m_zmax = zmax;
const kf_code kfc = m_flavs[2].Kfcode();
  if ( kfc % 2 == 1 && kfc / 10000 == 1 ) { // 1P1 state
    m_ctheta = 1.0;
    m_stheta = 0.0;
  } else if ( kfc % 2 == 1 && kfc / 10000 == 2 ) { // 3P1 state
    m_ctheta = 0.0;
    m_stheta = 1.0;
  } else {
    throw std::runtime_error("LF_FF2Sp1P1_Quarkonia_FF::OverIntegrated: unrecognized kf_code for quarkonium state.");
  } // this is a bit of a hack, but it works as long as mixing is not relevant
  const double prefactor = GetLDME(m_flavs[2].Kfcode()) * 4.0/27/pow(mj,5);
  return prefactor * sqr(p_cf->MaxCoupling(0))  * ri * cube(rij) / sqr(sqr(1-rij)) * (m_zmax - m_zmin);
}

double LF_F2Sp1P1F_Quarkonia_FF::OverEstimated(const double z, const double y) {
  double mj = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double ri = mj / (mj + mij);
  const double rij = mij / (mj + mij);
  const double prefactor = GetLDME(m_flavs[2].Kfcode()) * 4.0/27/pow(mj,5);
  return prefactor * sqr(p_cf->MaxCoupling(0))  * ri * cube(rij) / sqr(sqr(1-rij));
}

double LF_F2Sp1P1F_Quarkonia_FF::Z() {
  return m_zmin + (m_zmax - m_zmin) * ATOOLS::ran->Get();
}

double LF_F2Sp1P1F_Quarkonia_FI::operator()(const double z, const double y,
                                       const double eta, const double scale,
                                       const double Q2) {
  // morally y is x ij,a
  double mi  = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true); // works with the mapping c -> c J/psi
  double mj  = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  double ma  = ATOOLS::Flavour(m_flspec.Kfcode()).Mass(true);
  double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  double mui2 = sqr(mi) / Q2, muj2 = sqr(mj) / Q2, muk2 = sqr(ma) / Q2, muij2 = sqr(mij) / Q2;
  const double sij = (Q2 + sqr(mi)) * y / (1 - y) + sqr(mi) + sqr(mj);
  const double M = mj + mij;
  const double ri = mj / M;
  const double rij = mij / M;
  const double den = den;
  double value = 0;
  value  = sqr(m_ctheta) * D_FF1P1(1-z, sij, ri, rij, M) + sqr(m_stheta) * D_FF3P1(1-z, sij, ri, rij, M);
  value *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
  value *= 1. / (1 + sqr(1-z) * sqr(mi) / scale + sqr(z) * sqr(mj) / scale);
  double prefactor = GetLDME(m_flavs[1].Kfcode());
  prefactor *= 4.0 / 27 / pow(mj,5);
  return prefactor * sqr(p_cf->Coupling(scale, 0)) * value * JFI(y, eta, scale);
}

double LF_F2Sp1P1F_Quarkonia_FI::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  m_zmin = zmin;
  m_zmax = zmax;
  m_Jmax = 5.;
  double mj = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double ri = mj / (mj + mij);
  const double rij = mij / (mj + mij);
  const kf_code kfc = m_flavs[2].Kfcode();
  if ( kfc % 2 == 1 && kfc / 10000 == 1 ) { // 1P1 state
    m_ctheta = 1.0;
    m_stheta = 0.0;
  } else if ( kfc % 2 == 1 && kfc / 10000 == 2 ) { // 3P1 state
    m_ctheta = 0.0;
    m_stheta = 1.0;
  } else {
    throw std::runtime_error("LF_FF2Sp1P1_Quarkonia_FF::OverIntegrated: unrecognized kf_code for quarkonium state.");
  } // this is a bit of a hack, but it works as long as mixing is not relevant
  const double prefactor = GetLDME(m_flavs[2].Kfcode()) * 4.0/27/pow(mij,5);
  return prefactor * sqr(p_cf->MaxCoupling(0))  * ri * cube(rij) / sqr(sqr(1-rij)) * (m_zmax - m_zmin) * m_Jmax;
}

double LF_F2Sp1P1F_Quarkonia_FI::OverEstimated(const double z, const double y) {
  double mj = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double ri = mj / (mj + mij);
  const double rij = mij / (mj + mij);
  const double prefactor = GetLDME(m_flavs[2].Kfcode()) * 4.0/27/pow(mij,5);  
  return prefactor * sqr(p_cf->MaxCoupling(0))  * ri * cube(rij) / sqr(sqr(1-rij)) * m_Jmax;
}

double LF_F2Sp1P1F_Quarkonia_FI::Z() {
  return m_zmin + (m_zmax - m_zmin) * ATOOLS::ran->Get();
}


DECLARE_GETTER(LF_FF2Sp1P1_Quarkonia_FF, "FF2Sp1P1_Quarkonia", SF_Lorentz, SF_Key);

SF_Lorentz *ATOOLS::Getter<SF_Lorentz, SF_Key, LF_FF2Sp1P1_Quarkonia_FF>::operator()(
    const Parameter_Type &args) const {
  if (args.m_col < 0)
    return NULL;
  if ((args.m_mode == 0 && args.p_v->in[0].IntSpin() == 1 &&
       args.p_v->in[1].IntSpin() == 1 && args.p_v->in[2].IntSpin() == 2) ||
      (args.m_mode == 1 && args.p_v->in[0].IntSpin() == 1 &&
       args.p_v->in[1].IntSpin() == 2 && args.p_v->in[2].IntSpin() == 1)) {
    switch (args.m_type) {
    case cstp::FF: return new LF_FF2Sp1P1_Quarkonia_FF(args);
    case cstp::FI: return new LF_FF2Sp1P1_Quarkonia_FI(args);
    case cstp::IF: return new LF_FF2Sp1P1_Quarkonia_IF(args);
    case cstp::II: return new LF_FF2Sp1P1_Quarkonia_II(args);
    case cstp::none:
      break;
    }
  }
  if ((args.m_mode == 0 && args.p_v->in[0].IntSpin() == 1 &&
       args.p_v->in[1].IntSpin() == 2 && args.p_v->in[2].IntSpin() == 1) ||
      (args.m_mode == 1 && args.p_v->in[0].IntSpin() == 1 &&
       args.p_v->in[2].IntSpin() == 2 && args.p_v->in[1].IntSpin() == 1)) {
    switch (args.m_type) {
    case cstp::FF: return new LF_F2Sp1P1F_Quarkonia_FF(args);
    case cstp::FI: return new LF_F2Sp1P1F_Quarkonia_FI(args);
    case cstp::IF: NULL;
    case cstp::II: NULL;
    case cstp::none:
      break;
    }
  }
  return NULL;
}

void ATOOLS::Getter<SF_Lorentz, SF_Key, LF_FF2Sp1P1_Quarkonia_FF>::PrintInfo(
    std::ostream &str, const size_t width) const {
  str << "FF2Sp1P1 lorentz functions";
}



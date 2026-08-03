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
  const double mi = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double rij = mij / (mi + mij);
  const double ri = mi / (mi + mij);
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
  const double prefactor = GetLDME(m_flavs[2].Kfcode()) * 4.0/27/pow(mi,5);
  return prefactor * sqr(p_cf->MaxCoupling(0))  * ri * cube(rij) / sqr(sqr(1-rij)) * (m_zmax - m_zmin);
}

double LF_FF2Sp1P1_Quarkonia_FF::OverEstimated(const double z, const double y) {
  const double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double mi = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double rij = mij / (mi + mij);
  const double ri = mi / (mi + mij);
  return sqr(p_cf->MaxCoupling(0))  * ri * cube(rij) / sqr(sqr(1-rij));
}

double LF_FF2Sp1P1_Quarkonia_FF::Z() {
  return m_zmin + (m_zmax - m_zmin) * ATOOLS::ran->Get();
}

double LF_FF2Sp1P1_Quarkonia_FI::operator()(const double z, const double y,
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
  value  = sqr(m_ctheta) * D_FF1P1(z, sij, ri, rij, M) + sqr(m_stheta) * D_FF3P1(z, sij, ri, rij, M);
  value *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
  value *= 1. / (1 + sqr( 1 - z) * sqr(mi) / scale + sqr(z) * sqr(mj) / scale);
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *= 4.0 / 27 / pow(mi,5);
  return prefactor * sqr(p_cf->Coupling(scale, 0)) * value * JFI(y, eta, scale);
}

double LF_FF2Sp1P1_Quarkonia_FI::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
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
  const double mi = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double rij = mij / (mi + mij);
  const double ri = mi / (mi + mij);
  const double prefactor = GetLDME(m_flavs[2].Kfcode()) * 4.0/27/pow(mi,5);
  return prefactor * sqr(p_cf->MaxCoupling(0))  * ri * cube(rij) / sqr(sqr(1-rij)) * m_Jmax;
}

double LF_FF2Sp1P1_Quarkonia_FI::Z() {
  return m_zmin + (m_zmax - m_zmin) * ATOOLS::ran->Get();
}

double LF_FF2Sp1P1_Quarkonia_IF::operator()(const double z, const double y,
                                       const double eta, const double scale,
                                       const double Q2) {
  // c --> c (2S+1)P1 
  // morally y is xij,a
  double mai  = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true); // works with the mapping c -> c J/psi
  double mi  = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  double mk  = ATOOLS::Flavour(m_flspec.Kfcode()).Mass(true);
  double ma = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  double mui2 = sqr(mai) / Q2, muj2 = sqr(mi) / Q2, muk2 = sqr(mk) / Q2, muij2 = sqr(ma) / Q2;
  // const double sij = (Q2*(1-y) + sqr(mij) )/y;
  const double tai =  y * (Q2 - sqr(mai) - sqr(mi) - sqr(mk)) + (sqr(mai) + sqr(mi)); 
  const double M = mai + ma;
  const double ri = mai / M;
  const double rij = ma / M;
  const double den = tai - sqr(rij*M);
  double value = 0;
  value  = sqr(m_ctheta) * D_FF1P1(z, tai, ri, rij, M) + sqr(m_stheta) * D_FF3P1(z, tai, ri, rij, M);
  // value *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
  // value *= 1. / (1 + sqr( 1 - z) * sqr(mi) / scale + sqr(z) * sqr(mj) / scale);
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *=  4.0 / 27 / pow(mai,5);
  return prefactor * sqr(p_cf->Coupling(scale, 0)) * value * JFI(y, eta, scale);
}

double LF_FF2Sp1P1_Quarkonia_IF::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
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
  const double mai = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double ma = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double rij = ma / (mai + ma);
  const double ri = mai / (mai + ma);
  const double prefactor = GetLDME(m_flavs[2].Kfcode()) * 4.0/27/pow(mai,5);
  return prefactor * sqr(p_cf->MaxCoupling(0))  * ri * cube(rij) / sqr(sqr(1-rij)) * m_Jmax;
}

double LF_FF2Sp1P1_Quarkonia_IF::Z() {
  return m_zmin + (m_zmax - m_zmin) * ATOOLS::ran->Get();
}

double LF_FF2Sp1P1_Quarkonia_II::operator()(const double z, const double y,
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
  value  = sqr(m_ctheta) * D_FF1P1(z, tai, ri, rij, M) + sqr(m_stheta) * D_FF3P1(z, tai, ri, rij, M);
  // value *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
  // value *= 1. / (1 + sqr( 1 - z) * sqr(mi) / scale + sqr(z) * sqr(mj) / scale);
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *= 4.0 / 27 / pow(ma,5);
  return prefactor * sqr(p_cf->Coupling(scale, 0)) * value * JII(z, y, eta, scale);
}

double LF_FF2Sp1P1_Quarkonia_II::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
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
  const double mai = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double ma = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double rij = ma / (mai + ma);
  const double ri = mai / (mai + ma);
  const double prefactor = GetLDME(m_flavs[2].Kfcode()) * 4.0/27/pow(mai,5);
  return prefactor * sqr(p_cf->MaxCoupling(0))  * ri * cube(rij) / sqr(sqr(1-rij)) * m_Jmax;
}

double LF_FF2Sp1P1_Quarkonia_II::Z() {
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
  double prefactor = GetLDME(m_flavs[1].Kfcode());
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
  // const double sij = (Q2*(1-y) + sqr(mij) )/y;
  const double sij =  y * (Q2 - sqr(mi) - sqr(mj) - sqr(ma)) + (sqr(mi) + sqr(mj));
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
  double prefactor = GetLDME(m_flavs[1].Kfcode())*2;
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
  double prefactor = GetLDME(m_flavs[1].Kfcode())*2;
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



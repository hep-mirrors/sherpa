#include "CSSHOWER++/Showers/Splitting_Function_Base.H"
#include "ATOOLS/Phys/LDME.H"

namespace CSSHOWER {

class LF_FF3P0_Quarkonia_FF : public SF_Lorentz {
public:
  inline LF_FF3P0_Quarkonia_FF(const SF_Key &key) : SF_Lorentz(key) {}

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_FF3P0_Quarkonia_FI : public SF_Lorentz {

protected:
  double m_Jmax;

public:
  inline LF_FF3P0_Quarkonia_FI(const SF_Key &key) : SF_Lorentz(key) {}

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_FF3P0_Quarkonia_IF : public SF_Lorentz {

protected:
  double m_Jmax;

public:
  inline LF_FF3P0_Quarkonia_IF(const SF_Key &key) : SF_Lorentz(key) {}

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_FF3P0_Quarkonia_II : public SF_Lorentz {

protected:
  double m_Jmax;

public:
  inline LF_FF3P0_Quarkonia_II(const SF_Key &key) : SF_Lorentz(key) {}

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_F3P0F_Quarkonia_FF : public SF_Lorentz {
public:
  inline LF_F3P0F_Quarkonia_FF(const SF_Key &key) : SF_Lorentz(key) {}

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_F3P0F_Quarkonia_FI : public SF_Lorentz {

protected:
  double m_Jmax;

public:
  inline LF_F3P0F_Quarkonia_FI(const SF_Key &key) : SF_Lorentz(key) {}

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

double LF_FF3P0_Quarkonia_FF::operator()(const double z, const double y,
                                        const double eta, const double scale,
                                        const double Q2) {
  // based on Eq. 23 of arXiv:hep-ph/9405348v1
  // again this is b -> c 3P0 but holds for any 3P0 singlet state
  const double mi  = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true); 
  const double mj  = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  const double mk  = ATOOLS::Flavour(m_flspec.Kfcode()).Mass(true);
  const double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double mui2 = sqr(mi) / Q2, muj2 = sqr(mj) / Q2, muk2 = sqr(mk) / Q2, muij2 = sqr(mij) / Q2;
  const double M = mi + mij;
  const double ri = mi / M, rij = mij / M;
  const double sij = y * (Q2 - sqr(mi) - sqr(mj) - sqr(mk)) + (sqr(mi) + sqr(mj));
  const double den = sij - sqr(rij*M);
  double value = 0.0;
  value  = sqr(sqr(sqr(M)))/sqr(sqr(den))*( 64* sqr(ri) * cube(rij) * sqr(sqr(1 - rij*(1-z))) );
  value += sqr(cube(M))/cube(den) * ( 8*ri*rij*cube(1-rij*(1-z)) * ( 1 - 18*ri + 14*sqr(ri) -2*rij*(1-z)*(1-2*ri+7*sqr(ri)) + sqr(rij*(1-z))*(1+2*ri) ) );
  value += sqr(sqr(M))/sqr(den) * ( -sqr(1-rij*(1-z)) * (
        2*(1-4*ri)*(1+6*ri-4*sqr(ri))
        - (1-z)*(5+14*ri-8*sqr(ri))+80*cube(ri)-64*sqr(sqr(ri))
        + 2*rij*sqr(1-z)*(2+9*ri+18*sqr(ri)-28*cube(ri)-16*sqr(sqr(ri)))
        - sqr(rij)*cube(1-z)*(1+6*ri+16*sqr(ri)-32*cube(ri))));
  value += sqr(M)/(den) * z*sqr(1 - 4*ri - (1-z)*(1-4*ri)*(1-2*ri)-ri*rij*sqr(1-z)*(3-4*ri));
  value *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
  value *= 1. / (1 + sqr( 1 - z) * sqr(mi) / scale + sqr(z) * sqr(mj) / scale);
  value *= sqr(p_cf->Coupling(scale, 0)) * ri * cube(rij) / sqr(sqr(1-rij*(1-z))) * JFF(y, mui2, muj2, muk2, muij2);
  // BEWARE: this is missing a global constant factor in front (so does the overestimate)
  return value;
}

double LF_FF3P0_Quarkonia_FF::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  const double mi = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double rij = mij / (mi + mij);
  const double ri = mi / (mi + mij);
  m_zmin = zmin; 
  m_zmax = zmax;
  const double prefactor = 4./27*GetLDME(m_flavs[2].Kfcode())/pow(mi,5);
  return prefactor * 8./3 * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * (m_zmax - m_zmin); ;
  // extra factor 2 is spurios but necessary
}

double LF_FF3P0_Quarkonia_FF::OverEstimated(const double z, const double y) {
  const double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double mi = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double rij = mij / (mi + mij);
  const double ri = mi / (mi + mij);
  return sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij));
}

double LF_FF3P0_Quarkonia_FF::Z() {
  return m_zmin + (m_zmax - m_zmin) * ATOOLS::ran->Get();
  return 1. - (1. - m_zmin) * pow((1. - m_zmax) / (1. - m_zmin), ATOOLS::ran->Get());
}

double LF_FF3P0_Quarkonia_FI::operator()(const double z, const double y,
                                       const double eta, const double scale,
                                       const double Q2) {
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
  const double den = sij - sqr(rij*M);
  double value = 0;
  value  = sqr(sqr(sqr(M)))/sqr(sqr(den))*( 64* sqr(ri) * cube(rij) * sqr(sqr(1 - rij*(1-z))) );
  value += sqr(cube(M))/cube(den) * ( 8*ri*rij*cube(1-rij*(1-z)) * ( 1 - 18*ri + 14*sqr(ri) -2*rij*(1-z)*(1-2*ri+7*sqr(ri)) + sqr(rij*(1-z))*(1+2*ri) ) );
  value += sqr(sqr(M))/sqr(den) * ( -sqr(1-rij*(1-z)) * (
        2*(1-4*ri)*(1+6*ri-4*sqr(ri))
        - (1-z)*(5+14*ri-8*sqr(ri))+80*cube(ri)-64*sqr(sqr(ri))
        + 2*rij*sqr(1-z)*(2+9*ri+18*sqr(ri)-28*cube(ri)-16*sqr(sqr(ri)))
        - sqr(rij)*cube(1-z)*(1+6*ri+16*sqr(ri)-32*cube(ri))));
  value += sqr(M)/(den) * z*sqr(1 - 4*ri - (1-z)*(1-4*ri)*(1-2*ri)-ri*rij*sqr(1-z)*(3-4*ri));
  value *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
  value *= 1. / (1 + sqr( 1 - z) * sqr(mi) / scale + sqr(z) * sqr(mj) / scale);
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *= 4.0 / 27 / pow(mi,5);
  return prefactor * sqr(p_cf->Coupling(scale, 0)) * value * JFI(y, eta, scale);
}

double LF_FF3P0_Quarkonia_FI::OverIntegrated(const double zmin, const double zmax,
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
  return prefactor * 8./3 * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * (m_zmax - m_zmin) * m_Jmax;
  // extra factor 2 is spurios but necessary
}

double LF_FF3P0_Quarkonia_FI::OverEstimated(const double z, const double y) {
  const double mi = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double rij = mij / (mi + mij);
  const double ri = mi / (mi + mij);
  const double prefactor = 4./27*GetLDME(m_flavs[2].Kfcode())/pow(mi,5);
  return prefactor * 8./3 * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * m_Jmax;
}

double LF_FF3P0_Quarkonia_FI::Z() {
  return m_zmin + (m_zmax - m_zmin) * ATOOLS::ran->Get();
}

double LF_FF3P0_Quarkonia_IF::operator()(const double z, const double y,
                                       const double eta, const double scale,
                                       const double Q2) {
  // c --> c 3P0
  double mi  = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true); // works with the mapping c -> c J/psi
  double mj  = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  double mk  = ATOOLS::Flavour(m_flspec.Kfcode()).Mass(true);
  double ma = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  double muj2 = sqr(mj) / Q2, mui2 = sqr(mi) / Q2, muk2 = sqr(mk) / Q2, mua2 = sqr(ma) / Q2;
  // morally y is xij,a
  // const double sij = (Q2*(1-y) + sqr(mij) )/y;
  const double zt = (1.0 - y) / (z - y);
  const double saj = (sqr(ma) + sqr(mj)) * (1.0 - y / z) + (Q2 - sqr(mk)) * (y / z);
  const double M = mi + ma;
  const double ri = mi / M;
  const double rij = ma / M;
  const double den = saj - sqr(rij * M);
  double value = 0;
  value  = sqr(sqr(sqr(M)))/sqr(sqr(den))*( 64* sqr(ri) * cube(rij) * sqr(sqr(1 - rij*(1-zt))) );
  value += sqr(cube(M))/cube(den) * ( 8*ri*rij*cube(1-rij*(1-zt)) * ( 1 - 18*ri + 14*sqr(ri) -2*rij*(1-zt)*(1-2*ri+7*sqr(ri)) + sqr(rij*(1-zt))*(1+2*ri) ) );
  value += sqr(sqr(M))/sqr(den) * ( -sqr(1-rij*(1-zt)) * (
        2*(1-4*ri)*(1+6*ri-4*sqr(ri))
        - (1-zt)*(5+14*ri-8*sqr(ri))+80*cube(ri)-64*sqr(sqr(ri))
        + 2*rij*sqr(1-zt)*(2+9*ri+18*sqr(ri)-28*cube(ri)-16*sqr(sqr(ri)))
        - sqr(rij)*cube(1-zt)*(1+6*ri+16*sqr(ri)-32*cube(ri))));
  value += sqr(M)/(den) * zt*sqr(1 - 4*ri - (1-zt)*(1-4*ri)*(1-2*ri)-ri*rij*sqr(1-zt)*(3-4*ri));
  // value *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
  // value *= 1. / (1 + sqr( 1 - z) * sqr(mi) / scale + sqr(z) * sqr(mj) / scale);
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *=  4.0 / 27 / pow(rij * M,5);
  return prefactor * sqr(p_cf->Coupling(scale, 0)) * value * JIF(z, y, eta, scale);
}

double LF_FF3P0_Quarkonia_IF::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  m_zmin = zmin;
  m_zmax = zmax;
  m_Jmax = m_flavs[0].Kfcode() < 3 ? 5. : 1.;
  const double ma = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double mai = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double ra = ma / (mai + ma);
  const double rai = mai / (mai + ma);
  const double prefactor = 4./27*GetLDME(m_flavs[2].Kfcode())/pow(mai,5);
  return prefactor * 8./3 * sqr(p_cf->MaxCoupling(0)) * rai * cube(ra) / sqr(sqr(1-ra)) * (m_zmax - m_zmin) * m_Jmax;
  // extra factor 2 is spurios but necessary
}

double LF_FF3P0_Quarkonia_IF::OverEstimated(const double z, const double y) {
  const double mai = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double ma = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double ra = ma / (mai + ma);
  const double rai = mai / (mai + ma);
  const double prefactor = 4./27*GetLDME(m_flavs[2].Kfcode())/pow(mai,5);
  return prefactor * 8./3 * sqr(p_cf->MaxCoupling(0)) * rai * cube(ra) / sqr(sqr(1-ra)) * m_Jmax;
}

double LF_FF3P0_Quarkonia_IF::Z() {
  return m_zmin + (m_zmax - m_zmin) * ATOOLS::ran->Get();
}

double LF_FF3P0_Quarkonia_II::operator()(const double z, const double y,
                                       const double eta, const double scale,
                                       const double Q2) {
  double ma  = p_ms->Mass(m_flavs[0]); 
  double mi  = p_ms->Mass(m_flavs[2]);
  double mb  = p_ms->Mass(m_flspec);
  double mai = p_ms->Mass(m_flavs[1]);
  double mua2 = sqr(ma) / Q2,  mub2 = sqr(mb) / Q2, muai2 = sqr(mai) / Q2;
  const double zt = 1.0 / (z + y);
  const double sab = (Q2 - sqr(mi) - (1.0 - z) * (sqr(ma) + sqr(mb))) / z;
  const double saj = sqr(ma) + sqr(mi) - y * (sab - sqr(ma) - sqr(mb));
  const double M = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true) + ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double ri = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true) / M;
  const double rij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true) / M;
  const double den = saj - sqr(rij * M);
  double value = 0;
  value  = sqr(sqr(sqr(M)))/sqr(sqr(den))*( 64* sqr(ri) * cube(rij) * sqr(sqr(1 - rij*(1-zt))) );
  value += sqr(cube(M))/cube(den) * ( 8*ri*rij*cube(1-rij*(1-zt)) * ( 1 - 18*ri + 14*sqr(ri) -2*rij*(1-zt)*(1-2*ri+7*sqr(ri)) + sqr(rij*(1-zt))*(1+2*ri) ) );
  value += sqr(sqr(M))/sqr(den) * ( -sqr(1-rij*(1-zt)) * (
        2*(1-4*ri)*(1+6*ri-4*sqr(ri))
        - (1-zt)*(5+14*ri-8*sqr(ri))+80*cube(ri)-64*sqr(sqr(ri))
        + 2*rij*sqr(1-zt)*(2+9*ri+18*sqr(ri)-28*cube(ri)-16*sqr(sqr(ri)))
        - sqr(rij)*cube(1-zt)*(1+6*ri+16*sqr(ri)-32*cube(ri))));
  value += sqr(M)/(den) * zt*sqr(1 - 4*ri - (1-zt)*(1-4*ri)*(1-2*ri)-ri*rij*sqr(1-zt)*(3-4*ri));
  // value *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
  // value *= 1. / (1 + sqr( 1 - z) * sqr(mi) / scale + sqr(z) * sqr(mj) / scale);
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *= 4.0 / 27 / pow(ma,5);
  return prefactor * sqr(p_cf->Coupling(scale, 0)) * value * JII(z, y, eta, scale);
}

double LF_FF3P0_Quarkonia_II::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  m_zmin = zmin;
  m_zmax = zmax;
  m_Jmax = m_flavs[0].Kfcode() < 3 ? 5. : 1.;
  const double mai = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double ma = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double rai = mai / (mai + ma);
  const double ra = ma / (mai + ma);
  const double prefactor = 4./27*GetLDME(m_flavs[2].Kfcode())/pow(mai,5);
  return prefactor * 8./3 * sqr(p_cf->MaxCoupling(0)) * rai * cube(ra) / sqr(sqr(1-ra)) * (m_zmax - m_zmin) * m_Jmax;
  // extra factor 2 is spurious but necessary
}

double LF_FF3P0_Quarkonia_II::OverEstimated(const double z, const double y) {
  const double mai = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double ma = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double rij = ma / (mai + ma);
  const double ri = mai / (mai + ma);
  const double prefactor = 4./27*GetLDME(m_flavs[2].Kfcode())/pow(mai,5);
  return prefactor * 8./3 * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * m_Jmax;
  // extra factor 2 is spurios but necessary
}

double LF_FF3P0_Quarkonia_II::Z() {
  return m_zmin + (m_zmax - m_zmin) * ATOOLS::ran->Get();
}


double LF_F3P0F_Quarkonia_FF::operator()(const double z, const double y,
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
  const double den = sij - sqr(rij*M);
  double value = 0;
  value  = sqr(sqr(sqr(M)))/sqr(sqr(den))*( 64* sqr(ri) * cube(rij) * sqr(sqr(1 - rij*(z))) );
  value += sqr(cube(M))/cube(den) * ( 8*ri*rij*cube(1-rij*(z)) * ( 1 - 18*ri + 14*sqr(ri) -2*rij*(z)*(1-2*ri+7*sqr(ri)) + sqr(rij*(z))*(1+2*ri) ) );
  value += sqr(sqr(M))/sqr(den) * ( -sqr(1-rij*(z)) * (
        2*(1-4*ri)*(1+6*ri-4*sqr(ri))
        - (z)*(5+14*ri-8*sqr(ri))+80*cube(ri)-64*sqr(sqr(ri))
        + 2*rij*sqr(z)*(2+9*ri+18*sqr(ri)-28*cube(ri)-16*sqr(sqr(ri)))
        - sqr(rij)*cube(z)*(1+6*ri+16*sqr(ri)-32*cube(ri))));
  value += sqr(M)/(den) * (1-z)*sqr(1 - 4*ri - (z)*(1-4*ri)*(1-2*ri)-ri*rij*sqr(z)*(3-4*ri));
  value *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
  value *= 1. / (1 + sqr(1-z) * sqr(mi) / scale + sqr(z) * sqr(mj) / scale);
  double prefactor = GetLDME(m_flavs[1].Kfcode());
  prefactor *= 4.0/27 / pow(mj,5);
  return prefactor * sqr(p_cf->Coupling(scale, 0)) * value * JFF(y, mui2, muj2, muk2, muij2); //
}

double LF_F3P0F_Quarkonia_FF::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  double mj = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true); // works with the mapping c -> c J/psi
  double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double rj = mj / (mj + mij);
  const double rij = mij / (mj + mij);
  m_zmin = zmin; 
  m_zmax = zmax;
  double prefactor = GetLDME(m_flavs[1].Kfcode())*2;
  prefactor *=  4.0/27/pow(mj,5);
  return prefactor * 8./3 * sqr(p_cf->MaxCoupling(0)) * rj * cube(rij) / sqr(sqr(1-rij)) * (m_zmax - m_zmin);
  // extra factor 2 is spurios but necessary
}

double LF_F3P0F_Quarkonia_FF::OverEstimated(const double z, const double y) {
  double mj = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double rj = mj / (mj + mij);
  const double rij = mij / (mj + mij);
  double prefactor = GetLDME(m_flavs[1].Kfcode());
  prefactor *= 4.0/27/pow(mj,5);
  return prefactor * 8./3 * sqr(p_cf->MaxCoupling(0)) * rj * cube(rij) / sqr(sqr(1-rij));
  // extra factor 2 is spurios but necessary
}

double LF_F3P0F_Quarkonia_FF::Z() {
  return m_zmin + (m_zmax - m_zmin) * ATOOLS::ran->Get();
}

double LF_F3P0F_Quarkonia_FI::operator()(const double z, const double y,
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
  const double rj = mj / M;
  const double rij = mij / M;
  const double den = sij - sqr(rij*M);
  double value = 0;
  value  = sqr(sqr(sqr(M)))/sqr(sqr(den))*( 64* sqr(rj) * cube(rij) * sqr(sqr(1 - rij*(z))) );
  value += sqr(cube(M))/cube(den) * ( 8*rj*rij*cube(1-rij*(z)) * ( 1 - 18*rj + 14*sqr(rj) -2*rij*(z)*(1-2*rj+7*sqr(rj)) + sqr(rij*(z))*(1+2*rj) ) );
  value += sqr(sqr(M))/sqr(den) * ( -sqr(1-rij*(z)) * (
        2*(1-4*rj)*(1+6*rj-4*sqr(rj))
        - (z)*(5+14*rj-8*sqr(rj))+80*cube(rj)-64*sqr(sqr(rj))
        + 2*rij*sqr(z)*(2+9*rj+18*sqr(rj)-28*cube(rj)-16*sqr(sqr(rj)))
        - sqr(rij)*cube(z)*(1+6*rj+16*sqr(rj)-32*cube(rj))));
  value += sqr(M)/(den) * (1-z)*sqr(1 - 4*rj - (z)*(1-4*rj)*(1-2*rj)-rj*rij*sqr(z)*(3-4*rj));
  value *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
  value *= 1. / (1 + sqr(1-z) * sqr(mi) / scale + sqr(z) * sqr(mj) / scale);
  double prefactor = GetLDME(m_flavs[1].Kfcode());
  prefactor *= 4.0 / 27 / pow(mj,5);
  return prefactor * sqr(p_cf->Coupling(scale, 0)) * value * JFI(y, eta, scale);
}

double LF_F3P0F_Quarkonia_FI::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  m_zmin = zmin;
  m_zmax = zmax;
  m_Jmax = 5.;
  double mj = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double rj = mj / (mj + mij);
  const double rij = mij / (mj + mij);
  double prefactor = GetLDME(m_flavs[1].Kfcode())*2;
  prefactor *= 4.0/27/pow(mj,5);
  return prefactor * 8./3 * sqr(p_cf->MaxCoupling(0)) * rj * cube(rij) / sqr(sqr(1-rij)) * (m_zmax - m_zmin);
  // extra factor 2 is spurios but necessary
}

double LF_F3P0F_Quarkonia_FI::OverEstimated(const double z, const double y) {
  double mj = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double rj = mj / (mj + mij);
  const double rij = mij / (mj + mij);
  double prefactor = GetLDME(m_flavs[1].Kfcode())*2;
  prefactor *= 4.0/27/pow(mj,5);
  return prefactor * 8./3 * sqr(p_cf->MaxCoupling(0)) * rj * cube(rij) / sqr(sqr(1-rij));
}

double LF_F3P0F_Quarkonia_FI::Z() {
  return m_zmin + (m_zmax - m_zmin) * ATOOLS::ran->Get();
}

DECLARE_GETTER(LF_FF3P0_Quarkonia_FF, "FF3P0_Quarkonia", SF_Lorentz, SF_Key);

SF_Lorentz *ATOOLS::Getter<SF_Lorentz, SF_Key, LF_FF3P0_Quarkonia_FF>::operator()(
    const Parameter_Type &args) const {
  if (args.m_col < 0)
    return NULL;
  if ((args.m_mode == 0 && args.p_v->in[0].IntSpin() == 1 &&
       args.p_v->in[1].IntSpin() == 1 && args.p_v->in[2].IntSpin() == 0) ||
      (args.m_mode == 1 && args.p_v->in[0].IntSpin() == 1 &&
       args.p_v->in[1].IntSpin() == 0 && args.p_v->in[2].IntSpin() == 1)) {
    switch (args.m_type) {
    case cstp::FF: return new LF_FF3P0_Quarkonia_FF(args);
    case cstp::FI: return new LF_FF3P0_Quarkonia_FI(args);
    case cstp::IF: return new LF_FF3P0_Quarkonia_IF(args);
    case cstp::II: return new LF_FF3P0_Quarkonia_II(args);
    case cstp::none:
      break;
    }
  }
  if ((args.m_mode == 0 && args.p_v->in[0].IntSpin() == 1 &&
       args.p_v->in[1].IntSpin() == 2 && args.p_v->in[2].IntSpin() == 1) ||
      (args.m_mode == 1 && args.p_v->in[0].IntSpin() == 1 &&
       args.p_v->in[2].IntSpin() == 2 && args.p_v->in[1].IntSpin() == 1)) {
    switch (args.m_type) {
    case cstp::FF: return new LF_F3P0F_Quarkonia_FF(args);
    case cstp::FI: return new LF_F3P0F_Quarkonia_FI(args);
    case cstp::IF: NULL;
    case cstp::II: NULL;
    case cstp::none:
      break;
    }
  }
  return NULL;
}

void ATOOLS::Getter<SF_Lorentz, SF_Key, LF_FF3P0_Quarkonia_FF>::PrintInfo(
    std::ostream &str, const size_t width) const {
  str << "FF3P0 lorentz functions";
}



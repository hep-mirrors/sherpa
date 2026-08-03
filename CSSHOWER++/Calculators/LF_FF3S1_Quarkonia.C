#include "CSSHOWER++/Showers/Splitting_Function_Base.H"
#include "ATOOLS/Phys/LDME.H"

namespace CSSHOWER {

class LF_FF3S1_Quarkonia_FF : public SF_Lorentz {
public:
  LF_FF3S1_Quarkonia_FF(const SF_Key &key) : SF_Lorentz(key) {} ;

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};


class LF_FF3S1_Quarkonia_FI : public SF_Lorentz {

protected:
  double m_Jmax;

public:
  LF_FF3S1_Quarkonia_FI(const SF_Key &key) : SF_Lorentz(key) {};

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_FF3S1_Quarkonia_IF : public SF_Lorentz {

protected:
  double m_Jmax;

public:
  LF_FF3S1_Quarkonia_IF(const SF_Key &key) : SF_Lorentz(key) {};

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_FF3S1_Quarkonia_II : public SF_Lorentz {

protected:
  double m_Jmax;

public:
  LF_FF3S1_Quarkonia_II(const SF_Key &key) : SF_Lorentz(key) {};

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_F3S1F_Quarkonia_FF : public SF_Lorentz {
public:
  LF_F3S1F_Quarkonia_FF(const SF_Key &key) : SF_Lorentz(key) {};

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_F3S1F_Quarkonia_FI : public SF_Lorentz {

  protected:
    double m_Jmax;

  public:
    LF_F3S1F_Quarkonia_FI(const SF_Key &key) : SF_Lorentz(key) {};

    double operator()(const double, const double, const double, const double,
                      const double);
    double OverIntegrated(const double, const double, const double, const double);
    double OverEstimated(const double, const double);
    double Z();
};
} // namespace CSSHOWER

#include "ATOOLS/Math/Random.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include <unordered_map>
#include "ATOOLS/Phys/LDME.H"

using namespace CSSHOWER;
using namespace ATOOLS;


double LF_FF3S1_Quarkonia_FF::operator()(const double z, const double y,
                                       const double eta, const double scale,
                                       const double Q2) {
  double mi  = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true); // works with the mapping c -> c J/psi
  double mj  = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  double mk  = ATOOLS::Flavour(m_flspec.Kfcode()).Mass(true);
  double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  double mui2 = sqr(mi) / Q2, muj2 = sqr(mj) / Q2, muk2 = sqr(mk) / Q2, muij2 = sqr(mij) / Q2;
  const double sij = (sqr(mi) + sqr(mj))*(1-y) + (Q2 - sqr(mk)) * y;
  const double M = mi + mij;
  const double ri = mi / M;
  const double rij = mij / M;
  const double den = sij - sqr(rij*M);
  double value = 0;
  value += sqr(cube(M))/cube(den) * ( -12*ri*rij*sqr(1-rij*(1-z)) );
  value += sqr(sqr(M))/sqr(den) * ( -(1-rij*(1-z)))*( 2*(1+2*ri) - (1+12*ri-4*sqr(ri))*(1-z) - rij*(1+2*ri)*sqr(1-z) );
  value += sqr(M)/den * z * ( 1 + 2*ri*(1-z) + ( 2 + sqr(ri) )*sqr(1-z) );
  value *= ri*cube(rij)/sqr(1-rij*(1-z));
  // msg_Out() << "LF_FF3S1_Quarkonia_FF::operator() z=" << z << " y=" << y << " ri=" << ri << " rij=" << rij << " M=" << M << " scale=" << scale << " Q2=" << Q2 << " value=" << value << std::endl; 
  value *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
  value *= 1. / (1 + sqr( 1 - z) * sqr(mi) / scale + sqr(z) * sqr(mj) / scale);
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *= m_flavs[2].StrongCharge() == 0 ? 4.0 / 27 / cube(mi) : 4.0 / 27 / cube(mi) * 1.0 / 48;
  return prefactor * sqr(p_cf->Coupling(scale, 0)) * value * JFF(y, mui2, muj2, muk2, muij2);
  
  // same flavour case OLD
  // if ( m_flavs[0].Kfcode() == m_flavs[1].Kfcode()) {
  //   const double mij2 = mij*mij;
  //   const double mi2 = mi*mi;
  //   const double mj2 = mj*mj;
  //   double sameflavour = 1./sqr(sij - mij2)*(
  //     sqr(sij) - 2*mi2*sij - 47*sqr(mi2)
  //     - (1-z)*(sij - mi2)*(sij - (mi2 + 2*sqrt(mi2*mj2) +  mj2))
  //     + 4*sij*(sij-mi2)*z*(1-z)/(1+z)
  //     + 12*sqr(sij-mi2)*sqr(1-z)*z/sqr(1+z)
  //     - 4.*mi2*(sij-mi2)*(8.-7.*(1.-z)-5.*(1-z)*(1.-z))/(1.+z)
  //    );
  //   sameflavour /= sij - mij2;
  //   sameflavour *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
  //   sameflavour *= 1. / (1 + sqr( 1 - z) * mi2 / scale + sqr(z) * mj2 / scale);
  //   return sqr(p_cf->Coupling(scale, 0)) *  sameflavour * JFF(y, mui2, muj2, muk2, muij2);
  //   // prefactor should be CF^2/NC/mi*LDME, drops in the ratio with the overestimate,
  // }  
}

double LF_FF3S1_Quarkonia_FF::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double Q2) {
  const double mi  = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true); // works with the mapping c -> c J/psi
  const double mj  = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  const double mk  = ATOOLS::Flavour(m_flspec.Kfcode()).Mass(true);
  const double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double ri = mi / (mi + mij);
  const double rij = mij / (mi + mij);
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *= m_flavs[2].StrongCharge() == 0 ? 4.0 / 27 / cube(mi) : 4.0 / 27 / cube(mi) * 1.0 / 48;
  m_zmin = zmin;
  m_zmax = zmax;
  return prefactor * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.4 * (zmax - zmin);
  // return sqr(p_cf->MaxCoupling(0)) * 2.1 / (8*sqr(mij)) * (zmax - zmin); 
}

double LF_FF3S1_Quarkonia_FF::OverEstimated(const double z, const double y) {
  const double mij2 = p_ms->Mass2(m_flavs[0]);
  const double mi  = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true); // works with the mapping c -> c J/psi
  const double mj  = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  const double mk  = ATOOLS::Flavour(m_flspec.Kfcode()).Mass(true);
  const double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double ri = mi / (mi + mij);
  const double rij = mij / (mi + mij);
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *= m_flavs[2].StrongCharge() == 0 ? 4.0 / 27 / cube(mi) : 4.0 / 27 / cube(mi) * 1.0 / 48;
  return prefactor * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.4;
  // return sqr(p_cf->MaxCoupling(0)) * 2.1 / (8*mij2);
}

double LF_FF3S1_Quarkonia_FF::Z() {
  return m_zmin + (m_zmax - m_zmin) * ATOOLS::ran->Get();
}

double LF_FF3S1_Quarkonia_FI::operator()(const double z, const double y,
                                       const double eta, const double scale,
                                       const double Q2) {
  // morally y is x ij,a
  // works with the mapping c -> c J/psi
  double mi  = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true); 
  double mj  = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  double ma  = ATOOLS::Flavour(m_flspec.Kfcode()).Mass(true);
  double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  double mui2 = sqr(mi) / Q2, muj2 = sqr(mj) / Q2, muk2 = sqr(ma) / Q2, muij2 = sqr(mij) / Q2;
  // const double sij = (Q2*(1-y) + sqr(mij) )/y;
  const double sij = (sqr(mi) + sqr(mj))*(1-y) + (Q2 - sqr(ma)) * y;
  const double M = mi + mij;
  const double ri = mi / M;
  const double rij = mij / M;
  const double den = sij - sqr(rij*M);
  double value = 0;
  value += sqr(cube(M))/cube(den) * ( -12*ri*rij*sqr(1-rij*(1-z)) );
  value += sqr(sqr(M))/sqr(den) * ( -(1-rij*(1-z)))*( 2*(1+2*ri) - (1+12*ri-4*sqr(ri))*(1-z) - rij*(1+2*ri)*sqr(1-z) );
  value += sqr(M)/den * z * ( 1 + 2*ri*(1-z) + ( 2 + sqr(ri) )*sqr(1-z) );
  value *= ri*cube(rij)/sqr(1-rij*(1-z));
  value *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
  value *= 1. / (1 + sqr( 1 - z) * sqr(mi) / scale + sqr(z) * sqr(mj) / scale);
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *= m_flavs[2].StrongCharge() == 0 ? 4.0 / 27 / cube(mi) : 4.0 / 27 / cube(mi) * 1.0 / 48;
  return prefactor * sqr(p_cf->Coupling(scale, 0)) * value * JFI(y, eta, scale);
}

double LF_FF3S1_Quarkonia_FI::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  m_zmin = zmin;
  m_zmax = zmax;
  m_Jmax = 5.;
  const double mi  = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double mj  = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  const double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double ri = mi / (mi + mij);
  const double rij = mij / (mi + mij);
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *= m_flavs[2].StrongCharge() == 0 ? 4.0 / 27 / cube(mi) : 4.0 / 27 / cube(mi) * 1.0 / 48;
  return prefactor * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.4 * (m_zmax - m_zmin) * m_Jmax;
}

double LF_FF3S1_Quarkonia_FI::OverEstimated(const double z, const double y) {
  const double mi  = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true); // works with the mapping c -> c J/psi
  const double mj  = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  const double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double ri = mi / (mi + mij);
  const double rij = mij / (mi + mij);
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *= m_flavs[2].StrongCharge() == 0 ? 4.0 / 27 / cube(mi) : 4.0 / 27 / cube(mi) * 1.0 / 48;
  return prefactor * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.4 * m_Jmax;
}

double LF_FF3S1_Quarkonia_FI::Z() {
  return m_zmin + (m_zmax - m_zmin) * ATOOLS::ran->Get();
}

double LF_FF3S1_Quarkonia_IF::operator()(const double z, const double y,
                                       const double eta, const double scale,
                                       const double Q2) {
  // c --> c J/psi
  // (aj) (k) --> a  j  k
  // z = x_{jk,a};  y = u
  const double ma  = sqr(ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true));
  const double mi  = sqr(ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true));
  const double mai = sqr(ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true));
  const double mk = sqr(ATOOLS::Flavour(m_flspec.Kfcode()).Mass(true));
  // const double pkpa =  Q2 / 2 / z * (1 - y);
  // const double pkpj =  0.5 * (1 - z) / z * (mj2 + mk2 + ma2 - Q2);
  // const double taj  =  y * Q2 + ma2 * (1 - y) + mj2 + y * (mj2 + mk2 + 2 * pkpj); // (pa - pj)^2
  // const double z    =  pkpj / pkpa; // this is the momentum fraction of J/psi w.r.t. parent quark
  // const double tai = sqr(ma) + sqr(mi) - y / z * (Q2 - sqr(mb) - sqr(ma) - sqr(mi));
  // const double tai = sqr(ma) + sqr(mi) - y / z * (Q2 - sqr(mi) - sqr(ma) - sqr(mai));
  const double tai = (Q2 - sqr(mk))*(1-y) + (sqr(ma) + sqr(mi))*y;
  const double M = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true) + ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double ri = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true) / M;
  const double rij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true) / M;
  const double den = tai - sqr(rij*M);
  double value = 0;
  value += sqr(cube(M))/cube(den) * ( -12*ri*rij*sqr(1-rij*(1-z)) );
  value += sqr(sqr(M))/sqr(den) * ( -(1-rij*(1-z)))*( 2*(1+2*ri) - (1+12*ri-4*sqr(ri))*(1-z) - rij*(1+2*ri)*sqr(1-z) );
  value += sqr(M)/den * z * ( 1 + 2*ri*(1-z) + ( 2 + sqr(ri) )*sqr(1-z) );
  value *= ri*cube(rij)/sqr(1-rij*(1-z));
  // msg_Out() << "LF_FF3S1_Quarkonia_FF::operator() z=" << z << " y=" << y << " ri=" << ri << " rij=" << rij << " M=" << M << " scale=" << scale << " Q2=" << Q2 << " value=" << value << std::endl; 
  // value *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
  // value *= 1. / (1 + sqr( 1 - z) * sqr(mi) / scale + sqr(z) * sqr(mj) / scale);
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *= m_flavs[2].StrongCharge() == 0 ? 4.0 / 27 / cube(mi) : 4.0 / 27 / cube(mi) * 1.0 / 48;
  // return prefactor * sqr(p_cf->Coupling(scale, 0)) * value * JII(z, y, eta, scale);
  return sqr(p_cf->Coupling(scale, 0)) * value * JII(z, y, eta, scale);
  // const double newscale = sqrt(scale);
  // double value =
  //       1. / sqr(sqr(1 - muaj2)) *
  //       ((1 - 2 * muaj2 - 47 * sqr(muaj2)) -
  //        z * (1 - muaj2) * (1 - sqr(sqrt(fabs(mua2)) + sqrt(fabs(muj2)))) +
  //        4 * ( z * (1 - z)) / (2 - z) * (1 - mua2) -
  //        4 * (8 - 7 * z - 5 * z * z) / (2 - z) * muaj2 * (1 - muaj2) +
  //        12 * (z * z * (1 - z)) / sqr(2 - z) * sqr(1 - muaj2));
  //   const double LDME = 1. / sqrt(ma2) / sqr(taj) *
  //                       (m_flavs[2].IsOctetMeson() ? (1.5E-02) : pow(0.82, 3));
  //   return 16. / 27 * LDME * sqr(p_cf->Coupling(newscale, 0)) * value /
  //          sqr(taj) * JFI(y, eta, scale);
}

double LF_FF3S1_Quarkonia_IF::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  m_zmin = zmin;
  m_zmax = zmax;
  m_Jmax = m_flavs[0].Kfcode() < 3 ? 5. : 1.;
  const double mai  = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double ma = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double ri = mai / (ma + mai);
  const double rij = ma / (ma + mai);
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *= m_flavs[2].StrongCharge() == 0 ? 4.0 / 27 / cube(ma) : 4.0 / 27 / cube(ma) * 1.0 / 48;
  return prefactor*sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.4 * m_Jmax * (m_zmax - m_zmin);
}

double LF_FF3S1_Quarkonia_IF::OverEstimated(const double z, const double y) {
  const double mai  = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double ma = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double ri = mai / (ma + mai);
  const double rij = ma / (ma + mai);
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *= m_flavs[2].StrongCharge() == 0 ? 4.0 / 27 / cube(ma) : 4.0 / 27 / cube(ma) * 1.0 / 48;
  return prefactor*sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.4 * (m_zmax - m_zmin) * m_Jmax;
}

double LF_FF3S1_Quarkonia_IF::Z() {
  return m_zmin + (m_zmax - m_zmin) * ATOOLS::ran->Get();
}

double LF_FF3S1_Quarkonia_II::operator()(const double z, const double y,
                                       const double eta, const double scale,
                                       const double Q2) {
  // works with the mapping [c aj] ->[c a] [J/psi j]
  // z is effectively xiab and y vtildei

  double ma  = p_ms->Mass(m_flavs[1]); 
  double mi  = p_ms->Mass(m_flavs[2]);
  double mb  = p_ms->Mass(m_flspec);
  double mai = p_ms->Mass(m_flavs[0]);
  // msg_Out() << METHOD << "  ma: " << ma << " mi: " << mi << " mb: " << mb << " mai: " << mai << "\n";
  double mua2 = sqr(ma) / Q2,  mub2 = sqr(mb) / Q2, muai2 = sqr(mai) / Q2;
  // const double tai = sqr(ma) + sqr(mi) - y / z * (Q2 - sqr(mb) - sqr(ma) - sqr(mi));
  const double tai = (Q2 - sqr(mb))*(1-y) + (sqr(ma) + sqr(mi))*y;
  const double M = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true) + ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double ri = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true) / M;
  const double rij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true) / M;
  const double den = tai - sqr(rij*M);
  double value = 0;
  value += sqr(cube(M))/cube(den) * ( -12*ri*rij*sqr(1-rij*(1-z)) );
  value += sqr(sqr(M))/sqr(den) * ( -(1-rij*(1-z)))*( 2*(1+2*ri) - (1+12*ri-4*sqr(ri))*(1-z) - rij*(1+2*ri)*sqr(1-z) );
  value += sqr(M)/den * z * ( 1 + 2*ri*(1-z) + ( 2 + sqr(ri) )*sqr(1-z) );
  value *= ri*cube(rij)/sqr(1-rij*(1-z));
  // msg_Out() << "LF_FF3S1_Quarkonia_FF::operator() z=" << z << " y=" << y << " ri=" << ri << " rij=" << rij << " M=" << M << " scale=" << scale << " Q2=" << Q2 << " value=" << value << std::endl; 
  // value *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
  // value *= 1. / (1 + sqr( 1 - z) * sqr(mi) / scale + sqr(z) * sqr(mj) / scale);
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *= m_flavs[2].StrongCharge() == 0 ? 4.0 / 27 / cube(ma) : 4.0 / 27 / cube(ma) * 1.0 / 48;
  return prefactor * sqr(p_cf->Coupling(scale, 0)) * value * JII(z, y, eta, scale);
  
  // const double sab = (Q2 - mj2) / zz - (1. - zz) / zz * (ma2 + mb2);
  // const double z = 1 - (zz + y) + (ma2 + mb2 + mj2)/(sab - ma2 - mb2); // this is mom. fraction of J/psi w.r.t. to progenitor c
  // const double taj =  ma2 + mj2 - 2 * y * (sab - ma2 - mb2); // this is (pa - pj)^2
}

double LF_FF3S1_Quarkonia_II::OverIntegrated(const double zmin, const double zmax,
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
  return prefactor * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.4 * m_Jmax * (m_zmax - m_zmin);
  // return (4.0*p_cf->MaxCoupling(0) + 0.5*p_cf->MaxCoupling(1)) * log((1.-zmin)/(1.-zmax)) * m_Jmax;
}

double LF_FF3S1_Quarkonia_II::OverEstimated(const double z, const double y) {
  const double mai  = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double ma = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double ri = mai / (ma + mai);
  const double rij = ma / (ma + mai);
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *= m_flavs[2].StrongCharge() == 0 ? 4.0 / 27 / cube(ma) : 4.0 / 27 / cube(ma) * 1.0 / 48;
  return prefactor * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.4 * m_Jmax;
  // return 16. / 27 * pow(0.82, 3) / Flavour(kf_c).Mass() * sqr(p_cf->MaxCoupling(0)) * (0.5) / (1. - z) * m_Jmax;
}

double LF_FF3S1_Quarkonia_II::Z() {
  return m_zmin + (m_zmax - m_zmin) * ATOOLS::ran->Get();
}

double LF_F3S1F_Quarkonia_FF::operator()(const double z, const double y,
                                       const double eta, const double scale,
                                       const double Q2) {
  double mi  = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true); // works with the mapping c -> J/psi c
  double mj  = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  double mk  = ATOOLS::Flavour(m_flspec.Kfcode()).Mass(true);
  double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  double mui2 = sqr(mi) / Q2, muj2 = sqr(mj) / Q2, muk2 = sqr(mk) / Q2, muij2 = sqr(mij) / Q2;
  const double sij = (sqr(mi) + sqr(mj))*(1-y) + (Q2 - sqr(mk)) * y;
  const double M = mj + mij;
  const double ri = mj / M;
  const double rij = mij / M;
  const double den = sij - sqr(rij*M);
  double value = 0;
  value += sqr(cube(M))/cube(den) * ( -12*ri*rij*sqr(1-rij*(z)) );
  value += sqr(sqr(M))/sqr(den) * ( -(1-rij*(z)))*( 2*(1+2*ri) - (1+12*ri-4*sqr(ri))*(z) - rij*(1+2*ri)*sqr(z) );
  value += sqr(M)/den * (1-z) * ( 1 + 2*ri*(z) + ( 2 + sqr(ri) )*sqr(z) );
  value *= ri*cube(rij)/sqr(1-rij*(z));
  value *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
  value *= 1. / (1 + sqr(z) * sqr(mi) / scale + sqr(1-z) * sqr(mj) / scale);
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *= m_flavs[2].StrongCharge() == 0 ? 4.0 / 27 / cube(mj) : 4.0 / 27 / cube(mj) * 1.0 / 48;
  return prefactor * sqr(p_cf->Coupling(scale, 0)) * value * JFF(y, mui2, muj2, muk2, muij2);
}

double LF_F3S1F_Quarkonia_FF::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  const double mi  = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true); // works with the mapping c -> J/psi c
  const double mj  = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  const double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double ri = mj / (mj + mij);
  const double rij = mij / (mj + mij);
  double prefactor = GetLDME(m_flavs[1].Kfcode());
  prefactor *= m_flavs[1].StrongCharge() == 0 ? 4.0 / 27 / cube(mj) : 4.0 / 27 / cube(mj) * 1.0 / 48;
  m_zmin = zmin;
  m_zmax = zmax;
  return prefactor * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.4 * (zmax - zmin);
}

double LF_F3S1F_Quarkonia_FF::OverEstimated(const double z, const double y) {
  const double mi  = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true); 
  const double mj  = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  const double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double ri = mj / (mj + mij);
  const double rij = mij / (mj + mij);
  double prefactor = GetLDME(m_flavs[1].Kfcode());
  prefactor *= m_flavs[1].StrongCharge() == 0 ? 4.0 / 27 / cube(mj) : 4.0 / 27 / cube(mj) * 1.0 / 48;
  return prefactor * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.4;
}

double LF_F3S1F_Quarkonia_FF::Z() {
  return 1. - (m_zmin + (m_zmax - m_zmin) * ATOOLS::ran->Get());
}

double LF_F3S1F_Quarkonia_FI::operator()(const double z, const double y,
                                       const double eta, const double scale,
                                       const double Q2) {
  // morally y is x ij,a
  double mi  = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true); // works with the mapping c -> J/psi c
  double mj  = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  double ma  = ATOOLS::Flavour(m_flspec.Kfcode()).Mass(true);
  double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  double mui2 = sqr(mi) / Q2, muj2 = sqr(mj) / Q2, muk2 = sqr(ma) / Q2, muij2 = sqr(mij) / Q2;
  // const double sij = (Q2*(1-y) + sqr(mij) )/y;
  const double sij = (sqr(mi) + sqr(mj))*(1-y) + (Q2 - sqr(ma)) * y;
  const double M = mj + mij;
  const double ri = mj / M;
  const double rij = mij / M;
  const double den = sij - sqr(rij*M);
  double value = 0;
  value += sqr(cube(M))/cube(den) * ( -12*ri*rij*sqr(1-rij*(z)) );
  value += sqr(sqr(M))/sqr(den) * ( -(1-rij*(z)))*( 2*(1+2*ri) - (1+12*ri-4*sqr(ri))*(z) - rij*(1+2*ri)*sqr(z) );
  value += sqr(M)/den * (1-z) * ( 1 + 2*ri*(z) + ( 2 + sqr(ri) )*sqr(z) );
  value *= ri*cube(rij)/sqr(1-rij*(z));
  value *= 1. / ( (1 - mui2 - muj2 - muk2) + 1./ y * ( mui2 + muj2 - muij2 ) );
  value *= 1. / (1 + sqr(z) * sqr(mi) / scale + sqr(1-z) * sqr(mj) / scale);
  double prefactor = GetLDME(m_flavs[1].Kfcode());
  prefactor *= m_flavs[2].StrongCharge() == 0 ? 4.0 / 27 / cube(mj) : 4.0 / 27 / cube(mj) * 1.0 / 48;
  return prefactor * sqr(p_cf->Coupling(scale, 0)) * value * JFI(y, eta, scale);
}

double LF_F3S1F_Quarkonia_FI::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  m_zmin = zmin;
  m_zmax = zmax;
  m_Jmax = 5.;
  const double mi  = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double mj  = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  const double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double ri = mj / (mj + mij);
  const double rij = mij / (mj + mij);
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *= m_flavs[2].StrongCharge() == 0 ? 4.0 / 27 / cube(mj) : 4.0 / 27 / cube(mj) * 1.0 / 48;
  return prefactor * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.4 * (m_zmax - m_zmin) * m_Jmax;
}

double LF_F3S1F_Quarkonia_FI::OverEstimated(const double z, const double y) {
  const double mi  = ATOOLS::Flavour(m_flavs[1].Kfcode()).Mass(true);
  const double mj  = ATOOLS::Flavour(m_flavs[2].Kfcode()).Mass(true);
  const double mij = ATOOLS::Flavour(m_flavs[0].Kfcode()).Mass(true);
  const double ri = mj / (mj + mij);
  const double rij = mij / (mj + mij);
  double prefactor = GetLDME(m_flavs[2].Kfcode());
  prefactor *= m_flavs[2].StrongCharge() == 0 ? 4.0 / 27 / cube(mi) : 4.0 / 27 / cube(mi) * 1.0 / 48;
  return prefactor * sqr(p_cf->MaxCoupling(0)) * ri * cube(rij) / sqr(sqr(1-rij)) * 0.4 * m_Jmax;
}

double LF_F3S1F_Quarkonia_FI::Z() {
  return m_zmin + (m_zmax - m_zmin) * ATOOLS::ran->Get();
}


DECLARE_GETTER(LF_FF3S1_Quarkonia_FF, "F3S1F_Quarkonia", SF_Lorentz, SF_Key);

SF_Lorentz *ATOOLS::Getter<SF_Lorentz, SF_Key, LF_FF3S1_Quarkonia_FF>::operator()(
    const Parameter_Type &args) const {
  if (args.m_col < 0)
    return NULL;
  if ((args.m_mode == 0 &&
       args.p_v->in[0].IntSpin() == 1 &&
       args.p_v->in[1].IntSpin() == 1 &&
       args.p_v->in[2].IntSpin() == 2) ||
      (args.m_mode == 1 &&
       args.p_v->in[0].IntSpin() == 1 &&
       args.p_v->in[1].IntSpin() == 2 &&
       args.p_v->in[2].IntSpin() == 1)) {
    switch (args.m_type) {
    case cstp::FF:
      return new LF_FF3S1_Quarkonia_FF(args);
    case cstp::FI:
      return new LF_FF3S1_Quarkonia_FI(args);
    case cstp::IF:
      return new LF_FF3S1_Quarkonia_IF(args);
    case cstp::II:
      return new LF_FF3S1_Quarkonia_II(args);
    case cstp::none:
      break;
    }
  }
  if ((args.m_mode == 0 &&
       args.p_v->in[0].IntSpin() == 1 &&
       args.p_v->in[1].IntSpin() == 2 && 
       args.p_v->in[2].IntSpin() == 1) ||
      (args.m_mode == 1 && 
       args.p_v->in[0].IntSpin() == 1 &&
       args.p_v->in[1].IntSpin() == 1 &&
       args.p_v->in[2].IntSpin() == 2)) {
    switch (args.m_type) {
    case cstp::FF: return new LF_F3S1F_Quarkonia_FF(args);
    case cstp::FI: return new LF_F3S1F_Quarkonia_FI(args);
    case cstp::IF: return NULL;
    case cstp::II: return NULL;
    case cstp::none:
      break;
    }
  }
  return NULL;
}

void ATOOLS::Getter<SF_Lorentz, SF_Key, LF_FF3S1_Quarkonia_FF>::PrintInfo(
    std::ostream &str, const size_t width) const {
  str << "FF3S1 lorentz functions";
}



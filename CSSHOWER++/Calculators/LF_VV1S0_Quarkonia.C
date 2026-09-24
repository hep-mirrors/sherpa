#include "CSSHOWER++/Showers/Splitting_Function_Base.H"
#include "ATOOLS/Math/Random.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Phys/LDME.H"

#define ZMIN 1e-9
#define ZMAX 1.

namespace CSSHOWER {

class LF_VV1S0_Quarkonia_FF : public SF_Lorentz {
public:
  inline LF_VV1S0_Quarkonia_FF(const SF_Key &key) : SF_Lorentz(key) {}

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_VV1S0_Quarkonia_FI : public SF_Lorentz {

protected:
  double m_Jmax;
public:
  inline LF_VV1S0_Quarkonia_FI(const SF_Key &key) : SF_Lorentz(key) {}

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_VV1S0_Quarkonia_IF : public SF_Lorentz {

protected:
  double m_Jmax;
public:
  inline LF_VV1S0_Quarkonia_IF(const SF_Key &key) : SF_Lorentz(key) {}

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_VV1S0_Quarkonia_II : public SF_Lorentz {

protected:
  double m_Jmax;

public:
  inline LF_VV1S0_Quarkonia_II(const SF_Key &key) : SF_Lorentz(key) {}

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_V1S0V_Quarkonia_FF : public SF_Lorentz {
public:
  inline LF_V1S0V_Quarkonia_FF(const SF_Key &key) : SF_Lorentz(key) {}

  double operator()(const double, const double, const double, const double,
                    const double);
  double OverIntegrated(const double, const double, const double, const double);
  double OverEstimated(const double, const double);
  double Z();
};

class LF_V1S0V_Quarkonia_FI : public SF_Lorentz {

protected:
  double m_Jmax;
public:
  inline LF_V1S0V_Quarkonia_FI(const SF_Key &key) : SF_Lorentz(key) {}

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

double LF_VV1S0_Quarkonia_FF::operator()(
    const double z, const double y, const double eta, const double scale,
    const double Q2) { 
    const double mj  =  m_flavs[2].Mass();
    const double mk  =  m_flspec.Mass();
    const double muj2 = sqr(m_flavs[2].Mass(true))/Q2;
    const double muk2 = sqr(m_flspec.Mass(true))/Q2;
    const double sij = y * (Q2 - sqr(mk)) + (1-y)*(sqr(mj));
    const double M = mj; //2*Flavour((m_flavs[2].Kfcode() / 10) % 10).Mass(true);
    double value = 0;
    value  = sqr(sij) + sqr(sqr(M)) - 2 * (1 - z) * (sij + sqr(M))*sij + 2*sqr(sij*(1-z));
    value /= sij*sqr(sij-sqr(M));
    value *= 1. / ( (1 - muj2 - muk2) + 1./ y * ( muj2 ) );
    value *= 1. / (1 + sqr(z) * sqr(mj) / scale);
    return sqr(p_cf->Coupling(scale, 0)) * value * JFF(y, 0, muj2, muk2, 0);
}

double LF_VV1S0_Quarkonia_FF::OverIntegrated(const double zmin, const double zmax,
                                           const double scale,
                                           const double xbj) {
  m_zmin = std::max(zmin,ZMIN);
  m_zmax = std::min(zmax,1.);
  const double m = Flavour((m_flavs[2].Kfcode() / 10) % 10).Mass(true);
  const double prefactor = 4./27*GetLDME(m_flavs[2].Kfcode())/m;
  return prefactor * sqr(p_cf->MaxCoupling(0)) / sqr(2*m) * (zmax - zmin);
}

double LF_VV1S0_Quarkonia_FF::OverEstimated(const double z, const double y) {
  const double m = Flavour((m_flavs[2].Kfcode() / 10) % 10).Mass(true);
  return sqr(p_cf->MaxCoupling(0)) /sqr(2*m);
}

double LF_VV1S0_Quarkonia_FF::Z() {
  return m_zmin + (m_zmax - m_zmin) * ran->Get();
}

double LF_VV1S0_Quarkonia_FI::operator()(
  double z, const double y, const double eta, const double scale,
  const double Q2) {
  const double mi  =  0.;
  const double mj  =  m_flavs[2].Mass();
  const double ma  =  p_ms->Mass(m_flspec.Kfcode());
  const double q2 = -Q2;
  double xija = scale/(z*(1-z)) + sqr(mj)*z/(1-z);
  xija /= (q2 - sqr(ma) - sqr(mj));
  xija = 1./(1-xija);
  // msg_Out() << "xija = " << xija << "  " << (1-y)*(q2-sqr(ma)-sqr(mj))/(q2-sqr(ma)) << std::endl; 
  const double yt = (( q2 - sqr(ma)) / (q2 - sqr(ma) - sqr(mj)) - (1.0 - y)) / (1-y);
  // msg_Out() << "yt = " << yt << "  " << (1-xija)/xija  << std::endl;
  // msg_Out() << "xija = " << xija << "  " << 1.0/(1.0-(scale/(z*(1.0-z))+sqr(mi)*(1.0-z)/z+sqr(mj)*z/(1.0-z))/(q2-ma*ma-mi*mi-mj*mj)) << std::endl;
  // msg_Out() << "1-xija / xija = " << (1-xija)/xija << "  " << -(scale/(z*(1.0-z))+sqr(mi)*(1.0-z)/z+sqr(mj)*z/(1.0-z))/(q2-ma*ma-mi*mi-mj*mj) << std::endl;
  // msg_Out() << "t = " << -(1-xija)/xija*(q2-mj*mj)*(z*(1.0-z)) - sqr(mj*z) << "  " << scale << std::endl;
  const double sij = sqr(mj) * (1.0 + yt) - yt * (q2 - sqr(ma));
  const double M = 2*Flavour((m_flavs[2].Kfcode() / 10) % 10).Mass(true);
  const double muj2 = sqr(mj)/Q2;
  const double mu2 = sqr(M)/Q2;
  // if ( z > ((-1 + z) * (1 + muj2))/(-1 + z - muj2)) return 0;
  if ( sij < sqr(mj)/(1-z)) return 0;
  double value = 0;
  value  = sqr(sij) + sqr(sqr(M)) - 2 * (1 - z) * (sij + sqr(M))*sij + 2*sqr(sij*(1-z));
  value /= sij*sqr(sij-sqr(M));
  value *= (1+muj2)/(1 + muj2 - xija) * ( 1 - xija - xija * z * muj2 / ( 1 - z ));
  // msg_Out() << z << "  " << scale << "  " << "  " << q2 << "  " << mj << " --> " << yt << "  " << value << std::endl;
  return sqr(p_cf->Coupling(scale, 0)) * value; // * JFI(y, eta, scale);
}

double LF_VV1S0_Quarkonia_FI::OverIntegrated(const double zmin, const double zmax, const double scale, const double xbj) {
  m_zmin = std::max(zmin,ZMIN);
  m_zmax = std::min(zmax,ZMAX);
  m_Jmax = 5.;
  const double m = Flavour((m_flavs[2].Kfcode() / 10) % 10).Mass(true);
  const double prefactor = 4./27*GetLDME(m_flavs[2].Kfcode())/m;
  return prefactor * sqr(p_cf->MaxCoupling(0)) / sqr(2*m) * (m_zmax - m_zmin) * m_Jmax;
}

double LF_VV1S0_Quarkonia_FI::OverEstimated(const double z, const double y) {
  const double m = Flavour((m_flavs[2].Kfcode() / 10) % 10).Mass(true);
  return sqr(p_cf->MaxCoupling(0)) /sqr(2*m) * m_Jmax;
}

double LF_VV1S0_Quarkonia_FI::Z() {
  return m_zmin + (m_zmax - m_zmin) * ran->Get();
}

double LF_VV1S0_Quarkonia_IF::operator()(const double z, const double y, const double eta, const double scale, const double Q2) {
  const double mj  =  m_flavs[2].Mass(true);
  const double mk  =  m_flspec.Mass(true);
  const double muj2 = sqr(m_flavs[2].Mass(true))/Q2;
  const double muk2 = sqr(m_flspec.Mass(true))/Q2;
  const double sai = sqr(mj) * (1.0 - y / z) + (Q2 - sqr(mk)) * (y / z);
  const double M = 2*Flavour((m_flavs[2].Kfcode() / 10) % 10).Mass(true);
  double value = 0;
  value  = sqr(sai) + sqr(sqr(M)) - 2 * (1 - z) * (sai + sqr(M))*sai + 2*sqr(sai*(1-z));
  value /= sai*sqr(sai-sqr(M));

  value *= 1. / ( 1. + muj2 / (1 - muk2 - muj2) * z / y ) * ( 1 - y ) / ( 1 - 2*y);

  return sqr(p_cf->Coupling(scale, 0)) * value; //* JIF(z, y, eta, scale);
}

double LF_VV1S0_Quarkonia_IF::OverIntegrated(const double zmin, const double zmax, const double scale, const double xbj) {
  m_zmin = std::max(zmin,ZMIN);
  m_zmax = std::min(zmax,1.);
  m_Jmax = m_flavs[0].Kfcode() < 3 ? 5. : 1.;
  const double m = Flavour((m_flavs[2].Kfcode() / 10) % 10).Mass(true);
  const double prefactor = 4./27*GetLDME(m_flavs[2].Kfcode())/m;
  return prefactor * sqr(p_cf->MaxCoupling(0)) / sqr(2*m) * (zmax - zmin) * m_Jmax;
}

double LF_VV1S0_Quarkonia_IF::OverEstimated(const double z, const double y) {
  const double m = Flavour((m_flavs[2].Kfcode() / 10) % 10).Mass(true);
  return sqr(p_cf->MaxCoupling(0)) /sqr(2*m) * m_Jmax;
}

double LF_VV1S0_Quarkonia_IF::Z() {
  return m_zmin + (m_zmax - m_zmin) * ran->Get();
}

double LF_VV1S0_Quarkonia_II::operator()(const double z, const double y, const double eta, const double scale, const double Q2) {
  const double mj  =  m_flavs[2].Mass(true);
  const double mb  =  m_flspec.Mass(true);
  const double muj2 = sqr(m_flavs[2].Mass(true))/Q2;
  const double mub2 = sqr(m_flspec.Mass(true))/Q2;
  const double sab = Q2*((1 - muj2)/z - mub2*(1-z)/z);
  const double saj = sqr(mj) - y * (sab - sqr(mb));
  const double M = mj; //2*Flavour((m_flavs[2].Kfcode() / 10) % 10).Mass(true);
  double value = 0;
  value  = sqr(saj) + sqr(sqr(M)) - 2 * (1 - z) * (saj + sqr(M))*saj + 2*sqr(saj*(1-z));
  value /= (-saj)*sqr(saj-sqr(M));
  value = (1 - 3*z + 4 * sqr(z) - 2 * cube(z))/(scale);
  // saj is negative, so we need to flip the sign of the denominator to get a positive value.
  // This is because another 1/saj is factored out as part of the Sudakov
  // msg_Out() << "z = " << z << ",  y = " << y << ",  scale = " << scale << ",  saj = " << saj << ", M = " << mj << ", v = " << value << std::endl;
  // value *=  fabs(1. / ( sqr(mj) / (sab - sqr(mb)) / y - 1)); // / ( 1. + sqr(M)/scale); //(1 - z - y) / (1 - z - 2*y); 
  // PRINT_VAR(value);
  return 1.;
  return sqr(p_cf->Coupling(scale, 0)) * value; // * JII(z, y, eta, scale);
}

double LF_VV1S0_Quarkonia_II::OverIntegrated(const double zmin, const double zmax, const double scale, const double xbj) {
  m_zmin = std::max(zmin,ZMIN);
  m_zmax = std::min(zmax,ZMAX);
  // II needs a larger safety envelope in the small-z/small-y corner.
  m_Jmax = 10 * (m_flavs[0].Kfcode() < 3 ? 5. : 1.);
  const double m = m_flavs[2].Mass(true);
  //Flavour((m_flavs[2].Kfcode() / 10) % 10).Mass(true);
  const double prefactor = 4./27*GetLDME(m_flavs[2].Kfcode())/m;
  return prefactor * sqr(p_cf->MaxCoupling(0)) / sqr(m) * (m_zmax - m_zmin) * m_Jmax;
}

double LF_VV1S0_Quarkonia_II::OverEstimated(const double z, const double y) {
  const double m = m_flavs[2].Mass(true); //Flavour((m_flavs[2].Kfcode() / 10) % 10).Mass(true);
  return 1.;
  return sqr(p_cf->MaxCoupling(0)) * ( 5. ) / sqr(m)  * m_Jmax;
}

double LF_VV1S0_Quarkonia_II::Z() {
  return (m_zmax - m_zmin) * ran->Get() + m_zmin;
}

double LF_V1S0V_Quarkonia_FF::operator()(const double z, const double y, const double eta, const double scale, const double Q2) { 
    const double mi  =  m_flavs[1].Mass();
    const double mk  =  m_flspec.Mass();
    const double mui2 = mi/Q2;
    const double muk2 = sqr(m_flspec.Mass(true))/Q2;
    const double sij = y * (Q2 - sqr(mk)) + (1-y)*(sqr(mi));
    const double M = 2*Flavour((m_flavs[1].Kfcode() / 10) % 10).Mass(true);
    double value = 0;
    value  = sqr(sij) + sqr(sqr(M)) - 2 * (z) * (sij + sqr(M))*sij + 2*sqr(sij*(z));
    value /= sij*sqr(sij-sqr(M));
    value *= 1. / ( (1 - mui2 - muk2) + 1./ y * ( mui2 ) );
    value *= 1. / (1 + sqr(1-z) * sqr(mi) / scale);
    return sqr(p_cf->Coupling(scale, 0)) * value * JFF(y, 0, mui2, muk2, 0);
}

double LF_V1S0V_Quarkonia_FF::OverIntegrated(const double zmin, const double zmax, const double scale, const double xbj) {
  m_zmin = std::max(zmin,0.);
  m_zmax = std::min(zmax,1.-ZMIN);
  const double m = Flavour((m_flavs[2].Kfcode() / 10) % 10).Mass(true);
  const double prefactor = 4./27*GetLDME(m_flavs[1].Kfcode())/m;
  return prefactor * sqr(p_cf->MaxCoupling(0)) / sqr(2*m) * (m_zmax - m_zmin);
}

double LF_V1S0V_Quarkonia_FF::OverEstimated(const double z, const double y) {
  const double m = Flavour((m_flavs[1].Kfcode() / 10) % 10).Mass(true);
  return sqr(p_cf->MaxCoupling(0)) /sqr(2*m);
}

double LF_V1S0V_Quarkonia_FF::Z() {
  return m_zmin + (m_zmax - m_zmin) * ran->Get();
}

double LF_V1S0V_Quarkonia_FI::operator()(const double z, const double y, const double eta, const double scale, const double Q2) {
  const double mi  =  m_flavs[1].Mass();
  const double mk  =  m_flspec.Mass();
  const double mui2 = mi/Q2;
  const double muk2 = sqr(m_flspec.Mass(true))/Q2;
  const double sij = (Q2 + sqr(mi)) * y / (1 - y) + sqr(mi);
  const double M = 2*Flavour((m_flavs[1].Kfcode() / 10) % 10).Mass(true);
  double value = 0;
  value  = sqr(sij) + sqr(sqr(M)) - 2 * (z) * (sij + sqr(M))*sij + 2*sqr(sij*(z));
  value /= sij*sqr(sij-sqr(M));
  value *= 1. / ( (1 - mui2 - muk2) + 1./ y * ( mui2 ) );
  value *= 1. / (1 + sqr(1-z) * sqr(mi) / scale);
  return sqr(p_cf->Coupling(scale, 0)) * value * JFI(y, eta, scale);
}

double LF_V1S0V_Quarkonia_FI::OverIntegrated(const double zmin, const double zmax, const double scale, const double xbj) {
  m_zmin = std::max(zmin,ZMIN);
  m_zmax = std::min(zmax,1.);
  m_Jmax = 5.;//?
  const double m = Flavour((m_flavs[1].Kfcode() / 10) % 10).Mass(true);
  const double prefactor = 4./27*GetLDME(m_flavs[1].Kfcode())/m;
  return prefactor * sqr(p_cf->MaxCoupling(0)) / sqr(2*m) * (zmax - zmin) * m_Jmax;
}

double LF_V1S0V_Quarkonia_FI::OverEstimated(const double z, const double y) {
  const double m = Flavour((m_flavs[1].Kfcode() / 10) % 10).Mass(true);
  return sqr(p_cf->MaxCoupling(0)) /sqr(2*m) * m_Jmax;
}

double LF_V1S0V_Quarkonia_FI::Z() {
  return m_zmin + (m_zmax - m_zmin) * ran->Get();
}

DECLARE_GETTER(LF_VV1S0_Quarkonia_FF, "VV1S0_Quarkonia", SF_Lorentz, SF_Key);

SF_Lorentz *ATOOLS::Getter<SF_Lorentz, SF_Key, LF_VV1S0_Quarkonia_FF>::operator()(
    const Parameter_Type &args) const {
  if (args.m_col < 0)
    return NULL;
  if ((args.m_mode == 0 && args.p_v->in[0].IntSpin() == 2 &&
       args.p_v->in[1].IntSpin() == 2 && args.p_v->in[2].IntSpin() == 0) ||
      (args.m_mode == 1 && args.p_v->in[0].IntSpin() == 2 &&
       args.p_v->in[1].IntSpin() == 0 && args.p_v->in[2].IntSpin() == 2)) {
    switch (args.m_type) {
    case cstp::FF: return new LF_VV1S0_Quarkonia_FF(args);
    case cstp::FI: return new LF_VV1S0_Quarkonia_FI(args);
    case cstp::IF: return new LF_VV1S0_Quarkonia_IF(args);
    case cstp::II: return new LF_VV1S0_Quarkonia_II(args);
    case cstp::none:
      break;
    }
  }
if ((args.m_mode == 0 && args.p_v->in[0].IntSpin() == 2 &&
       args.p_v->in[1].IntSpin() == 0 && args.p_v->in[2].IntSpin() == 2) ||
      (args.m_mode == 1 && args.p_v->in[0].IntSpin() == 2 &&
       args.p_v->in[1].IntSpin() == 2 && args.p_v->in[2].IntSpin() == 0)) {
    switch (args.m_type) {
    case cstp::FF: return new LF_V1S0V_Quarkonia_FF(args);
    case cstp::FI: return new LF_V1S0V_Quarkonia_FI(args);
    case cstp::IF: NULL;
    case cstp::II: NULL;
    case cstp::none:
      break;
    }
  }
  return NULL;
}

void ATOOLS::Getter<SF_Lorentz, SF_Key, LF_VV1S0_Quarkonia_FF>::PrintInfo(
    std::ostream &str, const size_t width) const {
  str << "VV1S0 lorentz functions";
}
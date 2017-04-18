#include "AHADIC++/Tools/Z_Selector.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;


Z_Selector::Z_Selector() {}

void Z_Selector::Init()
{
  m_alpha[0] = hadpars->Get("alphaL");
  m_beta[0]  = hadpars->Get("betaL");
  m_gamma[0] = hadpars->Get("gammaL");
  m_alpha[1] = hadpars->Get("alphaH");
  m_beta[1]  = hadpars->Get("betaH");
  m_gamma[1] = hadpars->Get("gammaH");
  m_alpha[2] = hadpars->Get("alphaD");
  m_beta[2]  = hadpars->Get("betaD");
  m_gamma[2] = hadpars->Get("gammaD");
}

Z_Selector::~Z_Selector() {}

double Z_Selector::operator()(const double & zmin,const double & zmax,
			      const Flavour & flav) {
  m_flav  = flav;
  m_flcnt = 0;
  if (m_flav==Flavour(kf_b) || m_flav==Flavour(kf_b).Bar() ||
      m_flav==Flavour(kf_c) || m_flav==Flavour(kf_c).Bar()) {
    m_flcnt = 1;
  }
  else if (m_flav.IsDiQuark()) {
    m_flcnt = 2;
  }
  double m(m_flav.HadMass()), m2(m*m);
  // Maximum xmax of function  f(x) = x^alpha * (1-x)^beta   in [0,1]
  // did not think (yet) about making it mass-dependent - use m_gamma for it
  double xmax = m_alpha[m_flcnt]/(m_alpha[m_flcnt]+m_beta[m_flcnt]);
  double norm = WeightFunction(xmax);
  double ztest(-1.), weight;
  do {
    ztest  = zmin+(zmax-zmin)*ran->Get();
    weight = WeightFunction(ztest)*MassFunction(m2,ztest)/norm;
  } while (weight<ran->Get());
  //msg_Out()<<"   "<<METHOD<<"("<<flav<<", zmax = "<<zmax<<") cnt = "<<m_flcnt
  //	   <<", alpha = "<<m_alpha[m_flcnt]<<", beta = "<<m_beta[m_flcnt]
  //	   <<" --> z = "<<ztest<<".\n";
  return ztest;
}

double Z_Selector::WeightFunction(const double & z) {
  return (pow(z,m_alpha[m_flcnt])*pow(1.-z,m_beta[m_flcnt]));
}

double Z_Selector::MassFunction(const double & m2, const double & z) {
  return exp(-m2*m_gamma[m_flcnt]/z)/exp(-m2*m_gamma[m_flcnt]);
}


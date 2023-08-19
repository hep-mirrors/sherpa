#include "AHADIC++/Tools/KT_Selector.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;

KT_Selector::KT_Selector() {}

KT_Selector::~KT_Selector() {}

void KT_Selector::Init(const bool & isgluon) {
  m_sigma = hadpars->GetVec("kT_0");
  //m_isgluon = isgluon;
  // m_sigma.push_back(hadpars->Get("kT_0"));
  // m_sigma.push_back(hadpars->Get("kT_0")*1.5);
  // m_sigma.push_back(hadpars->Get("kT_0")*0.75);
  // m_sigma.push_back(hadpars->Get("kT_0"));
  //m_sigma2  = sqr(m_sigma);
}

double KT_Selector::operator()(const double & ktmax,const double & M2) {
  double kttest(-1.);
  do {
    kttest = dabs(m_sigma[0] * ran->GetGaussian());
  } while (kttest>ktmax);

  const double gaussian = Gaussian(kttest, m_sigma[0]);
  const double norm     = Erf(ktmax, m_sigma[0]);
  const double p0       = gaussian / norm;

  variation_weights.resize(m_sigma.size());
  tmp_variation_weights.resize(m_sigma.size());
  for(int i{0}; i<m_sigma.size(); ++i) {
    double g = Gaussian(kttest, m_sigma[i]);
    double n = Erf(ktmax, m_sigma[i]);
    double p = g / n;
    tmp_variation_weights[i] = p / p0;
  }
  return kttest;
}

double KT_Selector::Gaussian(const double kt, const double s) {
  const double s2 = s*s;
  return 2. / (std::sqrt(2*M_PI*s2)) * std::exp(-0.5 * kt*kt/ s2);
}

double KT_Selector::Inv_Gaussian(const double kt) {
  // Inverse gaussian for positive results, mean = 0, sigma = 1;
  return std::sqrt(std::abs(std::log(kt * std::sqrt(2*M_PI))));
}

double KT_Selector::Erf(const double kt, const double s) {
  // compute the part of the gaussian we chop of with ktmax
  return std::erf( kt / std::sqrt(2*s*s));
  //return 0.5 * ( 1.+std::erf( kt / std::sqrt(2) ));
}

double KT_Selector::WeightFunction(const double & kt) { return 1.; }

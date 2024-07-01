#include "AHADIC++/Tools/KT_Selector.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;

KT_Selector::KT_Selector() {
  std::cout << "Init KT_Selector" << std::endl;
}

KT_Selector::~KT_Selector() {}

void KT_Selector::Init() {
  const auto _tmp = hadpars->GetVec("kT_0");
  //m_sigma = hadpars->GetVec("kT_0");
  for (auto kt : _tmp)
    //m_sigma.push_back(1.0);
    m_sigma.push_back(kt);
}

double KT_Selector::Select_kt(const double ktmax) {
  double kt, kt_range {ktmax};
  int it {0};
  do {
    kt = ran->Get()*kt_range;
    auto sel_wgt = Gaussian(kt, m_sigma[0]) / 100;
    if(ran->Get() < sel_wgt) {
      kt_accepted(kt);
      break;
    }
    ++it;
    if(it % 10000 == 0)
      std::cout << "Z selection requires many iterations: " << it << " already." << std::endl;
    kt_rejected(kt);
  } while (true);
  return kt;
}

void KT_Selector::kt_accepted(double kt) {
  const auto wgt_old = Gaussian(kt, m_sigma[0]);
  for (int i{0}; i<m_sigma.size(); i++) {
    const auto wgt_new {Gaussian(kt, m_sigma[i])};
    tmp_variation_weights[i] *= wgt_new / wgt_old;
  }
}

void KT_Selector::kt_rejected(double kt) {
  const auto wgt_old = Gaussian(kt, m_sigma[0]) / 100.;
  for (int i{0}; i<m_sigma.size(); i++) {
    const auto wgt_new {Gaussian(kt, m_sigma[i]) / 100};
    tmp_variation_weights[i] *= (1. - wgt_new) / (1. - wgt_old);
  }
}


double KT_Selector::operator()(const double & ktmax) {
  double kttest(-1.);
  // do {
  //   kttest = dabs(m_sigma[0] * ran->GetGaussian());
  // } while (kttest>ktmax);
  std::fill(tmp_variation_weights.begin(), tmp_variation_weights.end(), 1);
  variation_weights.resize(m_sigma.size());
  tmp_variation_weights.resize(m_sigma.size());

  kttest = Select_kt(ktmax);
  // const double gaussian = Gaussian(kttest, m_sigma[0]);
  // const double norm     = Erf(ktmax, m_sigma[0]);
  // const double p0       = gaussian / norm;

  // for(int i{0}; i<m_sigma.size(); ++i) {
  //   double g = Gaussian(kttest, m_sigma[i]);
  //   double n = Erf(ktmax, m_sigma[i]);
  //   double p = g / n;
  //   variation_weights[i] *= p / p0;
  //   tmp_variation_weights[i] = p / p0;
  // }
  std::cout << "DEBUG: KT: " << kttest << " " << tmp_variation_weights << std::endl;
  // TODO needs fixing
  // accepted();
  return kttest;
}

double KT_Selector::Gaussian(const double kt, const double s) {
  const double s2 = s*s;
  return 2. / (std::sqrt(2*M_PI*s2)) * std::exp(-0.5 * kt*kt/ s2);
}


double KT_Selector::Erf(const double kt, const double s) {
  // compute the part of the gaussian we chop of with ktmax
  return std::erf( kt / std::sqrt(2*s*s));
}

double KT_Selector::WeightFunction(const double & kt) { return 1.; }

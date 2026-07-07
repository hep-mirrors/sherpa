#include "AHADIC++/Tools/KT_Selector.H"
#include "AHADIC++/Tools/Ahadic_Reweighting.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;

KT_Selector::KT_Selector(Ahadic_Reweighting * reweighting) :
  p_reweighting(reweighting) {}

KT_Selector::~KT_Selector() {}

void KT_Selector::Init() {
  m_sigma = p_reweighting->GetVariationVector("kT_0");
}

double KT_Selector::SelectKT(const double ktmax) {
  static const int max_iterations = 100000;
  double kt_range {ktmax};
  for (int it{0}; it<max_iterations; ++it) {
    double kt = ran->Get()*kt_range;
    auto sel_wgt = Gaussian(kt, m_sigma[0]) / 100;
    if(ran->Get() < sel_wgt) {
      KTAccepted(kt);
      return kt;
    }
    KTRejected(kt);
  }
  msg_Error() << METHOD << ": kt selection failed after " << max_iterations
              << " iterations (ktmax=" << ktmax << ")\n";
  return -1.;
}

void KT_Selector::KTAccepted(const double kt) {
  if (!p_reweighting->Active()) return;
  std::vector<double> probs(m_sigma.size());
  for (size_t i{0}; i<m_sigma.size(); i++) {
    probs[i] = Gaussian(kt, m_sigma[i]);
  }
  p_reweighting->KTSelectionReweighting(true, probs);
}

void KT_Selector::KTRejected(const double kt) {
  if (!p_reweighting->Active()) return;
  std::vector<double> probs(m_sigma.size());
  for (size_t i{0}; i<m_sigma.size(); i++) {
    probs[i] = Gaussian(kt, m_sigma[i]) / 100.;
  }
  p_reweighting->KTSelectionReweighting(false, probs);
}


double KT_Selector::operator()(const double & ktmax) {
  double kttest(-1.);
  p_reweighting->ResetKTSelectionWeights();
  kttest = SelectKT(ktmax);
  // The kt weights are never committed, so KT_0 variations remain incomplete.
  // An exact alternative to the trial-by-trial reweighting above would be the
  // erf-normalised density ratio of the truncated Gaussian:
  //   p_i = Gaussian(kttest, m_sigma[i]) / Erf(ktmax, m_sigma[i]);
  //   weight_i = p_i / p_0;
  // TODO needs fixing
  // p_reweighting->AcceptKTSelectionWeights();
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

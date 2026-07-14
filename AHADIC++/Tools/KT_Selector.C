#include "AHADIC++/Tools/KT_Selector.H"
#include "AHADIC++/Tools/Hadronisation_Reweighting.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;

KT_Selector::KT_Selector(Hadronisation_Reweighting * reweighting) :
  p_reweighting(reweighting) {}

KT_Selector::~KT_Selector() {}

void KT_Selector::Init() {
  m_sigma = p_reweighting->GetVariationVector("kT_0");
  m_ktmax = p_reweighting->GetVariationVector("kT_max");
  m_n_variations = p_reweighting->NumberOfVariations();
}

double KT_Selector::SelectKT(const double ktmax) {
  double kt;
  do { kt = ran->Get()*ktmax;
  } while (ran->Get() > Gaussian(kt, m_sigma[0]) / Gaussian(0., m_sigma[0]));
  return kt;
}

double KT_Selector::operator()(const double & ktmax, const bool vary_ptmax) {
  const double kt = SelectKT(ktmax);
  if (p_reweighting->Active()) {
    std::vector<double> probs(m_n_variations);
    probs[0] = Gaussian(kt, m_sigma[0]) / Erf(ktmax, m_sigma[0]);
    for (size_t ivar=1; ivar<m_n_variations; ++ivar) {
      const double L = vary_ptmax ? Min(m_ktmax[ivar], ktmax) : ktmax;
      probs[ivar] = kt > L ? 0. : Gaussian(kt, m_sigma[ivar]) / Erf(L, m_sigma[ivar]);
    }
    p_reweighting->KTSelectionReweighting(probs);
  }
  return kt;
}

double KT_Selector::Gaussian(const double kt, const double s) {
  const double s2 = s*s;
  return 2. / (std::sqrt(2*M_PI*s2)) * std::exp(-0.5 * kt*kt/ s2);
}


double KT_Selector::Erf(const double kt, const double s) {
  // compute the part of the gaussian we chop of with ktmax
  return std::erf( kt / std::sqrt(2*s*s));
}

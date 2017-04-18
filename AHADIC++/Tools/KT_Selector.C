#include "AHADIC++/Tools/KT_Selector.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;

KT_Selector::KT_Selector() {}

KT_Selector::~KT_Selector() {}

void KT_Selector::Init() {
  m_sigma  = hadpars->Get("kt_o");
  m_sigma2 = sqr(m_sigma);
}

double KT_Selector::operator()(const double & Emax) {
  double kttest(-1.);
  do {
    kttest = Emax*ran->Get();
  } while (WeightFunction(kttest)<ran->Get());
  return kttest;
}

double KT_Selector::WeightFunction(const double & kt) {
  return exp(-kt*kt/m_sigma2);
}

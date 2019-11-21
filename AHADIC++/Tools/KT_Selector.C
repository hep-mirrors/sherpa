#include "AHADIC++/Tools/KT_Selector.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;

KT_Selector::KT_Selector() {}

KT_Selector::~KT_Selector() {}

void KT_Selector::Init(const bool & isgluon) {
  m_isgluon = isgluon;
  m_sigma   = hadpars->Get("kt_o");
  m_sigma2  = sqr(m_sigma);
}

double KT_Selector::operator()(const double & ktmax,const double & M2) {
  double kttest(-1.);
  m_sig2 = M2*m_sigma2;
  do {
    kttest = ktmax*ran->Get();
  } while (ran->Get() > WeightFunction(kttest)); 
  return kttest;
}

double KT_Selector::WeightFunction(const double & kt) {
  return exp(-sqr(kt)/m_sig2);
}

#include "PHOTONS++/MEs/PHOTONS_ME_Base.H"
#include "MODEL/Main/Model_Base.H"

using namespace PHOTONS;
using namespace ATOOLS;
using namespace std;

PHOTONS_ME_Base::PHOTONS_ME_Base() {
  m_alpha = MODEL::s_model->ScalarConstant("alpha_QED(0)");
  m_e     = sqrt(4*M_PI*m_alpha);
  m_sW    = sqrt(MODEL::s_model->ScalarConstant("sin2_thetaW"));
  m_cW    = sqrt(MODEL::s_model->ScalarConstant("cos2_thetaW"));
  m_i     = Complex(0.,1.);
  p_boost = NULL;
}

PHOTONS_ME_Base::~PHOTONS_ME_Base() {
  if (p_boost != NULL) {
    delete p_boost, delete p_rot;
  }
}




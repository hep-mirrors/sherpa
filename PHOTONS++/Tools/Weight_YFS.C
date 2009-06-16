#include "PHOTONS++/Tools/Weight_YFS.H"

using namespace ATOOLS;
using namespace PHOTONS;
using namespace std;

Weight_YFS::Weight_YFS(const Particle_Vector newdip, const Particle_Vector olddip, const double omegaB, const double nbar) {
  m_weight    = 1;
  m_maxweight = 1;
  m_Y         = YFS_Form_Factor(newdip,omegaB).Get();     // in dipole-CMS
  m_Ymax      = YFS_Form_Factor(olddip,omegaB).Get();     // in dipole-CMS
  m_nbar      = nbar;
  CalculateWeight();
  CalculateMax();
}

Weight_YFS::~Weight_YFS() {
}

void Weight_YFS::CalculateWeight() {
  m_weight = exp( m_Y + m_nbar );
}

void Weight_YFS::CalculateMax() {
  m_maxweight = exp( m_Ymax + m_nbar );
}

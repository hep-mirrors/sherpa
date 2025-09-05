#include "PHOTONS++/Tools/Weight_YFS.H"
#include "PHOTONS++/Tools/YFS_Form_Factor.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Org/Message.H"
#include "PHOTONS++/Main/Photons.H"

using namespace ATOOLS;
using namespace PHOTONS;
using namespace std;

Weight_YFS::Weight_YFS
(const Particle_Vector& newdip, const Particle_Vector& olddip,
 const double& omegaB, const double& nbar) {
  DEBUG_FUNC("\\Omega_min="<<omegaB);
  m_Y         = YFS_Form_Factor(newdip,omegaB).Get();     // in dipole-CMS
  m_Ymax      = YFS_Form_Factor(olddip,omegaB).Get();     // in dipole-CMS
  m_nbar      = nbar;
  msg_Debugging()<<"Y="<<m_Y<<", Ymax="<<m_Ymax<<", nbar="<<nbar<<std::endl;
  CalculateWeight();
  CalculateMax();
}

Weight_YFS::~Weight_YFS() {
}

void Weight_YFS::CalculateWeight() {
  if      (Photons::s_expandtoorder==0) m_weight = exp( m_nbar );
  else if (Photons::s_expandtoorder==1) m_weight = (1.+m_Y)*exp(m_nbar);
  else if (Photons::s_expandtoorder==2) m_weight = (1.+m_Y+0.5*m_Y*m_Y)*exp(m_nbar);
  else                                  m_weight = exp( m_Y + m_nbar );
}

void Weight_YFS::CalculateMax() {
  if      (Photons::s_expandtoorder==0) m_maxweight = exp( m_nbar );
  else if (Photons::s_expandtoorder==1) m_maxweight = (1.+m_Ymax)*exp(m_nbar);
  else if (Photons::s_expandtoorder==2) m_maxweight = (1.+m_Ymax+0.5*m_Ymax*m_Ymax)*exp(m_nbar);
  else                                  m_maxweight = exp( m_Ymax + m_nbar );
}

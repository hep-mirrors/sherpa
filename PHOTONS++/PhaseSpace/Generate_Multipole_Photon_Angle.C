#include "PHOTONS++/PhaseSpace/Generate_Multipole_Photon_Angle.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Particle.H"
#include "PHOTONS++/PhaseSpace/Generate_Dipole_Photon_Angle.H"


using namespace PHOTONS;
using namespace ATOOLS;
using namespace std;

Generate_Multipole_Photon_Angle::Generate_Multipole_Photon_Angle
(const Particle_Vector& dip, const std::vector<double>& nbars) :
  m_dipole(dip), m_nbar(0.), m_nbars(nbars), m_theta(0.), m_phi(0.)
{
  for (unsigned int i=0; i<nbars.size(); i++) m_nbar = m_nbar + abs(nbars[i]);
  GenerateMultipoleAngle();
}

Generate_Multipole_Photon_Angle::~Generate_Multipole_Photon_Angle() {
}

double Generate_Multipole_Photon_Angle::CalculateBeta(const Vec4D& p) {
  return Vec3D(p).Abs()/p[0];
}

void Generate_Multipole_Photon_Angle::GenerateMultipoleAngle() {
  // choose according to which eikonal factor angle should be generated
  // p_i and p_j

  double rdm     = ran.Get();
  double s       = 0.;
  unsigned int k = 0;
  for (unsigned int i=0; i<m_nbars.size(); i++) {
    s = s + abs(m_nbars[i])/m_nbar;
    if (rdm < s) { k = i; break; }
  }
  IndexLookup(k);

  // generate this angle
  // generate massless unit vector in p_i-p_j rest frame then transform back
  // can be later streched to be massless vector of energy E

  Vec4D p1 = m_dipole[m_i]->Momentum();
  Vec4D p2 = m_dipole[m_j]->Momentum();

  Generate_Dipole_Photon_Angle gdpa(p1,p2);
  m_theta = gdpa.GetTheta();
  m_phi   = gdpa.GetPhi();

#ifdef PHOTONS_DEBUG
  msg_Info()<<"dipole generated: "<<k<<" -> "<<m_i<<" "<<m_j<<endl;
  msg_Info()<<"theta:        "<<m_theta<<" phi: "<<m_phi<<endl;
#endif
}

void Generate_Multipole_Photon_Angle::IndexLookup(unsigned int k) {
  switch (k) {
    case  0 : m_i = 0;  m_j = 1;  break;
    case  1 : m_i = 0;  m_j = 2;  break;
    case  2 : m_i = 1;  m_j = 2;  break;
    case  3 : m_i = 0;  m_j = 3;  break;
    case  4 : m_i = 1;  m_j = 3;  break;
    case  5 : m_i = 2;  m_j = 3;  break;
    case  6 : m_i = 0;  m_j = 4;  break;
    case  7 : m_i = 1;  m_j = 4;  break;
    case  8 : m_i = 2;  m_j = 4;  break;
    case  9 : m_i = 3;  m_j = 4;  break;
    case 10 : m_i = 0;  m_j = 5;  break;
    case 11 : m_i = 1;  m_j = 5;  break;
    case 12 : m_i = 2;  m_j = 5;  break;
    case 13 : m_i = 3;  m_j = 5;  break;
    case 14 : m_i = 4;  m_j = 5;  break;
    case 15 : m_i = 0;  m_j = 6;  break;
    case 16 : m_i = 1;  m_j = 6;  break;
    case 17 : m_i = 2;  m_j = 6;  break;
    case 18 : m_i = 3;  m_j = 6;  break;
    case 19 : m_i = 4;  m_j = 6;  break;
    case 20 : m_i = 5;  m_j = 6;  break;
    case 21 : m_i = 0;  m_j = 7;  break;
    case 22 : m_i = 1;  m_j = 7;  break;
    case 23 : m_i = 2;  m_j = 7;  break;
    case 24 : m_i = 3;  m_j = 7;  break;
    case 25 : m_i = 4;  m_j = 7;  break;
    case 26 : m_i = 5;  m_j = 7;  break;
    case 27 : m_i = 6;  m_j = 7;  break;
    case 28 : m_i = 0;  m_j = 8;  break;
    case 29 : m_i = 1;  m_j = 8;  break;
    case 30 : m_i = 2;  m_j = 8;  break;
    case 31 : m_i = 3;  m_j = 8;  break;
    case 32 : m_i = 4;  m_j = 8;  break;
    case 33 : m_i = 5;  m_j = 8;  break;
    case 34 : m_i = 6;  m_j = 8;  break;
    case 35 : m_i = 7;  m_j = 8;  break;
    default : m_i = 0;  m_j = 0;  break;
  }
}

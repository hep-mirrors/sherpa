# include "PHOTONS++/Tools/Generate_One_Photon.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Particle.H"
#include "PHOTONS++/PhaseSpace/Generate_Dipole_Photon_Angle.H"
#include "PHOTONS++/PhaseSpace/Generate_Multipole_Photon_Angle.H"

using namespace PHOTONS;
using namespace std;

// public members
Generate_One_Photon::Generate_One_Photon
(const double& b1, const double& b2, const double& E, Dipole_Type::code dtype) :
  m_w(E), m_dtype(dtype)
{
  // dipole radiation
  if (m_dtype == Dipole_Type::ff) {
    m_beta1 = b1;
    m_beta2 = b2;
  }
  else if (m_dtype == Dipole_Type::fi) {
    // define in different order, because b1 in -z direction in Dipole_FI
    m_beta1 = b2;
    m_beta2 = b1;
  }
  else if (m_dtype == Dipole_Type::ii) {
    m_beta1 = b1;
    m_beta2 = b2;
  }
  Generate_Dipole_Photon_Angle gdpa(m_beta1,m_beta2);
  m_theta = gdpa.GetTheta();
  m_phi   = gdpa.GetPhi();
  GeneratePhoton();
}

Generate_One_Photon::Generate_One_Photon
(const Particle_Vector& dip, const std::vector<double>& nbars,
 double E, Dipole_Type::code dtype) :
  m_w(E), m_dtype(dtype)
{
  // multipole radiation
  Generate_Multipole_Photon_Angle gmpa(dip,nbars);
  m_phi   = gmpa.GetPhi();
  m_theta = gmpa.GetTheta();
  GeneratePhoton();
}

Generate_One_Photon::~Generate_One_Photon() {
}

// private members
void Generate_One_Photon::GeneratePhotonAngleMassless() {
// Generation of theta for two massless particles
  double num = ran.Get();
  m_theta = acos(sqrt(1.-(sqr(sin(m_delta))/((1.-num)*sqr(sin(m_delta))+num))));
  if (ran.Get()>=0.5) m_theta = M_PI - m_theta;
}

void Generate_One_Photon::GeneratePhoton() {
  p_photon = new Particle(ATOOLS::Particle::Counter()+1,kf_photon,
                          Vec4D(0.,0.,0.,0.),'S');
  p_photon->SetMomentum(Vec4D(  m_w ,
                                m_w * sin(m_theta) * cos(m_phi) ,
                                m_w * sin(m_theta) * sin(m_phi) ,
                                m_w * cos(m_theta)  ));
}

void Generate_One_Photon::SetMinimalPhotonAngle(const double& del) {
  m_delta = del;
  delete p_photon;
  GeneratePhoton();
}

# include "PHOTONS++/PhaseSpace/Generate_Dipole_Photon_Angle.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Particle.H"

using namespace PHOTONS;
using namespace ATOOLS;
using namespace std;

Generate_Dipole_Photon_Angle::Generate_Dipole_Photon_Angle(Vec4D p1, Vec4D p2) {
  // determine location and axes of (p1+p2)-CMS
  Poincare boost(p1+p2);
  boost.Boost(p1);
  boost.Boost(p2);
  Poincare rotate(p1,Vec4D(0.,0.,0.,1.));
  // generate null vector of unit length in that frame
  m_b1 = CalculateBeta(p1);
  m_b2 = CalculateBeta(p2);
  if (!IsEqual(p1,p2)) GenerateDipoleAngle();
  else {
    msg_Error()<<METHOD<<"(): \\beta_1 = 0 / \\beta_2 = 0 !"<<std::endl;
    m_theta = 0.0;
    m_phi   = 2.*M_PI*ran->Get();
  }
  GenerateNullVector();
  // transform to original frame
  rotate.RotateBack(m_dir);
  boost.BoostBack(m_dir);
  // determine angles in that frame
  m_theta = acos(m_dir[3]/m_dir[0]);
  m_phi   = acos(m_dir[1]/(m_dir[0]*sin(m_theta)));
}

Generate_Dipole_Photon_Angle::Generate_Dipole_Photon_Angle
(const double& b1, const double& b2) : m_b1(b1), m_b2(b2) {
  // assume p1 and p2 are in their CMS and aligned along z-axis
  GenerateDipoleAngle();
}

Generate_Dipole_Photon_Angle::~Generate_Dipole_Photon_Angle() {
}

double Generate_Dipole_Photon_Angle::CalculateBeta(const Vec4D& p) {
  return Vec3D(p).Abs()/p[0];
}

void Generate_Dipole_Photon_Angle::GenerateDipoleAngle() {
  // Generation of theta for two massive particles
  double P  = log((1.+m_b1)/(1.-m_b1))
                /(log((1.+m_b1)/(1.-m_b1))+log((1.+m_b2)/(1.-m_b2)));
  while (true) {
    m_c = 0.;
    if (ran->Get() < P) {
      double rnd = ran->Get();
      double a   = 1./m_b1*log((1.+m_b1)/(1.-m_b1));
      m_c        = 1./m_b1*(1.-(1.+m_b1)*exp(-a*m_b1*rnd));
    }
    else {
      double rnd = ran->Get();
      double a   = 1./m_b2*log((1.+m_b2)/(1.-m_b2));
      m_c        = 1./m_b2*((1.-m_b2)*exp(a*m_b2*rnd)-1.);
    }
    double weight = 1.-((1.-m_b1*m_b1)/((1.-m_b1*m_c)*(1.-m_b1*m_c))
                        +(1.-m_b2*m_b2)/((1.+m_b2*m_c)*(1.+m_b2*m_c)))
                       /(2.*(1.+m_b1*m_b2)/((1.-m_b1*m_c)*(1.+m_b2*m_c)));
    if (ran->Get() < weight) break;
  }
  m_theta = acos(m_c);
  m_phi   = 2.*M_PI*ran->Get();
}

void Generate_Dipole_Photon_Angle::GenerateNullVector() {
  m_dir = Vec4D(1., sin(m_theta)*cos(m_phi), sin(m_theta)*sin(m_phi), cos(m_theta));
}

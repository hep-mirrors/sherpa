#include "BEAM/Main/RelicDensity_Kinematics.H"

using namespace BEAM;
using namespace ATOOLS;
using namespace std;

RelicDensity_Kinematics::RelicDensity_Kinematics(std::array<Beam_Base *, 2> beams) :
  Kinematics_Base(beams) {
  InitIntegration();
}

void RelicDensity_Kinematics::InitIntegration() {
  Beam_Parameters parameters;
  double Emin = m_m[0]+m_m[1];   // Equal to mass energy
  double Emax = Emin + parameters("RELIC_DENSITY_EMAX");
  m_smin = sqr(Emin);
  m_smax = sqr(Emax);
  m_S    = sqr(Emin);
  m_on   = true;
  m_exponent[0] = .5;
  m_exponent[1] = 2.;
}

void RelicDensity_Kinematics::AssignKeys(Integration_Info *const info) {
  m_sprimekey.Assign(m_keyid+string("s'"),5,0,info);
  SetLimits();
}

void RelicDensity_Kinematics::SetLimits() {
  m_sprimekey[0] = Max(m_smin, m_sminPS);
  m_sprimekey[1] = m_sprimekey[2] = m_smax;
  m_sprimekey[3] = m_S;
}

bool RelicDensity_Kinematics::operator()(ATOOLS::Vec4D_Vector& moms) {
  m_S = m_sprimekey[3];
  double Eprime = sqrt(m_S);
  if ( m_S<m_sprimekey[0] || m_S>m_sprimekey[1] ||
       m_sprimekey[0]==m_sprimekey[1] ||
       Eprime<m_m[0]+m_m[1]) return false;
  double x  = (m_S+m_m2[0]-m_m2[1])/(2.*m_S);
  double E1 = x*Eprime, E2 = Eprime-E1;
  moms[0] = Vec4D(E1,0.,0.,sqrt(sqr(E1)-m_m2[0]));
  moms[1] = Vec4D(E2,(-1.)*Vec3D(moms[0]));
  p_beams[0]->SetX(1.);
  p_beams[1]->SetX(1.);
  return true;
}


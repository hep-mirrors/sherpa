#include "BEAM/Main/DM_Annihilation_Kinematics.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace BEAM;
using namespace ATOOLS;

DM_Annihilation_Kinematics::DM_Annihilation_Kinematics(std::array<Beam_Base*, 2> beams) :
  Kinematics_Base(beams) {
  InitIntegration();
}

void DM_Annihilation_Kinematics::InitIntegration() {
  Beam_Parameters parameters;
  double Emin = m_m[0]+m_m[1];   // Equal to mass energy
  double Emax = Emin + parameters("RELIC_DENSITY_EMAX");
  m_smin = sqr(Emin);
  m_smax = sqr(Emax);
  m_S    = sqr(Emin);
  m_x[0] = 0.5;
  m_x[1] = 1 - m_x[0];
  // m_cosxi = 0.;
  m_on   = true;
  m_exponent[0] = .5;
  m_exponent[1] = 2.;
}

void DM_Annihilation_Kinematics::AssignKeys(Integration_Info *const info) {
  m_sprimekey.Assign(m_keyid + std::string("s'"), 5, 0, info);
  // TODO: Unsure about correct implementation of the chosen convention, see Collider_Kinematics.C
  m_xkey.Assign(m_keyid + std::string("x"), 3, 0, info);
  m_cosxikey.Assign(m_keyid + std::string("cosXi"), 3, 0, info);
  /////////////////////////////////////////////////////////////////////
  SetLimits();
}

void DM_Annihilation_Kinematics::SetLimits() {
  m_sprimekey[0] = Max(m_smin, m_sminPS);
  m_sprimekey[1] = m_sprimekey[2] = m_smax;
  m_sprimekey[3] = m_S;

  m_cosxikey[0] = -1;
  m_cosxikey[1] = 1;
  m_cosxikey[2] = m_cosxi;

  m_xkey[0] = (m_m[0]>m_m[1]) ? 0.5 - (m_m2[0]-m_m2[1])/(2.*m_S) : 0.5 + (m_m2[0]-m_m2[1])/(2*m_S);
  m_xkey[1] = 1 - m_xkey[0];
  m_xkey[2] = m_x[0]; // define to be beam 1
}

bool DM_Annihilation_Kinematics::operator()(Vec4D_Vector& p) {

  m_S = m_sprimekey[2] = m_sprimekey[3];
  double Eprime = sqrt(m_S);
  if ( m_S<m_sprimekey[0] || m_S>m_sprimekey[1] ||
       m_sprimekey[0]==m_sprimekey[1] ||
       Eprime<m_m[0]+m_m[1]) return false;

  m_cosxi = m_cosxikey[2];
  double sinxi = sqrt(1-sqr(m_cosxi));
  double x  = m_xkey[2];
  double E1 = x*Eprime;
  double E2 = (1-x)*Eprime;
  double p1 = sqrt(sqr(E1)-m_m2[0]);
  double p2 = sqrt(sqr(E2)-m_m2[1]);

  p[0] = Vec4D(E1,0.,0.,p1);
  p[1] = Vec4D(E2,p2*sinxi,0.,p2*m_cosxi);

  m_CMSBoost = Poincare(p[0]+p[1]);

  // initial state momenta
  BoostInCMS(p[0]);
  BoostInCMS(p[1]);

  // cms frame. No bunches so in/out momenta the same
  p_beams[0]->SetInMomentum(p[0]);
  p_beams[1]->SetInMomentum(p[1]);
  p_beams[0]->SetOutMomentum(p[0]);
  p_beams[1]->SetOutMomentum(p[1]);

  rpa->gen.SetEcms(Eprime);

  return true;
}

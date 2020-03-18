#include "BEAM/Main/DM_Annihilation_Kinematics.H"

using namespace BEAM;
using namespace ATOOLS;
using namespace std;

DM_Annihilation_Kinematics::DM_Annihilation_Kinematics(Beam_Base ** beams) :
  Kinematics_Base(beams) {
  InitIntegration();
}

DM_Annihilation_Kinematics::~DM_Annihilation_Kinematics() {}

void DM_Annihilation_Kinematics::InitIntegration() {
  Beam_Parameters parameters;
  double Emin = m_m[0]+m_m[1];   // Equal to mass energy
  double Emax = Emin + parameters("DM_ANNIHILATION_EMAX");
  m_smin = sqr(Emin);
  m_smax = sqr(Emax);
  m_S    = sqr(Emin);
  m_x[0] = 0.5;
  m_x[1] = 1 - m_x[0];
  m_cosxi = 0;
  m_on   = true;
  m_exponent[0] = .5;
  m_exponent[1] = 2.;
}

void DM_Annihilation_Kinematics::AssignKeys(Integration_Info *const info) {
  m_sprimekey.Assign(m_keyid+string("s'"),5,0,info);
  m_xkey.Assign(m_keyid+string("x"),3,0,info);
  m_cosxikey.Assign(m_keyid+string("cos(xi)"),3,0,info);
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

  m_xkey[0] = (m_m[0]>m_m[1]) ? 0.5 + (m_m2[0]-m_m2[1])/(2.*m_S) : 0.5 - (m_m2[0]-m_m2[1])/(2*m_S);
  m_xkey[1] = 1 - m_xkey[0];
  m_xkey[2] = m_x[0]; // define to be beam 1
}

bool DM_Annihilation_Kinematics::operator()(Vec4D * p) {
  m_S = m_sprimekey[3];
  double Eprime = sqrt(m_S);
  if ( m_S<m_sprimekey[0] || m_S>m_sprimekey[1] ||
       m_sprimekey[0]==m_sprimekey[1] ||
       Eprime<m_m[0]+m_m[1]) return false;

  double E1 = m_xkey[2]*Eprime;
  double E2 = Eprime-E1;
  double p1 = sqrt(sqr(E1)-m_m2[0]);
  double p2 = sqrt(sqr(E2)-m_m2[1]);
  double sinxi = sqrt(1-sqr(m_cosxi));  // take out semicolon for deliberate syntax error, debugging

  p[0] = Vec4D(E1,0.,0.,p1);
  p[1] = Vec4D(E2,p2*sinxi,0,p2*m_cosxi);
  //msg_Out()<<METHOD<<" --> "<<p[0]<<" & "<<p[1]<<"\n";

  p_beams[0]->SetX(m_x[0]);
  p_beams[1]->SetX(m_x[1]);
  return true;
}

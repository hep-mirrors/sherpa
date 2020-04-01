#include "BEAM/Main/DM_Annihilation_Weight.H"
#include "BEAM/Main/Weight_Base.H"
#include "ATOOLS/Org/Message.H"

using namespace BEAM;
using namespace ATOOLS;
using namespace std;

DM_Annihilation_Weight::DM_Annihilation_Weight(Kinematics_Base * kinematics) :
  Weight_Base(kinematics), m_relativistic(true)
{
  Beam_Parameters parameters;
  m_relativistic = parameters.On("DM_RELATIVISTIC");
  m_temperature  = parameters("DM_TEMPERATURE");
  for (size_t i=0;i<2;i++) {
    m_m[i]        = p_kinematics->m(i);
    m_m2[i]       = p_kinematics->m2(i);
    m_BesselK2[i] = cyl_bessel_2(m_m[i]/m_temperature);
    m_w[i]        = m_relativistic? 1./m_m2[i] : 1./(8.*m_m[i]+15.*m_temperature);
  }
  m_norm = (m_w[0]*m_w[1]);
  // Is this normalisation still correct?
  // Might add a check for whether we want weighted or not
  // if (m_relativistic) {
  //   m_norm /= (8.*m_temperature*m_BesselK2[0]*m_BesselK2[1]);
  // }
  // else {
  //   m_norm /= (sqrt(2.*M_PI*pow(m_temperature,3.)*sqrt(m_m[0]*m_m[1])));
  // }
}

DM_Annihilation_Weight::~DM_Annihilation_Weight() {}

void DM_Annihilation_Weight::AssignKeys(Integration_Info *const info) {
  m_sprimekey.Assign(m_keyid+std::string("s'"),5,0,info);
  m_xkey.Assign(m_keyid+string("x"),3,0,info);
  m_cosxikey.Assign(m_keyid+string("cosXi"),3,0,info);
}

bool DM_Annihilation_Weight::Calculate(const double & scale) {
  double s = m_sprimekey[3];
  if (s <= sqr(m_m[0]+m_m[1])) {
    m_weight = 0.;
    return true;
  }

  // Weight for each beam
  double x[2];
  x[0] = m_xkey[2];
  x[1] = 1-x[0];
  for (size_t i=0;i<2;i++) {
    m_w[i] = p_beams[i]->CalculateWeight(x[i],scale);
  }

  // Relative velocity - this is the Lorentz invariant one, may not be correct for XS
  // double lambda = (scale-sqr(m_m[0]+m_m[1]))*(scale-sqr(m_m[0]-m_m[1]));
  // double vrel = sqrt(lambda)/(scale-m_m2[0]-m_m2[1]);
  //
  // ////////////////////////////////////////////////
  // // alternatively:
  // double E1 = m_xkey[2]*sqrt(scale);
  // double E2 = sqrt(scale)-E1;
  // double p1 = sqrt(sqr(E1)-m_m2[0]);
  // double p2 = sqrt(sqr(E2)-m_m2[1]);
  // double sinxi = sqrt(1-sqr(m_cosxikey[2]));
  // Vec3D v1, v2;
  // v1 = Vec3D(0., 0., p1/m_m[0]);
  // v2 = Vec3D(p2*sinxi/m_m[1], 0, p2*m_cosxikey[2]/m_m[1]);
  // double vrel_alt = sqrt((v1-v2).Sqr() - ATOOLS::cross(v1,v2).Sqr()); // Sqr() is to Vec3D what Abs2() is to Vec4D
  // ////////////////////////////////////////////////

  m_weight      = m_norm * m_w[0] * m_w[1];
  //msg_Out()<<METHOD<<"(E = "<<sqrt(s)<<") = "<<m_weight<<"\n";
  return true;
}

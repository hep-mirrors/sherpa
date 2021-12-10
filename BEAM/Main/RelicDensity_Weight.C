#include "BEAM/Main/RelicDensity_Weight.H"
#include "BEAM/Main/Weight_Base.H"
#include "ATOOLS/Org/Message.H"

using namespace BEAM;
using namespace ATOOLS;

RelicDensity_Weight::RelicDensity_Weight(Kinematics_Base * kinematics) :
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
  
  if (m_relativistic) {
    m_norm /= (8.*m_temperature*m_BesselK2[0]*m_BesselK2[1]);
  }
  else {
    m_norm *= sqrt(2/(M_PI*pow(m_temperature,3.)*(m_m[0]+m_m[1])));
  }
}

RelicDensity_Weight::~RelicDensity_Weight() {}

void RelicDensity_Weight::AssignKeys(Integration_Info *const info) {
  m_sprimekey.Assign(m_keyid+std::string("s'"),5,0,info);
}

bool RelicDensity_Weight::Calculate(const double & s) {
  if (s <= sqr(m_m[0]+m_m[1])) {
    m_weight = 0.;
    return true;
  }
  double lambda = (s-sqr(m_m[0]+m_m[1]))*(s-sqr(m_m[0]-m_m[1]));
  double E      = sqrt(s);
  m_weight      = m_norm*lambda;
  if (m_relativistic) {
    m_weight *= cyl_bessel_1(E/m_temperature) / E;
  }
  else {
    long double arg = (E-m_m[0]-m_m[1])/m_temperature;
    m_weight *= pow(s,-3./4.)*(1.+3.*m_temperature/(8.*E))*expl(-arg);
  }
  // msg_Out()<<"s="<<s<<", weight="<<m_weight<<"\n"; //debugging
  return true;
}

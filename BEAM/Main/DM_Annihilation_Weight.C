#include "BEAM/Main/DM_Annihilation_Weight.H"
#include "BEAM/Main/Weight_Base.H"
#include "ATOOLS/Org/Message.H"

using namespace BEAM;
using namespace ATOOLS;

DM_Annihilation_Weight::DM_Annihilation_Weight(Kinematics_Base * kinematics) :
  Weight_Base(kinematics), m_relativistic(true)
{
  auto& s = Settings::GetMainSettings();
  m_relativistic = s["DM_RELATIVISTIC"].Get<bool>();
  m_temperature  = s["DM_TEMPERATURE"].Get<double>();
  for (size_t i=0;i<2;i++) {
    m_m[i]        = p_kinematics->m(i);
    m_m2[i]       = p_kinematics->m2(i);
    m_w[i]        = 1;
  }
  m_norm = 8*sqr(M_PI) * (m_w[0]*m_w[1]);
}

DM_Annihilation_Weight::~DM_Annihilation_Weight() {}

void DM_Annihilation_Weight::AssignKeys(Integration_Info *const info) {
  m_sprimekey.Assign(m_keyid+std::string("s'"),5,0,info);
  // TODO: Unsure about correct implementation of the chosen convention, see Collider_Kinematics.C
  m_xkey.Assign(m_keyid + std::string("xDM"), 3, 0, info);
  m_cosxikey.Assign(m_keyid + std::string("cosXi"), 3, 0, info);
}

bool DM_Annihilation_Weight::Calculate(const double & scale) {
  if (scale <= sqr(m_m[0]+m_m[1])) {
    m_weight = 0.;
    return true;
  }

  // Weight for each beam
  double x[2];
  x[0] = m_xkey[2];
  x[1] = 1-x[0];

  for (size_t i=0;i<2;i++) {
    p_beams[i]->CalculateWeight(x[i],scale);
    m_w[i] = p_beams[i]->Weight();
  }

  m_weight      = m_norm * m_w[0] * m_w[1];
  // m_weight = 1; //TEST
  return true;
}

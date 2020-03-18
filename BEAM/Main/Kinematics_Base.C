#include "BEAM/Main/Kinematics_Base.H"
#include "ATOOLS/Org/Message.H"

using namespace BEAM;
using namespace ATOOLS;

Kinematics_Base::Kinematics_Base(Beam_Base * beams[2]) :
  m_on(false), m_keyid("BEAM::")
{
  for (size_t i=0;i<2;i++) {
    p_beams[i] = beams[i];
    m_m[i]     = p_beams[i]->Bunch().Mass(); m_m2[i] = sqr(m_m[i]);
  }
}

Kinematics_Base::~Kinematics_Base() {}

void Kinematics_Base::SetLimits()
{
  double           m_splimits[3], m_ylimits[2], m_exponent[2];
  ATOOLS::Info_Key m_sprimekey, m_ykey, m_costhetakey;

}


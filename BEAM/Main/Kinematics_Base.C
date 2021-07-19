#include "BEAM/Main/Kinematics_Base.H"
#include "ATOOLS/Org/Message.H"

using namespace BEAM;
using namespace ATOOLS;

Kinematics_Base::Kinematics_Base(Beam_Base * beams[2]) :
  m_on(false), m_keyid("BEAM::"), m_Plab(Vec4D(0.,0.,0.,0.))
{
  for (size_t i=0;i<2;i++) {
    p_beams[i] = beams[i];
    m_m[i]     = p_beams[i]->Bunch().Mass(); m_m2[i] = sqr(m_m[i]);
    m_Plab    += p_beams[i]->InMomentum();
  }
  m_S = m_Plab.Abs2();
}

Kinematics_Base::~Kinematics_Base() {}



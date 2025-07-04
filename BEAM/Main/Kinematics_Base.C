#include "BEAM/Main/Kinematics_Base.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace BEAM;
using namespace ATOOLS;

Kinematics_Base::Kinematics_Base(std::array<Beam_Base *, 2> beams) :
  m_on(false), m_keyid("BEAM::"), m_smin(0.),
  m_smax(sqr(rpa->gen.Ecms())), m_sminPS(0.), m_smaxPS(sqr(rpa->gen.Ecms())),
  m_Plab(Vec4D(0.,0.,0.,0.))
{
  for (size_t i=0;i<2;i++) {
    p_beams[i] = beams[i];
    m_m[i] = p_beams[i]->Bunch().Mass();
    m_m2[i] = sqr(m_m[i]);
    m_Plab += p_beams[i]->InMomentum();
  }
  m_S = m_Plab.Abs2();
}

Kinematics_Base::~Kinematics_Base() = default;

#include "SHERPA/SoftPhysics/Hadron_Rescattering_Handler.H"

using namespace SHERPA;
using namespace HADRON_RESCATTERING;
using namespace ATOOLS;

Hadron_Rescattering_Handler::Hadron_Rescattering_Handler() :
  m_name("Hadron_Rescattering"), m_on(false), 
  m_rescattering(Hadron_Rescatterings(true))
{
  if (!m_on) return;
  m_rescattering.Initialize();
}

Hadron_Rescattering_Handler::~Hadron_Rescattering_Handler() {}


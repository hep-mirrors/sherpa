#include "SHERPA/SoftPhysics/Hadron_Rescattering_Handler.H"

using namespace SHERPA;
using namespace HADRON_RESCATTERING;
using namespace ATOOLS;

Hadron_Rescattering_Handler::Hadron_Rescattering_Handler() :
  p_rescattering(NULL)
{}

Hadron_Rescattering_Handler::~Hadron_Rescattering_Handler() {
  if (p_rescattering) { delete p_rescattering; p_rescattering = NULL; }
}


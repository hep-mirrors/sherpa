#include "SHERPA/SoftPhysics/Colour_Reconnection_Handler.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace SHERPA;
using namespace RECONNECTIONS;
using namespace ATOOLS;
using namespace std;

Colour_Reconnection_Handler::Colour_Reconnection_Handler() :
  p_reconnections(NULL)
{
  Settings& s = Settings::GetMainSettings();
  m_on = s["COLOUR_RECONNECTIONS"].SetDefault(false).Get<bool>();
  p_reconnections = new Reconnection_Handler(m_on);
  p_reconnections->Initialize();
}

Colour_Reconnection_Handler::~Colour_Reconnection_Handler() {
  if (p_reconnections) delete p_reconnections;
}

Return_Value::code Colour_Reconnection_Handler::operator()(Blob_List *const blobs) {
  if (!m_on) return Return_Value::Nothing;
  return (*p_reconnections)(blobs);
}

void Colour_Reconnection_Handler::CleanUp(const size_t & mode) {
  p_reconnections->Reset();
}

void Colour_Reconnection_Handler::Output() { }

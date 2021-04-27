#include "SHERPA/SoftPhysics/Colour_Reconnection_Handler.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace SHERPA;
using namespace RECONNECTIONS;
using namespace ATOOLS;
using namespace std;

Colour_Reconnection_Handler::Colour_Reconnection_Handler() :
  p_reconnections(NULL)
{
  const Scoped_Settings& s =
    Settings::GetMainSettings()["COLOUR_RECONNECTIONS"];
  string mode = s["Mode"].SetDefault(string("Off")).Get<string>(); 
  m_on = (mode!=string("Off"));
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

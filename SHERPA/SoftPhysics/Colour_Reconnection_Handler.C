#include "SHERPA/SoftPhysics/Colour_Reconnection_Handler.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace SHERPA;
using namespace RECONNECTIONS;
using namespace ATOOLS;
using namespace std;

Colour_Reconnection_Handler::Colour_Reconnection_Handler() :
  p_reconnections(NULL)
{
  auto s = Settings::GetMainSettings()["COLOUR_RECONNECTIONS"];
  m_on   = s["MODE"].SetDefault(false).Get<bool>();
  p_reconnections = new Reconnection_Handler(m_on);
  p_reconnections->Initialize();
}

Colour_Reconnection_Handler::~Colour_Reconnection_Handler() {
  if (p_reconnections) delete p_reconnections;
}

Return_Value::code
Colour_Reconnection_Handler::operator()(Blob_List *const blobs) {
  if (!m_on) {
    UpdateStatus(blobs);
    return Return_Value::Nothing;
  }
  return (*p_reconnections)(blobs);
}

void Colour_Reconnection_Handler::UpdateStatus(Blob_List *const blobs) {
  for (Blob_List::iterator blit=blobs->begin();blit!=blobs->end();++blit) {
    if ((*blit)->Has(blob_status::needs_reconnections)) {
      (*blit)->UnsetStatus(blob_status::needs_reconnections);
      (*blit)->SetStatus(blob_status::needs_hadronization);
    }
  }
}

void Colour_Reconnection_Handler::CleanUp(const size_t & mode) {
  p_reconnections->Reset();
}

void Colour_Reconnection_Handler::Output() { }

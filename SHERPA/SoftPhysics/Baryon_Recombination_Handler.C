#include "SHERPA/SoftPhysics/Baryon_Recombination_Handler.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace SHERPA;
using namespace RECOMBINATION;
using namespace ATOOLS;
using namespace std;

Baryon_Recombination_Handler::Baryon_Recombination_Handler() :
  p_recombination(NULL)
{
  // TODO: Tidy up/update input structures at the very end.
  // auto s = Settings::GetMainSettings()["BARYON_RECOMBINATION"];
  // m_on = s["MODE"].SetDefault(false).Get<bool>();
  m_on = false;
  p_recombination = new Recombination_Handler(m_on);
  p_recombination->Initialize();
}

Baryon_Recombination_Handler::~Baryon_Recombination_Handler() {
  if (p_recombination) delete p_recombination;
}

Return_Value::code
Baryon_Recombination_Handler::operator()(Blob_List *const blobs) {
  if (!m_on) return Return_Value::Nothing;
  return (*p_recombination)(blobs);
}

void Baryon_Recombination_Handler::CleanUp(const size_t & mode) {
  p_recombination->Reset();
}

void Baryon_Recombination_Handler::Output() { }

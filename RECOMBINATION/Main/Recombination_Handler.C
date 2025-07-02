#include "RECOMBINATION/Main/Recombination_Handler.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace RECOMBINATION;
using namespace ATOOLS;
using namespace std;

Recombination_Handler::Recombination_Handler(const bool & on) :
  m_on(on)
{}

Recombination_Handler::~Recombination_Handler() {
  if (m_on) {
    msg_Info()<<METHOD<<".\n";
  }
}

void Recombination_Handler::Initialize() {
  if (m_on) {}
}

void Recombination_Handler::Reset() {
  if (m_on) {}
}

Return_Value::code Recombination_Handler::operator()(Blob_List *const blobs,
						    Particle_List *const parts) {
  if (!m_on) return Return_Value::Nothing;
  return Return_Value::Success; 
}

/*
Blob * Recombination_Handler::MakeRecombinationBlob() {
  Blob * blob = new Blob();
  blob->SetType(btp::Fragmentation);
  blob->SetId();
  return blob;
}
*/

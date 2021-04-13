#include "RECONNECTIONS/Main/Reconnection_Handler.H"
#include "RECONNECTIONS/Main/Reconnect_By_Singlet.H"
#include "RECONNECTIONS/Main/Reconnect_Statistical.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace RECONNECTIONS;
using namespace ATOOLS;
using namespace std;

Reconnection_Handler::Reconnection_Handler() :
  m_on(false),
  p_reconnector(new Reconnect_Statistical())
{}

Reconnection_Handler::~Reconnection_Handler() {
  delete p_reconnector;
}

void Reconnection_Handler::Initialize(Default_Reader *const reader) {
  string onoff = reader->GetStringNormalisingNoneLikeValues("COLOUR_RECONNECTIONS",
							    string("Off"));
  m_on = bool(onoff==string("On"));
  if (m_on) p_reconnector->Initialize(reader);
}

void Reconnection_Handler::Reset() {
  if (m_on) p_reconnector->Reset();
}

Return_Value::code Reconnection_Handler::operator()(Blob_List *const blobs,
						    Particle_List *const parts) {
  if (!m_on) return Return_Value::Nothing;
  if ((*p_reconnector)(blobs)) {
    AddReconnectionBlob(blobs);
    p_reconnector->Reset();
    return Return_Value::Success; 
  }
  return Return_Value::New_Event;
}
  

void Reconnection_Handler::AddReconnectionBlob(Blob_List *const blobs) {
  Blob * blob = new Blob();
  blob->AddStatus(blob_status::needs_hadronization);
  blob->SetType(btp::Fragmentation);
  blob->SetId();
  Part_List * particles = p_reconnector->GetParticles();
  while (!particles->empty()) {
    Particle * part = particles->front();
    part->DecayBlob()->AddToOutParticles(part);
    part->SetDecayBlob(NULL);
    blob->AddToInParticles(part);
    particles->pop_front();
  }
  blobs->push_back(blob);
}


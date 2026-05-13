#include "SHERPA/Single_Events/Hadron_Rescattering.H"
#include "SHERPA/SoftPhysics/Hadron_Rescattering_Handler.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/Message.H"

using namespace SHERPA;
using namespace ATOOLS;

Hadron_Rescattering::
Hadron_Rescattering(Hadron_Rescattering_Handler * reschandler) :
  m_on(true),
  p_reschandler(reschandler)
{
  msg_Out()<<METHOD<<"\n";
  m_name = std::string("Hadron Rescattering");
  m_type = eph::Hadronization;
}

Hadron_Rescattering::~Hadron_Rescattering() {}

ATOOLS::Return_Value::code Hadron_Rescattering::Treat(Blob_List* blobs) {
  if (!m_on) return Return_Value::Nothing;
  if (blobs->empty()) {
    msg_Error()<<METHOD<<"("<<blobs<<"):\n"
	       <<"   Blob list contains "<<blobs->size()<<" entries.\n"
	       <<"   Continue and hope for the best.\n";
    return Return_Value::Error;
  }
  msg_Out()<<"=== "<<METHOD<<" for "<<blobs->size()<<" blobs.\n";
  bool treated = false;
  for (size_t blit(0);blit<blobs->size();++blit) {
    Blob* blob=(*blobs)[blit];
    if (blob->Has(blob_status::needs_hadronRescatter)) {
      if (m_on && p_reschandler) {
	p_reschandler->HarvestParticles(blob);
	treated = true;
      }
      blob->UnsetStatus(blob_status::needs_hadronRescatter);
    }
  }
  if (!treated) return Return_Value::Nothing;
  Blob * rblob = nullptr;
  bool resc    = false;
  while ((*p_reschandler)()) {
    blobs->push_back(rblob = p_reschandler->GetBlob());
    msg_Out()<<"   * rescattering blob added: "<<(*rblob)<<"\n";
  }
  return resc ? Return_Value::Success : Return_Value::Nothing;
}

void Hadron_Rescattering::CleanUp(const size_t & mode) {
  p_reschandler->CleanUp(mode);
}

void Hadron_Rescattering::Finish(const std::string &) {}

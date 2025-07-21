#include "SHERPA/Single_Events/Hadron_Rescattering.H"
#include "SHERPA/SoftPhysics/Hadron_Rescattering_Handler.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/Message.H"

using namespace SHERPA;
using namespace ATOOLS;

Hadron_Rescattering::Hadron_Rescattering(Hadron_Rescattering_Handler * reschandler) :
  m_on(true),
  p_reschandler(reschandler)
{}

Hadron_Rescattering::~Hadron_Rescattering() {}

ATOOLS::Return_Value::code Hadron_Rescattering::Treat(Blob_List* blobs) {
  msg_Out()<<METHOD<<" ==============================================\n";
  if (blobs->empty()) {
    msg_Error()<<METHOD<<"("<<blobs<<"):\n"
	       <<"   Blob list contains "<<blobs->size()<<" entries.\n"
	       <<"   Continue and hope for the best.\n";
    return Return_Value::Error;
  }
  if (!m_on) return Return_Value::Nothing;
  for (size_t blit(0);blit<blobs->size();++blit) {
    Blob* blob=(*blobs)[blit];
    if (p_reschandler && blob->Has(blob_status::needs_hadronRescatter)) {
      p_reschandler->HarvestParticles(blob);
      if ((*p_reschandler)()) {
	msg_Out()<<METHOD<<" yields kaboom.\n";
	exit(1);
      }
    }
  }
  return Return_Value::Nothing;
}

void Hadron_Rescattering::CleanUp(const size_t & mode) {
  p_reschandler->CleanUp(mode);
}

void Hadron_Rescattering::Finish(const std::string &) {}

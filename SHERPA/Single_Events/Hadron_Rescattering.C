#include "SHERPA/Single_Events/Hadron_Rescattering.H"
#include "SHERPA/SoftPhysics/Hadron_Rescattering_Handler.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/Message.H"

using namespace SHERPA;
using namespace ATOOLS;

Hadron_Rescattering::Hadron_Rescattering(Hadron_Rescattering_Handler * reschandler) :
  p_reschandler(reschandler)
{}

Hadron_Rescattering::~Hadron_Rescattering() {}

ATOOLS::Return_Value::code Hadron_Rescattering::Treat(Blob_List* blobs) {
  return Return_Value::Nothing;
}
void Hadron_Rescattering::CleanUp(const size_t & mode) {}
void Hadron_Rescattering::Finish(const std::string &) {}

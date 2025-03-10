#include "SHERPA/Single_Events/Rescattering.H"
#include "ATOOLS/Org/Message.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Rescattering::Rescattering(Rescattering_Handler * reschandler) :
  m_on(reschandler->Name()!="None"),
  p_reschandler(reschandler)
{
  m_name = std::string("Rescattering: ")+p_reschandler->Name();
  m_type = eph::Perturbative;
}

Rescattering::~Rescattering() {}

Return_Value::code Rescattering::Treat(ATOOLS::Blob_List* bloblist)
{
  if (bloblist->empty()) {
    msg_Error()<<"Rescattering::Treat("<<bloblist<<"): "<<endl
	       <<"   Blob list contains "<<bloblist->size()<<" entries."<<endl
	       <<"   Continue and hope for the best."<<endl;
    return Return_Value::Error;
  }
  if (!m_on) return Return_Value::Nothing;
  Return_Value::code ret = (*p_reschandler)(bloblist);
  if (ret!=Return_Value::Success && ret!=Return_Value::Nothing)
    THROW(fatal_error,"unexpected (and undefined) result.");
  return ret;
}

void Rescattering::CleanUp(const size_t & mode) {}

void Rescattering::Finish(const std::string &) {}


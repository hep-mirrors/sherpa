#include "SHERPA/Single_Events/NP_Colour_Reconnections.H"
#include "ATOOLS/Org/Message.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

NP_Colour_Reconnections::
NP_Colour_Reconnections(Colour_Reconnection_Handler * reconnections) :
  m_on(reconnections->On()),
  p_reconnectionhandler(reconnections)
{
  m_name = string("NP_Colour_Reconnections: ")+string(m_on ? "On" : "Off");
  m_type = eph::Hadronization;
}

NP_Colour_Reconnections::~NP_Colour_Reconnections() {}

Return_Value::code NP_Colour_Reconnections::Treat(Blob_List* bloblist)
{
  if (bloblist->empty()) {
    msg_Error()<<METHOD<<"("<<bloblist<<"):\n"
	       <<"   Blob list contains "<<bloblist->size()<<" entries.\n"
	       <<"   Continue and hope for the best.\n";
    return Return_Value::Error;
  }
  switch (int(m_singlets(bloblist))) {
  case int(Return_Value::Success)   : break;
  case int(Return_Value::New_Event) : return Return_Value::New_Event;
  case int(Return_Value::Nothing)   : return Return_Value::Nothing;
  case int(Return_Value::Error)     : return Return_Value::Error;
  default :
    msg_Error()<<"ERROR in "<<METHOD<<":\n"
	       <<"   ExtractSinglets yields unknown return value.\n"
	       <<"   Return 'Retry_Event' and hope for the best.\n";
    return Return_Value::Retry_Event;
  }
  return ( m_on ?
	   (*p_reconnectionhandler)(bloblist) :
	   Return_Value::Nothing );
}

void NP_Colour_Reconnections::CleanUp(const size_t & mode) {}

void NP_Colour_Reconnections::Finish(const std::string &) {}


#include "SHERPA/Single_Events/Colour_Reconnections.H"
#include "ATOOLS/Org/Message.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Colour_Reconnections::
Colour_Reconnections(Colour_Reconnection_Handler * reconnections) :
  p_reconnectionhandler(reconnections)
{
  m_name = ( string("Colour_Reconnections: ") +
	     p_reconnectionhandler->Name() );
  m_type = eph::Hadronization;
  m_on = Settings::GetMainSettings()["FRAGMENTATION"].Get<bool>();
}

Colour_Reconnections::~Colour_Reconnections() {}

Return_Value::code Colour_Reconnections::Treat(ATOOLS::Blob_List* bloblist)
{
  if (!m_on && !p_reconnectionhandler->On()) return Return_Value::Nothing;
  if (bloblist->empty()) {
    msg_Error()<<"Colour_Reconnections::Treat("<<bloblist<<"): "<<endl
	       <<"   Blob list contains "<<bloblist->size()<<" entries."<<endl
	       <<"   Continue and hope for the best."<<endl;
    return Return_Value::Error;
  }
  int res = m_singlets(bloblist);
  switch (int(res)) {
  case int(Return_Value::Success)   : break;
  case int(Return_Value::New_Event) : return Return_Value::New_Event;
  case int(Return_Value::Nothing)   : return Return_Value::Nothing;
  case int(Return_Value::Error)     : return Return_Value::Error;
  default :
    msg_Error()<<"ERROR in "<<METHOD<<":"<<std::endl
	       <<"   ExtractSinglets yields unknown return value."<<std::endl
	       <<"   Return 'Retry_Event' and hope for the best."<<std::endl;
    return Return_Value::Retry_Event;
  }
  Return_Value::code ret = (*p_reconnectionhandler)(bloblist);
  if (ret!=Return_Value::Success && ret!=Return_Value::Nothing)
    THROW(fatal_error,"unexpected (and undefined) result.");
  return ret;
}

void Colour_Reconnections::CleanUp(const size_t & mode) {}

void Colour_Reconnections::Finish(const std::string &) {}


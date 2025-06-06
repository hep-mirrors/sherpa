#include "SHERPA/Single_Events/Nuclear_Formation.H"
#include "ATOOLS/Org/Message.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Nuclear_Formation::Nuclear_Formation(Colour_Reconnection_Handler * reconnections,
			     Fragmentation_Base * fragmentation) :
  m_on(fragmentation->Name()!="None"),
  p_reconnectionhandler(reconnections),
  p_fragmentationhandler(fragmentation)
{
  m_name = std::string("Nuclear_Formation: ")+p_fragmentationhandler->Name();
  m_type = eph::Nuclear_Formation;
}

Nuclear_Formation::~Nuclear_Formation() {}

Return_Value::code Nuclear_Formation::Treat(ATOOLS::Blob_List* bloblist)
{
  if (bloblist->empty()) {
    msg_Error()<<"Nuclear_Formation::Treat("<<bloblist<<"): "<<endl
	       <<"   Blob list contains "<<bloblist->size()<<" entries."<<endl
	       <<"   Continue and hope for the best."<<endl;
    return Return_Value::Error;
  }
  if (!m_on) return Return_Value::Nothing;
  switch (int(m_singlets(bloblist))) {
  case int(Return_Value::Success) : break;
  case int(Return_Value::New_Event) : return Return_Value::New_Event;
  case int(Return_Value::Nothing) : return Return_Value::Nothing;
  case int(Return_Value::Error)   : return Return_Value::Error;
  default :
    msg_Error()<<"ERROR in "<<METHOD<<":"<<std::endl
	       <<"   ExtractSinglets yields unknown return value."<<std::endl
	       <<"   Return 'Retry_Event' and hope for the best."<<std::endl;
    return Return_Value::Retry_Event;
  }
  Return_Value::code ret = (*p_reconnectionhandler)(bloblist);
  if (ret!=Return_Value::Success && ret!=Return_Value::Nothing)
    THROW(fatal_error,"unexpected (and undefined) result.")
  return p_fragmentationhandler->Hadronize(bloblist);
}

void Nuclear_Formation::CleanUp(const size_t & mode) {}

void Nuclear_Formation::Finish(const std::string &) {}


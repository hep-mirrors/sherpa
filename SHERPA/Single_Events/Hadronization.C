#include "SHERPA/Single_Events/Hadronization.H"
#include "ATOOLS/Org/Message.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Hadronization::Hadronization(Colour_Reconnection_Handler * reconnections,
			     Fragmentation_Handler * fragmentation) :
  m_on(fragmentation->Mode()>0),
  p_reconnectionhandler(reconnections),
  p_fragmentationhandler(fragmentation)
{
  m_name = std::string("Hadronization:")+p_fragmentationhandler->FragmentationModel();
  m_type = eph::Hadronization;
}

Hadronization::~Hadronization() {}

Return_Value::code Hadronization::Treat(ATOOLS::Blob_List* bloblist)
{
  if (bloblist->empty()) {
    msg_Error()<<"Hadronization::Treat("<<bloblist<<"): "<<endl
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
  if (ret!=Return_Value::Success && ret!=Return_Value::Nothing) exit(1);
  return (*p_fragmentationhandler)(bloblist);
}

void Hadronization::CleanUp(const size_t & mode) {}

void Hadronization::Finish(const std::string &) {}


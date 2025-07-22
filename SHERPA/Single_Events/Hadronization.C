#include "SHERPA/Single_Events/Hadronization.H"
#include "ATOOLS/Org/Message.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Hadronization::Hadronization(Fragmentation_Base * fragmentation) :
  m_on(fragmentation->Name()!="None"),
  p_fragmentationhandler(fragmentation)
{
  m_name = std::string("Hadronization: ")+p_fragmentationhandler->Name();
  m_type = eph::Hadronization;
}

Hadronization::~Hadronization() {}

Return_Value::code Hadronization::Treat(Blob_List* bloblist)
{
  if (bloblist->empty()) {
    msg_Error()<<"Hadronization::Treat("<<bloblist<<"): "<<endl
	       <<"   Blob list contains "<<bloblist->size()<<" entries."<<endl
	       <<"   Continue and hope for the best."<<endl;
    return Return_Value::Error;
  }
  return (m_on ?
	  p_fragmentationhandler->Hadronize(bloblist) :
	  Return_Value::Nothing );
}

void Hadronization::CleanUp(const size_t & mode) {}

void Hadronization::Finish(const std::string &) {}


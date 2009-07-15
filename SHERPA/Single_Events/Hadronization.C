#include "SHERPA/Single_Events/Hadronization.H"
#include "ATOOLS/Org/Message.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

#ifdef PROFILE__all
#define PROFILE__Hadronization
#endif
#ifdef PROFILE__Hadronization
#include "prof.hh"
#else 
#define PROFILE_HERE {}
#define PROFILE_LOCAL(LOCALNAME) {}
#endif


Hadronization::Hadronization(Fragmentation_Handler * fragmentation) :
  p_fragmentationhandler(fragmentation)
{
  m_name = std::string("Hadronization:")+p_fragmentationhandler->FragmentationModel();
  m_type = eph::Hadronization;
}

Hadronization::~Hadronization() {}

Return_Value::code Hadronization::Treat(ATOOLS::Blob_List *bloblist,double &weight) 
{
  PROFILE_LOCAL("Hadronization::Treat");
  if (bloblist->empty()) {
    msg_Error()<<"Hadronization::Treat("<<bloblist<<","<<weight<<"): "<<endl
	       <<"   Blob list contains "<<bloblist->size()<<" entries."<<endl
	       <<"   Continue and hope for the best."<<endl;
    return Return_Value::Error;
  }
  return p_fragmentationhandler->PerformFragmentation(bloblist);
}


void Hadronization::CleanUp() {}

void Hadronization::Finish(const std::string &) {}


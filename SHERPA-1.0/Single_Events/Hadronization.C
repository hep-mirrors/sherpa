#include "Hadronization.H"

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


Hadronization::Hadronization(Beam_Remnant_Handler * _beamremnant,Fragmentation_Handler * _fragmentation) :
  p_beamremnanthandler(_beamremnant), p_fragmentationhandler(_fragmentation)
{
  m_name = std::string("Hadronization:")+
    p_fragmentationhandler->FragmentationModel();
  m_type = eph::Hadronization;
}

Hadronization::~Hadronization() {}

bool Hadronization::Treat(ATOOLS::Blob_List *bloblist,double &weight) 
{
  PROFILE_LOCAL("Hadronization::Treat");
  if (bloblist->empty()) {
    msg.Error()<<"Hadronization::Treat("<<bloblist<<","<<weight<<"): "<<endl
	       <<"   Blob list contains "<<bloblist->size()<<" entries."<<endl
	       <<"   Continue and hope for the best."<<endl;
    return false;
  }
  p_beamremnanthandler->FillBeamBlobs(bloblist);
  p_beamremnanthandler->FillBunchBlobs(bloblist);
  p_fragmentationhandler->PerformFragmentation(bloblist);
  return false;
}


void Hadronization::CleanUp() {}

void Hadronization::Finish(const std::string &) {}


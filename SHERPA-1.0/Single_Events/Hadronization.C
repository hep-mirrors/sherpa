#include "Hadronization.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

#ifdef PROFILE__Hadronization
#include "prof.hh"
#else 
#define PROFILE_HERE {}
#define PROFILE_LOCAL(LOCALNAME) {}
#endif


Hadronization::Hadronization(Beam_Remnant_Handler * _beamremnant,Fragmentation_Handler * _fragmentation) :
  p_beamremnanthandler(_beamremnant), p_fragmentationhandler(_fragmentation)
{
  m_name = string("Hadronization : ")+p_fragmentationhandler->FragmentationModel();
  m_type = string("Hadronization");
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
  bool showered=true;
  if (p_beamremnanthandler->BeamParticle(0)->Type()==Remnant_Base::Hadron ||
      p_beamremnanthandler->BeamParticle(1)->Type()==Remnant_Base::Hadron) {
    showered=false;
    for (ATOOLS::Blob_List::iterator blit=bloblist->begin();blit!=bloblist->end();++blit) {
      if ((*blit)->Type()==btp::IS_Shower) showered=true;
    }
  }
  if (showered) {
    p_beamremnanthandler->FillBeamBlobs(bloblist);
    p_beamremnanthandler->FillBunchBlobs(bloblist);
    p_fragmentationhandler->PerformFragmentation(bloblist);
  }
  return false;
}


void Hadronization::CleanUp() {}

void Hadronization::Finish(const std::string &) {}


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
  m_name      = string("Hadronization : ")+p_fragmentationhandler->FragmentationModel();
  m_type      = string("Hadronization");
}

Hadronization::~Hadronization() {}

bool Hadronization::Treat(ATOOLS::Blob_List * _bloblist, double &) 
{
  PROFILE_LOCAL("Hadronization::Treat");
  if (_bloblist->empty()) {
    msg.Error()<<"Potential error in Hadronization::Treat."<<endl
	       <<"   Incoming blob list contains "<<_bloblist->size()<<" entries."<<endl
	       <<"   Continue and hope for the best."<<endl;
    return false;
  }
  p_beamremnanthandler->FillBeamBlobs(_bloblist);
  p_beamremnanthandler->FillBunchBlobs(_bloblist);

  // Here we should boost -> It's only here we have the full beam information.

  p_fragmentationhandler->PerformFragmentation(_bloblist);

  return false;
}


void Hadronization::CleanUp() {}

void Hadronization::Finish(const std::string &) {}


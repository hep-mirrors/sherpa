#include "SHERPA/Single_Events/Beam_Remnants.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

#ifdef PROFILE__all
#define PROFILE__Beam_Remnants
#endif
#ifdef PROFILE__Beam_Remnants
#include "prof.hh"
#else 
#define PROFILE_HERE {}
#define PROFILE_LOCAL(LOCALNAME) {}
#endif


Beam_Remnants::Beam_Remnants(Beam_Remnant_Handler * _beamremnant) :
  p_beamremnanthandler(_beamremnant)
{
  m_name = std::string("Beam_Remnants");
  m_type = eph::Hadronization;
}

Beam_Remnants::~Beam_Remnants() {}

Return_Value::code Beam_Remnants::Treat(ATOOLS::Blob_List *bloblist,double &weight) 
{
  PROFILE_LOCAL("Beam_Remnants::Treat");
  if (bloblist->empty()) {
    msg_Error()<<"Beam_Remnants::Treat("<<bloblist<<","<<weight<<"): "<<endl
	       <<"   Blob list contains "<<bloblist->size()<<" entries."<<endl
	       <<"   Continue and hope for the best."<<endl;
    return Return_Value::Error;
  }
  Blob *signal(bloblist->FindFirst(btp::Signal_Process));
  if (!signal || signal->Has(blob_status::needs_signal)) return Return_Value::Nothing;
  return p_beamremnanthandler->FillBeamAndBunchBlobs(bloblist);
}


void Beam_Remnants::CleanUp() 
{
  p_beamremnanthandler->CleanUp();
}

void Beam_Remnants::Finish(const std::string &) {}


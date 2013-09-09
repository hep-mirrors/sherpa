#include "SHERPA/Single_Events/Beam_Remnants.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Beam_Remnants::Beam_Remnants(Beam_Remnant_Handler * _beamremnant) :
  p_beamremnanthandler(_beamremnant)
{
  m_name = std::string("Beam_Remnants");
  m_type = eph::Hadronization;
}

Beam_Remnants::~Beam_Remnants() {}

Return_Value::code Beam_Remnants::Treat(ATOOLS::Blob_List *bloblist,double &weight) 
{
  if (bloblist->empty()) {
    msg_Error()<<"Beam_Remnants::Treat("<<bloblist<<","<<weight<<"): "<<endl
	       <<"   Blob list contains "<<bloblist->size()<<" entries."<<endl
	       <<"   Continue and hope for the best."<<endl;
    return Return_Value::Error;
  }
  Blob *signal(bloblist->FindFirst(btp::Signal_Process));
  if (signal && signal->NInP()<2) return Return_Value::Nothing;
  if (!signal || signal->Has(blob_status::needs_signal)) {
    Blob * hard  = bloblist->FindFirst(btp::Hard_Collision);
    Blob * qelas = bloblist->FindFirst(btp::QElastic_Collision);
    if (!hard && !qelas) {
      return Return_Value::Nothing;
    }
  }
  Blob *beam(bloblist->FindFirst(btp::Beam));
  if (beam && !beam->Has(blob_status::needs_beams)) 
    return Return_Value::Nothing;
  return p_beamremnanthandler->FillBeamAndBunchBlobs(bloblist);
}


void Beam_Remnants::CleanUp(const size_t & mode) 
{
  p_beamremnanthandler->CleanUp(mode);
}

void Beam_Remnants::Finish(const std::string &) {}


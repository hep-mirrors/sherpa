#include "Hadronization.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Hadronization::Hadronization(Beam_Remnant_Handler * _beamremnant,Fragmentation_Handler * _fragmentation) :
  p_beamremnanthandler(_beamremnant), p_fragmentationhandler(_fragmentation)
{
  m_name      = string("Hadronization : ")+p_fragmentationhandler->FragmentationModel();
  m_type      = string("Hadronization");
}

Hadronization::~Hadronization() {
  // are deleted in the Initialization_Handler
//   if (p_beamremnanthandler)   { delete p_beamremnanthandler;   p_beamremnanthandler   = NULL; }
//   if (p_fragmentationhandler) { delete p_fragmentationhandler; p_fragmentationhandler = NULL; }
// 

}

bool Hadronization::Treat(ATOOLS::Blob_List * _bloblist, double &) {
  if (_bloblist->empty()) {
    msg.Error()<<"Potential error in Hadronization::Treat."<<endl
	       <<"   Incoming blob list contains "<<_bloblist->size()<<" entries."<<endl
	       <<"   Continue and hope for the best."<<endl;
    return 0;
  }
  p_beamremnanthandler->FillBeamBlobs(_bloblist);
  p_beamremnanthandler->FillBunchBlobs(_bloblist);

  // Here we should boost -> It's only here we have the full beam information.

  p_fragmentationhandler->PerformFragmentation(_bloblist);
}


void Hadronization::CleanUp() {}


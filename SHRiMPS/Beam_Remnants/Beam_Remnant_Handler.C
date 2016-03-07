#include "SHRiMPS/Beam_Remnants/Beam_Remnant_Handler.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace SHRIMPS;
using namespace ATOOLS;
using namespace std;

Beam_Remnant_Handler::
Beam_Remnant_Handler(BEAM::Beam_Spectra_Handler * beamspectra,
		     vector<Continued_PDF> & pdfs) :
  p_softblob(NULL)
{
  for (int beam=0;beam<2;beam++) {
    m_hadrons.push_back(new Hadron_Dissociation(beamspectra->GetBeam(beam),
						&pdfs[beam]));
  }
}

Beam_Remnant_Handler::~Beam_Remnant_Handler() {
  for (int beam=0;beam<2;beam++) {
    delete m_hadrons[beam];
  }
  m_hadrons.clear();
}

bool Beam_Remnant_Handler::InitialiseCollision() {
  Reset();
}

void Beam_Remnant_Handler::Reset() {
  for (int beam=0;beam<2;beam++) m_hadrons[beam]->Reset();
}

Return_Value::code Beam_Remnant_Handler::FillBeamBlobs(Blob_List * blobs) {
  AddBeamBlobs(blobs);
  return Return_Value::Success;
}

void Beam_Remnant_Handler::AddBeamBlobs(Blob_List * blobs) {
  for (size_t beam=0;beam<2;beam++) {
    m_hadrons[beam]->FillBeamBlob(blobs);
    blobs->push_front(m_hadrons[beam]->GetBeamBlob());
  }
}


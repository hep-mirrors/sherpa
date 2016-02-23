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
  p_blob(NULL), p_pdfs(&pdfs)
{
  m_outmoms.push_back(Vec4D(0,0,0,0)); 
  m_outmoms.push_back(Vec4D(0,0,0,0)); 
  for (int beam=0;beam<2;beam++) {
    m_beams.push_back(beamspectra->GetBeam(beam));
    double E(beamspectra->GetBeam(beam)->OutMomentum()[0]);
    int dir(beamspectra->GetBeam(beam)->OutMomentum()[3]>0?1:-1);
    m_beamvecs.push_back(Vec4D(E,0,0,dir*E));
  }
}

bool Beam_Remnant_Handler::InitialiseCollision() {
  Reset();
  p_blob = new Blob();
  p_blob->SetType(btp::Soft_Collision);
  p_blob->SetTypeSpec("Four_Momentum_Compensation");
  p_blob->SetId();
  p_blob->SetStatus(blob_status::needs_hadronization |
		    blob_status::needs_beams);
}

void Beam_Remnant_Handler::Reset() {
  m_outmoms[0] = m_outmoms[1] = Vec4D(0,0,0,0);
}

Return_Value::code Beam_Remnant_Handler::FillBeamBlobs(Blob_List * blobs) {
  AddBeamBlobs(blobs);
  UnsetNeedsBeamsStatus(blobs);
  return Return_Value::Success;
}

void Beam_Remnant_Handler::AddBeamBlobs(Blob_List * blobs) {
  for (size_t beam=0;beam<2;beam++) {
    //m_hadrons[beam]->FillBeamBlob();
    //blobs->push_front(m_hadrons[beam]->GetBeamBlob());
  }
}

void Beam_Remnant_Handler::UnsetNeedsBeamsStatus(Blob_List * blobs) {
  Blob * blob;
  for (Blob_List::iterator biter=blobs->begin();biter!=blobs->end();biter++) {
    blob = (*biter);
    if (blob->Has(blob_status::needs_beams) &&
	(blob->Type()==btp::Beam || blob->Type()==btp::Shower)) {
      blob->UnsetStatus(blob_status::needs_beams);
    }
  }
}

#include "SHRiMPS/Beam_Remnants/Remnant_Handler.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "PDF/Main/ISR_Handler.H"
#include "REMNANTS/Main/Remnant_Base.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace SHRIMPS;
using namespace ATOOLS;
using namespace std;

Remnant_Handler::Remnant_Handler(PDF::ISR_Handler *const isr) :
  p_colourgenerator(NULL)
{
  for (int beam=0;beam<2;beam++) {
    m_hadrons.
      push_back(new Hadron_Dissociation(beam,
					isr->GetRemnant(beam)->InMomentum(),
					isr->Flav(beam),
					new Continued_PDF(isr->PDF(beam),
							  isr->Flav(beam))));
  }
}

Remnant_Handler::~Remnant_Handler() {
  for (int beam=0;beam<2;beam++) delete m_hadrons[beam];
  m_hadrons.clear();
}

void Remnant_Handler::Reset() {
  for (int beam=0;beam<2;beam++) m_hadrons[beam]->Reset();
}

Return_Value::code Remnant_Handler::
FillBeamBlobs(Blob_List * blobs,const double & B) {
  msg_Out()<<METHOD<<"(B = "<<B<<")\n";
  InitialiseCollision(blobs);
  for (size_t beam=0;beam<2;beam++) {
    if (!m_hadrons[beam]->FillBeamBlob(blobs, B)) {
      return Return_Value::New_Event;
    }
  }
  for (Blob_List::iterator bit=blobs->begin();bit!=blobs->end();bit++)
    (*bit)->UnsetStatus(blob_status::needs_beams);
  msg_Out()<<(*(m_hadrons[0]->GetSoftBlob()))<<"\n";
  return Return_Value::Success;
}
 
void Remnant_Handler::InitialiseCollision(Blob_List * blobs) {
  Blob * softblob = blobs->FindFirst(btp::Soft_Collision);
  if (!softblob || softblob->NInP()>0 || softblob->NOutP()>0) {
    softblob = new Blob();
    softblob->SetType(btp::Soft_Collision);
    softblob->SetId();
    blobs->push_front(softblob);
  }
  softblob->SetTypeSpec("Four_Momentum_Compensation");
  softblob->UnsetStatus(blob_status::needs_minBias);
  softblob->SetStatus(blob_status::needs_hadronization);
  for (int beam=0;beam<2;beam++) m_hadrons[beam]->SetSoftBlob(softblob);
}

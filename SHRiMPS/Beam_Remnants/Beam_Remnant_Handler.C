#include "SHRiMPS/Beam_Remnants/Beam_Remnant_Handler.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "ATOOLS/Phys/Momenta_Stretcher.H"
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

void Beam_Remnant_Handler::SetEikonal(Omega_ik * eikonal) {
  m_hadrons[0]->SetFormFactor(eikonal->FF1());
  m_hadrons[1]->SetFormFactor(eikonal->FF2());
  m_hadrons[0]->SetQTMap(&m_qtmap);
  m_hadrons[1]->SetQTMap(&m_qtmap);
}

void Beam_Remnant_Handler::Reset() {
  for (int beam=0;beam<2;beam++) m_hadrons[beam]->Reset();
  m_qtmap.clear();
}

void Beam_Remnant_Handler::
SetBeamBlob(ATOOLS::Blob *const beamblob,const int & beam) {
  m_hadrons[beam]->SetBeamBlob(beamblob);
} 

Return_Value::code Beam_Remnant_Handler::FillBeamBlobs(Blob_List * blobs) {
  AddBeamBlobs(blobs);
  AddTransverseMomentaToSpectators(blobs);
  for (Blob_List::iterator biter=blobs->begin();biter!=blobs->end();biter++) {
    (*biter)->UnsetStatus(blob_status::needs_beams);
    //msg_Out()<<METHOD<<": "<<(*biter)->Type()<<" "
    //	     <<(*biter)->CheckMomentumConservation()<<"\n";
  }
  return Return_Value::Success;
}

void Beam_Remnant_Handler::AddBeamBlobs(Blob_List * blobs) {
  for (size_t beam=0;beam<2;beam++) {
    m_hadrons[beam]->FillBeamBlob(blobs);
    blobs->push_front(m_hadrons[beam]->GetBeamBlob());
  }
}

void Beam_Remnant_Handler::AddTransverseMomentaToSpectators(Blob_List * blobs) {
  return;
  msg_Out()<<(*blobs)<<"\n";
  exit(0);
  Blob * softblob(blobs->FindFirst(btp::Soft_Collision));
  softblob->BoostInCMS();
  msg_Out()<<METHOD<<" for\n"<<(*softblob)<<"\n";
  std::vector<Particle *> parts;
  std::vector<Vec4D>      moms;
  std::vector<double>     masses;
  for (size_t i=0;i<softblob->NOutP();i++) {
    Particle * part(softblob->OutParticle(i));
    parts.push_back(part);
    moms.push_back(part->Momentum());
    masses.push_back(part->Flav().HadMass());
  }
  Momenta_Stretcher stretcher;
  stretcher.ZeroThem(moms.size(),moms);
  msg_Out()<<METHOD<<" for\n"<<(*softblob)<<"\n";
}


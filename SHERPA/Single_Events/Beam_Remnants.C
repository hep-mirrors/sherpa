#include "SHERPA/Single_Events/Beam_Remnants.H"
#include "ATOOLS/Org/Message.H"

using namespace SHERPA;
using namespace ATOOLS;

Beam_Remnants::Beam_Remnants(Beam_Remnant_Handler* _beamremnant)
    : m_ana(false), p_beamremnanthandler(_beamremnant)
{
  m_name = "Beam_Remnants:" + (p_beamremnanthandler->Fill() == 1
                                       ? p_beamremnanthandler->Name()
                                       : std::string("None"));
  m_type = eph::Hadronization;
  if (m_ana) InitHistos();
}

Beam_Remnants::~Beam_Remnants() {
  if (m_ana && !m_histos.empty()) {
    for (std::map<std::string, Histogram*>::iterator hit = m_histos.begin();
         hit != m_histos.end(); hit++) {
      hit->second->Finalize();
      hit->second->Output(std::string("ImpactParameter_Analysis/") +
                          hit->first + std::string(".dat"));
      delete hit->second;
    }
  }
}

Return_Value::code Beam_Remnants::Treat(Blob_List* bloblist)
{
  switch (EstablishNeed(bloblist)) {
  case 10: return DealWithRescattering(bloblist);
  case 3:  return DealWithShowerFromBeams(bloblist);
  case 2:  return StandardTreatment(bloblist,false);
  case 1:  return StandardTreatment(bloblist,true);
  default:
    break;
  }
  return Return_Value::Nothing;
}

Return_Value::code Beam_Remnants::StandardTreatment(Blob_List*  bloblist,
                                                    const bool& onlyBunch)
{
  Return_Value::code rv =
          p_beamremnanthandler->FillBeamAndBunchBlobs(bloblist, onlyBunch);
  if (m_ana) Analyse(bloblist);
  return rv;
}

Return_Value::code Beam_Remnants::DealWithShowerFromBeams(Blob_List* bloblist)
{
  msg_Out()<<"**** "<<METHOD<<"\n";
  Return_Value::code rv = p_beamremnanthandler->FillBunchBlobsFromShower(bloblist);
  if (m_ana) Analyse(bloblist);
  //msg_Out()<<"-----------------------------------------------------\n"
  //	   <<METHOD<<"\n"<<(*bloblist)<<"\n"
  //	   <<"-----------------------------------------------------\n";
  return rv;
}

Return_Value::code Beam_Remnants::DealWithRescattering(Blob_List* bloblist)
{
  Blob * shower = bloblist->FindLast(btp::Shower);
  if (shower) {
    if (p_beamremnanthandler->NeedsToDealWithRescattering()) {
      return p_beamremnanthandler->FillRescatterBeamBlobs(bloblist);
    }
  }
  //////////////////////////////////////////////////////////////////////////////
  // No blobs to be created or filled for bunch-rescattering - there is no
  // hard process in need of it
  // TODO: Extend this for elastic, diffractive, soft etc. scattering
  //////////////////////////////////////////////////////////////////////////////
  for (Blob_List::iterator bit=bloblist->begin();
       bit!=bloblist->end();bit++) {
    if ((*bit)->Type() == btp::Bunch)
      (*bit)->UnsetStatus(blob_status::needs_beamRescatter);
  }
  return Return_Value::Nothing;
}


int Beam_Remnants::EstablishNeed(Blob_List * bloblist) {
  if (bloblist->empty()) THROW(fatal_error,"Blob list is empty.");
  //////////////////////////////////////////////////////////////////////////////
  // Nothing to be done for 1->N processes or if we already have filled bunch
  // blobs without rescattering
  //////////////////////////////////////////////////////////////////////////////
  Blob *signal(bloblist->FindFirst(btp::Signal_Process));
  if (signal && signal->NInP()<2)                             return 0;
  //////////////////////////////////////////////////////////////////////////////
  // Some open bunch rescattering business to be taken care off.
  //////////////////////////////////////////////////////////////////////////////
  Blob * bunch(bloblist->FindFirst(btp::Bunch));
  if (bunch && bunch->Has(blob_status::needs_beamRescatter)) return 10;
  //////////////////////////////////////////////////////////////////////////////
  // Beam blobs already exists, and they do not need beams any more -
  // so nothing to be done here.
  //////////////////////////////////////////////////////////////////////////////
  Blob * beam(bloblist->FindFirst(btp::Beam));
  if (beam && !beam->Has(blob_status::needs_beams))           return 0;
  // Funny or absent signal blob - suggest we have a case of soft or semi-soft
  // scattering - check if we need beam blobs at all.
  if (!signal || signal->Has(blob_status::needs_signal)) {
    Blob * hard  = bloblist->FindFirst(btp::Hard_Collision);
    Blob * qelas = bloblist->FindFirst(btp::Elastic_Collision);
    if (!qelas) qelas = bloblist->FindFirst(btp::Soft_Diffractive_Collision);
    if (!qelas) qelas = bloblist->FindFirst(btp::Quasi_Elastic_Collision);
    if (!hard && !qelas) return 0;
    if (qelas)           return 1;
    if (CountBunchCandidates(bloblist->FindFirst(btp::Shower))==3) return 3;
  }
  else {
    btp::code signal_type = signal->Type();
    if (signal_type == btp::Elastic_Collision ||
        signal_type == btp::Soft_Diffractive_Collision ||
        signal_type == btp::Quasi_Elastic_Collision)
      return 1;
  }
  // Standard case, fill beam blobs in full beam remnant treatment
  return 2;
}

size_t Beam_Remnants::CountBunchCandidates(Blob * blob) {
  // Counting the number of bunch particles that are entering a shower blob.
  // This should only happen for diffractive events where we have hadron
  // splitting into their constituents and showering.
  size_t bunch_candidates = 0;
  if (blob) {
    for (size_t i=0;i<blob->NInP();i++) {
      Particle * part = blob->InParticle(i);
      if (part->Info()=='I' && part->Beam()>-1) {
	REMNANTS::Remnant_Base * remnant =
	  p_beamremnanthandler->GetRemnants()->GetRemnant(part->Beam());
	if (part->Flav()==remnant->InFlav() &&
	    part->Momentum()==remnant->InMomentum()) bunch_candidates+=part->Beam()+1;
      }
    }
  }
  return bunch_candidates;
}

void Beam_Remnants::CleanUp(const size_t& mode)
{
  p_beamremnanthandler->CleanUp(mode);
}

void Beam_Remnants::Finish(const std::string&) {}

void Beam_Remnants::InitHistos() {
  m_histos["phi_1"] = new Histogram(0,0.,360.,36);
  m_histos["phi_2"] = new Histogram(0,0.,360.,36);
  m_histos["b_1"]   = new Histogram(0,0.,10.,100);
  m_histos["b_2"]   = new Histogram(0,0.,10.,100);
  m_histos["B"]     = new Histogram(0,0.,20.,200);
}

void Beam_Remnants::Analyse(Blob_List * blobs) {
  Vec4D B1 = p_beamremnanthandler->GetRemnants()->GetRemnant(0)->Position();
  Vec4D B2 = -B1;
  m_histos["B"]->Insert(2.*B1.PPerp());
  for (auto blob : *blobs) {
    if (blob->Type()==btp::Hard_Collision) {
      Vec4D pos   = blob->Position();
      double b1   = (pos-B1).PPerp(), b2 = (pos-B2).PPerp();
      double phi1 = acos((pos[1]-B1[1])/b1) * 360./(2.*M_PI);
      double phi2 = acos((pos[1]-B2[1])/b2) * 360./(2.*M_PI);
      m_histos["b_1"]->Insert(b1);
      m_histos["b_2"]->Insert(b2);
      m_histos["phi_1"]->Insert(phi1);
      m_histos["phi_2"]->Insert(phi2);
    }
  }
}


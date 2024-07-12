#include "SHERPA/Single_Events/Multiple_Interactions.H"

#include "ATOOLS/Org/My_Limits.H"
#include "REMNANTS/Main/Remnant_Base.H"
#include "BEAM/Main/Beam_Base.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "SHERPA/PerturbativePhysics/Matrix_Element_Handler.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "MODEL/Main/Running_AlphaS.H"

using namespace SHERPA;
using namespace ATOOLS;

Multiple_Interactions::Multiple_Interactions(MI_Handler_Map *mihandlers):
  p_mihandlers(mihandlers), m_result(Return_Value::Nothing)
{
  m_type   = eph::Perturbative;
  m_name   = std::string("Multiple_Interactions: ")+MakeNameSpec();
  if (!CheckMIHandlers()) THROW(fatal_error,"No beam remnant handler found.");
  Settings& s = Settings::GetMainSettings();
  m_hardveto  = s["MPI_PT_MAX"].SetDefault(1e12).Get<double>();
  m_ptmax_fac = s["MPI_PT_Max_Fac"].SetDefault(1.).Get<double>();
  m_bbr_veto  = s["BBR_VETO"].Get<bool>();
}

std::string Multiple_Interactions::MakeNameSpec() {
  std::string spec = std::string("");
  bool add = false;
  if (p_mihandlers->find(PDF::isr::hard_subprocess)!=p_mihandlers->end()) {
    spec += (*p_mihandlers)[PDF::isr::hard_subprocess]->Name();
    add   = true;
  }
  if (p_mihandlers->find(PDF::isr::bunch_rescatter)!=p_mihandlers->end()) {
    if (add) spec += std::string(" + ");
    spec          += ( (*p_mihandlers)[PDF::isr::bunch_rescatter]->Name()+
		       std::string(" (rescatter)") );
  }
  return spec;
}

bool Multiple_Interactions::CheckMIHandlers() {
  for (MI_Handler_Map::iterator mihit=p_mihandlers->begin();
       mihit!=p_mihandlers->end();mihit++) {
    if (mihit->second->Generator()!=MI_Handler::genID::none &&
	mihit->second->Generator()!=MI_Handler::genID::unknown) {
      for (size_t i=0;i<2;i++) {
	if (mihit->second->Remnants()->GetRemnant(i)==NULL) return false;
      }
    }
  }
  return true;
}

bool Multiple_Interactions::CheckForMinBias() {
  ////////////////////////////////////////////////////////////////////////////
  // Min Bias events are initialised by having an empty Soft_Collsion blob
  // which needs_minBias.
  ////////////////////////////////////////////////////////////////////////////
  if (p_bloblist->empty()) THROW(fatal_error, "Empty bloblist - this is odd.");
  if (p_bloblist->size()==1) {
    Blob * signal = (*p_bloblist)[0];
    if (signal->Has(ATOOLS::blob_status::needs_minBias) &&
	signal->Type()==ATOOLS::btp::Soft_Collision) return true;
  }
  return false;
}

bool Multiple_Interactions::CheckForMPIs() {
  ////////////////////////////////////////////////////////////////////////////
  // In MPI events we need a signal event that has already been filled.
  // We check for stray blobs that need still need showering later.
  ////////////////////////////////////////////////////////////////////////////
  if (p_bloblist->empty()) THROW(fatal_error,"Empty bloblist - this is odd.");
  if (!m_newevent[0]) return false;
  Blob * blob = p_bloblist->FindFirst(btp::Signal_Process);
  return (blob && !blob->Has(blob_status::needs_signal));
}

Return_Value::code Multiple_Interactions::InitMinBias() {
  ////////////////////////////////////////////////////////////////////////////
  // No new_event flag or no meaningful MI_Handler - this is somewhat odd.
  ////////////////////////////////////////////////////////////////////////////
  if (!m_newevent[0] ||
      p_mihandlers->find(PDF::isr::hard_subprocess)==p_mihandlers->end() ||
      !(*p_mihandlers)[PDF::isr::hard_subprocess]->IsMinBias() ||
      !(*p_mihandlers)[PDF::isr::hard_subprocess]->On()) {
    return Return_Value::Nothing;
  }
  ////////////////////////////////////////////////////////////////////////////
  // Setup the MI_Handler: select the right one, fix initial hadron energies
  // and switch to the right perturbative parameters (including x-rescaling
  // of the PDFs in the MPIs, which depend on the hadron energies)
  ////////////////////////////////////////////////////////////////////////////
  p_activeMI = (*p_mihandlers)[PDF::isr::hard_subprocess];
  FixMaxEnergies();
  SwitchPerturbativeInputsToMIs();
  ////////////////////////////////////////////////////////////////////////////
  // Generate a first hard scatter - if successful move particles from
  // scatter blob into signal blob, add required information to the latter,
  // delete the former, and inform the remnants about the impact parameter
  // of the collision
  ////////////////////////////////////////////////////////////////////////////
  Blob * signal = (*p_bloblist)[0];
  p_lastblob = p_activeMI->GenerateHardProcess(MI_Handler::typeID::minbias,signal);
  if (p_lastblob) {
    Blob * signal         = (*p_bloblist)[0];
    Particle_Vector * ins = p_lastblob->InParticles();
    while (!ins->empty()) {
      signal->AddToInParticles(p_lastblob->RemoveInParticle(ins->back()));
    }
    Particle_Vector * outs = p_lastblob->OutParticles();
    while (!outs->empty()) {
      signal->AddToOutParticles(p_lastblob->RemoveOutParticle(outs->back()));
    }
    signal->SetStatus(blob_status::code(p_lastblob->Status()));
    signal->SetType(p_lastblob->Type());
    signal->SetTypeSpec(p_lastblob->TypeSpec());
    signal->SetPosition(p_lastblob->Position());
    signal->AddData("WeightsMap",
                    new Blob_Data<Weights_Map>((*p_lastblob)["WeightsMap"]->Get<Weights_Map>()));
    signal->AddData("Renormalization_Scale",
		    new Blob_Data<double>((*p_lastblob)["Renormalization_Scale"]
					  ->Get<double>()));
    signal->AddData("Factorization_Scale",
		    new Blob_Data<double>((*p_lastblob)["Factorization_Scale"]
					  ->Get<double>()));
    signal->AddData("Resummation_Scale",
		    new Blob_Data<double>((*p_lastblob)["Resummation_Scale"]
					  ->Get<double>()));
    delete p_lastblob;
    p_activeMI->Remnants()->SetImpactParameter(p_activeMI->ImpactParameter());
    m_newevent[0] = false;
    return Return_Value::Success;
  }
  ////////////////////////////////////////////////////////////////////////////
  // No meaningful first scatter in MinBias event produced - ask for new event
  ////////////////////////////////////////////////////////////////////////////
  return Return_Value::New_Event;
}

bool Multiple_Interactions::InitMPIs() {
  ////////////////////////////////////////////////////////////////////////////
  // No new_event flag or no meaningful MI_Handler - this is somewhat odd.
  ////////////////////////////////////////////////////////////////////////////
  if (!m_newevent[0] ||
      p_mihandlers->find(PDF::isr::hard_subprocess)==p_mihandlers->end() ||
      (*p_mihandlers)[PDF::isr::hard_subprocess]->IsMinBias() ||
      !(*p_mihandlers)[PDF::isr::hard_subprocess]->On()) {
    m_newevent[0] = false;
    return false;
  }
  ////////////////////////////////////////////////////////////////////////////
  // Setup the MI_Handler: select the right one, fix initial hadron energies
  // and switch to the right perturbative parameters (including x-rescaling
  // of the PDFs in the MPIs, which depend on the hadron energies)
  ////////////////////////////////////////////////////////////////////////////
  p_activeMI = (*p_mihandlers)[PDF::isr::hard_subprocess];
  FixMaxEnergies();
  SwitchPerturbativeInputsToMIs();
  Blob * signal = p_bloblist->FindFirst(btp::Signal_Process);
  if (p_activeMI->GenerateHardProcess(MI_Handler::typeID::MPI,signal)) {
    p_activeMI->Remnants()->SetImpactParameter(p_activeMI->ImpactParameter());
    Blob * shower = signal->OutParticle(0)->DecayBlob();
    if (shower) shower->SetPosition(signal->Position());
    m_newevent[0] = false;
    return true;
  }
  return false;
}

bool Multiple_Interactions::CheckForRescatter() {
  ////////////////////////////////////////////////////////////////////////////
  // In rescatter events - similar to MinBias - we look for a soft collision
  // blob that still needs to be filled with MinBias methods.
  ////////////////////////////////////////////////////////////////////////////
  if (p_bloblist->empty()) THROW(fatal_error, "Empty bloblist - this is odd.");
  if (!m_newevent[1]) return false;
  Blob * blob = p_bloblist->FindLast(ATOOLS::btp::Soft_Collision);
  return (blob && blob->Has(ATOOLS::blob_status::needs_beamRescatter));
}

ATOOLS::Return_Value::code Multiple_Interactions::InitRescatter() {
  if (!m_newevent[1] ||
      p_mihandlers->find(PDF::isr::bunch_rescatter) == p_mihandlers->end() ||
      !(*p_mihandlers)[PDF::isr::bunch_rescatter]->IsMinBias() ||
      !(*p_mihandlers)[PDF::isr::bunch_rescatter]->On()) {
    return Return_Value::Nothing;
  }
  ////////////////////////////////////////////////////////////////////////////
  // Setup the MI_Handler: select the right one, fix initial hadron energies
  // and switch to the right perturbative parameters (including x-rescaling
  // of the PDFs in the MPIs, which depend on the hadron energies)
  ////////////////////////////////////////////////////////////////////////////
  p_activeMI = (*p_mihandlers)[PDF::isr::bunch_rescatter];
  FixMaxEnergies(true);
  SwitchPerturbativeInputsToMIs();
  ////////////////////////////////////////////////////////////////////////////
  // Generate a first hard scatter - if successful move particles from
  // scatter blob into soft collision blob, add required information to the
  // latter, delete the former.  If no scatter found, delete soft blob.
  ////////////////////////////////////////////////////////////////////////////
  Blob * soft = p_bloblist->FindLast(ATOOLS::btp::Soft_Collision);
  if (p_activeMI->GenerateHardProcess(MI_Handler::typeID::rescatter,soft)) {
    if (m_bbr_veto) return Return_Value::New_Event;
    m_newevent[0] = m_newevent[1] = false;
    return Return_Value::Success;
  }
  ////////////////////////////////////////////////////////////////////////////
  // No meaningful first scatter in MinBias event produced - ask for new event
  ////////////////////////////////////////////////////////////////////////////
  p_bloblist->Delete(soft);
  return Return_Value::Nothing;
}

void Multiple_Interactions::FixMaxEnergies(const bool & updateResidualE) {
  ////////////////////////////////////////////////////////////////////////////
  // The emax has to be set here (instead of e.g. the CleanUp()) to ensure
  // that the correct energy is taken in case of EPA-approximated beams.
  ////////////////////////////////////////////////////////////////////////////
  for (short unsigned int i = 0; i < 2; ++i) {
    m_emax[i] =
      ( (p_activeMI->Remnants()->Id()==PDF::isr::bunch_rescatter) ?
	(p_activeMI->Remnants()->GetRemnant(i)->GetBeam()->InMomentum()-
	 p_activeMI->Remnants()->GetRemnant(i)->GetBeam()->OutMomentum())[0] :
	p_activeMI->Remnants()->GetRemnant(i)->GetBeam()->OutMomentum()[0]);
    if (updateResidualE)
      p_activeMI->Remnants()->GetRemnant(i)->SetResidualEnergy();
  }
  p_activeMI->SetMaxEnergies(m_emax[0],m_emax[1]);
}

void Multiple_Interactions::SwitchPerturbativeInputsToMIs() {
  ////////////////////////////////////////////////////////////////////////////
  // Use the right alphaS and initialise the x-rescaling of the beams, i.e. by
  // multiplying them with a ratio of the energies.
  ////////////////////////////////////////////////////////////////////////////
  MODEL::as->SetActiveAs(PDF::isr::hard_subprocess);
  for (size_t i=0;i<2;i++) {
    double x_resc =
      ((p_activeMI->Remnants()->Id()!=PDF::isr::bunch_rescatter) ?
       m_emax[i]/
       p_activeMI->Remnants()->GetRemnant(i)->GetBeam()->OutMomentum()[0] :
       m_emax[i]/
       (p_activeMI->Remnants()->GetRemnant(i)->GetBeam()->InMomentum()[0]-
	p_activeMI->Remnants()->GetRemnant(i)->GetBeam()->OutMomentum()[0]));
    p_activeMI->ISRHandler()->SetRescaleFactor(x_resc,i);
  }
}

Return_Value::code Multiple_Interactions::Treat(Blob_List *bloblist) {
  ////////////////////////////////////////////////////////////////////////////
  // Check for - and if necessary initialise - Minimum Bias, bunch
  // rescattering, and MPIs.  Return nothing if there is no suitable
  // MI_Handler left.  The initialisation will either lead to a first blob
  // to be copied into the signal (MinBias) or the definition of suitable
  // starting conditions, impact parameters etc. (MPIs and Rescattering).
  ////////////////////////////////////////////////////////////////////////////
  p_bloblist = bloblist;
  if (CheckForMinBias())                 return InitMinBias();
  if (CheckForRescatter())               return InitRescatter();
  if (CheckForMPIs() && !InitMPIs())     return Return_Value::Nothing;
  if (!p_activeMI || p_activeMI->Done()) return Return_Value::Nothing;
  ////////////////////////////////////////////////////////////////////////////
  // Sanity checks on blob_list: four-momentum is conserved, no blob in there
  // that needs to parton shower, beams are viable.
  ////////////////////////////////////////////////////////////////////////////
  m_result = Return_Value::Nothing;
  if (!p_bloblist->FourMomentumConservation()) return Return_Value::Retry_Event;
  if (!CheckBlobList())                        return Return_Value::Nothing;
  if (!BeamsViable())                          return m_result;
  ////////////////////////////////////////////////////////////////////////////
  // Try to colour-connect the last interaction with the remnants, if necessary:
  // this is the case only,  if the shower_blob has the "needs_beams" tag still
  // attached.
  ////////////////////////////////////////////////////////////////////////////
  Blob * lastShower = p_bloblist->FindLast(btp::Shower);
  if (lastShower && lastShower->Has(blob_status::needs_beams)) {
    if (!p_activeMI->ConnectColours(p_bloblist->FindLast(btp::Shower)))
      return Return_Value::New_Event;
  }
  ////////////////////////////////////////////////////////////////////////////
  // Try to produce a hard scattering blob and test if this works.
  ////////////////////////////////////////////////////////////////////////////
  //<<Blob::Counter()<<"/"<<Particle::Counter()<<".\n";
  Blob * blob = new Blob();
  if (p_activeMI->GenerateHardProcess(MI_Handler::typeID::MPI,blob)) {
    //////////////////////////////////////////////////////////////////////////
    // Check that the partons can be extracted from remnant - mainly a
    // confirmation that the remnant has enough energy to accommodate
    // the extra parton.
    //////////////////////////////////////////////////////////////////////////
    if (!TestHardScatter(blob)) {
      delete blob;
      return Return_Value::Success;
      return Return_Value::Retry_Event;
    }
    if (p_activeMI->Id()==PDF::isr::bunch_rescatter)
      blob->AddStatus(blob_status::needs_beamRescatter);
    p_bloblist->push_back(blob);
    return Return_Value::Success;
  }
  //////////////////////////////////////////////////////////////////////////
  // If we have reached the end of MPI production with a meaningful event,
  // we can stop here.
  //////////////////////////////////////////////////////////////////////////
  else if (p_activeMI->Done()) {
    delete blob;
    return Return_Value::Nothing;
  }
  //////////////////////////////////////////////////////////////////////////
  // If it is a MinBias event where the event handler didn't manage to
  // produce a first scatter (i.e. the first blob still needs a signal).
  //////////////////////////////////////////////////////////////////////////
  return Return_Value::New_Event;
}

bool Multiple_Interactions::TestHardScatter(Blob * blob) {
  for (size_t i=0;i<(size_t)blob->NInP();++i) {
    if (!p_activeMI->Remnants()->GetRemnant(i)->
	TestExtract(blob->InParticle(i))) {
      return false;
    }
  }
  return true;
}

bool Multiple_Interactions::CheckBlobList()
{
  //////////////////////////////////////////////////////////////////////////
  // check if there is a blob that must shower first.
  //////////////////////////////////////////////////////////////////////////
  for (Blob_List::const_iterator bit=p_bloblist->begin();
       bit!=p_bloblist->end();++bit) {
    if (((*bit)->Type()==btp::Hard_Collision ||
         (*bit)->Type() == btp::Signal_Process) &&
        (*bit)->Has(blob_status::needs_showers)) {
      return false;
    }
  }
  return true;
}

bool Multiple_Interactions::BeamsViable() {
  //////////////////////////////////////////////////////////////////////////
  // Checking if the total energy in shower initiators exceeds the
  // energies of the incoming beams.  If yes, undo showering and
  // start with the signal process (if this happens after first shower)
  // or undo last MI interaction + shower.
  // This method uses some implicit knowledge -- it knows the sequence
  // of shower blobs is such that the last Hard Collision (i.e. last
  // MI interaction) gives rise to last shower blob, while the Signal
  // Process comes first.  It also knows that the shower initiators
  // in the In-state of the shower blobs are sorted ....
  //////////////////////////////////////////////////////////////////////////
  Blob_List isr=p_bloblist->Find(btp::Shower);
  for (Blob_List::iterator bit=isr.begin();bit!=isr.end();++bit) {
    if (m_ISblobs.find((*bit))!=m_ISblobs.end()) continue;
    if (!ExtractISInfo((*bit))) return false;
    m_ISblobs.insert((*bit));
  }
  return true;
}

bool Multiple_Interactions::ExtractISInfo(Blob * blob) {
  for (size_t i=0;i<blob->NInP();++i) {
    Particle *particle(blob->InParticle(i));
    if (particle->ProductionBlob()) continue;
    size_t beam = particle->Beam();
    if (!p_activeMI->Remnants()->GetRemnant(beam)->TestExtract(particle)) {
      if (!blob->IsConnectedTo(btp::Signal_Process)) {
        p_bloblist->DeleteConnected(blob);
        m_result = Return_Value::Retry_Phase;
      }
      else {
        m_result = Return_Value::Retry_Event;
      }
      return false;
    }
    m_emax[beam] -= particle->Momentum()[0];
  }
  return true;
}

void Multiple_Interactions::Reset() {
  m_newevent[0] = m_newevent[1] = true;
}

void Multiple_Interactions::Finish(const std::string &resultpath) {}

void Multiple_Interactions::CleanUp(const size_t& mode)
{
  for (MI_Handler_Map::iterator mihit=p_mihandlers->begin();
       mihit!=p_mihandlers->end();mihit++) {
    mihit->second->CleanUp();
  }
  m_vetoed      = false;
  m_newevent[0] = m_newevent[1] = true;
  p_activeMI    = NULL;
}

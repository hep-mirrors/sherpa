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
  p_mihandlers(mihandlers), m_result(Return_Value::Nothing),
  m_newevent(true)
{
  m_name   = std::string("Multiple_Interactions: ");
  bool add = false;
  if (p_mihandlers->find(PDF::isr::hard_subprocess)!=p_mihandlers->end()) {
    m_name += (*p_mihandlers)[PDF::isr::hard_subprocess]->Name();
    add     = true;
  }
  if (p_mihandlers->find(PDF::isr::bunch_rescatter)!=p_mihandlers->end()) {
    if (add)  m_name += std::string(" + "); 
    m_name += ( (*p_mihandlers)[PDF::isr::bunch_rescatter]->Name()+
		std::string(" (rescatter)") );
  }
  m_type = eph::Perturbative;
  for (MI_Handler_Map::iterator mihit=p_mihandlers->begin();
       mihit!=p_mihandlers->end();mihit++) {
    //msg_Out()<<METHOD<<"("<<int(mihit->first)<<"): "<<mihit->second->Type()<<".\n";
    if (mihit->second->Type()!=0) {
      for (size_t i=0;i<2;i++) {
	if (mihit->second->Remnants()->GetRemnant(i)==NULL) {
	  THROW(fatal_error,"No beam remnant handler found.");
	}
      }
    }
  }
  Settings& s = Settings::GetMainSettings();
  m_hardveto  = s["MPI_PT_MAX"].SetDefault(1e12).Get<double>();
  m_ptmax_fac = s["MPI_PT_Max_Fac"].SetDefault(1.).Get<double>();
  if (p_mihandlers->find(PDF::isr::hard_subprocess)!=p_mihandlers->end())
    p_activeMI = (*p_mihandlers)[PDF::isr::hard_subprocess];
  else
    THROW(fatal_error,"No suitable MI handler found.");
}

Multiple_Interactions::~Multiple_Interactions() { }

void Multiple_Interactions::Reset() {
  if (p_mihandlers->find(PDF::isr::hard_subprocess)!=p_mihandlers->end() &&
      (*p_mihandlers)[PDF::isr::hard_subprocess]->Type()!=MI_Handler::none)
    p_activeMI = (*p_mihandlers)[PDF::isr::hard_subprocess];
  else SwapToRescatter();
  ResetMIHandler();
}

void Multiple_Interactions::SwapToRescatter() {
  if (p_mihandlers->find(PDF::isr::bunch_rescatter)!=p_mihandlers->end() &&
      (*p_mihandlers)[PDF::isr::bunch_rescatter]->Type()!=MI_Handler::none)
    p_activeMI = (*p_mihandlers)[PDF::isr::bunch_rescatter];
  else p_activeMI = NULL;
}

void Multiple_Interactions::ResetMIHandler() {
  if (p_activeMI!=NULL) {
    p_lastblob = NULL;
    m_ISblobs.clear();
    if (p_activeMI->Id()==PDF::isr::bunch_rescatter)
      msg_Out()<<"################################################# "<<METHOD<<"\n";
    if (p_activeMI->Type()!=0) p_activeMI->Reset();
  }
  msg_Out()<<METHOD<<": "<<p_activeMI<<"\n";
}


Return_Value::code Multiple_Interactions::Treat(Blob_List *bloblist)
{
  msg_Out()<<"#################################################\n"
  	   <<"Entering "<<METHOD<<": ";
  if (p_activeMI)
    msg_Out()<<"Id = "<<p_activeMI->Id()<<", done = "<<p_activeMI->Done()<<", "
  	     <<"MI MinBias/Type = "<<p_activeMI->IsMinBias()<<"/"<<int(p_activeMI->Type())<<", "
  	     <<"done = "<<p_activeMI->Done()<<".\n";
  else
    msg_Out()<<"No MI Handler.\n";
  /////////////////////////////////////////////////////////////////////////////////
  // The "regular" multi-parton interactions are done - do we need to do some
  // beam rescattering (which will be some MinBias collisions?
  /////////////////////////////////////////////////////////////////////////////////
  if (p_activeMI && p_activeMI->Done() && p_activeMI->Id()==PDF::isr::hard_subprocess) SwapToRescatter();
  if (!p_activeMI || p_activeMI->Type()==MI_Handler::none) return Return_Value::Nothing; 
  m_result   = Return_Value::Nothing;  
  p_bloblist = bloblist;
  /////////////////////////////////////////////////////////////////////////////////
  // Try to colour-connect the last interaction with the remnants, if necessary:
  // this is the case if the shower_blob has the "needs_beams" tag still attached.
  /////////////////////////////////////////////////////////////////////////////////
  Blob * lastShower = p_bloblist->FindLast(btp::Shower);
  if (lastShower->Has(blob_status::needs_beams)) {
    p_activeMI->ConnectColours(p_bloblist->FindLast(btp::Shower));
  }
  m_isfirstMB = ( (p_bloblist->size()==1 && p_activeMI->IsMinBias()) ||
		  p_activeMI->IsFirstRescatter() );
  /////////////////////////////////////////////////////////////////////////////////
  // CheckBlobList makes sure a new interaction can be added.
  // If its the first then a completely new chain of 2->2 scatters 
  // must be initialised.  This is steered by a flag m_newevent, which is 
  // set to true in the CleanUp() method.
  /////////////////////////////////////////////////////////////////////////////////
  if (m_newevent) {
    msg_Out()<<"   * resetting energies and remnants.\n";
    ///////////////////////////////////////////////////////////////////////////////
    // The flag m_newevent is then set to false after the bloblist and the
    // energies have been checked, i.e. in the InitNewEvent function.
    // The emax has to be set here (instead of e.g. the CleanUp()) to ensure
    // that the correct energy is taken in case of EPA-approximated beams.
    ///////////////////////////////////////////////////////////////////////////////
    for (short unsigned int i = 0; i < 2; ++i) {
      m_emax[i] = ( (p_activeMI->Remnants()->Id()==PDF::isr::bunch_rescatter) ?
		    (p_activeMI->Remnants()->GetRemnant(i)->GetBeam()->InMomentum()-
		     p_activeMI->Remnants()->GetRemnant(i)->GetBeam()->OutMomentum())[0] :
		    p_activeMI->Remnants()->GetRemnant(i)->GetBeam()->OutMomentum()[0]); 
    }
  }
  if (!CheckBlobList() || !InitNewEvent() || !MIKinematics()) {
    return m_result;
  }
  /////////////////////////////////////////////////////////////////////////////////
  // Possibly switch to new PDF and alphaS.
  // TODO: will have to check that this happens.
  /////////////////////////////////////////////////////////////////////////////////
  SwitchPerturbativeInputsToMIs();
  p_lastblob = p_activeMI->GenerateHardProcess();
  if (p_lastblob) {
    ///////////////////////////////////////////////////////////////////////////////
    // This assumes that the scatters are ordered in transverse momentum.
    // Then maximal scale of subsequent scatters is given by the pT of the
    // previous ones.
    ///////////////////////////////////////////////////////////////////////////////
    m_ptmax = p_lastblob->OutParticle(0)->Momentum().PPerp();
    ///////////////////////////////////////////////////////////////////////////////
    // Check that the partons can be extracted from remnant - mainly a
    // confirmation that the remnant has enough energy to accommodate
    // the extra parton.
    ///////////////////////////////////////////////////////////////////////////////
    for (size_t i=0;i<(size_t)p_lastblob->NInP();++i) {
      if (!p_activeMI->Remnants()->GetRemnant(i)->TestExtract(p_lastblob->InParticle(i))) {
        delete p_lastblob;
        return Return_Value::Retry_Event;
      }
    }
    ///////////////////////////////////////////////////////////////////////////////
    // If it is a MinBias event, the first blob is a dummy soft collision blob.
    // We have to fill it with the content of the actual blob created by
    // the MI_Handler
    ///////////////////////////////////////////////////////////////////////////////
    if (m_isfirstMB) InitMinBiasEvent();
    else p_bloblist->push_back(p_lastblob);
    if (m_ptmax > m_hardveto) return Return_Value::New_Event;
    return Return_Value::Success;
  }
  /////////////////////////////////////////////////////////////////////////////////
  // If we have reached the end of MPI production with a meaningful event,
  // we can stop here.
  /////////////////////////////////////////////////////////////////////////////////
  if (p_activeMI->Done()) {
    if (!(p_activeMI->IsMinBias() &&
	  bloblist->size()==1 &&
	  ((*p_bloblist)[0]->Has(blob_status::needs_signal) ||
	   (*p_bloblist)[0]->Has(blob_status::needs_minBias))))
	return Return_Value::Nothing;
  }
  /////////////////////////////////////////////////////////////////////////////////
  // If it is a MinBias event where the event handler didn't manage to produce a
  // first scatter (i.e. the first blob still needs a signal) then we have to
  // produce a new event.
  /////////////////////////////////////////////////////////////////////////////////
  return Return_Value::New_Event;
}

bool Multiple_Interactions::CheckBlobList() 
{
  msg_Out()<<"   * "<<METHOD<<"\n";
  /////////////////////////////////////////////////////////////////////////////////
  // don't need to check trivial first MB blob
  /////////////////////////////////////////////////////////////////////////////////
  if (m_isfirstMB) return true;
  /////////////////////////////////////////////////////////////////////////////////
  // naive checks on blob list - does it exist and conserve momentum.
  /////////////////////////////////////////////////////////////////////////////////
  if (p_bloblist->empty()) {
    msg_Error()<<METHOD<<": incoming blob list is empty.\n";
    m_result = Return_Value::Error;
    return false;
  }
  if (!p_bloblist->FourMomentumConservation()) {
    msg_Tracking()<<METHOD<<" does not conserve four-momentum.\n";
    m_result = Return_Value::Retry_Event;
    return false;
  }
  /////////////////////////////////////////////////////////////////////////////////
  // check if there is a blob that must shower first.
  /////////////////////////////////////////////////////////////////////////////////
  for (Blob_List::const_iterator bit=p_bloblist->begin();
       bit!=p_bloblist->end();++bit) {
    if (((*bit)->Type()==btp::Hard_Collision ||
	 (*bit)->Type()==btp::Signal_Process) && 
	(*bit)->Has(blob_status::needs_showers)) {
      m_result = Return_Value::Nothing;
      return false;
    }
  }
  return BeamsViable();
}

bool Multiple_Interactions::BeamsViable() {
  /////////////////////////////////////////////////////////////////////////////////
  // Checking if the total energy in shower initiators exceeds the
  // energies of the incoming beams.  If yes, undo showering and
  // start with the signal process (if this happens after first shower)
  // or undo last MI interaction + shower.
  // This method uses some implicit knowledge -- it knows the sequence
  // of shower blobs is such that the last Hard Collision (i.e. last
  // MI interaction) gives rise to last shower blob, while the Signal
  // Process comes first.  It also knows that the shower initiators
  // in the In-state of the shower blobs are sorted ....
  /////////////////////////////////////////////////////////////////////////////////
  msg_Out()<<"   * "<<METHOD<<"\n";
  Blob_List isr=p_bloblist->Find(btp::Shower);
  for (Blob_List::iterator bit=isr.begin();bit!=isr.end();++bit) {
    if (m_ISblobs.find((*bit))!=m_ISblobs.end()) continue;
    if (!ExtractISInfo((*bit))) return false;
    m_ISblobs.insert((*bit));
  }
  m_result = Return_Value::Success;
  return true;
}

bool Multiple_Interactions::ExtractISInfo(Blob * blob) {
  msg_Out()<<"   * "<<METHOD<<"\n";
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

bool Multiple_Interactions::InitNewEvent() {
  msg_Out()<<"   * "<<METHOD<<"(new = "<<m_newevent<<", firstMB = "<<m_isfirstMB<<")\n";
  if (!m_newevent || m_isfirstMB) return true;
  p_lastblob = p_bloblist->FindFirst(btp::Signal_Process);
  if (p_lastblob->Has(blob_status::needs_signal)) return false;
  Blob_Data_Base * ptinfo=(*p_lastblob)["MI_Scale"];
  if (ptinfo==NULL) THROW(fatal_error,"No starting scale info in signal blob");
  m_ptmax=ptinfo->Get<double>();
  if (m_ptmax==std::numeric_limits<double>::max()) return false;
  double ptfac=sqrt((*p_lastblob)["Factorisation_Scale"]->Get<double>());
  double ptren=sqrt((*p_lastblob)["Renormalization_Scale"]->Get<double>());
  m_ptmax = ptfac;
  if (!IsZero(ptfac-ptren)) m_ptmax += ptren;
  p_activeMI->InitialiseMPIs(m_ptmax_fac*m_ptmax);
  p_lastblob->SetPosition(p_activeMI->SelectPositionForScatter());
  Blob * showerblob = p_lastblob->OutParticle(0)->DecayBlob();
  if (showerblob) showerblob->SetPosition(p_lastblob->Position());
  p_activeMI->Remnants()->SetImpactParameter(p_activeMI->ImpactParameter());
  m_newevent = false;
  msg_Out()<<"     maximal pt = "<<m_ptmax<<" from factors = "<<ptfac<<"/"<<ptren<<"\n";
  return true;
}  


void Multiple_Interactions::InitMinBiasEvent() {
  /*
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
  signal->AddData("Renormalization_Scale",
		  new Blob_Data<double>((*p_lastblob)["Renormalization_Scale"]->Get<double>()));
  signal->AddData("Factorization_Scale",
		  new Blob_Data<double>((*p_lastblob)["Factorization_Scale"]->Get<double>()));
  signal->AddData("Resummation_Scale",
		  new Blob_Data<double>((*p_lastblob)["Resummation_Scale"]->Get<double>()));
  delete p_lastblob;
  p_activeMI->Remnants()->SetImpactParameter(p_activeMI->ImpactParameter());
  m_newevent = m_isfirstMB = false;
  */
}


void Multiple_Interactions::SwitchPerturbativeInputsToMIs() {
  MODEL::as->SetActiveAs(PDF::isr::hard_subprocess);
  for (size_t i=0;i<2;i++) {
    double x_resc = m_emax[i]/p_activeMI->Remnants()->GetRemnant(i)->GetBeam()->OutMomentum()[0];
    p_activeMI->ISRHandler()->SetRescaleFactor(x_resc,i);
  }
}

bool Multiple_Interactions::MIKinematics() {
  p_activeMI->SetMaxEnergies(m_emax[0],m_emax[1]);
  return true; 
}

void Multiple_Interactions::Finish(const std::string &resultpath) {}

void Multiple_Interactions::CleanUp(const size_t & mode) 
{
  for (MI_Handler_Map::iterator mihit=p_mihandlers->begin();
       mihit!=p_mihandlers->end();mihit++) {
    mihit->second->CleanUp();
  }
  m_vetoed   = false;
  m_newevent = true;
  if (p_mihandlers->find(PDF::isr::hard_subprocess)!=p_mihandlers->end()) 
    p_activeMI = (*p_mihandlers)[PDF::isr::hard_subprocess];
}

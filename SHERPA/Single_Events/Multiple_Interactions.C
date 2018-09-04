#include "SHERPA/Single_Events/Multiple_Interactions.H"

#include "ATOOLS/Org/My_Limits.H"
#include "REMNANTS/Main/Remnant_Base.H"
#include "BEAM/Main/Beam_Base.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "SHERPA/PerturbativePhysics/Matrix_Element_Handler.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Default_Reader.H"
#include "MODEL/Main/Running_AlphaS.H"

using namespace SHERPA;
using namespace ATOOLS;

Multiple_Interactions::Multiple_Interactions(MI_Handler *mihandler):
  p_mihandler(mihandler), m_result(Return_Value::Nothing),
  m_newevent(true)
{
  m_name = std::string("Multiple_Interactions:")+p_mihandler->Name();
  m_type = eph::Perturbative;
  if (p_mihandler->Type()!=0) {
    m_ecms = sqrt(p_mihandler->ISRHandler()->Pole());
    for (size_t i=0;i<2;i++) p_remnants[i] = mihandler->Remnants()->GetRemnant(i);
    if (p_remnants[0]==NULL || p_remnants[1]==NULL) {
      THROW(fatal_error,"No beam remnant handler found.");
    }
  }
  Default_Reader read;
  read.SetInputPath(rpa->GetPath());
  read.SetInputFile(rpa->gen.Variable("RUN_DATA_FILE"));
  m_hardveto=read.GetValue<double>("MPI_PT_MAX",1.0e12);
  ResetIS();
}

Multiple_Interactions::~Multiple_Interactions() { }

Return_Value::code Multiple_Interactions::Treat(Blob_List *bloblist,double &weight)
{
  m_result   = Return_Value::Nothing; 
  if (p_mihandler->Type()==MI_Handler::None || p_mihandler->Done()) return m_result;
  p_bloblist = bloblist;
  // Try to colour-connect the last interaction with the remnants
  p_mihandler->ConnectColours(p_bloblist->FindLast(btp::Shower));
  // CheckBlobList makes sure a new interaction can be added.
  // If its the first then a completely new chain of 2->2 scatters 
  // must be initialised.  This is steered by a flag m_newevent, which is 
  // set to true in the CleanUp() method. 
  if (!CheckBlobList() || !InitNewEvent() || !MIKinematics()) return m_result;
  // Possibly switch to new PDF and alphaS.
  // TODO: will have to check that this happens.
  SwitchPerturbativeInputsToMIs();
  p_lastblob = p_mihandler->GenerateHardProcess();
  if (p_lastblob) {
    // This assumes that the scatters are ordered in transverse momentum.
    // Then maximal scale of subsequent scatters is given by the pT of the
    // previous ones.
    m_ptmax = p_lastblob->OutParticle(0)->Momentum().PPerp();
    // Check that the partons can be extracted from remnant - mainly a
    // confirmation that the remnant has enough energy to accommodate
    // the extra parton.
    for (size_t i=0;i<(size_t)p_lastblob->NInP();++i) {
      if (!p_remnants[i]->TestExtract(p_lastblob->InParticle(i))) {
	delete p_lastblob;
	return Return_Value::Retry_Event;
      }
    }
    bloblist->push_back(p_lastblob);
    if (m_ptmax > m_hardveto) return Return_Value::New_Event;
    return Return_Value::Success;
  }
  return Return_Value::Nothing;
}

bool Multiple_Interactions::CheckBlobList() 
{
  // naive checks on blob list - does it exist and conserve momentum.
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
  // check if there is a blob that must shower first.
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

void Multiple_Interactions::ResetIS() {
  if (p_mihandler->Type()!=0) {
    for (short unsigned int i=0;i<2;++i) {
      m_emax[i] = p_remnants[i]->GetBeam()->Energy();
      p_remnants[i]->Reset();
      p_mihandler->ISRHandler()->ResetRescaleFactor(i);
      p_mihandler->ISRHandler()->Reset(i);
    }
  }
  p_lastblob = NULL;
  m_ISblobs.clear();
}

bool Multiple_Interactions::BeamsViable() {
  // Checking if the total energy in shower initiators exceeds the
  // energies of the incoming beams.  If yes, undo showering and
  // start with the signal process (if this happens after first shower)
  // or undo last MI interaction + shower.
  // This method uses some implicit knowledge -- it knows the sequence
  // of shower blobs is such that the last Hard Collision (i.e. last
  // MI interaction) gives rise to last shower blob, while the Signal
  // Process comes first.  It also knows that the shower initiators
  // in the In-state of the shower blobs are sorted ....
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
  for (size_t i=0;i<blob->NInP();++i) {
    Particle *particle(blob->InParticle(i));
    if (particle->ProductionBlob()) continue;
    size_t beam = particle->Beam();
    if (!p_remnants[beam]->TestExtract(particle)) {
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
  if (!m_newevent) return true;
  p_lastblob = p_bloblist->FindFirst(btp::Signal_Process);
  if (p_lastblob->Has(blob_status::needs_signal)) return false;
  Blob_Data_Base * ptinfo=(*p_lastblob)["MI_Scale"];
  if (ptinfo==NULL) THROW(fatal_error,"No starting scale info in signal blob");
  m_ptmax=ptinfo->Get<double>();
  if (m_ptmax!=std::numeric_limits<double>::max()) {
    p_mihandler->InitialiseMPIs(4.*m_ptmax);
    p_lastblob->SetPosition(p_mihandler->SelectPositionForScatter());
    Blob * showerblob = p_lastblob->OutParticle(0)->DecayBlob();
    if (showerblob) showerblob->SetPosition(p_lastblob->Position());
    //msg_Out()<<METHOD<<" set x = "<<p_lastblob->Position()<<".\n";
    m_newevent = false;
    return true;
  }
  return false;
}
  
void Multiple_Interactions::SwitchPerturbativeInputsToMIs() {
  MODEL::as->SetActiveAs(PDF::isr::hard_subprocess);
  for (size_t i=0;i<2;i++) {
    double x_resc = m_emax[i]/p_remnants[i]->GetBeam()->Energy();
    //msg_Out()<<METHOD<<" sets PDF("<<i<<") x-rescaling = "<<x_resc<<" for E = "<<m_emax[i]<<"\n";
    p_mihandler->ISRHandler()->SetRescaleFactor(x_resc,i);
  }
}

bool Multiple_Interactions::MIKinematics() {
  p_mihandler->SetMaxEnergies(m_emax[0],m_emax[1]);
  return true; 
}

void Multiple_Interactions::Finish(const std::string &resultpath) {}

void Multiple_Interactions::CleanUp(const size_t & mode) 
{
  p_mihandler->CleanUp();
  ResetIS();
  m_vetoed   = false;
  m_newevent = true;
}

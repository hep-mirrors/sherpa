#include "SHERPA/Single_Events/Signal_Processes.H"

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/MCatNLO_Process.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "METOOLS/SpinCorrelations/Amplitude2_Tensor.H"
#include "METOOLS/SpinCorrelations/Decay_Matrix.H"
#include "METOOLS/SpinCorrelations/Spin_Density.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Random.H"
#include "MODEL/Main/Running_AlphaS.H"

using namespace SHERPA;
using namespace METOOLS;
using namespace ATOOLS;
using namespace PHASIC;

Signal_Processes::Signal_Processes(Matrix_Element_Handler * mehandler) :
  p_mehandler(mehandler), p_atensor(NULL), m_overweight(0.0)
{
  m_name="Signal_Processes";
  m_type=eph::Perturbative;
  p_remnants[0]=mehandler->GetISR()->GetRemnant(0);
  p_remnants[1]=mehandler->GetISR()->GetRemnant(1);
  if (p_remnants[0]==NULL || p_remnants[1]==NULL)
    THROW(critical_error,"No beam remnant handler found.");
}

Signal_Processes::~Signal_Processes()
{
  if (p_atensor) delete p_atensor;
}

Return_Value::code Signal_Processes::Treat(Blob_List * bloblist, double & weight)
{
  Blob *blob(bloblist->FindFirst(btp::Signal_Process));
  if (blob && blob->Has(blob_status::needs_signal)) {
    MODEL::as->SetActiveAs(PDF::isr::hard_process);
    while (true) {
      if (m_overweight>0.0) {
	if (m_overweight<ran->Get()) {
	  m_overweight=0.0;
	  CleanUp();
	  continue;
	}
	double overweight(m_overweight-1.0);
	if (!FillBlob(bloblist,blob))
	  THROW(fatal_error,"Internal error");
	(*blob)["Trials"]->Set(0.0);
	m_overweight=Max(overweight,0.0);
	weight = p_mehandler->WeightInfo().m_weight;
	return Return_Value::Success; 
      }
      if (p_mehandler->GenerateOneEvent() &&
	  FillBlob(bloblist,blob)) {
	weight = p_mehandler->WeightInfo().m_weight;
	return Return_Value::Success; 
      }
      else return Return_Value::New_Event;
    }
  }
  return Return_Value::Nothing;
}

bool Signal_Processes::FillBlob(Blob_List *const bloblist,Blob *const blob)
{
  DEBUG_FUNC(blob->Id());
  PHASIC::Process_Base *proc(p_mehandler->Process());
  blob->SetPosition(Vec4D(0.,0.,0.,0.));
  blob->SetTypeSpec(proc->Parent()->Name());
  if (p_mehandler->NLOMode()==3 && proc->Parent()->Info().m_fi.NLOType()!=nlo_type::lo) {
    MCatNLO_Process* powproc=dynamic_cast<MCatNLO_Process*>(proc->Parent());
    if (powproc) {
      if (powproc->WasSEvent()) blob->SetTypeSpec(proc->Parent()->Name()+"+S");
      else blob->SetTypeSpec(proc->Parent()->Name()+"+H");
    }
  }
  Vec4D cms = Vec4D(0.,0.,0.,0.);
  for (size_t i=0;i<proc->NIn();i++) 
    cms += proc->Integrator()->Momenta()[i];
  blob->SetCMS(cms);
  blob->DeleteOwnedParticles();
  blob->ClearAllData();
  bool success(true);
  Particle *particle(NULL);
  if (proc->Info().m_nlomode!=1)
    blob->SetStatus(blob_status::needs_showers | blob_status::needs_harddecays);
  else blob->SetStatus();
  for (unsigned int i=0;i<proc->NIn();i++) {
    particle = new Particle(0,proc->Flavours()[i],
			    proc->Integrator()->Momenta()[i]);
    particle->SetNumber(0);
    particle->SetStatus(part_status::decayed);
    particle->SetInfo('G');
    blob->AddToInParticles(particle);
    if (p_remnants[i]!=NULL) {
      p_remnants[i]->QuickClear();
      if (!p_remnants[i]->TestExtract(particle)) success=false;
    }
    else THROW(fatal_error,"No remnant found.");
  }
  for (unsigned int i=proc->NIn();
       i<proc->NIn()+proc->NOut();i++) {
    particle = new Particle(0,proc->Flavours()[i],
			    proc->Integrator()->Momenta()[i]);
    particle->SetNumber(0);
    particle->SetStatus(part_status::active);
    particle->SetInfo('H');
    blob->AddToOutParticles(particle);
  }
  DEBUG_VAR(*blob);
  PHASIC::Weight_Info winfo(p_mehandler->WeightInfo());
  double weight(winfo.m_weight);
  if (p_mehandler->EventGenerationMode()==1) {
    m_overweight=p_mehandler->WeightFactor()-1.0;
    if (m_overweight<0.0) m_overweight=0.0;
    else weight/=m_overweight+1.0;
  }
  blob->AddData("Weight",new Blob_Data<double>(weight));
  blob->AddData("Trials",new Blob_Data<double>(winfo.m_ntrial));
  blob->AddData("Enhance",new Blob_Data<double>
		(p_mehandler->Process()->Integrator()->EnhanceFactor()));
  blob->AddData("Factorisation_Scale",new Blob_Data<double>
		(sqrt(winfo.m_muf12*winfo.m_muf22)));
  blob->AddData("XF1",new Blob_Data<double>(winfo.m_xf1));
  blob->AddData("XF2",new Blob_Data<double>(winfo.m_xf2));

  ME_wgtinfo* wgtinfo=proc->GetMEwgtinfo();
  if (wgtinfo) {
    blob->AddData("ME_wgtinfo",new Blob_Data<ME_wgtinfo*>(wgtinfo));
    blob->AddData("Renormalization_Scale",new Blob_Data<double>(wgtinfo->m_mur2));
  }
  NLO_subevtlist* nlos=proc->GetSubevtList();
  if (nlos) blob->AddData("NLO_subeventlist",new Blob_Data<NLO_subevtlist*>(nlos));

  if (rpa->gen.HardSC() || rpa->gen.SoftSC()) {
    DEBUG_INFO("Filling amplitude tensor for spin correlations.");
    std::vector<Spin_Amplitudes> amps;
    std::vector<std::vector<Complex> > cols;
    proc->FillAmplitudes(amps, cols);
    DEBUG_VAR(amps[0]);
    Particle_Vector outparts=blob->GetOutParticles();
    Particle_Vector inparts=blob->GetInParticles();
    std::deque<Particle*> parts;
    parts.insert(parts.end(), inparts.begin(), inparts.end());
    parts.insert(parts.end(), outparts.begin(), outparts.end());
    if (p_atensor) delete p_atensor;
    p_atensor=new Amplitude2_Tensor(parts, amps, cols);
    DEBUG_VAR(*p_atensor);
    for (size_t i=0; i<inparts.size(); ++i) {
      Particle* part=inparts[i];
      DEBUG_INFO("contracting incoming particle "<<part->Flav());
      Spin_Density sigma(part);
      p_atensor->Contract(&sigma);
      DEBUG_VAR(*p_atensor);
    }
    for (size_t i=0; i<outparts.size(); ++i) {
      Particle* part=outparts[i];
      if (part->Flav().IsStable()) {
        DEBUG_INFO("contracting stable particle "<<part->Flav());
        Decay_Matrix dm(part);
        p_atensor->Contract(&dm);
        DEBUG_VAR(*p_atensor);
      }
    }
    DEBUG_VAR(*p_atensor);
    blob->AddData("ATensor",new Blob_Data<Amplitude2_Tensor*>(p_atensor));
  }
  return success;
}

void Signal_Processes::CleanUp() 
{ 
  if (m_overweight>0.0) return;
  if (p_mehandler)
    if (p_mehandler->Process())
      if (p_mehandler->Process()->Parent())
        if (p_mehandler->Process()->Parent()->Integrator())
          p_mehandler->Process()->Parent()->Integrator()->RestoreInOrder();
}

void Signal_Processes::Finish(const std::string &) 
{
}

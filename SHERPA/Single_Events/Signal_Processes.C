#include "SHERPA/Single_Events/Signal_Processes.H"

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace SHERPA;
using namespace ATOOLS;

Signal_Processes::Signal_Processes(Matrix_Element_Handler * mehandler) :
  p_mehandler(mehandler)
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
}

Return_Value::code Signal_Processes::Treat(Blob_List * bloblist, double & weight)
{
  Blob *blob(bloblist->FindFirst(btp::Signal_Process));
  if (blob==NULL) THROW(fatal_error,"Internal error");
  if (blob->Has(blob_status::needs_signal))
    while (true) {
      if (p_mehandler->GenerateOneEvent() &&
	  FillBlob(bloblist,blob)) {
	weight = p_mehandler->WeightInfo().m_weight;
	return Return_Value::Success; 
      }
    }
  return Return_Value::Nothing;
}

bool Signal_Processes::FillBlob(Blob_List *const bloblist,Blob *const blob)
{
  PHASIC::Process_Base *proc(p_mehandler->Process());
  blob->SetPosition(Vec4D(0.,0.,0.,0.));
  blob->SetTypeSpec(proc->Name());
  Vec4D cms = Vec4D(0.,0.,0.,0.);
  for (size_t i=0;i<proc->NIn();i++) 
    cms += proc->Integrator()->Momenta()[i];
  blob->SetCMS(cms);
  blob->SetBeam(-1);
  blob->DeleteOwnedParticles();
  blob->ClearAllData();
  bool success(true);
  Particle *particle(NULL);
  if (proc->Info().m_nlomode!=1)
    blob->SetStatus(blob_status::needs_showers |
		    blob_status::needs_harddecays);
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
  PHASIC::Weight_Info winfo(p_mehandler->WeightInfo());
  blob->AddData("Weight",new Blob_Data<double>(winfo.m_weight));
  blob->AddData("Trials",new Blob_Data<double>(winfo.m_ntrial));
  blob->AddData("Enhance",new Blob_Data<double>
		(p_mehandler->Process()->Integrator()->EnhanceFactor()));

  PHASIC::NLO_subevtlist* nlos=proc->GetSubevtList();
  if (nlos) blob->AddData("NLO_subeventlist",new Blob_Data<PHASIC::NLO_subevtlist*>(nlos));

  return success;
}

void Signal_Processes::CleanUp() 
{ 
}

void Signal_Processes::Finish(const std::string &) 
{
}

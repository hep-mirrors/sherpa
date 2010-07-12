#include "SHERPA/Single_Events/Signal_Processes.H"

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "METOOLS/Main/Spin_Structure.H"
#include "AMEGIC++/Main/Single_Process.H"

using namespace SHERPA;
using namespace ATOOLS;

Signal_Processes::Signal_Processes(Matrix_Element_Handler * mehandler) :
  p_mehandler(mehandler), p_amps(NULL)
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
  if (p_amps) delete p_amps; p_amps=NULL;
}

Return_Value::code Signal_Processes::Treat(Blob_List * bloblist, double & weight)
{
  Blob *blob(bloblist->FindFirst(btp::Signal_Process));
  if (blob && blob->Has(blob_status::needs_signal))
    while (true) {
      if (p_mehandler->GenerateOneEvent() &&
	  FillBlob(bloblist,blob)) {
	weight = p_mehandler->WeightInfo().m_weight;
	return Return_Value::Success; 
      }
      else return Return_Value::New_Event;
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
  blob->AddData("Factorisation_Scale",new Blob_Data<double>
		(sqrt(winfo.m_muf12*winfo.m_muf22)));
  blob->AddData("XF1",new Blob_Data<double>(winfo.m_xf1));
  blob->AddData("XF2",new Blob_Data<double>(winfo.m_xf2));

  PHASIC::ME_wgtinfo* wgtinfo=proc->GetMEwgtinfo();
  if (wgtinfo) {
    blob->AddData("ME_wgtinfo",new Blob_Data<PHASIC::ME_wgtinfo*>(wgtinfo));
    blob->AddData("Renormalization_Scale",new Blob_Data<double>
		(wgtinfo->m_renscale));
  }
  PHASIC::NLO_subevtlist* nlos=proc->GetSubevtList();
  if (nlos) blob->AddData("NLO_subeventlist",new Blob_Data<PHASIC::NLO_subevtlist*>(nlos));
  
  if (rpa.gen.SpinCorrelation()) {
    Particle_Vector inparticles = blob->GetInParticles();
    Particle_Vector outparticles = blob->GetOutParticles();
    Particle_Vector particles(inparticles.begin(),inparticles.end());
    particles.insert(particles.end(),outparticles.begin(),outparticles.end());
    if (p_amps) delete p_amps;
    p_amps = new METOOLS::Amplitude_Tensor(particles);
    
    AMEGIC::Single_Process* aproc = dynamic_cast<AMEGIC::Single_Process*>(proc);
    if (!aproc) {
      THROW(fatal_error, "Spin correlations can only be enabled when using "+
            std::string("Amegic as matrix element generator."));
    }
    aproc->FillAmplitudes(p_amps);
    blob->AddData("amps",new Blob_Data<METOOLS::Amplitude_Tensor*>(p_amps));
  }

  return success;
}

void Signal_Processes::CleanUp() 
{ 
}

void Signal_Processes::Finish(const std::string &) 
{
}

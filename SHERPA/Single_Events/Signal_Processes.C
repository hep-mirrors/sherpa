#include "SHERPA/Single_Events/Signal_Processes.H"

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/MCatNLO_Process.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "METOOLS/SpinCorrelations/Amplitude2_Tensor.H"
#include "METOOLS/SpinCorrelations/Decay_Matrix.H"
#include "METOOLS/SpinCorrelations/Spin_Density.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Math/Random.H"
#include "MODEL/Main/Running_AlphaS.H"

using namespace SHERPA;
using namespace METOOLS;
using namespace ATOOLS;
using namespace PHASIC;
using namespace std;

Signal_Processes::Signal_Processes(Matrix_Element_Handler * mehandler) :
  p_mehandler(mehandler), p_atensor(NULL), m_overweight(0.0)
{
  m_name="Signal_Processes";
  m_type=eph::Perturbative;
  p_remnants[0]=mehandler->GetISR()->GetRemnant(0);
  p_remnants[1]=mehandler->GetISR()->GetRemnant(1);
  if (p_remnants[0]==NULL || p_remnants[1]==NULL)
    THROW(critical_error,"No beam remnant handler found.");
  Data_Reader read(" ",";","!","=");
  m_setcolors=read.GetValue<int>("SP_SET_COLORS",0);
  m_cmode=ToType<int>(rpa->gen.Variable("METS_CLUSTER_MODE"));
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
  Cluster_Amplitude *ampl(NULL);
  if (p_mehandler->HasNLO()==3 &&
      proc->Parent()->Info().m_fi.NLOType()!=nlo_type::lo) {
    MCatNLO_Process* mcatnloproc=dynamic_cast<MCatNLO_Process*>(proc->Parent());
    if (mcatnloproc) {
      if (mcatnloproc->WasSEvent())
        blob->SetTypeSpec(proc->Parent()->Name()+"+S");
      else
        blob->SetTypeSpec(proc->Parent()->Name()+"+H");
      if (m_setcolors) ampl=mcatnloproc->GetAmplitude();
    }
  }
  else {
    if (m_setcolors) ampl=proc->Get<Single_Process>()->
		       Cluster(proc->Integrator()->Momenta(),m_cmode);
  }
  Vec4D cms = Vec4D(0.,0.,0.,0.);
  for (size_t i=0;i<proc->NIn();i++) 
    cms += proc->Integrator()->Momenta()[i];
  blob->SetCMS(cms);
  blob->DeleteOwnedParticles();
  blob->ClearAllData();
  bool success(true);
  Particle *particle(NULL);
  blob->SetStatus(blob_status::needs_harddecays);
  if (proc->Info().m_nlomode!=1)
    blob->AddStatus(blob_status::needs_showers);
  const DecayInfo_Vector &decs(proc->DecayInfos());
  blob->AddData("Decay_Info",new Blob_Data<DecayInfo_Vector>(decs));
  for (unsigned int i=0;i<proc->NIn();i++) {
    particle = new Particle(0,proc->Flavours()[i],
			    proc->Integrator()->Momenta()[i]);
    particle->SetNumber(0);
    particle->SetStatus(part_status::decayed);
    particle->SetInfo('G');
    blob->AddToInParticles(particle);
    if (ampl) {
      particle->SetFlow(1,ampl->Leg(i)->Col().m_j);
      particle->SetFlow(2,ampl->Leg(i)->Col().m_i);
    }
    if (p_remnants[i]!=NULL) {
      if (proc->NIn()>1) {
      p_remnants[i]->QuickClear();
      if (!p_remnants[i]->TestExtract(particle)) success=false;
      }
    }
    else THROW(fatal_error,"No remnant found.");
  }
  for (unsigned int i=proc->NIn();
       i<proc->NIn()+proc->NOut();i++) {
    particle = new Particle(0,proc->Flavours()[i],
			    proc->Integrator()->Momenta()[i]);
    particle->SetNumber(0);
    for (size_t j(0);j<decs.size();++j)
      if (decs[j]->m_id&(1<<i)) particle->SetMEId(1<<i);
    particle->SetStatus(part_status::active);
    particle->SetInfo('H');
    blob->AddToOutParticles(particle);
    if (ampl) {
      particle->SetFlow(1,ampl->Leg(i)->Col().m_i);
      particle->SetFlow(2,ampl->Leg(i)->Col().m_j);
    }
  }
  if (ampl && p_mehandler->HasNLO()==3 &&
      proc->Parent()->Info().m_fi.NLOType()!=nlo_type::lo) {
    if (ampl->NLO()&4) {
      blob->AddData("MC@NLO_KT2_Stop",new Blob_Data<double>(0.0));
      blob->AddData("MC@NLO_KT2_Start",new Blob_Data<double>(ampl->MuQ2()));
    }
    else if (ampl->Next()) {
      DEBUG_VAR(*ampl->Next());
      blob->AddData("MC@NLO_KT2_Stop",new Blob_Data<double>(ampl->KT2()));
      blob->AddData("MC@NLO_KT2_Start",new Blob_Data<double>
                    (ampl->Next()->Next()?ampl->Next()->KT2():ampl->MuQ2()));
    }
    blob->AddData("Resummation_Scale",new Blob_Data<double>(ampl->MuQ2()));
  }
  if (ampl) ampl->Delete();
  PHASIC::Weight_Info winfo(p_mehandler->WeightInfo());
  double weight(winfo.m_weight);
  if (p_mehandler->EventGenerationMode()==1) {
    m_overweight=p_mehandler->WeightFactor()-1.0;
    if (m_overweight<0.0) m_overweight=0.0;
    else weight/=m_overweight+1.0;
  }
  blob->AddData("Weight",new Blob_Data<double>(weight));
  blob->AddData("MEWeight",new Blob_Data<double>(winfo.m_dxs));
  blob->AddData("Weight_Norm",new Blob_Data<double>
		(p_mehandler->Sum()*rpa->Picobarn()));
  blob->AddData("Trials",new Blob_Data<double>(winfo.m_ntrial));
  blob->AddData("Enhance",new Blob_Data<double>
		(p_mehandler->Process()->Integrator()->EnhanceFactor()));
  blob->AddData("Factorisation_Scale",new Blob_Data<double>
		(sqrt(winfo.m_muf12*winfo.m_muf22)));
  blob->AddData("XF1",new Blob_Data<double>(winfo.m_xf1));
  blob->AddData("XF2",new Blob_Data<double>(winfo.m_xf2));
  blob->AddData("OQCD",new Blob_Data<int>
		(p_mehandler->Process()->OrderQCD()));
  blob->AddData("OEW",new Blob_Data<int>
		(p_mehandler->Process()->OrderEW()));
  blob->AddData("NLOQCDType",new Blob_Data<std::string>
		(ToString(p_mehandler->Process()->Info().m_fi.m_nloqcdtype)));
  blob->AddData("NLOEWType",new Blob_Data<std::string>
		(ToString(p_mehandler->Process()->Info().m_fi.m_nloewtype)));

  ME_wgtinfo* wgtinfo=proc->GetMEwgtinfo();
  if (wgtinfo) {
    blob->AddData("ME_wgtinfo",new Blob_Data<ME_wgtinfo*>(wgtinfo));
    blob->AddData("Renormalization_Scale",new Blob_Data<double>(wgtinfo->m_mur2));
  }
  NLO_subevtlist* nlos=proc->GetSubevtList();
  if (nlos) blob->AddData("NLO_subeventlist",new Blob_Data<NLO_subevtlist*>(nlos));

  if (rpa->gen.HardSC() || (rpa->gen.SoftSC() && !Flavour(kf_tau).IsStable())) {
    DEBUG_INFO("Filling amplitude tensor for spin correlations.");
    std::vector<Spin_Amplitudes> amps;
    std::vector<std::vector<Complex> > cols;
    proc->FillAmplitudes(amps, cols);
    DEBUG_VAR(amps[0]);
    Particle_Vector inparts=blob->GetInParticles();
    Particle_Vector outparts=blob->GetOutParticles();
    vector<pair<Particle*, size_t> > parts(inparts.size()+outparts.size());
    for (size_t i=0; i<inparts.size(); ++i)
      parts[i]=make_pair(inparts[i], i);
    for (size_t i=inparts.size(); i<inparts.size()+outparts.size(); ++i)
      parts[i]=make_pair(outparts[i-inparts.size()], i);

    DEBUG_INFO("particles before stability sorting:");
    for (size_t i=0; i<parts.size(); ++i) DEBUG_INFO(parts[i].first->RefFlav());
    stable_sort(parts.begin(), parts.end(), Amplitude2_Tensor::SortCrit);
    DEBUG_INFO("particles after stability sorting:");
    for (size_t i=0; i<parts.size(); ++i) DEBUG_INFO(parts[i].first->RefFlav());

    vector<int> permutation(parts.size(), -1);
    for (size_t i=0; i<parts.size(); ++i) permutation[parts[i].second]=i;
    DEBUG_INFO("permutation:");
    for (size_t i=0; i<parts.size(); ++i) DEBUG_INFO(permutation[i]);

    vector<int> spin_i(parts.size(), -1), spin_j(parts.size(), -1);
    vector<Particle*> partsonly(parts.size());
    for (size_t i=0; i<parts.size(); ++i) partsonly[i]=parts[i].first;
    if (p_atensor) delete p_atensor;
    p_atensor=new Amplitude2_Tensor(partsonly, permutation, 0, amps, cols,
                                    spin_i, spin_j);
    DEBUG_VAR(*p_atensor);
    blob->AddData("ATensor",new Blob_Data<Amplitude2_Tensor*>(p_atensor));
  }
  return success;
}

void Signal_Processes::CleanUp(const size_t & mode) 
{ 
  if (m_overweight>0.0) return;
}

void Signal_Processes::Finish(const std::string &) 
{
}

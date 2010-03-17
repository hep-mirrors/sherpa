#include "SHERPA/PerturbativePhysics/Perturbative_Interface.H"

#include "SHERPA/PerturbativePhysics/Shower_Handler.H"
#include "SHERPA/PerturbativePhysics/Matrix_Element_Handler.H"
#include "SHERPA/PerturbativePhysics/MI_Handler.H"
#include "SHERPA/SoftPhysics/Hadron_Decay_Handler.H"
#include "PDF/Main/Shower_Base.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Math/Random.H"

using namespace SHERPA;
using namespace PHASIC;
using namespace ATOOLS;

Perturbative_Interface::Perturbative_Interface
(Matrix_Element_Handler *const meh,Shower_Handler *const psh):
  p_me(meh), p_mi(NULL), p_hd(NULL), p_shower(psh), p_ampl(NULL) {}

Perturbative_Interface::Perturbative_Interface
(MI_Handler *const mi,Shower_Handler *const psh):
  p_me(NULL), p_mi(mi), p_hd(NULL), p_shower(psh), p_ampl(NULL) {}

Perturbative_Interface::Perturbative_Interface
(Hadron_Decay_Handler *const hdh,Shower_Handler *const psh):
  p_me(NULL), p_mi(NULL), p_hd(hdh), p_shower(psh), p_ampl(NULL) {}

Perturbative_Interface::~Perturbative_Interface() 
{
  if (p_ampl) p_ampl->Delete();
}

bool Perturbative_Interface::SetColours
(Cluster_Amplitude *ampl,Blob *const bl)
{
  msg_Debugging()<<METHOD<<"(): {\n";
  while (ampl->Prev()) ampl=ampl->Prev();
  for (size_t i(0);i<ampl->Legs().size();++i) {
    Cluster_Leg *cl(ampl->Leg(i));
    if (i<ampl->NIn()) {
      bl->InParticle(i)->SetFlow(1,Max(0,cl->Col().m_j-100));
      bl->InParticle(i)->SetFlow(2,Max(0,cl->Col().m_i-100));
    }
    else {
      bl->OutParticle(i-ampl->NIn())->SetFlow(1,Max(0,cl->Col().m_i-100));
      bl->OutParticle(i-ampl->NIn())->SetFlow(2,Max(0,cl->Col().m_j-100));
    }
  }
  msg_Debugging()<<*bl<<"\n";
  msg_Debugging()<<"}\n";
  return true;
}

Return_Value::code Perturbative_Interface::
DefineInitialConditions(ATOOLS::Blob *blob) 
{
  if (blob==NULL) {
    msg_Error()<<METHOD<<"(): Signal process not found."<<std::endl;
    return Return_Value::Error;
  }
  p_hard=blob;
  if (!p_shower->On()) {
    m_weight=1.0;
    return Return_Value::Success;
  }
  if (p_ampl) {
    p_ampl->Delete();
    p_ampl=NULL;
  }
  p_shower->CleanUp();
  msg_Indent();
  if (p_mi) {
    p_ampl=p_mi->ClusterConfiguration();
    if (!SetColours(p_ampl,blob)) return Return_Value::New_Event;
    if (!p_shower->GetShower()->PrepareShower(p_ampl)) return Return_Value::New_Event;
    return Return_Value::Success;
  }
  if (p_hd) {
    p_ampl=p_hd->ClusterConfiguration(blob);
    if (!SetColours(p_ampl,blob)) return Return_Value::New_Event;
    if (!p_shower->GetShower()->PrepareShower(p_ampl)) return Return_Value::New_Event;
    return Return_Value::Success;
  }
  p_me->Process()->Generator()->SetClusterDefinitions
    (p_shower->GetShower()->GetClusterDefinitions());
  p_ampl=p_me->Process()->Generator()->ClusterConfiguration(p_me->Process());
  if (p_ampl==NULL) return Return_Value::New_Event;
  if (p_me->Process()->Parent()->Info().m_fi.NLOType()&nlo_type::born) 
    p_ampl->SetRBMax(p_me->Process()->Integrator()->RBMax());
  if (!SetColours(p_ampl,blob)) return Return_Value::New_Event;
  if (!(p_me->Process()->Info().m_ckkw&1)) {
    m_weight=1.0;
  }
  else {
    m_weight=p_shower->GetShower()->CalculateWeight(p_ampl);
    blob->AddData("Sud_Weight",new Blob_Data<double>(m_weight));
    if (p_me->EventGenerationMode()==1) {
      if (m_weight>=ran.Get()) m_weight=1.0;
      else return Return_Value::New_Event;
    }
    Blob_Data_Base *winfo((*blob)["Weight"]);
    if (!winfo) THROW(fatal_error,"No weight information in signal blob");
    double meweight(winfo->Get<double>());
    blob->AddData("Weight",new Blob_Data<double>(meweight*m_weight));
  }
  if (!p_shower->GetShower()->PrepareShower(p_ampl)) return Return_Value::New_Event;
  return Return_Value::Success;
}

bool Perturbative_Interface::FillBlobs(ATOOLS::Blob_List *blobs)
{
  if (p_hard==NULL) return false;
  Blob *sblob = new Blob();
  sblob->SetType(btp::Shower);
  sblob->SetStatus(blob_status::needs_showers);
  sblob->SetId();
  if (p_shower->On()) {
    if (!p_hd)
      for (int i(0);i<p_hard->NInP();++i)
	sblob->AddToOutParticles(p_hard->InParticle(i));
    for (int i(0);i<p_hard->NOutP();++i)
      sblob->AddToInParticles(p_hard->OutParticle(i));
  }
  blobs->push_back(sblob);
  p_shower->FillBlobs(blobs); 
  return true;
}

int Perturbative_Interface::PerformShowers()
{
  PDF::Shower_Base *csh(p_shower->GetShower());
  if (p_me && 
      (p_me->Process()->Parent()->Info().
       m_fi.NLOType()&nlo_type::born)) {
    while (csh->TrialEmission()) {
      Cluster_Amplitude *ampl(csh->GetRealEmissionAmplitude());
      double ps(csh->TrialWeight(ampl)), me(ps);
      msg_Debugging()<<METHOD<<"():  me / ps = "
		     <<me<<" / "<<ps<<" = "<<me/ps;
      double weight(me/ps);
      if (weight>1.0) {
	msg_Info()<<METHOD<<"(): '"<<p_me->Process()->Name()
		  <<"' w_{me}/w_{ps} = "<<weight<<std::endl;
	Process_Integrator *pint(p_me->Process()->Integrator());
	pint->AddRBPoint(weight*pint->RBMax());
	pint->WriteOutRB(p_me->RPath()+"/"+p_me->Process()->
			 Generator()->Name()+"_"+csh->Name());
      }
      if (weight>ran.Get()) {
	msg_Debugging()<<" -> accept\n";
	return csh->PerformShowers();
      }
      msg_Debugging()<<" -> reject\n";
      for (size_t i(0);i<p_ampl->Legs().size();++i)
	p_ampl->Leg(i)->SetQ2Shower(ampl->Leg(i)->Q2Shower());
      csh->CleanUp();
      csh->PrepareShower(p_ampl);
    }
    return true;
  }
  int stat=csh->PerformShowers();
  double weight=csh->Weight();
  p_hard->AddData("Shower_Weight",new Blob_Data<double>(weight));
  Blob_Data_Base *winfo((*p_hard)["Weight"]);
  if (!winfo) THROW(fatal_error,"No weight information in signal blob");
  double meweight(winfo->Get<double>());
  p_hard->AddData("Weight",new Blob_Data<double>(meweight*weight));
  return stat;
}

int Perturbative_Interface::PerformDecayShowers()
{ 
  return p_shower->GetShower()->PerformDecayShowers(); 
}

void Perturbative_Interface::CleanUp()
{
}


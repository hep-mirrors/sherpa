#include "SHERPA/Single_Events/Multiple_Interactions.H"

#include "ATOOLS/Org/My_Limits.H"
#include "PDF/Remnant/Remnant_Base.H"
#include "BEAM/Main/Beam_Base.H"
#include "PDF/Main/ISR_Handler.H"
#include "PHASIC++/Selectors/KT_Finder.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "SHERPA/PerturbativePhysics/Matrix_Element_Handler.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "MODEL/Main/Running_AlphaS.H"

#ifdef USING__Amisic
#include "AMISIC++/Main/Amisic.H"
#endif

using namespace SHERPA;
using namespace ATOOLS;
#ifdef USING__Amisic
using namespace AMISIC;
#endif

Multiple_Interactions::Multiple_Interactions(MI_Handler *mihandler):
  p_mihandler(mihandler), p_jetfinder(NULL)
{
  m_name = std::string("Multiple_Interactions:")+p_mihandler->Name();
  m_type = eph::Perturbative;
  if (p_mihandler->Type()!=0) {
    m_ecms = sqrt(p_mihandler->ISRHandler()->Pole());
    p_remnants[0]=mihandler->ISRHandler()->GetRemnant(0);
    p_remnants[1]=mihandler->ISRHandler()->GetRemnant(1);
    if (p_remnants[0]==NULL || p_remnants[1]==NULL) {
      THROW(fatal_error,"No beam remnant handler found.");
    }
  }
}

Multiple_Interactions::~Multiple_Interactions() 
{
  if (p_jetfinder==NULL) delete p_jetfinder;
}

Return_Value::code Multiple_Interactions::CheckBlobList(ATOOLS::Blob_List *const bloblist) 
{
  p_bloblist=bloblist;
  if (m_vetoed) return Return_Value::Nothing;
  if (!p_bloblist->FourMomentumConservation()) {
    msg_Error()<<"Multiple_Interactions::CheckBlobList(..): "
	       <<"Retry event "<<rpa->gen.NumberOfGeneratedEvents()<<std::endl;
    return Return_Value::Retry_Event;
  }
  for (Blob_List::const_iterator bit=bloblist->begin();
       bit!=bloblist->end();++bit) {
    if ((*bit)->Type()==btp::Hard_Collision ||
	(*bit)->Type()==btp::Signal_Process) 
      if ((*bit)->Has(blob_status::needs_showers)) 
	return Return_Value::Nothing;
  }
  for (short unsigned int i=0;i<2;++i) {
    m_emax[i]=p_remnants[i]->GetBeam()->Energy();
    p_mihandler->ISRHandler()->Reset(i);
    p_remnants[i]->QuickClear();
  }
  Blob_List isr=bloblist->Find(btp::Shower);
  for (Blob_List::iterator iit=isr.begin();iit!=isr.end();++iit) {
    for (int beam(0), i(0);i<(*iit)->NInP();++i) {
      Particle *cp((*iit)->InParticle(i));
      if (cp->ProductionBlob()) continue;
      m_emax[beam]-=cp->Momentum()[0];
      p_mihandler->ISRHandler()->
	Extract(cp->Flav(),cp->Momentum()[0],beam);
      if (!p_remnants[beam]->Extract(cp)) {
	msg_Tracking()<<METHOD<<"(): Cannot extract parton from hadron. \n"
		      <<*cp<<std::endl;
	if (!(*iit)->IsConnectedTo(btp::Signal_Process))
	  p_bloblist->DeleteConnected(*iit);
	else return Return_Value::Retry_Event;
	return Return_Value::Retry_Phase;
      }
      ++beam;
    } 
  }
  if (m_generated) return Return_Value::Success;
  Blob * signal=bloblist->FindFirst(btp::Signal_Process);
  if (signal->Has(blob_status::needs_signal)) return Return_Value::Nothing;
  Blob_Data_Base *ptinfo=(*signal)["MI_Scale"];
  if (ptinfo==NULL)
    THROW(fatal_error,"No starting scale info in signal blob");
  m_ptmax=ptinfo->Get<double>();
  if (m_ptmax!=std::numeric_limits<double>::max()) return Return_Value::Success;
  return Return_Value::Nothing;
}

Return_Value::code Multiple_Interactions::Treat(ATOOLS::Blob_List *bloblist,double &weight)
{
#ifdef USING__Amisic
  if (p_mihandler->Type()==MI_Handler::None ||
      MI_Base::StopGeneration()) return Return_Value::Nothing;
  if (bloblist->empty()) {
    msg_Error()<<"Multiple_Interactions::Treat(): "
		       <<"Incoming blob list is empty!"<<std::endl;
    return Return_Value::Error;
  }
  Return_Value::code cbc(CheckBlobList(bloblist));
  if (cbc!=Return_Value::Success) return cbc;
  MODEL::as->SetActiveAs(PDF::isr::hard_subprocess);
  p_mihandler->SetScaleMax(m_emax[0],2);
  p_mihandler->SetScaleMax(m_emax[1],3);
  if (!m_generated) {
    p_mihandler->SetScaleMax(m_ptmax,0);
    p_mihandler->SetScaleMin(p_mihandler->ScaleMin(0),0);
    p_mihandler->Reset();
    m_generated=true;
  }
  Blob *blob(NULL);
  bool success=false;
  if (!m_vetoed) {
    if (m_ptmax<=p_mihandler->ScaleMin(0)) {
      return Return_Value::Nothing;
    }
    else {
      blob = new Blob();
      success=p_mihandler->GenerateHardProcess(blob);
    }
  }
  else if (m_vetoed) {
    blob = new Blob();
    success=p_mihandler->GenerateSoftProcess(blob);
    // dummy settings for analysis
    blob->SetType(btp::Soft_Collision);
    blob->SetTypeSpec("Soft UE");
    blob->SetStatus(blob_status::needs_showers);
    blob->AddData("ME_Weight",new Blob_Data<double>(m_weight));
    blob->AddData("ME_NumberOfTrials",new Blob_Data<int>((int)m_ntrials));
  }
  if (success) {
    blob->SetId(bloblist->size());
    m_ptmax=blob->OutParticle(0)->Momentum().PPerp();
    for (size_t i=0;i<(size_t)blob->NInP();++i) {
      if (!p_remnants[i]->Extract(blob->InParticle(i))) {
	msg_Tracking()<<"Multiple_Interactions::Treat(..): "
		      <<"Cannot extract parton from hadron. \n"
		      <<*blob->InParticle(0)<<std::endl;
	delete blob;
	return Return_Value::Retry_Phase;
      }
    }
    blob->SetStatus(blob_status::needs_showers);
    blob->AddData("Weight",new Blob_Data<double>(1.0));
    bloblist->push_back(blob);
    return Return_Value::Success;
  }
  delete blob;
  p_mihandler->ISRHandler()->Reset(0);
  p_mihandler->ISRHandler()->Reset(1);
  if (!MI_Base::StopGeneration()) return Return_Value::Retry_Phase;
#endif
  return Return_Value::Nothing;
}

bool Multiple_Interactions::VetoHardProcess(ATOOLS::Blob *const blob)
{
  if (p_mihandler->VetoHardProcess(blob)) {
    m_weight=(*blob)["ME_Weight"]->Get<double>();
    m_ntrials=(*blob)["ME_NumberOfTrials"]->Get<int>();
    p_bloblist->DeleteConnected(blob);
    p_bloblist->AddBlob(btp::Signal_Process);
    return m_vetoed=true;
  }
  return false;
}

void Multiple_Interactions::Finish(const std::string &resultpath) 
{
}

void Multiple_Interactions::CleanUp(const size_t & mode) 
{
  p_mihandler->CleanUp();
  m_ptmax=std::numeric_limits<double>::max();
  m_vetoed=false;
  m_generated=false;
}

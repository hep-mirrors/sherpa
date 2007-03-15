#include "Multiple_Interactions.H"

#include "My_Limits.H"
#include "Remnant_Base.H"
#include "Beam_Base.H"
#include "ISR_Handler.H"
#include "Jet_Finder.H"
#include "Kt_Algorithm.H"
#include "Run_Parameter.H"
#include "Exception.H"
#include "Matrix_Element_Handler.H"
#include "XS_Base.H"
#include "MyStrStream.H"

#ifdef USING__Amisic
#include "Amisic.H"
#endif

#ifdef PROFILE__all
#define PROFILE__Multiple_Interactions
#endif
#ifdef PROFILE__Multiple_Interactions
#include "prof.hh"
#else 
#define PROFILE_HERE 
#endif

using namespace SHERPA;
using namespace EXTRAXS;
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
    msg.Error()<<"Multiple_Interactions::CheckBlobList(..): "
	       <<"Retry event "<<rpa.gen.NumberOfDicedEvents()<<std::endl;
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
  Blob_List isr=bloblist->Find(btp::IS_Shower);
  for (Blob_List::reverse_iterator iit=isr.rbegin();
       iit!=isr.rend();++iit) {
    if ((*iit)->InParticle(0)->Momentum().Nan()) {
      msg.Error()<<METHOD<<"(): Nan momentum for "
		 <<*(*iit)->InParticle(0)<<"\n  Kill subprocess."<<std::endl;
      if (!(*iit)->IsConnectedTo(btp::Signal_Process))
	p_bloblist->DeleteConnected(*iit);
      else return Return_Value::Retry_Event;
      return Return_Value::Retry_Phase;
    }
    m_emax[(*iit)->Beam()]-=(*iit)->InParticle(0)->Momentum()[0];
    p_mihandler->ISRHandler()->
      Extract((*iit)->InParticle(0)->Flav(),
	      (*iit)->InParticle(0)->Momentum()[0],(*iit)->Beam());
    if (!p_remnants[(*iit)->Beam()]->Extract((*iit)->InParticle(0))) {
      msg_Tracking()<<METHOD<<"(): Cannot extract parton from hadron. \n"
		    <<*(*iit)->InParticle(0)<<std::endl;
      if (!(*iit)->IsConnectedTo(btp::Signal_Process))
	p_bloblist->DeleteConnected(*iit);
      else return Return_Value::Retry_Event;
      return Return_Value::Retry_Phase;
    } 
  }
  if (m_diced) return Return_Value::Success;
  Blob * signal=bloblist->FindFirst(btp::Hard_Collision);
  if (signal==NULL) signal=bloblist->FindFirst(btp::Signal_Process);
  if (signal->Has(blob_status::needs_signal)) return Return_Value::Nothing;
  if (!m_diced) {
    m_ptmax=ATOOLS::rpa.gen.Ecms()/2.0;
    if (VetoHardProcess(signal)) {
      m_ptmax=0.0;
      return Return_Value::Success;
    }
  }
  switch (p_mihandler->ScaleScheme()) {
    // min p_{T, out}
  case 2: {
    bool construct=false;
    Blob_Data_Base *xsinfo=
      (*bloblist->FindLast(btp::ME_PS_Interface_FS))["Core_Process"];
    if (xsinfo==NULL) {
      construct=true;
    }
    else {
      XS_Base *xs=xsinfo->Get<XS_Base*>();
      if (xs==NULL) {
	construct=true;
      }
      else {
	bool one=false;
	for (short unsigned int i=2;i<4;++i) 
	  if (xs->Flavours()[i].Strong()) {
	    m_ptmax=Min(m_ptmax,xs->Momenta()[i].PPerp());
	    one=true;
	  }
	if (!one) construct=true;
      }
    }
    if (construct) {
      Blob_List shower=bloblist->
	FindConnected(m_diced?bloblist->FindLast(btp::Hard_Collision):signal);
      Particle_List jets=shower.ExtractLooseParticles(1);
      Is_Parton isparton;
      jets.Keep(&isparton);
      if (p_jetfinder==NULL) 
	p_jetfinder = new Jet_Finder(ToString(p_mihandler->YCut()),4);
      p_jetfinder->ConstructJets(&jets,1,true);
      if (jets.size()>0) {
	for (size_t i=0;i<jets.size();++i)
	  m_ptmax=Min(m_ptmax,jets[i]->Momentum().PPerp());
      }
      jets.Clear();
    }
    break;
  }
  case 1: {
    Blob_List shower=bloblist->
      FindConnected(m_diced?bloblist->FindLast(btp::Hard_Collision):signal);
    int njets=2;
    if (!m_diced) {
      Blob_Data_Base *xsinfo=
	(*bloblist->FindFirst(btp::ME_PS_Interface_FS))["Core_Process"];
      if (xsinfo==NULL) {
	njets=1;
      }
      else {
	XS_Base *xs=xsinfo->Get<XS_Base*>();
	if (xs!=NULL) if (xs->NStrong()<4) njets=1;
      }
    }
    Particle_List jets=shower.ExtractLooseParticles(1);
    if (p_jetfinder==NULL) 
      p_jetfinder = new Jet_Finder(ToString(p_mihandler->YCut()),4);
    p_jetfinder->ConstructJets(&jets,njets,true);
    double ptmax=0.0;
    for (Particle_List::const_iterator pit=jets.begin();
	 pit!=jets.end();++pit)
      ptmax=Max(ptmax,(*pit)->Momentum().PPerp());
    jets.Clear();
    m_ptmax=Min(m_ptmax,ptmax);
    break;
  }
    // mean p_{T, out}
  case 0: {
    Blob_Data_Base *xsinfo=
      (*bloblist->FindLast(btp::ME_PS_Interface_FS))["Core_Process"];
    if (xsinfo==NULL) THROW(critical_error,"No merging information.");
    XS_Base *xs=xsinfo->Get<XS_Base*>();
    if (xs==NULL) THROW(critical_error,"No merging information.");
    m_ptmax=0.0;
    double nstrong=0.0;
    for (short unsigned int i=2;i<4;++i) 
      if (xs->Flavours()[i].Strong()) {
	m_ptmax+=xs->Momenta()[i].PPerp();
	++nstrong;
      }
    m_ptmax/=nstrong;
    break;
  }
    // abort
  default: THROW(not_implemented,"Wrong mi scale scheme.");
  }
  if (!m_diced) {
    signal->AddData("MI_Scale",new Blob_Data<double>(m_ptmax));
  }
  if (m_ptmax!=std::numeric_limits<double>::max()) return Return_Value::Success;
  return Return_Value::Nothing;
}

Return_Value::code Multiple_Interactions::Treat(ATOOLS::Blob_List *bloblist,double &weight)
{
#ifdef USING__Amisic
  PROFILE_HERE;
  if (p_mihandler->Type()==MI_Handler::None ||
      MI_Base::StopGeneration()) return Return_Value::Nothing;
  if (bloblist->empty()) {
    msg.Error()<<"Multiple_Interactions::Treat(): "
		       <<"Incoming blob list is empty!"<<std::endl;
    return Return_Value::Error;
  }
  Return_Value::code cbc(CheckBlobList(bloblist));
  if (cbc!=Return_Value::Success) return cbc;
  p_mihandler->SetScaleMax(m_emax[0],2);
  p_mihandler->SetScaleMax(m_emax[1],3);
  if (!m_diced) {
    p_mihandler->SetScaleMax(m_ptmax,0);
    p_mihandler->SetScaleMin(p_mihandler->ScaleMin(0),0);
    p_mihandler->Reset();
    m_diced=true;
  }
  Blob *blob(NULL);
  bool success=false;
  if (!m_vetoed) {
    if (m_ptmax<=p_mihandler->ScaleMin(0)) {
      return Return_Value::Nothing;
    }
    else {
      blob = new Blob();
      blob->AddData("MI_Scale",new Blob_Data<double>(m_ptmax));
      success=p_mihandler->GenerateHardProcess(blob);
    }
  }
  else if (m_vetoed) {
    blob = new Blob();
    blob->AddData("MI_Scale",new Blob_Data<double>(m_ptmax));
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
    bloblist->push_back(blob);
    return Return_Value::Success;
  }
  delete blob;
  p_mihandler->ISRHandler()->Reset(0);
  p_mihandler->ISRHandler()->Reset(1);
#endif
  if (!MI_Base::StopGeneration()) return Return_Value::Retry_Phase;
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

void Multiple_Interactions::CleanUp() 
{
  p_mihandler->CleanUp();
  m_ptmax=std::numeric_limits<double>::max();
  m_vetoed=false;
  m_diced=false;
}

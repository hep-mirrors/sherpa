#include "Multiple_Interactions.H"

#include "Amisic.H"
#include "My_Limits.H"
#include "Remnant_Base.H"
#include "ISR_Handler.H"
#include "Jet_Finder.H"
#include "Kt_Algorithm.H"
#include "Run_Parameter.H"

#ifdef PROFILE__all
#define PROFILE__Multiple_Interactions
#endif
#ifdef PROFILE__Multiple_Interactions
#include "prof.hh"
#else 
#define PROFILE_HERE 
#endif

using namespace SHERPA;
using namespace AMISIC;
using namespace EXTRAXS;
using namespace ATOOLS;

Multiple_Interactions::Multiple_Interactions(MI_Handler *mihandler):
  p_mihandler(mihandler)
{
  m_name = std::string("Multiple_Interactions:")+p_mihandler->Name();
  m_type = eph::Perturbative;
  m_ecms = sqrt(p_mihandler->ISRHandler()->Pole());
  p_remnants[0]=GET_OBJECT(Remnant_Base,"Remnant_Base_0");
  p_remnants[1]=GET_OBJECT(Remnant_Base,"Remnant_Base_1");
  if (p_remnants[0]==NULL || p_remnants[1]==NULL) {
    THROW(fatal_error,"No beam remnant handler found.");
  }
}

Multiple_Interactions::~Multiple_Interactions() 
{
}

bool Multiple_Interactions::CheckBlobList(ATOOLS::Blob_List *const bloblist) 
{
  p_bloblist=bloblist;
  if (m_vetoed) return false;
  if (!p_bloblist->FourMomentumConservation()) {
    msg.Error()<<"Multiple_Interactions::CheckBlobList(..): "
	       <<"Retry event "<<rpa.gen.NumberOfDicedEvents()<<std::endl;
    p_bloblist->Clear();
    return false;
  }
  for (Blob_List::const_iterator bit=bloblist->begin();
       bit!=bloblist->end();++bit) {
    if ((*bit)->Type()==btp::Hard_Collision ||
	(*bit)->Type()==btp::Signal_Process) 
      if ((*bit)->Status()!=0) return false;
  }
  for (short unsigned int i=0;i<2;++i) {
    m_emax[i]=p_remnants[i]->BeamEnergy();
    p_mihandler->ISRHandler()->Reset(i);
    p_remnants[i]->QuickClear();
  }
  Blob_List isr=bloblist->Find(btp::IS_Shower);
  static double ntot=0.0;
  ++ntot;
  for (Blob_List::reverse_iterator iit=isr.rbegin();
       iit!=isr.rend();++iit) {
    m_emax[(*iit)->Beam()]-=(*iit)->InParticle(0)->Momentum()[0];
    p_mihandler->ISRHandler()->
      Extract((*iit)->InParticle(0)->Flav(),
	      (*iit)->InParticle(0)->Momentum()[0],(*iit)->Beam());
    if (!p_remnants[(*iit)->Beam()]->Extract((*iit)->InParticle(0))) {
      msg_Tracking()<<"Multiple_Interactions::CheckBlobList(..): "
		    <<"Cannot extract parton from hadron. \n"
		    <<*(*iit)->InParticle(0)<<std::endl;
      p_bloblist->DeleteConnected(*iit);
      if (bloblist->empty()) {
	Blob *blob = new Blob();
	blob->SetType(btp::Signal_Process);
	blob->SetStatus(-1);
	blob->SetId();
	blob->SetStatus(2);
	bloblist->push_back(blob);	  
      }
      static double nrej=0.0;
      if (10*++nrej>rpa.gen.NumberOfDicedEvents())
	ATOOLS::msg.Error()<<"Multiple_Interactions::CheckBlobList(..): "
			   <<"Shower rejection rate is "
			   <<nrej/ntot<<"."<<std::endl;
      return false;
    } 
  }
  if (m_diced) return true;
  Blob *signal=bloblist->FindFirst(btp::Signal_Process);
  if (!m_diced) {
    m_ptmax=ATOOLS::rpa.gen.Ecms()/2.0;
    if (VetoHardProcess(signal)) {
      m_ptmax=0.0;
      return true;
    }
  }
  switch (p_mihandler->ScaleScheme()) {
    // min p_{T, out}
  case 1: {
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
      Jet_Finder finder(p_mihandler->YCut(),4);
      finder.ConstructJets(&jets,1,true);
      if (jets.size()>0) {
	for (size_t i=0;i<jets.size();++i)
	  m_ptmax=Min(m_ptmax,jets[i]->Momentum().PPerp());
      }
      jets.Clear();
    }
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
  return m_ptmax!=std::numeric_limits<double>::max();
}

bool Multiple_Interactions::Treat(ATOOLS::Blob_List *bloblist,double &weight)
{
  PROFILE_HERE;
  if (p_mihandler->Type()==MI_Handler::None ||
      MI_Base::StopGeneration()) return false;
  if (bloblist->empty()) {
    msg.Error()<<"Multiple_Interactions::Treat(): "
		       <<"Incoming blob list is empty!"<<std::endl;
    return false;
  }
  if (!CheckBlobList(bloblist)) return false;
  p_mihandler->SetScaleMax(m_emax[0],2);
  p_mihandler->SetScaleMax(m_emax[1],3);
  if (!m_diced) {
    p_mihandler->SetScaleMax(m_ptmax,0);
    p_mihandler->SetScaleMin(p_mihandler->ScaleMin(0),0);
    p_mihandler->Reset();
    m_diced=true;
  }
  Blob *blob = new Blob();
  blob->AddData("MI_Scale",new Blob_Data<double>(m_ptmax));
  bool success=false;
  if (!m_vetoed && m_ptmax>p_mihandler->ScaleMin(0)) {
    success=p_mihandler->GenerateHardProcess(blob);
  }
  else if (m_vetoed) {
    success=p_mihandler->GenerateSoftProcess(blob);
    // dummy settings for analysis
    blob->SetType(btp::Signal_Process);
    blob->SetTypeSpec("Soft UE");
    blob->SetStatus(5);
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
	return true;
      }
    }
    bloblist->push_back(blob);
    return true;
  }
  delete blob;
  p_mihandler->ISRHandler()->Reset(0);
  p_mihandler->ISRHandler()->Reset(1);
  return false;
}

bool Multiple_Interactions::VetoHardProcess(ATOOLS::Blob *const blob)
{
  if (p_mihandler->VetoHardProcess(blob)) {
    m_weight=(*blob)["ME_Weight"]->Get<double>();
    m_ntrials=(*blob)["ME_NumberOfTrials"]->Get<int>();
    p_bloblist->DeleteConnected(blob);
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

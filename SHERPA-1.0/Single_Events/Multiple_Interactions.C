#include "Multiple_Interactions.H"

#include "Amisic.H"
#include "My_Limits.H"
#include "Remnant_Base.H"
#include "ISR_Handler.H"

#ifdef PROFILE__all
#define PROFILE__Multiple_Interactions
#endif
#ifdef PROFILE__Multiple_Interactions
#include "prof.hh"
#else 
#define PROFILE_HERE 
#endif

using namespace SHERPA;

Multiple_Interactions::Multiple_Interactions(MI_Handler *mihandler):
  p_mihandler(mihandler)
{
  m_name = std::string("Multiple_Interactions:")+p_mihandler->Name();
  m_type = eph::Perturbative;
  m_ecms = sqrt(p_mihandler->ISRHandler()->Pole());
  p_remnants[0]=GET_OBJECT(Remnant_Base,"Remnant_Base_0");
  p_remnants[1]=GET_OBJECT(Remnant_Base,"Remnant_Base_1");
  if (p_remnants[0]==NULL || p_remnants[1]==NULL) {
    throw(ATOOLS::Exception(ATOOLS::ex::fatal_error,"No beam remnant handler found.",
			    "Multiple_Interactions","Multiple_Interactions"));
  }
}

Multiple_Interactions::~Multiple_Interactions() 
{
}

bool Multiple_Interactions::CheckBlobList(ATOOLS::Blob_List *const bloblist) 
{
  p_bloblist=bloblist;
  for (size_t i=0;i<2;++i) {
    m_xmax[i]=1.0;
    p_remnants[i]->QuickClear();
  }
  for (ATOOLS::Blob_List::const_iterator bit=bloblist->begin();
       bit!=bloblist->end();++bit) {
    if ((*bit)->Type()==ATOOLS::btp::Hard_Collision) m_diced=true;
    if ((*bit)->Type()==ATOOLS::btp::Hard_Collision ||
	(*bit)->Type()==ATOOLS::btp::Signal_Process) 
      if ((*bit)->Status()!=0) return false;
  }
  for (ATOOLS::Blob_List::const_iterator bit=bloblist->end();bit!=bloblist->begin();) {
    --bit;
    if ((*bit)->Type()==ATOOLS::btp::ME_PS_Interface_FS) {
      if (!m_diced) {
	ATOOLS::Blob_Data_Base *info=(*(*bit))["MI_Scale"];
	if (info==NULL) {
	  ATOOLS::msg.Error()<<"Multiple_Interactions::CheckBlobList(..): "
			     <<"No scale information in merging blob. "
			     <<"Taking P_\\perp."<<std::endl;
	  for (int i=0;i<(*bit)->NInP();++i) 
	    m_ptmax=ATOOLS::Min(m_ptmax,(*bit)->InParticle(i)->Momentum().PPerp());
	}
	else {
	  m_ptmax=sqrt(ATOOLS::Min(m_ptmax,info->Get<double>()));
	}
	if (VetoHardProcess(*bit)) {
	  m_ptmax=0.0;
	  break;
	}
      }
    }
    else if ((*bit)->Type()==ATOOLS::btp::IS_Shower) {
      m_xmax[(*bit)->Beam()]-=2.0*(*bit)->InParticle(0)->Momentum()[0]/m_ecms;
      p_mihandler->ISRHandler()->
	Extract((*bit)->InParticle(0)->Flav(),
		(*bit)->InParticle(0)->Momentum()[0],(*bit)->Beam());
      if (!p_remnants[(*bit)->Beam()]->Extract((*bit)->InParticle(0))) {
	msg_Tracking()<<"Multiple_Interactions::CheckBlobList(..): "
		      <<"Cannot extract parton from hadron. \n"
		      <<*(*bit)->InParticle(0)<<std::endl;
	m_deleted.clear();
	DeleteConnected(*bit);
	if (bloblist->empty()) {
	  ATOOLS::Blob *blob = new ATOOLS::Blob();
	  blob->SetType(ATOOLS::btp::Signal_Process);
	  blob->SetStatus(-1);
	  blob->SetId();
	  blob->SetStatus(2);
	  bloblist->push_back(blob);	  
	}
	return true; // blob-list changed
      } 
    }
  }
  if (m_ptmax==std::numeric_limits<double>::max()) {
    for (ATOOLS::Blob_List::const_iterator bit=bloblist->begin();
	 bit!=bloblist->end();++bit) {
      if ((*bit)->Type()==ATOOLS::btp::Hard_Collision || 
	  (*bit)->Type()==ATOOLS::btp::Signal_Process) {
	for (int i=0;i<(*bit)->NOutP();++i) 
	  m_ptmax=ATOOLS::Min(m_ptmax,(*bit)->OutParticle(i)->Momentum().PPerp());
	if (VetoHardProcess(*bit)) {
	  m_ptmax=0.0;
	  break;
	}
	for (int i=0;i<(*bit)->NInP();++i) {
	  p_mihandler->ISRHandler()->
	    Extract((*bit)->InParticle(i)->Flav(),
		    (*bit)->InParticle(i)->Momentum()[0],i);
	  p_remnants[i]->Extract((*bit)->InParticle(i));
	}
      }
    }
  }
  return m_ptmax!=std::numeric_limits<double>::max();
}

bool Multiple_Interactions::Treat(ATOOLS::Blob_List *bloblist,double &weight)
{
  PROFILE_HERE;
  if (p_mihandler->Type()==MI_Handler::None ||
      AMISIC::MI_Base::StopGeneration()) return false;
  if (bloblist->empty()) {
    ATOOLS::msg.Error()<<"Multiple_Interactions::Treat(): "
		       <<"Incoming blob list is empty!"<<std::endl;
    return false;
  }
  if (!CheckBlobList(bloblist)) return false;
  p_mihandler->SetScaleMax(m_ptmax,0);
  p_mihandler->SetScaleMin(p_mihandler->ScaleMin(0),0);
  p_mihandler->SetScaleMax(m_xmax[0],2);
  p_mihandler->SetScaleMax(m_xmax[1],3);
  ATOOLS::Blob *blob = new ATOOLS::Blob();
  p_mihandler->Reset();
  bool success=false;
  if (!m_vetoed && m_ptmax>p_mihandler->ScaleMin(0)) {
    success=p_mihandler->GenerateHardProcess(blob);
  }
  else if (m_vetoed) {
    success=p_mihandler->GenerateSoftProcess(blob);
    // dummy settings for unweighted analysis
    blob->SetType(ATOOLS::btp::Signal_Process);
    blob->SetTypeSpec("Soft UE");
    blob->SetStatus(5);
    blob->AddData("ME_Weight",new ATOOLS::Blob_Data<double>(1.0));
    blob->AddData("ME_NumberOfTrials",new ATOOLS::Blob_Data<int>(1));
  }
  if (success) {
    blob->SetId(bloblist->size());
    m_ptmax=blob->OutParticle(0)->Momentum().PPerp();
    for (size_t i=0;i<(size_t)blob->NInP();++i) {
      p_mihandler->ISRHandler()->Reset(i);
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

void Multiple_Interactions::DeleteConnected(ATOOLS::Blob *const blob)
{
  if (m_deleted.find(blob)!=m_deleted.end()) return;
  m_deleted.insert(blob);
  for (int i=blob->NOutP()-1;i>=0;i=ATOOLS::Min(blob->NOutP()-1,i-1)) {
    ATOOLS::Blob *dblob=blob->OutParticle(i)->DecayBlob();
    if (dblob!=NULL) DeleteConnected(dblob);
  }
  for (int i=blob->NInP()-1;i>=0;i=ATOOLS::Min(blob->NInP()-1,i-1)) {
    ATOOLS::Blob *pblob=blob->InParticle(i)->ProductionBlob();
    if (pblob!=NULL) DeleteConnected(pblob);
  }
  delete blob;
  for (ATOOLS::Blob_List::iterator bit=p_bloblist->begin();
       bit!=p_bloblist->end();++bit)
    if (*bit==blob) {
      p_bloblist->erase(bit);
      break;
    }
}

bool Multiple_Interactions::VetoHardProcess(ATOOLS::Blob *const blob)
{
  if (p_mihandler->VetoHardProcess(blob)) {
    m_deleted.clear();
    DeleteConnected(blob);
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

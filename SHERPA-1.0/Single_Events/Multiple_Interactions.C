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

bool Multiple_Interactions::CheckBlobList(const ATOOLS::Blob_List *bloblist) 
{
  double pperpmax=m_pperpmax=std::numeric_limits<double>::max();
  for (size_t i=0;i<2;++i) {
    m_xmax[i]=1.0;
    p_remnants[i]->QuickClear();
  }
  for (ATOOLS::Blob_List::const_iterator bit=bloblist->begin();
       bit!=bloblist->end();++bit) {
    if ((*bit)->Type()==ATOOLS::btp::Hard_Collision || 
	(*bit)->Type()==ATOOLS::btp::Signal_Process) {
      if ((*bit)->Status()!=0) return false;
    }
    else if ((*bit)->Type()==ATOOLS::btp::ME_PS_Interface_FS) {
      for (int i=0;i<(*bit)->NInP();++i) {
	pperpmax=ATOOLS::Min(pperpmax,
			     (*bit)->InParticle(i)->Momentum().PPerp());
      }
    }
    else if ((*bit)->Type()==ATOOLS::btp::IS_Shower) {
      m_xmax[(*bit)->Beam()]-=2.0*(*bit)->InParticle(0)->Momentum()[0]/m_ecms;
      p_mihandler->ISRHandler()->
	Extract((*bit)->InParticle(0)->Flav(),
		(*bit)->InParticle(0)->Momentum()[0],(*bit)->Beam());
      if (!p_remnants[(*bit)->Beam()]->Extract((*bit)->InParticle(0))) 
	VetoHardProcess(*bit);
    }
  }
  if (pperpmax>=m_pperpmax) {
    for (ATOOLS::Blob_List::const_iterator bit=bloblist->begin();
	 bit!=bloblist->end();++bit) {
      if ((*bit)->Type()==ATOOLS::btp::Hard_Collision || 
	  (*bit)->Type()==ATOOLS::btp::Signal_Process) {
	for (int i=0;i<(*bit)->NOutP();++i) {
	  pperpmax=ATOOLS::Min(pperpmax,
			       (*bit)->OutParticle(i)->Momentum().PPerp());
	}
	for (int i=0;i<(*bit)->NInP();++i) {
	  p_mihandler->ISRHandler()->
	    Extract((*bit)->InParticle(i)->Flav(),
		    (*bit)->InParticle(i)->Momentum()[0],i);
	  if (!p_remnants[i]->Extract((*bit)->InParticle(i))) VetoHardProcess(*bit);
	}
      }
    }
  }
  return (m_pperpmax=pperpmax)!=std::numeric_limits<double>::max();
}

bool Multiple_Interactions::Treat(ATOOLS::Blob_List *bloblist,double &weight)
{
  PROFILE_HERE;
  if (p_mihandler->Type()==MI_Handler::None) return false;
  if (bloblist->empty()) {
    ATOOLS::msg.Error()<<"Multiple_Interactions::Treat(): "
		       <<"Incoming blob list is empty!"<<std::endl;
    return false;
  }
  if (!CheckBlobList(bloblist)) return false;
  p_mihandler->SetScaleMax(m_pperpmax,0);
  p_mihandler->SetScaleMax(m_xmax[0],2);
  p_mihandler->SetScaleMax(m_xmax[1],3);
  ATOOLS::Blob *blob = new ATOOLS::Blob();
  blob->SetType(ATOOLS::btp::Hard_Collision);
  p_mihandler->Reset();
  if (p_mihandler->GenerateHardProcess(blob)) {
    blob->SetId(bloblist->size());
    blob->SetStatus(1);
    for (size_t i=0;i<(size_t)blob->NInP();++i) {
      p_mihandler->ISRHandler()->Reset(i);
      if (!p_remnants[i]->Extract(blob->InParticle(i))) {
	AMISIC::MI_Base::SetStopGeneration();
	delete blob;
	return false;
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

void Multiple_Interactions::VetoHardProcess(ATOOLS::Blob *const blob)
{
  ATOOLS::msg.Error()<<"Veto "<<*blob<<std::endl;
}

void Multiple_Interactions::Finish(const std::string &resultpath) 
{
}

void Multiple_Interactions::CleanUp() 
{
  p_mihandler->CleanUp();
}

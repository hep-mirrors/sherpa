#include "Multiple_Interactions.H"

#include "Amisic.H"
#include "My_Limits.H"

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
}

Multiple_Interactions::~Multiple_Interactions() {}


bool Multiple_Interactions::CheckBlobList(const ATOOLS::Blob_List *bloblist) 
{
  m_xmax[1]=m_xmax[0]=1.0;
  m_pperpmax=std::numeric_limits<double>::max();
  for (Blob_List::const_iterator bit=bloblist->begin();bit!=bloblist->end();++bit) {
    if ((*bit)->Type()==btp::Hard_Collision || 
	(*bit)->Type()==btp::Signal_Process) {
      if ((*bit)->Status()!=0) return false;
    }
    else if ((*bit)->Type()==btp::ME_PS_Interface_FS) {
      for (int i=0;i<(*bit)->NInP();++i) {
	m_pperpmax=ATOOLS::Min(m_pperpmax,(*bit)->InParticle(i)->Momentum().PPerp());
      }
    }
    else if ((*bit)->Type()==btp::IS_Shower) {
      m_xmax[(*bit)->Beam()]-=2.0*(*bit)->InParticle(0)->Momentum()[0]/m_ecms;
//       p_mihandler->ISRHandler()->Extract((*bit)->InParticle(0)->Flav(),
// 					 (*bit)->InParticle(0)->Momentum()[0],(*bit)->Beam());
    }
  }
  if (m_pperpmax==std::numeric_limits<double>::max()) {
    for (Blob_List::const_iterator bit=bloblist->begin();bit!=bloblist->end();++bit) {
      if ((*bit)->Type()==btp::Hard_Collision || 
	  (*bit)->Type()==btp::Signal_Process) {
	for (int i=0;i<(*bit)->NOutP();++i) {
	  m_pperpmax=ATOOLS::Min(m_pperpmax,(*bit)->OutParticle(i)->Momentum().PPerp());
	}
      }
    }
  }
  return m_pperpmax!=std::numeric_limits<double>::max();
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
    bloblist->push_back(blob);
    blob->SetStatus(1);
//     p_mihandler->ISRHandler()->Reset(0);
//     p_mihandler->ISRHandler()->Reset(1);
    return true;
  }
  delete blob;
//   p_mihandler->ISRHandler()->Reset(0);
//   p_mihandler->ISRHandler()->Reset(1);
  return false;
}

void Multiple_Interactions::Finish(const std::string &resultpath) 
{
}

void Multiple_Interactions::CleanUp() 
{
  p_mihandler->CleanUp();
}

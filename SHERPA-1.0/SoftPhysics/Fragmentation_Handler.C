#include "Fragmentation_Handler.H"

#include "Data_Read.H"
#include "Run_Parameter.H"
#include "Exception.H"

#ifdef PROFILE__all
#define PROFILE__Fragmentation_Handler
#endif
#ifdef PROFILE__Fragmentation_Handler
#include "prof.hh" 
#else
#define PROFILE_HERE
#endif

using namespace SHERPA;

Fragmentation_Handler::Fragmentation_Handler(std::string _dir,std::string _file):
  m_dir(_dir), 
  m_file(_file),
  m_mode(0),
  m_maxtrials(1),
  p_lund(NULL)
{
  ATOOLS::Data_Read dr(m_dir+m_file);
  m_fragmentationmodel=dr.GetValue<std::string>("FRAGMENTATION",std::string("Pythiav6.214"));
  if (m_fragmentationmodel==std::string("Lund")) {
    std::string lundfile=dr.GetValue<std::string>("LUND_FILE",std::string("Lund.dat"));
    p_lund = new Lund_Interface(m_dir,lundfile,true);
    m_mode=1;
    return;
  }
  if (m_fragmentationmodel==std::string("Off")) {
    return;
  }
  throw(ATOOLS::Exception(ATOOLS::ex::critical_error,"Fragmentation model not implemented.",
			  "Fragmentation_Handler","Fragmentation_Handler"));
}
   
Fragmentation_Handler::~Fragmentation_Handler() 
{
  if (p_lund!=NULL) delete p_lund;
}

bool Fragmentation_Handler::PerformFragmentation(ATOOLS::Blob_List *bloblist,
						 ATOOLS::Particle_List *particlelist) 
{
  PROFILE_HERE;
  if (m_mode==0 || bloblist->size()==0) return 1;
  p_blob = new ATOOLS::Blob();
  p_blob->SetId();
  p_blob->SetType(ATOOLS::btp::Fragmentation);
  bloblist->push_back(p_blob);
  std::vector<ATOOLS::Particle*> startpoints;
  for (ATOOLS::Blob_Iterator bit=bloblist->begin();bit!=bloblist->end();++bit) {
    for (size_t i=0;i<(size_t)(*bit)->NOutP();++i) {
      ATOOLS::Particle *cur=(*bit)->OutParticle(i);
      if (cur->DecayBlob()==NULL && cur->Status()==1 && 
	  (cur->Info()=='F' || cur->Info()=='H') &&
	  cur->GetFlow(1)!=0 && cur->GetFlow(2)==0) {
	startpoints.push_back(cur);
      }
    }
  }
  m_used.clear();
  for (std::vector<ATOOLS::Particle*>::iterator sit=startpoints.begin();
       sit!=startpoints.end();++sit) {
    ATOOLS::Particle *cur=*sit, *comp;  
    p_blob->AddToInParticles(cur);
    m_used.insert(cur);
    do {
      bool found=false;
      for (ATOOLS::Blob_Iterator bit=bloblist->begin();bit!=bloblist->end();++bit) {
	for (size_t i=0;i<(size_t)(*bit)->NOutP();++i) {
	  comp=(*bit)->OutParticle(i);
	  bool test=false;
	  if (comp->DecayBlob()==NULL) test=true;
	  else if (comp->DecayBlob()->Type()==ATOOLS::btp::Fragmentation) test=true;
	  if (test && comp->Status()==1 && comp->GetFlow(2)==cur->GetFlow(1) &&
	      (comp->Info()=='F' || comp->Info()=='H')) {
	    if (m_used.find(comp)==m_used.end()) {
	      p_blob->AddToInParticles(comp);
	      m_used.insert(comp);
	      found=true;
	      break;
	    }
	  }
	}
	if (found) break;
      }
      if (!found) {
	ATOOLS::msg.Error()<<"Fragmentation_Handler::PerformFragmentation(..): "
			   <<"Cannot find connected parton for parton ("
			   <<cur->Number()<<") in event ["
			   <<ATOOLS::rpa.gen.NumberOfDicedEvents()<<"]"<<std::endl;
	msg_Tracking()<<"   Empty blob list and retry event. {\n"<<*bloblist<<"   }"<<std::endl;
	while (bloblist->size()>0) {
	  delete *bloblist->begin();
	  bloblist->erase(bloblist->begin());
	}
	return false;
      }
    } while ((cur=comp)->Flav().IsGluon());
  }
  for (size_t trials=0;trials<m_maxtrials;++trials) {
    if (p_lund->Hadronize(p_blob,bloblist,particlelist)) return true;
    if (m_maxtrials>1) ATOOLS::msg.Error()<<"Fragmentation_Handler::PerformFragmentation(..): "
					  <<"Hadronization failed. Retry."<<std::endl;
  }
  ATOOLS::msg.Error()<<"Fragmentation_Handler::PerformFragmentation(..): "
		     <<"Hadronization failed."<<std::endl;
  while (bloblist->size()>0) {
    delete *bloblist->begin();
    bloblist->erase(bloblist->begin());
  }
  return false;
}


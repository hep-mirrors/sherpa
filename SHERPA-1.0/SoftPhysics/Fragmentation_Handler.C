#include "Fragmentation_Handler.H"

#include "Data_Read.H"
#include "Exception.H"

using namespace SHERPA;

Fragmentation_Handler::Fragmentation_Handler(std::string _dir,std::string _file):
  m_dir(_dir), 
  m_file(_file),
  m_mode(0),
  p_lund(NULL)
{
  ATOOLS::Data_Read dr(m_dir+m_file);
  m_fragmentationmodel=dr.GetValue<std::string>("FRAGMENTATION",std::string("Lund"));
  if (m_fragmentationmodel==std::string("Lund")) {
    std::string lundfile=dr.GetValue<std::string>("LUND_FILE",std::string("Lund.dat"));
    ATOOLS::msg.Events()<<"Fragmentation_Handler::Fragmentation_Handler(..): "
			<<"Initialize Lund Fragmentation according to "<<lundfile<<std::endl;
    p_lund = new Lund_Interface(m_dir,lundfile);
    m_mode=1;
    return;
  }
  if (m_fragmentationmodel==std::string("Off")) {
    ATOOLS::msg.Out()<<"Fragmentation_Handler::Fragmentation_Handler(..): "
		     <<"WARNING: The Fragmentation is switched off "<<std::endl;
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
  if (m_mode==0) return 1;
  p_blob = new ATOOLS::Blob();
  bloblist->push_back(p_blob);
  p_blob->SetId(bloblist->size());
  p_blob->SetType(ATOOLS::btp::Fragmentation);
  std::set<ATOOLS::Particle*> startpoints;
  for (ATOOLS::Blob_Iterator bit=bloblist->begin();bit!=bloblist->end();++bit) {
    for (size_t i=0;i<(size_t)(*bit)->NOutP();++i) {
      ATOOLS::Particle *cur=(*bit)->OutParticle(i);
      if (cur->DecayBlob()==NULL && cur->Status()==1 && 
	  (cur->Info()=='F' || cur->Info()=='H') &&
	  cur->GetFlow(1)!=0 && cur->GetFlow(2)==0) {
	startpoints.insert(cur);
      }
    }
  }
  for (std::set<ATOOLS::Particle*>::iterator sit=startpoints.begin();
       sit!=startpoints.end();++sit) {
    ATOOLS::Particle *cur=*sit, *comp;  
    p_blob->AddToInParticles(cur);
    m_used.clear();
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
			   <<cur->Number()<<")"<<std::endl;
	return false;
      }
    } while ((cur=comp)->Flav().IsGluon());
  }
  return p_lund->Hadronize(p_blob,bloblist,particlelist);
}


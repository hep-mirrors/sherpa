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
  for (ATOOLS::Blob_Iterator bit=bloblist->begin();
       bit!=bloblist->end();++bit) {
    for (size_t i=0;i<(size_t)(*bit)->NOutP();++i) {
      ATOOLS::Particle *cur=(*bit)->OutParticle(i);
      if (cur->DecayBlob()==NULL && cur->Status()==1 && 
	  (cur->Info()=='F' || cur->Info()=='H')) {
	p_blob->AddToInParticles(cur);
      }
    }
  }
  for (size_t i=0;i<(size_t)p_blob->NInP();++i) {
    ATOOLS::Particle *comp, *cur=p_blob->InParticle(i);  
    bool anti=cur->Flav().IsAnti();
    for (size_t step=0;step<2;++step) {
      comp=cur;
      if (step==anti || cur->Flav().IsGluon()) {
	if (FindConnected(comp,cur->GetFlow(step+1),anti,true,1)) {
	  bool add=true;
	  for (size_t j=0;j<(size_t)p_blob->NInP();++j) if (comp==p_blob->InParticle(j)) add=false;	
	  if (add) p_blob->AddToInParticles(comp);
	}
	else {
	  ATOOLS::msg.Error()<<"Fragmentation_Handler::PerformFragmentation(..): "
			     <<"Cannot find connected parton for parton ("<<cur->Number()<<")"<<std::endl;
	  return false;
	}
      }
    }
  }
  return p_lund->Hadronize(p_blob,bloblist,particlelist);
}

bool Fragmentation_Handler::FindConnected(ATOOLS::Particle *particle,unsigned int color,
					  bool anti,bool forward,unsigned int catcher)
{
  if (++catcher>100) {
    ATOOLS::msg.Tracking()<<"Fragmentation_Handler::FindConnected(..): "
			  <<"Colour nesting is too deep."<<std::endl;
    return false;
  }
  ATOOLS::Blob *cur;
  if (forward) cur=particle->DecayBlob();
  else cur=particle->ProductionBlob();
  if (cur!=NULL) {
    if (cur->Type()==ATOOLS::btp::Fragmentation) return true;
    for (int i=0;i<cur->NOutP();++i) {
      for (int j=1;j<3;++j) {
	ATOOLS::Particle *help=cur->OutParticle(i);
	if (help->GetFlow(j)==(int)color && 
	    (help->Flav().IsAnti()!=anti || help->Flav().IsGluon())) {
	  return FindConnected(help,color,anti,true,catcher);
	}
      }
    }
    for (int i=0;i<cur->NInP();++i) {
      for (int j=1;j<3;++j) {
	ATOOLS::Particle *help=cur->InParticle(i);
	if (help->GetFlow(j)==(int)color && 
	    (help->Flav().IsAnti()!=anti || help->Flav().IsGluon())) {
	  return FindConnected(help,color,anti,true,catcher);
	}
      }
    }
  }
  else {
    return true;
  }
  return false;
}



#include "Fragmentation_Handler.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "Random.H"
#include "Vector.H"
#include "Data_Read.H"

#include <stdio.h>

using namespace SHERPA;
using namespace ATOOLS;


Fragmentation_Handler::Fragmentation_Handler(std::string _dir,std::string _file) :
  m_dir(_dir), m_file(_file)
{
  Data_Read dr(m_dir+m_file);
  m_fragmentationmodel = dr.GetValue<std::string>("FRAGMENTATION",std::string("Lund"));
  
  if (m_fragmentationmodel==std::string("Lund")) {
    m_lund_a     = dr.GetValue<double>("LUND_A",0.4);
    m_lund_b     = dr.GetValue<double>("LUND_B",0.85);
    m_lund_sigma = dr.GetValue<double>("LUND_SIGMA",0.36);
    msg.Events()<<"Initialize Lund Fragmentation : "<<std::endl
		<<"  LUND_A = "<<m_lund_a<<", LUND_B = "<<m_lund_b
		<<", LUND_SIGMA = "<<m_lund_sigma<<std::endl;
    p_lund       = new Lund_Fortran_Interface(m_lund_a,m_lund_b,m_lund_sigma);
    m_mode       = 1;
    return;
  }
  if (m_fragmentationmodel==std::string("Off") ||
      m_fragmentationmodel==std::string("off") ) {
    msg.Out()<<"WARNING : The Fragmentation is switched off "<<std::endl;
    p_lund       = 0;
    m_mode       = 0;
    return;
  }

  msg.Error()<<"ERROR in Fragmentation_Handler::Fragmentation_Handler."<<std::endl
	     <<"    please choose between <Lund> and <Off> as Fragmentation model"<<std::endl
	     <<"    the Fragmentation model <"<<m_fragmentationmodel
	     <<"> is not implemented yet. Abort."<<std::endl;
  abort();
}
   
Fragmentation_Handler::~Fragmentation_Handler() {
  if (p_lund)      delete p_lund;
}




bool Fragmentation_Handler::PerformFragmentation(ATOOLS::Blob_List * bl,
						 ATOOLS::Particle_List * pl) 
{
  if (m_mode==0) return 1;

  if (!ExtractSinglets(bl,pl)) return 0;

  bool okay = 1;
  for (Blob_Iterator biter=bl->begin();biter!=bl->end();++biter) {
    if ( ((*biter)->Type()==std::string("Fragmentation")) && 
	 ((*biter)->Status()==1) ) {
      //(*biter)->BoostInCMS();
      (*biter)->SetCMS();
       okay = okay && p_lund->Hadronize((*biter),bl,pl);
      //(*biter)->BoostInLab();
      (*biter)->SetStatus(0);
    }
  }
  return okay;
}

bool Fragmentation_Handler::ExtractSinglets(Blob_List * _bloblist,Particle_List * pl) 
{
  Blob       * newb = NULL;
  Particle     * part;
  bool use_one_blob = 1;

  bool foundatall   = 0;
  bool found        = 1;
  bool active;
  while (found) {
    found = 0;
    for (Blob_Iterator blit=_bloblist->begin();blit!=_bloblist->end();++blit) {
      if (((*blit)->Status()==1) || ((*blit)->Status()==2)) {
	for (int i=0;i<(*blit)->NOutP();i++) {
	  part = (*blit)->OutParticle(i);
	  if ( (part->Info()=='F' || part->Info()=='H') && part->Status()==1) {
	    if (( (part->Flav().IsQuark()   && !part->Flav().IsAnti() ) ||
		  (part->Flav().IsDiQuark() && part->Flav().IsAnti())) && part->GetFlow(1)>0 ) {
	      if (use_one_blob==0 || newb==0) {
		newb = new Blob();
		newb->SetId(_bloblist->size());
		newb->SetStatus(1);
		newb->SetType(std::string("Fragmentation"));
		_bloblist->push_back(newb);
	      }
	      part->SetStatus(2);
	      newb->AddToInParticles(part);
	      foundatall = found = 1;
	      if (!(FindConnected(_bloblist,part,newb))) {
		msg.Error()<<"Fragmentation_Handler::ExtractSinglets :"
			   <<"Could not find connected parton for quark:"<<std::endl
			   <<part<<std::endl;
		return 0;
	      }
	    }
	  }
	}      
      }
    }
    for (Blob_Iterator blit=_bloblist->begin();blit!=_bloblist->end();++blit) {
      if (!((*blit)->Type()==std::string("Fragmentation"))) {
	active = 0;
	for (int i=0;i<(*blit)->NOutP();i++) {
	  if ((*blit)->OutParticle(i)->Status()==1) { 
	    active = 1; 
	    break; 
	  }
	}
	if (!active) (*blit)->SetStatus(0);
      }
    }
  }
  for (Blob_Iterator blit=_bloblist->begin();blit!=_bloblist->end();++blit) {
    if (!((*blit)->Type()==std::string("Fragmentation"))) {
      active = 0;
      for (int i=0;i<(*blit)->NOutP();i++) {
	if ((*blit)->OutParticle(i)->Status()==1) { 
	  active = 1; 
	  break; 
	}
      }
      if (!active) (*blit)->SetStatus(0);
    }
  }
  return foundatall;
}


bool Fragmentation_Handler::FindConnected(Blob_List * _bloblist,
					  Particle * compare,Blob * blob) {
  Particle * part;
  for (Blob_Iterator blit=_bloblist->begin();blit!=_bloblist->end();++blit) {
    for (int i=0;i<(*blit)->NOutP();i++) {
      part = (*blit)->OutParticle(i);
      if (part==compare || part->Status()!=1 || part->DecayBlob()!=NULL) continue;
      if (part->Info()=='F' || part->Info() == 'H') {
	if (part->GetFlow(2)==compare->GetFlow(1)) {
	  part->SetStatus(2);
	  blob->AddToInParticles(part);
	  if (part->GetFlow(1)==0) {
	    if ( part->Flav().IsQuark() && part->Flav().IsAnti()) {
	      return 1;
	    }
	    if ( part->Flav().IsDiQuark() && (!part->Flav().IsAnti()) ) {
	      return 1;
	    }
	  }
	  else { return FindConnected(_bloblist,part,blob); }
	}
      }
    }
  }
  msg.Error()<<"Fragmentation_Handler::FindConnected : No closed singlet line !"<<std::endl;
  return 0;
}

Lund_Fortran_Interface * Fragmentation_Handler::GetLundFortranInterface() 
{ 
  if (p_lund) return p_lund; 
  msg.Out()<<"WARNING: in Fragmentation_Handler::GetLundFortranInterface()."<<std::endl
	   <<"   Not yet initialized. This is an inconsistent option at the moment."<<std::endl;
    //	   <<"   Abort program. "<<std::endl;
    //  abort();
  return p_lund;
}


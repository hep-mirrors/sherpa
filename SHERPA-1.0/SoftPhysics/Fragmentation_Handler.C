#include "Fragmentation_Handler.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "Random.H"
#include "Vector.H"
#include "Data_Read.H"

using namespace SHERPA;
using namespace AMATOOLS;
using namespace APHYTOOLS;
using namespace AORGTOOLS;


Fragmentation_Handler::Fragmentation_Handler(string _dir,string _file) :
  m_dir(_dir), m_file(_file)
{
  Data_Read dr(m_dir+m_file);
  m_fragmentationmodel = dr.GetValue<string>("FRAGMENTATION",string("Lund"));
  
  if (m_fragmentationmodel==string("Lund")) {
    m_lund_a     = dr.GetValue<double>("LUND_A",0.4);
    m_lund_b     = dr.GetValue<double>("LUND_B",0.85);
    m_lund_sigma = dr.GetValue<double>("LUND_SIGMA",0.36);
    msg.Out()<<"Initialize Lund Fragmentation : "<<endl
	     <<"  LUND_A = "<<m_lund_a<<", LUND_B = "<<m_lund_b<<", LUND_SIGMA = "<<m_lund_sigma<<endl;
    p_lund       = new Lund_Fortran_Interface(m_lund_a,m_lund_b,m_lund_sigma);
    m_mode       = 1;
    return;
  }

  msg.Error()<<"Error in Fragmentation_Handler::Fragmentation_Handler."<<endl
	     <<"    Fragmentation model "<<m_fragmentationmodel<<" not implemented yet. Abort."<<endl;
  abort;
}
   
Fragmentation_Handler::~Fragmentation_Handler() {
  msg.Tracking()<<"+++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  msg.Tracking()<<"In Fragmentation_Handler::~Fragmentation_Handler :"<<std::endl;
  if (p_lund)      delete p_lund;
  
  msg.Tracking()<<"Out Fragmentation_Handler::~Fragmentation_Handler :"<<std::endl;
  msg.Tracking()<<"+++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
}; 




bool Fragmentation_Handler::PerformFragmentation(APHYTOOLS::Blob_List * bl,
						 APHYTOOLS::Parton_List * pl) 
{
  if (!ExtractSinglets(bl,pl)) return 0;
  bool okay = 1;
  for (Blob_Iterator biter=bl->begin();biter!=bl->end();++biter) {
    msg.Debugging()<<"Take "<<(*biter)->Type()<<" "<<(*biter)->Status()<<endl;
    if ( ((*biter)->Type()==std::string("Fragmentation")) && 
	 ((*biter)->Status()==1) ) {
      (*biter)->BoostInCMS();
      okay = okay && p_lund->Hadronize((*biter),bl,pl);
      (*biter)->BoostInLab();
      (*biter)->SetStatus(0);
    }
  }
  return okay;
}

bool Fragmentation_Handler::ExtractSinglets(Blob_List * _bloblist,Parton_List * pl) 
{
  msg.Debugging()<<"Fragmentation_Handler::ExtractSinglets :"<<std::endl<<std::endl;
  Blob       * newb = NULL;
  Parton     * part;
  bool use_one_blob = 1;

  bool foundatall   = 0;
  bool found        = 1;
  bool active;
  while (found) {
    found = 0;
    for (Blob_Iterator blit=_bloblist->begin();blit!=_bloblist->end();++blit) {
      if ((*blit)->Status()==1) {
	for (int i=0;i<(*blit)->NOutP();i++) {
	  part = (*blit)->OutParton(i);
	  if ( (part->Info()=='F' || part->Info()=='H') && part->Status()==1) {
	    if (( (part->Flav().IsQuark()    && !part->Flav().IsAnti() ) ||
		  (part->Flav().IsDiQuark()) && part->Flav().IsAnti()       ) && part->GetFlow(1)>0 ) {
	      msg.Debugging()<<"   Start new singlet chain for parton : "<<part<<std::endl;
	      if (use_one_blob==0 || newb==0) {
		newb = new Blob();
		newb->SetId(_bloblist->size());
		newb->SetStatus(1);
		newb->SetType(std::string("Fragmentation"));
		_bloblist->push_back(newb);
	      }
	      part->SetStatus(2);
	      part->SetDecayBlob(newb);
	      newb->AddToInPartons(part);
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
      if ((*blit)->Type()==string("Fragmentation")) continue;
      active = 0;
      for (int i=0;i<(*blit)->NOutP();i++) {
	if ((*blit)->OutParton(i)->Status()==1) { active = 1; break; }
      }
      if (!active) (*blit)->SetStatus(0);
    }
  }
  for (Blob_Iterator blit=_bloblist->begin();blit!=_bloblist->end();++blit) {
    if ((*blit)->Type()==string("Fragmentation")) continue;
    active = 0;
    for (int i=0;i<(*blit)->NOutP();i++) {
      if ((*blit)->OutParton(i)->Status()==1) { active = 1; break; }
    }
    if (!active) (*blit)->SetStatus(0);
  }
  msg.Debugging()<<"Out Extract Singlets "<<foundatall<<endl;
  return foundatall;
}


bool Fragmentation_Handler::FindConnected(Blob_List * _bloblist,
					  Parton * compare,Blob * blob) {
  msg.Debugging()<<"Fragmentation_Handler::FindConnected :"<<std::endl;
  Parton * part;
  for (Blob_Iterator blit=_bloblist->begin();blit!=_bloblist->end();++blit) {
    if ((*blit)->Status()==1) {
      for (int i=0;i<(*blit)->NOutP();i++) {
	part = (*blit)->OutParton(i);
	if (part==compare) continue;
	if ( part->Info()!='F' && part->Info() != 'H') continue;
	if (part->Status()!=1) continue; 
	if (part->GetFlow(2)==compare->GetFlow(1)) {
	  msg.Debugging()<<"Colour match for "<<compare->GetFlow(1)<<std::endl;
	  part->SetStatus(2);
	  part->SetDecayBlob(blob);
	  blob->AddToInPartons(part);
	  if (part->GetFlow(1)==0) {
	    if ( part->Flav().IsQuark() && part->Flav().IsAnti()) {
	      msg.Debugging()<<"Closed singlet list with an antiquark."<<std::endl;
	      return 1;
	    }
	    if ( part->Flav().IsDiQuark() && (!part->Flav().IsAnti()) ) {
	      msg.Debugging()<<"Closed singlet list with a diquark."<<std::endl;
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
  msg.Error()<<"Error in Fragmentation_Handler::GetLundFortranInterface()."<<endl
	     <<"   Not yet initialized. This is an inconsistent option at the moment."<<endl
	     <<"   Abort program. "<<endl;
  abort();
}


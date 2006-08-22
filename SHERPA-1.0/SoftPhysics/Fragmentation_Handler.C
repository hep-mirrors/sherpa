#include "Fragmentation_Handler.H"

#include "Data_Read.H"
#include "Run_Parameter.H"
#include "Exception.H"
#include "Return_Value.H"

#ifdef PROFILE__all
#define PROFILE__Fragmentation_Handler
#endif
#ifdef PROFILE__Fragmentation_Handler
#include "prof.hh" 
#else
#define PROFILE_HERE
#endif

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;
#ifdef USING__Ahadic
using namespace AHADIC;
#endif

Fragmentation_Handler::Fragmentation_Handler(string _dir,string _file):
  m_dir(_dir), m_file(_file), m_mode(0), 
#ifdef USING__Ahadic
  p_ahadic(NULL),
#endif
  p_lund(NULL)
{
  Data_Read dr(m_dir+m_file);
  m_fragmentationmodel=dr.GetValue<string>("FRAGMENTATION",string("Pythiav6.214"));
  if (m_fragmentationmodel==string("Lund")) {
    string lundfile=dr.GetValue<string>("LUND_FILE",string("Lund.dat"));
    p_lund = new Lund_Interface(m_dir,lundfile,true);
    m_mode=1;
    return;
  }
#ifdef USING__Ahadic
  else if (m_fragmentationmodel==string("Ahadic")) {
    string clusterfile=dr.GetValue<string>("AHADIC_FILE",string("Cluster.dat"));
    p_ahadic = new AHADIC::Ahadic(m_dir,clusterfile,true);
    m_mode=2;
    return;
  }
#endif
  else if (m_fragmentationmodel==string("Off")) return;
  THROW(critical_error,"Fragmentation model not implemented.");
}
   
Fragmentation_Handler::~Fragmentation_Handler() 
{
  if (p_lund!=NULL)   { delete p_lund;   p_lund   = NULL;   }
#ifdef USING__Ahadic
  if (p_ahadic!=NULL) { delete p_ahadic; p_ahadic = NULL;   }
#endif
}

Return_Value::code Fragmentation_Handler::PerformFragmentation(Blob_List *bloblist,
							       Particle_List *particlelist) 
{
  PROFILE_HERE;
  if (m_mode==0 || bloblist->size()==0) return Return_Value::Nothing;
  switch (int(ExtractSinglets(bloblist))) {
    case int(Return_Value::Success) : break;
    case int(Return_Value::Nothing) : return Return_Value::Nothing;
    case int(Return_Value::Error)   : return Return_Value::Error;
    default :
      msg.Error()<<"ERROR in "<<METHOD<<":"<<std::endl
		 <<"   ExtractSinglets yields unknown return value."<<std::endl
		 <<"   Return 'Retry_Event' and hope for the best."<<std::endl;
      return Return_Value::Retry_Event;
  }
  switch (m_mode) {
    case 1  : return p_lund->Hadronize(bloblist);
#ifdef USING__Ahadic
    case 2  : return p_ahadic->Hadronize(bloblist);
#endif
    default : 
      msg.Error()<<"ERROR in "<<METHOD<<":"<<std::endl
		 <<"   Unknown hadronization model in mode = "<<m_mode<<"."<<std::endl
		 <<"   Abort the run."<<std::endl;
      THROW(critical_error,"Fragmentation model not implemented.");
  }
}

Return_Value::code Fragmentation_Handler::ExtractSinglets(Blob_List * bloblist)
{
  Particle  * part(NULL);
  Blob      * blob(NULL);
  Part_List * plist = new Part_List;
  for (Blob_List::iterator blit=bloblist->begin();blit!=bloblist->end();++blit) {
    if ((*blit)->Has(blob_status::needs_hadronization)) {
      for (int i=0;i<(*blit)->NOutP();i++) {
	part = (*blit)->OutParticle(i); 
	if (part->Status()==part_status::active) {
	  if (part->GetFlow(1)!=0 || part->GetFlow(2)!=0) {
	    plist->push_back(part);
	    part->SetStatus(part_status::fragmented);
	  }
	  else if (part->Flav()==Flavour(kf::tau) ||
		   part->Flav()==Flavour(kf::tau).Bar()) {
	    blob = new Blob();
	    blob->SetId();
	    blob->SetType(btp::Fragmentation);
	    blob->SetStatus(blob_status::needs_hadrondecays);
	    blob->AddToInParticles(part);
	    blob->AddToOutParticles(new Particle((*part)));
	    blob->InParticle(0)->SetStatus(part_status::decayed);
	    blob->OutParticle(0)->SetStatus(part_status::active);
	  }
	}
      }
      (*blit)->UnsetStatus(blob_status::needs_hadronization);
    }
  }
  if (plist->empty()) {
    msg.Debugging()<<"WARNING in Lund_Interface::PrepareFragmentationBlob:"<<endl
		   <<"   No coloured particle found leaving shower blobs."<<endl;
    return Return_Value::Nothing;
  }


  int  col1, col2;
  bool hit1, hit2;
  Part_List * pli(NULL);
  vector<Part_List *> partlists; 
  do {
    hit1 = false;
    for (Part_Iterator pit=plist->begin();pit!=plist->end();++pit) {
      col1 = (*pit)->GetFlow(1);
      col2 = (*pit)->GetFlow(2);
      if (col1!=0 && col2==0) {
	hit1 = true;
	pli  = new Part_List;
	pli->push_back((*pit));
	pit  = plist->erase(pit);
	partlists.push_back(pli);
	do {
	  hit2 = false;
	  for (Part_Iterator pit1=plist->begin();pit1!=plist->end();++pit1) {
	    if ((int)((*pit1)->GetFlow(2))==col1) {
	      col1 = (*pit1)->GetFlow(1);
	      pli->push_back((*pit1));
	      pit1 = plist->erase(pit1);
	      hit2 = true;
	      break;
	    }
	  }
	} while (hit2 && col1!=0);
      }
      if (hit1) break;
    }
    if (!hit1) {
      for (Part_Iterator pit=plist->begin();pit!=plist->end();++pit) {
	col1 = (*pit)->GetFlow(1);
	col2 = (*pit)->GetFlow(2);
	if (col1!=0 && col2!=0) {
	  hit1 = true;
	  pli  = new Part_List;
	  pli->push_back((*pit));
	  pit  = plist->erase(pit);
	  partlists.push_back(pli);
	  do {
	    hit2 = false;
	    for (Part_Iterator pit1=plist->begin();pit1!=plist->end();++pit1) {
	      if ((int)((*pit1)->GetFlow(2))==col1) {
		col1 = (*pit1)->GetFlow(1);
		pli->push_back((*pit1));
		pit1 = plist->erase(pit1);
		hit2 = true;
		break;
	      }
	    }
	  } while (hit2 && col1!=col2);
	}
	if (hit1) break;
      }
    }
  } while(plist->size()>0);

  if (plist->empty()) {
    blob = new Blob();
    blob->SetId();
    blob->SetType(btp::Fragmentation);
    blob->SetStatus(blob_status::needs_hadronization);
    bloblist->push_back(blob);
    for (vector<Part_List *>::iterator pliter=partlists.begin();
	 pliter!=partlists.end();pliter++) {
      while (!(*pliter)->empty()) {
	blob->AddToInParticles((*pliter)->front());
	(*pliter)->pop_front();
      }
    }
    
    return Return_Value::Success;
  }
  return Return_Value::Error;
}


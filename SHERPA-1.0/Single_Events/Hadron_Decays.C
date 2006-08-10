#include "Hadron_Decays.H"
#include "Message.H"
#include "Random.H"
#include <utility>

#ifdef USING__Hadrons
#include "Hadrons.H"
#endif
#include "Spin_Density_Matrix.H"
#include "Spin_Correlation_Tensor.H"

#ifdef PROFILE__Hadron_Decays
#include "prof.hh"
#else 
#define PROFILE_HERE {}
#define PROFILE_LOCAL(LOCALNAME) {}
#endif
#include <algorithm>

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Hadron_Decays::Hadron_Decays(HDHandlersMap * _dechandlers) :
  p_dechandlers(_dechandlers)
{
  m_name      = std::string("Hadron_Decays");
  m_type      = eph::Hadronization;
}

Hadron_Decays::~Hadron_Decays() 
{
}

Return_Value::code Hadron_Decays::Treat(ATOOLS::Blob_List * bloblist, double & weight) 
{
  PROFILE_HERE;
  if(p_dechandlers->empty()) return Return_Value::Nothing;

  if (bloblist->empty()) {
    msg.Error()<<"Potential error in Hadron_Decays::Treat."<<endl
      <<"   Incoming blob list contains "<<bloblist->size()<<" entries."<<endl
      <<"   Continue and hope for the best."<<endl;
    return Return_Value::Error;
  }

  Hadron_Decay_Handler * hdhandler;		// pointer on considered HD Handler
  Particle * part;
  bool found(true), didit(false);
  while (found) {
    found = false;
    for (Blob_List::iterator blit=bloblist->begin();
	 blit!=bloblist->end();++blit) {
      //std::cout<<METHOD<<" : check = "<<(*blit)->Type()<<"/"<<(*blit)->Status()<<std::endl;
      if ((*blit)->Has(blob_status::needs_hadrondecays)) {
	//std::cout<<"   Iterate over "<<(*blit)->NOutP()<<" particles."<<std::endl;
	for (int i=0;i<(*blit)->NOutP();i++) {
	  part = (*blit)->OutParticle(i);
	  if (part->Status()==part_status::active && 
	      !part->Flav().IsStable() && part->DecayBlob()==NULL) {
	    //std::cout<<METHOD<<" : check for "<<part->Flav()<<std::endl;
	    hdhandler = NULL;
	    for (HDHandlersIter hd=p_dechandlers->begin();hd!=p_dechandlers->end();hd++) {
	      if (hd->second->CanDealWith(part->Flav().Kfcode())) { hdhandler = hd->second; break; }
	    }
	    /*
	      To do: Return_Value for hdhandler->FillHadronDecayBlob for the BW smearing.
	             Distance measure for decay blob: add it to actual position.
	    */
	    if (hdhandler) {
	      //std::cout<<METHOD<<" : "<<hdhandler->Name()<<" will decay "<<part->Flav()<<std::endl;
	      if (hdhandler->FillHadronDecayBlob(part,bloblist)) {
		found=true;
		//std::cout<<(*bloblist->back())<<std::endl;
	      }
	      else {
		msg.Error()<<"ERROR in "<<METHOD<<":"<<std::endl
			   <<"   Hadron_Decay_Handler "<<hdhandler->Name()<<" failed to decay "<<std::endl
			   <<part<<","<<std::endl
			   <<"   Will retry event."<<std::endl;
		return Return_Value::Retry_Event;
	      }
	    }
	    else {
	      msg.Error()<<"ERROR in "<<METHOD<<":"<<std::endl
			 <<"   Unstable particle found ("<<part->Flav()<<"), "
			 <<"   but no handler found to deal with it."<<std::endl
			 <<"   Will continue and hope for the best."<<std::endl;
	    }
	  }
	}


	(*blit)->UnsetStatus(blob_status::needs_hadrondecays);
      }
    }
  }
  msg_Tracking()<<"--------- Hadron_Decays::Treat - FINISH -------------"<<endl; 
  return (didit?Return_Value::Success:Return_Value::Nothing);
}

void Hadron_Decays::CleanUp() { 
  
}

void Hadron_Decays::Finish(const std::string &) {}

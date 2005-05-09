#include "Hadron_Decays.H"
#include "Message.H"
#include "Hadrons.H"

#ifdef PROFILE__Hadron_Decays
#include "prof.hh"
#else 
#define PROFILE_HERE {}
#define PROFILE_LOCAL(LOCALNAME) {}
#endif

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

bool Hadron_Decays::Treat(ATOOLS::Blob_List * _bloblist, double & weight) 
{
  PROFILE_HERE;
  
  if(p_dechandlers->empty()) return false;

  if (_bloblist->empty()) {
    msg.Error()<<"Potential error in Hadron_Decays::Treat."<<endl
	       <<"   Incoming blob list contains "<<_bloblist->size()<<" entries."<<endl
	       <<"   Continue and hope for the best."<<endl;
    return 0;
  }

  Blob     * myblob, * decblob;
  Particle * check;
  bool found = true;
  msg_Tracking()<<"Hadron_Decays::Treat"<<endl
	<<"Check if Sherpa can handle an occuring hadron decay"<<endl;
  Hadron_Decay_Handler * hdhandler;				// pointer on considered HD Handler
  while (found) {
	found = false;
	for (Blob_List::iterator blit=_bloblist->begin();blit!=_bloblist->end();++blit) {
	  msg.Debugging()<<"status: "<<(*blit)->Status()
		<<" type: "<<(*blit)->Type()
		<<" ID: "<<(*blit)->Id()<<endl;
	  if ( (*blit)->Status()==0 && (*blit)->Type()==btp::Fragmentation  ||
		   (*blit)->Status()==0 && (*blit)->Type()==btp::Hadron_Decay      ) {
		myblob = (*blit);
		found  = true;
		msg.Debugging()<<"Found a blob : "<<myblob->Id()<<"  "<<_bloblist->size()<<endl;
		bool decayed(false); 
		for (int i=0;i<myblob->NOutP();i++) {
		  decayed = false;
		  // check if sherpa can cope with OutParticle and pick implemented ones
		  if( p_dechandlers->find("Sherpa") != p_dechandlers->end() ) {
			hdhandler = (*p_dechandlers)["Sherpa"]; 
			if( hdhandler->GetHadrons()->FindDecay(myblob->OutParticle(i)->RefFlav()) ) 
			{
			  if (myblob->OutParticle(i)->DecayBlob()==NULL) {
				msg.Debugging()<<"Hadron_Decays::Treat (Sherpa): Decay for "
				  <<myblob->OutParticle(i)->Flav()<<std::endl;
				hdhandler->FillHadronDecayBlobs( myblob->OutParticle(i), _bloblist );
			  }
			  else {
				msg.Error()<<"Error in Hadron_Decays::Treat (Sherpa): "<<endl
				  <<"   Unstable particle has a decayblob."<<endl
				  <<"   Terminate run."<<endl;
				abort();
			  }
			  decayed = true;						
			}
		  }
		  // check if Lund can cope with OutParticle 
		  if( p_dechandlers->find("Lund") != p_dechandlers->end() &&
			  !decayed ) {
			hdhandler = (*p_dechandlers)["Lund"];
			if( hdhandler->GetLund()->FindDecay(myblob->OutParticle(i)) ) {
			  if (myblob->OutParticle(i)->DecayBlob()==NULL) {
				msg.Debugging()<<"Hadron_Decays::Treat (Lund): Decay for "
				  <<myblob->OutParticle(i)->Flav()
				  <<" "<<myblob->OutParticle(i)->Flav().Kfcode()<<std::endl;
				hdhandler->FillHadronDecayBlobs( myblob->OutParticle(i), _bloblist );
			  }
			  else {
				msg.Error()<<"Error in Hadron_Decays::Treat (Lund): "<<endl
				  <<"   Unstable particle has a decayblob."<<endl
				  <<"   Terminate run."<<endl;
				abort();
			  }
			}
		  }
		}
		myblob->SetStatus(1);
	  }
	}
  }
  msg_Tracking()<<"--------- Hadron_Decays::Treat - FINISH -------------"<<endl; 
  return false;							// can't do anything more
}

void Hadron_Decays::CleanUp() { 

}

//void Hadron_Decays::ConstructBlob(ATOOLS::Particle * part,ATOOLS::Blob_List * bloblist)
//{
//  msg_Tracking()<<"Hadron_Decays::ConstructBlob"<<endl;
//  p_dechandler->FillHadronDecayBlobs(part,bloblist);
//  msg_Tracking()<<"-------------------- Hadron_Decays::ConstructBlob: ready -------"<<endl;
//}

void Hadron_Decays::Finish(const std::string &) {}

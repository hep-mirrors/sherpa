#include "Hard_Decays.H"
#include "Message.H"

#ifdef PROFILE__Hard_Decays
#include "prof.hh"
#else 
#define PROFILE_HERE {}
#define PROFILE_LOCAL(LOCALNAME) {}
#endif

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Hard_Decays::Hard_Decays(Hard_Decay_Handler * _dechandler) :
  p_dechandler(_dechandler)
{
  m_name      = std::string("Hard_Decays:")+p_dechandler->Name();
  m_type      = eph::Perturbative;
}

Hard_Decays::~Hard_Decays() 
{
}

bool Hard_Decays::Treat(ATOOLS::Blob_List * _bloblist, double & weight) 
{
  PROFILE_HERE;

  if (_bloblist->empty()) {
    msg.Error()<<"Potential error in Hard_Decays::Treat."<<endl
	       <<"   Incoming blob list contains "<<_bloblist->size()<<" entries."<<endl
	       <<"   Continue and hope for the best."<<endl;
    return 0;
  }

  Blob     * myblob, * decblob;
  Particle * check;
  bool found = 1;
  bool hit   = 0;
  while (found) {
    found = 0;
    for (Blob_List::iterator blit=_bloblist->begin();blit!=_bloblist->end();++blit) {
      if ((*blit)->Status()==1 && (*blit)->Type()==btp::FS_Shower) {
	myblob = (*blit);
	//std::cout<<"Found a blob : "<<std::endl<<myblob;
	found  = 1;
	for (int i=0;i<myblob->NOutP();i++) {
	  if ((!myblob->OutParticle(i)->Flav().IsStable()) &&
	      (myblob->OutParticle(i)->DecayBlob()!=NULL)) {
	    //std::cout<<"Dec for : "<<myblob->OutParticle(i)<<std::endl;
	    decblob = myblob->OutParticle(i)->DecayBlob();
	    if (decblob->NInP()!=1) {
	      msg.Error()<<"Error in Hard_Decays::Treat : "<<endl
			 <<"   wrong number of incoming particles for decayblob."<<endl
			 <<"   Terminate run."<<endl;
	      abort();
	    }
	    check = decblob->RemoveInParticle(0,0);
	    if (check->Flav()!=myblob->OutParticle(i)->Flav()) {
	      msg.Error()<<"Error in Hard_Decays::Treat : "<<endl
			 <<"   wrong incoming particle for decayblob : "
			 <<check->Flav()<<" vs. "<<myblob->OutParticle(i)->Flav()<<endl
			 <<"   Terminate run."<<endl;
	      abort();
	    }
	    decblob->AddToInParticles(myblob->OutParticle(i));
	    FillBlob(check,decblob);
	    decblob->SetId();
	    _bloblist->push_back(decblob);
	  }
	}
	//std::cout<<std::endl<<std::endl<<std::endl;
	myblob->SetStatus(0);
	return 1;
      }
    }
  }
  return hit;
}

void Hard_Decays::CleanUp() { 
  p_dechandler->ResetTables(); 
  p_dechandler->ResetSelect(); 
}

void Hard_Decays::FillBlob(Particle * _orig,Blob * _blob)
{
  Particle * decayer = _blob->InParticle(0);
  double lifetime    = decayer->LifeTime();
  _blob->SetPosition(decayer->XProd()+Vec4D(lifetime,decayer->Distance(lifetime)));
  p_dechandler->PerformDecay(_blob);
}

void Hard_Decays::Finish(const std::string &) {}

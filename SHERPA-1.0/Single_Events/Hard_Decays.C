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

Return_Value::code Hard_Decays::Treat(ATOOLS::Blob_List * _bloblist, double & weight) 
{
  PROFILE_HERE;

  if (_bloblist->empty()) {
    msg.Error()<<"Potential error in Hard_Decays::Treat."<<endl
	       <<"   Incoming blob list contains "<<_bloblist->size()<<" entries."<<endl
	       <<"   Continue and hope for the best."<<endl;
    return Return_Value::Error;
  }

  Blob     * myblob, * decblob;
  Particle * check;
  bool found(true),didit(false),switchit(false);
  while (found) {
    found = false;
    for (Blob_List::iterator blit=_bloblist->begin();blit!=_bloblist->end();++blit) {
      if ((*blit)->Has(blob_status::needs_harddecays)) {
	myblob   = (*blit);
	switchit = true;
	for (int i=0;i<myblob->NOutP();i++) {
	  if ((!myblob->OutParticle(i)->Flav().IsStable()) &&
	      (myblob->OutParticle(i)->Status()==1) &&
	      (myblob->OutParticle(i)->DecayBlob()!=NULL)) {
	    decblob = myblob->OutParticle(i)->DecayBlob();
	    if (decblob->NInP()!=1) {
	      msg.Error()<<"Error in Hard_Decays::Treat : "<<endl
			 <<"   wrong number of incoming particles for decayblob."<<endl
			 <<"   Try a new event."<<endl;
	      rvalue.IncError(METHOD);
	      return Return_Value::Error;
	    }
	    check = decblob->RemoveInParticle(0,0);
	    if (check->Flav()!=myblob->OutParticle(i)->Flav()) {
	      msg.Error()<<"Error in Hard_Decays::Treat : "<<endl
			 <<"   wrong incoming particle for decayblob : "
			 <<check->Flav()<<" vs. "<<myblob->OutParticle(i)->Flav()<<endl
			 <<"   Try to repeat the event."<<endl;
	      rvalue.IncError(METHOD);
	      return Return_Value::Error;
	    }
	    found    = didit = true;
	    switchit = false;
	    decblob->AddToInParticles(myblob->OutParticle(i));
	    FillBlob(check,decblob);
	    decblob->SetId();
	    _bloblist->push_back(decblob);
	  }
	}
	if (switchit) myblob->UnsetStatus(blob_status::needs_harddecays);
	//std::cout<<std::endl<<std::endl<<std::endl;
	//myblob->SetStatus(0);
      }
    }
  }
  if (didit) return Return_Value::Success;
  return Return_Value::Nothing;
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

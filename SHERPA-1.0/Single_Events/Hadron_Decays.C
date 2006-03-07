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
   
  // get spin correlation tensor from Signal_Process blob
  bool spincorr (false); 
  Blob_Data_Base *info=(*_bloblist->FindFirst(btp::Signal_Process))
    ["Spin_Correlation_Tensor"];
  SP(Spin_Correlation_Tensor) sct = NULL;
  if (info) {
    spincorr = true;
    sct=info->Get<SP(Spin_Correlation_Tensor)>();
  }

  // treat blob
  Blob * myblob(NULL);
  bool found(true);
  Hadron_Decay_Handler * hdhandler;				// pointer on considered HD Handler
  Particle * myout(NULL);
  int index (0);
  Spin_Density_Matrix * sigma   = NULL;         // SDM for decaying particle
  Spin_Density_Matrix * decmatr = NULL;         // decay matrix
  while (found) {
    found = false;
    bool keeprunning (true);
    for (Blob_List::iterator blit=_bloblist->begin();
	 blit!=_bloblist->end() && keeprunning;++blit) {
      if ((*blit)->Type()==btp::Fragmentation ||
          (*blit)->Type()==btp::Hadron_Decay) {
        myblob = (*blit);
        switch (myblob->Status()) {
          case 2:
            if( p_dechandlers->find("Lund")!= p_dechandlers->end()) {
              hdhandler = (*p_dechandlers)["Lund"];
              hdhandler->PrepareDecays(myblob);
              found = true;
            }
          case 0: {
            bool decayed(false); 
             
            // get vector of out particles and shuffle randomly (if applicable)
            Particle_Vector outparticles = myblob->GetOutParticles();
            random_shuffle( outparticles.begin(), outparticles.end() );
             
            // go through all out particles
            for( Particle_Vector::iterator outit = outparticles.begin();
                outit!=outparticles.end(); ++outit ) {
              myout = (*outit);                         // pointer on out particle
              index = myout->Number();                  // index for this particle in SCT
              if( myout->Flav().IsStable() ) continue;
              decayed = false;
              // check if sherpa can cope with OutParticle and pick implemented ones
              int spin  = myout->Flav().IntSpin();
#ifdef USING__Hadrons
              if( p_dechandlers->find("Sherpa") != p_dechandlers->end() ) {
                hdhandler = (*p_dechandlers)["Sherpa"]; 
                if( hdhandler->GetHadrons()->FindDecay(myout->RefFlav()) ) {
                  if (myout->DecayBlob()==NULL) {
                    msg.Debugging()<<"Hadron_Decays::Treat (Sherpa): Decay for "
                      <<myout->Flav()<<std::endl;
                    sigma   = NULL; 
                    decmatr = NULL; 
                    if( spin && spincorr ) {
                      sigma = new Spin_Density_Matrix( sct->GetSigma(index) );
                      sigma->Normalise();
                      decmatr = new Spin_Density_Matrix();      // decay matrix
                    }
                    hdhandler->FillHadronDecayBlobs( myout, _bloblist, sigma, decmatr );
                    if( decmatr && spincorr ) {
                      decmatr->Normalise();
                      sct->Contract(index,decmatr);             // contract over decay matrix
                    }
                    found       = true;
                    keeprunning = false;		// start again !
                    if( sigma )   delete sigma; 
                    if( decmatr ) delete decmatr;
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
#endif
              // check if Lund can cope with OutParticle 
              if( p_dechandlers->find("Lund") != p_dechandlers->end() && !decayed ) {
                hdhandler = (*p_dechandlers)["Lund"];
                if(hdhandler->GetLund()->FindDecay(myout) ) {
                  if (myout->DecayBlob()==NULL) {
                    hdhandler->FillHadronDecayBlobs( myout, _bloblist );
                    // contract over unit decay matrix to reduce size of sct
                    if( spin && spincorr ) sct->Contract(index,NULL );
                    found       = true;
                    keeprunning = false;		// start again !
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
            if( spincorr )
            if( sct->GetDepth() ) {
              for ( int i=0;i<myblob->NOutP();i++) {
                PRINT_INFO(myblob->OutParticle(i)->Number()<<" "<<myblob->OutParticle(i)->Flav() );
              }
              PRINT_INFO("Not eth. was contracted. "<<sct );
              PRINT_INFO(om::red<<"Will abort."<<om::reset);
              abort();
            }
            myblob->SetStatus(1);
                  }
        }
        if ( (*blit)->Id() == (*(_bloblist->end()-1))->Id() ) {
          keeprunning = false;
        }
      }
    }
  }
//  PRINT_INFO("---------------");
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

#include "Signal_Processes.H"

#ifdef PROFILE__Signal_Processes
#include "prof.hh"
#endif

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Signal_Processes::Signal_Processes(Matrix_Element_Handler * _mehandler,
				   Hard_Decay_Handler * _hdhandler) :
  p_mehandler(_mehandler), p_hdhandler(_hdhandler)
{
  m_name      = string("Signal_Processes : ")+p_mehandler->Name();
  m_type      = string("Perturbative");
}

Signal_Processes::~Signal_Processes()
{
  //  if (p_mehandler) { delete p_mehandler; p_mehandler = NULL; }
  //  Matrix Element Handler will be deleted in Initialization Handler
}


bool Signal_Processes::Treat(Blob_List * _bloblist, double & weight)
{
#ifdef PROFILE__Signal_Processes
  PROFILE_HERE;
#endif
  //  if (_bloblist->size()>1) return 0;
  if (_bloblist->empty()) {
    msg.Error()<<"Potential error in Signal_Processes::Treat."<<endl
	       <<"   Incoming blob list contains "<<_bloblist->size()<<" entries."<<endl
	       <<"   Continue and hope for the best."<<endl;
    return 0;
  }
  
  Blob * myblob;
  bool found = 1;
  bool hit   = 0;
  
  while (found) {
    found = 0;
    for (Blob_Iterator blit=_bloblist->begin();blit!=_bloblist->end();++blit) {
      if (((*blit)->Type()==string("Signal Process : ")) && ((*blit)->Status()==0)) {
	myblob = (*blit);
	found  = 1;
	if (p_mehandler->GenerateOneEvent()) {
	  FillBlob(myblob);
	  weight=p_mehandler->Weight();
	  hit = 1;
	}
      }
      else if (((*blit)->Type()==string("Signal Process : ")) && ((*blit)->Status()==-1)) {
	myblob = (*blit);
	found  = 1;
	if (p_mehandler->GenerateSameEvent()) {
	  FillBlob(myblob);
	  weight=p_mehandler->Weight();
	  hit = 1;
	}
      }
    }
  }
  return hit;
}

void Signal_Processes::CleanUp() { return; }

void Signal_Processes::FillBlob(Blob * _blob)
{
#ifdef PROFILE__Signal_Processes
  PROFILE_HERE;
#endif
  _blob->SetPosition(Vec4D(0.,0.,0.,0.));
  _blob->SetType(_blob->Type()+p_mehandler->ProcessName());
  _blob->SetStatus(1);

  Vec4D cms = Vec4D(0.,0.,0.,0.);
  for (size_t i=0;i<p_mehandler->NIn();i++) cms += p_mehandler->Momenta()[i];
  _blob->SetCMS(cms);
  _blob->SetBeam(-1);

  // make sure that blob is empty
  _blob->DeleteOwnedParticles();

  Particle * particle;
  for (unsigned int i=0;i<p_mehandler->NIn();i++) {
    particle = new Particle(i,p_mehandler->Flavs()[i],p_mehandler->Momenta()[i]);
    particle->SetNumber((long int)particle);
    particle->SetStatus(2);
    particle->SetInfo('G');
    _blob->AddToInParticles(particle);
  }
  bool unstable = false; 
  for (unsigned int i=p_mehandler->NIn();i<p_mehandler->NIn()+p_mehandler->NOut();i++) {
    particle = new Particle(i,p_mehandler->Flavs()[i],p_mehandler->Momenta()[i]);
    if (!(particle->Flav().IsStable())) unstable = true;
    particle->SetStatus(1);
    particle->SetInfo('H');
    _blob->AddToOutParticles(particle);
  }
  if (unstable) {
    if (p_hdhandler->On()) {
      p_hdhandler->ResetTables();
      p_hdhandler->DefineSecondaryDecays(_blob);
      return;
    }
    else {
      msg.Error()<<"Error in Signal_Processes::FillBlob."<<endl
		 <<"   Should treat unstable particles without Hard_Decay_Handler = On."<<endl
		 <<"   Assume no reasonable decay tables. Will abort."<<endl;
      abort();
    }
  }
}


#include "Signal_Processes.H"

#ifdef PROFILE__Signal_Processes
#include "prof.hh"
#else 
#define PROFILE_HERE {}
#define PROFILE_LOCAL(LOCALNAME) {}
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


bool Signal_Processes::Treat(Blob_List * bloblist, double & weight)
{
  PROFILE_HERE;
  // if (bloblist->size()>1) return 0;
  if (bloblist->empty()) {
    msg.Error()<<"Potential error in Signal_Processes::Treat."<<endl
	       <<"   Incoming blob list contains "<<bloblist->size()<<" entries."<<endl
	       <<"   Continue and hope for the best."<<endl;
    return 0;
  }
  
  Blob * myblob;
  bool found = 1;
  bool hit   = 0;
  
  while (found) {
    found = 0;
    for (Blob_Iterator blit=bloblist->begin();blit!=bloblist->end();++blit) {
      if ((*blit)->Type()==btp::Signal_Process && (*blit)->Status()==2) {
	myblob = (*blit);
	found  = 1;
	msg.Tracking()<<" calling GenerateOneEvent() "<<endl;
	if (p_mehandler->GenerateOneEvent()) {
	  weight=p_mehandler->Weight();
	  int  ntrial =p_mehandler->NumberOfTrials();
	  FillBlob(myblob,weight,ntrial);
	  hit = 1;
	}
      }
      else if (((*blit)->Type()==btp::Signal_Process) && ((*blit)->Status()==-1)) {
	myblob = (*blit);
	found  = 1;
	msg.Tracking()<<" calling GenerateSameEvent() "<<endl;
	if (p_mehandler->GenerateSameEvent()) {
	  weight=p_mehandler->Weight();
	  int  ntrial =p_mehandler->NumberOfTrials();
	  FillBlob(myblob,weight,ntrial);
	  hit = 1;
	}
      }
    }
  }
  return hit;
}

void Signal_Processes::CleanUp() { return; }

void Signal_Processes::FillBlob(Blob * blob, const double, const int)
{
  PROFILE_HERE;
  blob->SetPosition(Vec4D(0.,0.,0.,0.));
  blob->SetTypeSpec(p_mehandler->ProcessName());
  blob->SetStatus(1);

  Vec4D cms = Vec4D(0.,0.,0.,0.);
  for (size_t i=0;i<p_mehandler->NIn();i++) cms += p_mehandler->Momenta()[i];
  blob->SetCMS(cms);
  blob->SetBeam(-1);
  
  // make sure that blob is empty
  blob->DeleteOwnedParticles();
  blob->ClearAllData();

  Particle * particle;
  for (unsigned int i=0;i<p_mehandler->NIn();i++) {
    particle = new Particle(i,p_mehandler->Flavours()[i],p_mehandler->Momenta()[i]);
    particle->SetNumber((long int)particle);
    particle->SetStatus(2);
    particle->SetInfo('G');
    blob->AddToInParticles(particle);
  }
  bool unstable = false; 
  for (unsigned int i=p_mehandler->NIn();i<p_mehandler->NIn()+p_mehandler->NOut();i++) {
    particle = new Particle(i,p_mehandler->Flavours()[i],p_mehandler->Momenta()[i]);
    if (!(particle->Flav().IsStable())) unstable = true;
    particle->SetStatus(1);
    particle->SetInfo('H');
    blob->AddToOutParticles(particle);
  }
  if (unstable) {
    if (p_hdhandler->On()) {
      p_hdhandler->ResetTables();
      p_hdhandler->DefineSecondaryDecays(blob);
      return;
    }
    else {
      msg.Error()<<"Error in Signal_Processes::FillBlob."<<endl
		 <<"   Should treat unstable particles without Hard_Decay_Handler = On."<<endl
		 <<"   Assume no reasonable decay tables. Will abort."<<endl;
      abort();
    }
  }
  msg.Tracking()<<" ME_Weight = "<<p_mehandler->Weight()<<std::endl;
  msg.Tracking()<<" Blob: "<<blob<<std::endl;

  // store some additional information
  blob->AddData("ME_Weight",new Blob_Data<double>(p_mehandler->Weight()));
  blob->AddData("ME_NumberOfTrials",new Blob_Data<int>(p_mehandler->NumberOfTrials()));
}

void Signal_Processes::Finish(const std::string &) {}

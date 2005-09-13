#include "Signal_Processes.H"

#include "XS_Base.H"

#ifdef PROFILE__all
#define PROFILE__Signal_Processes
#endif
#ifdef PROFILE__Signal_Processes
#include "prof.hh"
#else 
#define PROFILE_HERE {}
#define PROFILE_LOCAL(LOCALNAME) {}
#endif

using namespace SHERPA;
using namespace ATOOLS;

Signal_Processes::Signal_Processes(Matrix_Element_Handler * mehandler,
				   Hard_Decay_Handler * hdhandler) :
  p_mehandler(mehandler), p_hdhandler(hdhandler),
  m_addedxs(false)
{
  m_name      = std::string("Signal_Processes:")+p_mehandler->Name();
  m_type      = eph::Perturbative;
  p_remnants[0]=mehandler->GetISR_Handler()->GetRemnant(0);
  p_remnants[1]=mehandler->GetISR_Handler()->GetRemnant(1);
  if (p_remnants[0]==NULL || p_remnants[1]==NULL) {
    THROW(critical_error,"No beam remnant handler found.");
  }
}

Signal_Processes::~Signal_Processes()
{
}

bool Signal_Processes::Treat(Blob_List * bloblist, double & weight)
{
  PROFILE_HERE;
  if (bloblist->empty()) {
    msg.Error()<<"Potential error in Signal_Processes::Treat."<<std::endl
	       <<"   Incoming blob list contains no entries."<<std::endl
	       <<"   Continue and hope for the best."<<std::endl;
    return 0;
  }
  
  Blob * myblob;
  bool found = 1;
  bool hit   = 0;

  while (found) {
    found = 0;
    for (Blob_List::iterator blit=bloblist->begin();
	 blit!=bloblist->end();++blit) {
      if ((*blit)->Type()==btp::Signal_Process) {
	if ((*blit)->Status()==0 && !m_addedxs) {
 	  Blob_Data_Base * message = (*(*blit))["ME_Weight"];
 	  double lastweight = message->Get<double>();
 	  message = (*(*blit))["ME_Weight_One"];
	  double nljweight = lastweight;
	  if (message) nljweight = message->Get<double>();
	  p_mehandler->AddEvent(nljweight/rpa.Picobarn(),
				lastweight/rpa.Picobarn(),0);
	  m_addedxs=true;
	}
	else if ((*blit)->Status()==2) {
	  myblob = (*blit);
	  found  = 1;
	  bool success=false;
	  ATOOLS::Blob *isr[2]={NULL,NULL};
	  while (!success) {
	    success=true;
	    if (p_mehandler->GenerateOneEvent()) {
	      EXTRAXS::XS_Base *xs=p_mehandler->GetXS(1);
	      if (xs!=NULL && xs->NAddOut()!=0) {
		xs->SetColours(xs->Momenta());
		for (size_t stop=xs->NAddOut(), i=0;i<stop;++i) {
		  isr[i] = new ATOOLS::Blob();
		  isr[i]->SetType(ATOOLS::btp::IS_Shower);
		  isr[i]->SetTypeSpec("KMR DUPDF");
		  isr[i]->SetId();
		  isr[i]->SetStatus(1);
		  ATOOLS::Particle *parton1 = 
		    new ATOOLS::Particle(-1,xs->AddFlavours()[i],
					 xs->AddMomenta()[i]);
		  parton1->SetNumber();
		  parton1->SetStatus(1);
		  isr[i]->AddToOutParticles(parton1);
		  ATOOLS::Particle *
		    parton2 = new ATOOLS::Particle(-1,xs->AddFlavours()[i],
						   xs->AddMomenta()[i]);
		  parton2->SetNumber();
		  parton2->SetStatus(2);
		  parton2->SetMomentum(parton2->Momentum()
				      +xs->Momenta()[i]);
		  isr[i]->AddToInParticles(parton2);
		  isr[i]->SetBeam(i);
		  if (p_remnants[i]!=NULL) {
		    p_remnants[i]->QuickClear();
		    if (!p_remnants[i]->Extract(parton2)) success=false;
		  }
		  else THROW(fatal_error,"No remnant found.");
		  if (xs->Flavours()[i].IsQuark()) {
		    int anti=xs->Flavours()[i].IsAnti();
		    int newc=ATOOLS::Flow::Counter();
		    if (xs->AddFlavours()[i].IsQuark()) {
		      parton1->SetFlow(2-anti,newc);
		      parton2->SetFlow(2-anti,newc);
		      parton2->SetFlow(1+anti,xs->Colours()[i][anti]);
		    }
		    else {
		      parton1->SetFlow(2-anti,xs->Colours()[i][1-anti]);
		      parton1->SetFlow(1+anti,newc);
		      parton2->SetFlow(1+anti,newc);
		    }
		  }
		  else {
		    int anti=xs->AddFlavours()[i].IsAnti();
		    if (xs->AddFlavours()[i].IsQuark()) {
		      parton1->SetFlow(1+anti,xs->Colours()[i][1-anti]);
		      parton2->SetFlow(1+anti,xs->Colours()[i][anti]);
		    }
		    else {
		      int newc=ATOOLS::Flow::Counter();
		      parton1->SetFlow(1+anti,xs->Colours()[i][1-anti]);
		      parton1->SetFlow(2-anti,newc);
		      parton2->SetFlow(2-anti,newc);
		      parton2->SetFlow(1+anti,xs->Colours()[i][anti]);
		    }
		  }
		}
		blit=bloblist->begin();
	      }
	      if (success) {
		if (!FillBlob(myblob,false,false)) success=false;
		// moved to Beam_Remnant_Handler
		// else p_mehandler->ResetNumberOfTrials();
		weight = p_mehandler->Weight();
		if (isr[0]!=NULL && isr[1]!=NULL) {
		  for (short unsigned int i=0;i<2;++i) {
		    bloblist->push_front(isr[i]);
		    isr[i]->AddToOutParticles(myblob->InParticle(i));
		    ATOOLS::Vec4D sum=isr[i]->CheckMomentumConservation();
		    if (!(sum==ATOOLS::Vec4D()))
		      ATOOLS::msg.Error()
			<<"Signal_Processes::Treat(): "
			<<"4-momentum not conserved: sum = "
			<<sum<<"."<<std::endl;
		  }
		  myblob->SetStatus(0);
		  for (short unsigned int i=0;i<2;++i) {
		    ATOOLS::Blob *fsr = new ATOOLS::Blob();
		    fsr->SetType(ATOOLS::btp::FS_Shower);
		    fsr->SetTypeSpec("KMR DUPDF");
		    fsr->SetId();
		    fsr->SetStatus(1);
		    fsr->AddToInParticles(myblob->OutParticle(i));
		    ATOOLS::Particle *parton = 
		      new ATOOLS::Particle(*myblob->OutParticle(i));
		    parton->SetNumber();
		    parton->SetStatus(1);
		    parton->SetFlow(1,xs->Colours()[2+i][0]);
		    parton->SetFlow(2,xs->Colours()[2+i][1]);
		    fsr->AddToOutParticles(parton);
		    bloblist->push_back(fsr);
		  }
		}
		hit = 1;
	      }
	      else {
		for (short unsigned int i=0;i<2;++i) 
		  if (isr[i]!=NULL) delete isr[i];
	      }
	    }
	  }
	}
	else if (((*blit)->Type()==btp::Signal_Process) &&
		 ((*blit)->Status()==-1)) {
	  myblob = (*blit);
	  found  = 1;
	  bool success=false;
	  while (!success) {
	    success=true;
	    if (p_mehandler->GenerateSameEvent()) {
	      if (!FillBlob(myblob,true)) success=false;
	      weight = p_mehandler->Weight();
	      hit = 1;
	    }
	  }
	}
      }
    }
  }
  return hit;
}

void Signal_Processes::CleanUp() 
{ 
  m_addedxs=false;
}

bool Signal_Processes::FillBlob(Blob * blob,const bool sameevent,
				const bool extract)
{
  PROFILE_HERE; 

  double  weight = p_mehandler->Weight();
  double  procweight = p_mehandler->ProcessWeight();
  int  ntrial = p_mehandler->NumberOfTrials();

  double weight_one=0.;
  int ntrial_one=0;
  if (sameevent) {
    Blob_Data_Base * info=(*blob)["ME_Weight_One"];
    if (info) {
      weight_one = info->Get<double>();
      ntrial_one = (*blob)["ME_NumberOfTrials_One"]->Get<int>();
    }
    else {
      info=(*blob)["ME_Weight"];
      if (info) {
	weight_one = info->Get<double>();
	ntrial_one = (*blob)["ME_NumberOfTrials"]->Get<int>();
      }
      else {
	msg.Error()<<"Signal_Processes::FillBlob(..): "
		   <<"Missing call to OneEvent() before SameEvent() !"
		   <<std::endl;
      }
    }
  }
  blob->SetPosition(Vec4D(0.,0.,0.,0.));
  blob->SetTypeSpec(p_mehandler->ProcessName());
  blob->SetStatus(1);

  Vec4D cms = Vec4D(0.,0.,0.,0.);
  for (size_t i=0;i<p_mehandler->NIn();i++) 
    cms += p_mehandler->Momenta()[i];
  blob->SetCMS(cms);
  blob->SetBeam(-1);
  
  // make sure that blob is empty
  blob->DeleteOwnedParticles();
  blob->ClearAllData();

  bool success=true;

  Particle * particle;
  EXTRAXS::XS_Base *xs=p_mehandler->GetXS(1);
  for (unsigned int i=0;i<p_mehandler->NIn();i++) {
    particle = new Particle(i,p_mehandler->Flavours()[i],
			    p_mehandler->Momenta()[i]);
    particle->SetNumber(0);
    particle->SetStatus(2);
    particle->SetInfo('G');
    if (xs!=NULL) {
      particle->SetFlow(1,xs->Colours()[i][0]);
      particle->SetFlow(2,xs->Colours()[i][1]);
    }
    blob->AddToInParticles(particle);
    if (extract)
      if (p_remnants[i]!=NULL) {
	p_remnants[i]->QuickClear();
	if (!p_remnants[i]->Extract(particle)) success=false;
      }
      else THROW(fatal_error,"No remnant found.");
  }
  bool unstable = false; 
  for (unsigned int i=p_mehandler->NIn();
       i<p_mehandler->NIn()+p_mehandler->NOut();i++) {
    particle = new Particle(i,p_mehandler->Flavours()[i],
			    p_mehandler->Momenta()[i]);
    if (!(particle->Flav().IsStable())) unstable = true;
    particle->SetNumber(0);
    particle->SetStatus(1);
    particle->SetInfo('H');
    if (xs!=NULL) {
      particle->SetFlow(1,xs->Colours()[i][0]);
      particle->SetFlow(2,xs->Colours()[i][1]);
    }
    blob->AddToOutParticles(particle);
  }
  if (unstable) {
    if (p_hdhandler->On()) {
      p_hdhandler->ResetTables();
      p_hdhandler->DefineSecondaryDecays(blob);

      // consider rejection by remnants !!
      return success;
    }
    else {
      msg.Error()<<"Error in Signal_Processes::FillBlob."<<std::endl
		 <<"   No decay tables yet. Abort."<<std::endl;
      abort();
    }
  }

  if (!success && p_mehandler->Weight()!=1.) {
    p_mehandler->SaveNumberOfTrials();
  }
  //  moved to Beam_Remnant_Handler::FillBeamBlobs(..)
  //  else p_mehandler->ResetNumberOfTrials();

  // store some additional information
  if (sameevent) {
    blob->AddData("ME_Weight_One",new Blob_Data<double>(weight_one));
    blob->AddData("ME_NumberOfTrials_One",new Blob_Data<int>(ntrial_one));
  }
  blob->AddData("ME_Weight",new Blob_Data<double>(weight));
  blob->AddData("ME_NumberOfTrials",new Blob_Data<int>(ntrial));
  blob->AddData("Process_Weight",new Blob_Data<double>(procweight));
  return success;
}

void Signal_Processes::Finish(const std::string &) 
{
}

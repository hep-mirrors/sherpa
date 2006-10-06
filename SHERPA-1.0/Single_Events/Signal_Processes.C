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

Return_Value::code Signal_Processes::Treat(Blob_List * bloblist, double & weight)
{
  PROFILE_HERE;
  if (bloblist->empty()) {
    msg.Error()<<"Potential error in Signal_Processes::Treat."<<std::endl
	       <<"   Incoming blob list contains no entries."<<std::endl
	       <<"   Continue and hope for the best."<<std::endl;
    return Return_Value::Error;
  }

  bool found(true), hit(false);

  while (found) {
    found = false;
    for (Blob_List::iterator blit=bloblist->begin();
	 blit!=bloblist->end();++blit) {
      if ((*blit)->Type()==btp::Signal_Process) {
	if ((*blit)->Status()==blob_status::inactive && !m_addedxs) {
 	  Blob_Data_Base * message = (*(*blit))["ME_Weight"];
 	  double lastweight = message->Get<double>();
 	  message = (*(*blit))["ME_Weight_One"];
	  double nljweight = lastweight;
	  if (message) nljweight = message->Get<double>();
	  p_mehandler->AddEvent(nljweight/rpa.Picobarn(),
				lastweight/rpa.Picobarn(),0);
	  m_addedxs=true;
	}
	if ((*blit)->Has(blob_status::needs_signal)) {
	  found  = true;
	  hit    = false;
	  //cout<<"Check : "<<blob_status::internal_flag<<" && "<<(*blit)->Status()
	  //   <<" = "<<(blob_status::internal_flag&(*blit)->Status())
	  //   <<", "<<(0 & 0)<<endl;
	  for (;;) {
	    if ((*blit)->Has(blob_status::internal_flag)) {
	      //cout<<"Same : "<<(*blit)->Has(blob_status::internal_flag)<<endl;
	      if (p_mehandler->GenerateSameEvent() &&
		  FillBlob((*blit),true,true)) hit = true;
	    }
	    else {
	      //cout<<"New  : "<<(*blit)->Has(blob_status::internal_flag)<<endl;
	      if (p_mehandler->GenerateOneEvent() &&
		  FillBlob((*blit),false,true)) hit = true;
	    }
	    if (hit) { 
	      (*blit)->SetStatus(blob_status::needs_showers);
	      //cout<<"Out "<<METHOD<<":"<<std::endl<<(**blit)<<endl;
	      weight = p_mehandler->Weight(); 
	      return Return_Value::Success; 
	    }
	  }
	}
      }
    }
  }
  return Return_Value::Nothing;
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
  double  facscale = p_mehandler->FactorisationScale();
  double  xsecweight = p_mehandler->XSecWeight();
  int  xsecntrial = p_mehandler->NumberOfXSecTrials();
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
  blob->SetStatus(blob_status::needs_showers |
		  blob_status::needs_harddecays);

  Vec4D cms = Vec4D(0.,0.,0.,0.);
  for (size_t i=0;i<p_mehandler->NIn();i++) cms += p_mehandler->Momenta()[i];
  blob->SetCMS(cms);
  blob->SetBeam(-1);
  
  // make sure that blob is empty
  blob->DeleteOwnedParticles();
  blob->ClearAllData();

  bool success=true;

  Particle * particle(NULL);
  EXTRAXS::XS_Base *xs=p_mehandler->GetXS(1);
  for (unsigned int i=0;i<p_mehandler->NIn();i++) {
    particle = new Particle(i,p_mehandler->Flavours()[i],
			    p_mehandler->Momenta()[i]);
    particle->SetNumber(0);
    particle->SetStatus(part_status::decayed);
    particle->SetInfo('G');
    if (xs!=NULL) {
      particle->SetFlow(1,xs->Colours()[i][0]);
      particle->SetFlow(2,xs->Colours()[i][1]);
    }
    blob->AddToInParticles(particle);
    if (extract)
      if (p_remnants[i]!=NULL) {
	p_remnants[i]->QuickClear();
	if (!p_remnants[i]->TestExtract(particle)) success=false;
      }
      else THROW(fatal_error,"No remnant found.");
  }
  bool unstable = false; 
  for (unsigned int i=p_mehandler->NIn();
       i<p_mehandler->NIn()+p_mehandler->NOut();i++) {
    particle = new Particle(i,p_mehandler->Flavours()[i],
			    p_mehandler->Momenta()[i]);
    if( particle->Flav().Kfcode() != kf::tau )
      if (!(particle->Flav().IsStable())) unstable = true;
    particle->SetNumber(0);
    particle->SetStatus(part_status::active);
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
		 <<"   No hard decay tables for "<<particle->Flav()<<". Abort."<<std::endl;
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
  blob->AddData("Factorisation_Scale",new Blob_Data<double>(facscale));
  blob->AddData("XS_Weight",new Blob_Data<double>(xsecweight));
  blob->AddData("XS_NumberOfTrials",new Blob_Data<int>(xsecntrial));
//   Spin_Correlation_Tensor* SCT = p_mehandler->GetSpinCorrelations();
//   if (SCT!=NULL)
//     blob->AddData("Spin_Correlation_Tensor",
// 		  new Blob_Data<SP(Spin_Correlation_Tensor) >(SCT));
  return success;
}

void Signal_Processes::Finish(const std::string &) 
{
}

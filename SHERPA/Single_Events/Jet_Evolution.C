#include "SHERPA/Single_Events/Jet_Evolution.H"

#include "SHERPA/PerturbativePhysics/Perturbative_Interface.H"
#include "REMNANTS/Main/Remnant_Handler.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Phys/NLO_Types.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_MPI.H"
#include "MODEL/Main/Running_AlphaS.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace PDF;
using namespace std;

Jet_Evolution::Jet_Evolution(Matrix_Element_Handler *_mehandler,
                             Hard_Decay_Handler * _dechandler,
                             Decay_Handler_Base *_hdhandler,
			     MI_Handler *_mihandler,
			     Soft_Collision_Handler *_schandler,
                             const Shower_Handler_Map& showers)
{
  Shower_Handler_Map::const_iterator shIter=showers.find(isr::hard_process);
  m_name      = string("Jet_Evolution:")+shIter->second->ShowerGenerator();
  m_type      = eph::Perturbative;

  Perturbative_Interface * interface;
  shIter=showers.find(isr::hard_process);
  interface = new Perturbative_Interface(_mehandler,
                                         _dechandler,
					 shIter->second);
  if (interface!=NULL) m_interfaces.insert(make_pair("SignalMEs",interface));

  shIter=showers.find(isr::hard_subprocess);
  interface = new Perturbative_Interface(_hdhandler, shIter->second);
  if (interface!=NULL) 
    m_interfaces.insert(make_pair("HadronDecays",interface));

  if (_mihandler) {
    interface = new Perturbative_Interface(_mihandler, shIter->second);
    if (interface!=NULL) m_interfaces.insert(make_pair("MPIs",interface));
  }
  if (_schandler) {
    interface = new Perturbative_Interface(_schandler, shIter->second);
    if (interface!=NULL) 
      m_interfaces.insert(make_pair("SoftCollisions",interface));
  }
  p_remnants = _mehandler->Remnants();
}

Jet_Evolution::~Jet_Evolution() 
{ 
  while (m_interfaces.size()>0) {
    if (m_interfaces.begin()->second!=NULL) delete m_interfaces.begin()->second;
    m_interfaces.erase(m_interfaces.begin());
  }
}


Return_Value::code Jet_Evolution::Treat(Blob_List * bloblist)
{
  if (bloblist->empty()) {
    msg_Error()<<"Potential error in Jet_Evolution::Treat."<<endl
	       <<"   Incoming blob list contains "<<bloblist->size()
	       <<" entries."<<endl
	       <<"   Continue and hope for the best."<<endl;
    return Return_Value::Error;
  }
  PertInterfaceIter piIter;
  string tag("SignalMEs");
  bool hit(false), found(true);
  while (found) {
    found = false;
    for (size_t i=0;i<bloblist->size();++i) {
      Blob * meblob = (*bloblist)[i];
      if (meblob->Has(blob_status::needs_showers) &&
          meblob->Type()!=btp::Hard_Decay) {
	switch (int(meblob->Type())) {
	  case (int(btp::Signal_Process)) :
            tag = string("SignalMEs");
            MODEL::as->SetActiveAs(PDF::isr::hard_process);
	    break;
	  case (int(btp::Hard_Collision)) : 
	    tag = string("MPIs"); 
	    if (meblob->TypeSpec()=="MinBias") 
	      tag = string("SoftCollisions"); 
            MODEL::as->SetActiveAs(PDF::isr::hard_subprocess);
	    break;
	  case (int(btp::Hadron_Decay))   : 
	    tag = string("HadronDecays"); 
            MODEL::as->SetActiveAs(PDF::isr::hard_subprocess);
	    break;
	  default:
	    msg_Error()<<"ERROR in "<<METHOD<<": "
		       <<"Do not have an interface for this type of blob.\n"
		       <<(*meblob)<<"\n   Will abort.\n";
	    THROW(fatal_error,"No perturbative interface found.");
	}
	piIter = m_interfaces.find(tag);
	if (piIter==m_interfaces.end()) {
	  msg_Error()<<"Error in Jet_Evolution::Treat: "
		     <<"No Perturbative_Interface found for type "<<tag<<"\n"
		     <<"   Abort the run.\n";
	  THROW(fatal_error,"No perturbative interface found.");
	}	
	switch (AttachShowers(meblob,bloblist,piIter->second)) {
	case Return_Value::Success:
	  found = hit = true;
	  break;
	case Return_Value::New_Event  : return Return_Value::New_Event;
	case Return_Value::Retry_Event: return Return_Value::Retry_Event;
	case Return_Value::Nothing    : return Return_Value::Nothing;
	case Return_Value::Error      : return Return_Value::Error;
	default:
	  msg_Error()<<"ERROR in "<<METHOD<<": Unexpected status of AttachShowers for\n"
		     <<(*meblob)<<"   Return 'Error' and hope for the best.\n";
	  return Return_Value::Error;
	}
      }
    }
    if (found) hit = true;
    Reset();
  }
  if (hit) {
    // enable shower generator independent FS QED correction to ME
    // TODO: check first, whether shower did FS QED
    if (!bloblist->FourMomentumConservation()) {
      msg_Tracking()<<METHOD<<" found four momentum conservation error.\n";
      return Return_Value::New_Event;
    }
    Blob * showerblob = bloblist->FindLast(btp::Shower);
    showerblob->AddStatus(blob_status::needs_extraQED);
    return Return_Value::Success;
  }
  // Capture potential problem with empty remnants here.
  // This should only happen after retrying an event has been called.  In this case
  // we find the last (and hopefully only) shower blob and extract its initiators.
  Blob * showerblob = bloblist->FindLast(btp::Shower);
  if (showerblob!=NULL && !p_remnants->ExtractShowerInitiators(showerblob)) return Return_Value::New_Event;
  return Return_Value::Nothing;
}

Return_Value::code Jet_Evolution::
AttachShowers(Blob * blob,Blob_List * bloblist,
	      Perturbative_Interface * interface) 
{
  if (!interface->Shower()->On() ||
      (interface->MEHandler() && 
       interface->MEHandler()->Process()->Info().m_nlomode==nlo_mode::fixedorder)) {
    AftermathOfNoShower(blob,bloblist);
    return Return_Value::Nothing;
  }
  int shower(0);
  Return_Value::code stat(interface->DefineInitialConditions(blob,bloblist));
  if (stat==Return_Value::New_Event ||
      stat==Return_Value::Retry_Event) {
    interface->CleanUp();
    return stat;
  }
  if (blob->Type()!=::btp::Hadron_Decay) {
    msg_Debugging()<<METHOD<<"(): Setting scale for MI {\n";
    double scale(0.0);
    Cluster_Amplitude *ampl(interface->Amplitude());
    while (ampl->Next()) ampl=ampl->Next();
    msg_Debugging()<<*ampl<<"\n";
    scale=sqrt(ampl->MuQ2());
    blob->AddData("MI_Scale",new Blob_Data<double>(scale));
    msg_Debugging()<<"} -> p_T = "<<scale<<"\n";
  }
  switch (stat) {
  case Return_Value::Success:
    if (blob->Type()!=::btp::Hadron_Decay) 
      DefineInitialConditions(blob,bloblist, interface);
    if (blob->NInP()==1) shower = interface->PerformDecayShowers();
    if (blob->NInP()==2) shower = interface->PerformShowers();
    switch (shower) {
    case 1:
      // No Sudakov rejection
      Reset();
      if (AftermathOfSuccessfulShower(blob,bloblist,interface)) {
	interface->CleanUp();
	return Return_Value::Success;
      }
      blob->SetStatus(blob_status::inactive);
      CleanUp();
      return Return_Value::New_Event;
    case 0:
      // Sudakov rejection
      Reset();
      CleanUp();
      return Return_Value::New_Event;
    default:
      THROW(fatal_error,"Invalid return value from shower");
    }
  case Return_Value::Nothing:
    if (AftermathOfNoShower(blob,bloblist)) {
      interface->CleanUp();
      return Return_Value::Success;
    }
    blob->SetStatus(blob_status::inactive);
    CleanUp();
    return Return_Value::New_Event;
  case Return_Value::Error:
    msg_Error()<<"ERROR in "<<METHOD<<":"<<std::endl
	       <<"   DefineInitialConditions yields an error for "<<std::endl<<(*blob)
	       <<"   Return 'Error' and hope for the best."<<std::endl;
    blob->SetStatus(blob_status::inactive);
    CleanUp();
    return Return_Value::Error;
  default :
    msg_Error()<<"ERROR in "<<METHOD<<":"<<std::endl
	       <<"   Unexpected status of DefineInitialConditions for "<<std::endl<<(*blob)
	       <<"   Return 'Error' and hope for the best."<<std::endl;
    blob->SetStatus(blob_status::inactive);
    CleanUp();
    return Return_Value::Error;    
  }
  return Return_Value::Error;    
}


bool Jet_Evolution::AftermathOfNoShower(Blob * blob,Blob_List * bloblist)
{
  Blob * noshowerblob = new Blob();
  noshowerblob->SetType(btp::Shower);
  for (size_t i=0; i<blob->GetInParticles().size();++i) {
    noshowerblob->AddToOutParticles(blob->InParticle(i));
    noshowerblob->AddToInParticles(new Particle(*blob->InParticle(i)));
    noshowerblob->InParticle(i)->SetBeam(i);
    blob->InParticle(i)->SetStatus(part_status::decayed);
  }
  for (size_t i=0; i<blob->GetOutParticles().size();++i) {
    if (blob->OutParticle(i)->DecayBlob()) continue;
    noshowerblob->AddToInParticles(blob->OutParticle(i));
    noshowerblob->AddToOutParticles(new Particle(*blob->OutParticle(i)));
    blob->OutParticle(i)->SetStatus(part_status::decayed);
  }
  noshowerblob->SetStatus(blob_status::needs_beams|blob_status::needs_hadronization);
  if (blob->Type()!=btp::Hadron_Decay) {
    noshowerblob->AddStatus(blob_status::needs_reconnections);
  }
  noshowerblob->SetId();
  noshowerblob->SetTypeSpec("No_Shower");
  bloblist->push_back(noshowerblob);
  blob->SetStatus(blob_status::inactive);
  return p_remnants->ExtractShowerInitiators(noshowerblob);
}

bool Jet_Evolution::
AftermathOfSuccessfulShower(Blob * blob,Blob_List * bloblist,
			    Perturbative_Interface * interface)
{
  if (blob->NInP()==1 && 
      blob->Type()!=btp::Hadron_Decay) blob->InParticle(0)->SetInfo('h');
  interface->FillBlobs();
  blob->UnsetStatus(blob_status::needs_showers);
  Blob * showerblob = (!interface->Shower()->On()?
		       CreateMockShowerBlobs(blob,bloblist):
		       bloblist->FindLast(btp::Shower));
  if (showerblob!=NULL) {
    if (blob->Type()!=btp::Hadron_Decay) {
      showerblob->AddStatus(blob_status::needs_reconnections);
    }
    return p_remnants->ExtractShowerInitiators(showerblob);
  }
  return true;
}

ATOOLS::Blob * Jet_Evolution::
CreateMockShowerBlobs(Blob * const meblob,Blob_List * const bloblist) {
  Blob * ISRblob = NULL;
  if (meblob->NInP()!=1) {
    for (int i=0;i<2;i++) {
      // new ISR Blob
      ISRblob = new Blob();
      ISRblob->SetType(btp::Shower);
      ISRblob->SetStatus(blob_status::needs_beams);
      Particle * part = new Particle(*meblob->InParticle(i));
      part->SetStatus(part_status::decayed);
      part->SetBeam(int(meblob->InParticle(1-i)->Momentum()[3] 
			 > meblob->InParticle(i)->Momentum()[3]));
      ISRblob->AddToInParticles(part);
      ISRblob->AddToOutParticles(meblob->InParticle(i));
      meblob->InParticle(i)->SetStatus(part_status::decayed);
      ISRblob->SetId();
      bloblist->insert(bloblist->begin(),ISRblob);
    }
  }
  for (int i=0;i<meblob->NOutP();i++) {
    Blob * FSRblob = new Blob();
    FSRblob->SetType(btp::Shower);
    FSRblob->SetStatus(blob_status::needs_hadronization);
    if (meblob->Type()!=btp::Hadron_Decay) {
      FSRblob->AddStatus(blob_status::needs_reconnections);
    }
    Particle * part = new Particle(*meblob->OutParticle(i));
    if (meblob->OutParticle(i)->DecayBlob()) {
      Blob * dec  = meblob->OutParticle(i)->DecayBlob();
      if (dec->Type()==btp::Hard_Decay) {
	dec->RemoveInParticle(meblob->OutParticle(i));
	dec->AddToInParticles(part);
      }
    }
    FSRblob->AddToInParticles(meblob->OutParticle(i));
    meblob->OutParticle(i)->SetStatus(part_status::decayed);
    FSRblob->AddToOutParticles(part);
    FSRblob->SetId();
    bloblist->push_back(FSRblob);
  }
  return ISRblob;
}

void Jet_Evolution::CleanUp(const size_t & mode) 
{ 
  for (PertInterfaceIter piIter=m_interfaces.begin();
       piIter!=m_interfaces.end(); ++piIter) {
    piIter->second->CleanUp();
  }
}

void Jet_Evolution::Reset()
{
  for (PertInterfaceIter piIter=m_interfaces.begin();
       piIter!=m_interfaces.end(); ++piIter) {
    piIter->second->Shower()->GetISRHandler()->Reset(0);
    piIter->second->Shower()->GetISRHandler()->Reset(1);
  }
}

bool Jet_Evolution::
DefineInitialConditions(const Blob *blob,const Blob_List *bloblist,
			Perturbative_Interface *interface)
{ 
  Reset();
  msg_Debugging()<<METHOD<<"(): {\n";
  for (::Blob_List::const_iterator blit=bloblist->begin();
       blit!=bloblist->end();++blit) {
    if ((*blit)->Type()==::btp::Shower) {
      //Update(*blit,0, interface);
      //Update(*blit,1, interface);
    }
  }
  msg_Debugging()<<"}\n";
  return true;
}

void Jet_Evolution::Update(Blob *blob,const size_t beam,
                           Perturbative_Interface *interface)
{ 
  size_t cbeam=0;
  for (int i=0;i<blob->NInP();++i) {
    Particle *cur=blob->InParticle(i);
    if (!cur->Flav().Strong() || cur->ProductionBlob()) continue;
    if (cbeam==beam) {
      msg_Debugging()<<"  "<<*cur<<", beam = "<<beam<<"\n";
      p_remnants->Extract(cur,beam);
      return;
    }
    ++cbeam;
  }
}

void Jet_Evolution::Finish(const string &) 
{
}

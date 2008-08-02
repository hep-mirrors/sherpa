#include "Jet_Evolution.H"
#include "HadronDecays_Apacic_Interface.H"
#include "SimpleXS_Apacic_Interface.H"
#include "Amegic_Apacic_Interface.H"

#ifdef PROFILE__all
#define PROFILE__Jet_Evolution
#endif
#ifdef PROFILE__Jet_Evolution
#include "prof.hh"
#else 
#define PROFILE_HERE {}
#define PROFILE_LOCAL(LOCALNAME) {}
#endif


using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Jet_Evolution::Jet_Evolution(MEHandlersMap *_mehandlers,Shower_Handler *_showerhandler) :
  p_showerhandler(_showerhandler)
{
  m_name      = string("Jet_Evolution:")+p_showerhandler->ShowerGenerator();
  m_type      = eph::Perturbative;

  Perturbative_Interface * interface;
  MEHandlerIter            meIter;
  for (meIter=_mehandlers->begin();meIter!=_mehandlers->end();++meIter) {
    interface=NULL;
    if (meIter->second->Name()==string("Amegic") &&
	p_showerhandler->ShowerGenerator()==string("Apacic")) 
      interface = new Amegic_Apacic_Interface(meIter->second,p_showerhandler);
    if (meIter->second->Name()==string("SimpleXS") &&
	p_showerhandler->ShowerGenerator()==string("Apacic")) 
      interface = new SimpleXS_Apacic_Interface(meIter->second,p_showerhandler);
    if (interface!=NULL) m_interfaces.insert(make_pair(meIter->first,interface));
  }
  interface = new HadronDecays_Apacic_Interface(NULL,p_showerhandler);
  if (interface!=NULL) m_interfaces["HadronDecays"]=interface;
}

Jet_Evolution::~Jet_Evolution() 
{ 
  while (m_interfaces.size()>0) {
    if (m_interfaces.begin()->second!=NULL) delete m_interfaces.begin()->second;
    m_interfaces.erase(m_interfaces.begin());
  }
}


Return_Value::code Jet_Evolution::Treat(Blob_List * bloblist, double & weight)
{
  PROFILE_LOCAL("Jet_Evolution::Treat");
  if (bloblist->empty()) {
    msg_Error()<<"Potential error in Jet_Evolution::Treat."<<endl
	       <<"   Incoming blob list contains "<<bloblist->size()<<" entries."<<endl
	       <<"   Continue and hope for the best."<<endl;
    return Return_Value::Error;
  }
  PertInterfaceIter piIter;
  string tag("SignalMEs");
  bool hit(false), found(true);
  Blob * blob;
  while (found) {
    found = false;
    for (size_t i=0;i<bloblist->size();++i) {
      blob = (*bloblist)[i];
      //std::cout<<METHOD<<" for "<<int(blob->Type())
      //	       <<"; check for status "<<int(blob->Status())<<endl;
      if (blob->Has(blob_status::needs_showers)) {
	FillDecayBlobMap(blob,bloblist);
	switch (int(blob->Type())) {
	  case (int(btp::Signal_Process)) : break;
	  case (int(btp::Hard_Collision)) : tag = string("MIMEs"); break;
	  case (int(btp::Hard_Decay))     : tag = string("HardDecays"); break;
	  case (int(btp::Hadron_Decay))   : tag = string("HadronDecays"); break;
	  default:
	    msg_Error()<<"ERROR in "<<METHOD<<":"<<std::endl
		       <<"   Do not have an interface for this type of blob:"<<std::endl
		       <<(*blob)<<std::endl
		       <<"   Will abort."<<std::endl;
	    abort();
	}
	piIter = m_interfaces.find(tag);
	if (piIter==m_interfaces.end()) {
	  msg_Error()<<"Error in Jet_Evolution::Treat :"<<endl
		     <<"   No Perturbative_Interface found for type : "<<tag<<endl
		     <<"   Abort the run."<<endl;
	  abort();
	}	
	switch (AttachShowers(blob,bloblist,piIter->second)) {
	case Return_Value::Success:
	  found = hit = true;
	  weight *= piIter->second->Weight();
	  break;
	case Return_Value::New_Event  : return Return_Value::New_Event;
	case Return_Value::Retry_Event: return Return_Value::Retry_Event;
	case Return_Value::Error      : return Return_Value::Error;
	default:
	  msg_Error()<<"ERROR in "<<METHOD<<":"<<std::endl
		     <<"   Unexpected status of AttachShowers for "<<std::endl<<(*blob)
		     <<"   Return 'Error' and hope for the best."<<std::endl;
	  return Return_Value::Error;
	}
      }
    }
    if (found) hit = true;
    Reset();
  }
  if (hit) return Return_Value::Success;
  return Return_Value::Nothing;
}

Return_Value::code Jet_Evolution::AttachShowers(Blob * blob,Blob_List * bloblist,
						Perturbative_Interface * interface) 
{
  if (blob->Type()==btp::Hard_Collision && !p_showerhandler->ShowerMI() && 
      blob->Has(blob_status::needs_showers)) {
    AftermathOfNoShower(blob,bloblist);
    return Return_Value::Success;
  }
  int shower(0);
  switch (interface->DefineInitialConditions(blob)) {
  case Return_Value::Success:
    DefineInitialConditions(blob,bloblist);
    if (blob->NInP()==1) shower = interface->PerformDecayShowers();
    if (blob->NInP()==2) shower = interface->PerformShowers();
    switch (shower) {
    case 1: 
      Reset();
      AftermathOfSuccessfulShower(blob,bloblist,interface);    
      return Return_Value::Success;
    case -1:
      Reset();
      p_showerhandler->CleanUp();
      if (blob->Type()==btp::Hadron_Decay) {
        Particle* inpart = blob->InParticle(0);
        inpart->SetStatus(part_status::active);
        inpart->ProductionBlob()->SetStatus(blob_status::needs_hadrondecays);
	return Return_Value::Success;
      }
      else if (blob->Type()!=btp::Signal_Process) 
	interface->CleanBlobList(bloblist,blob->Type());
      return Return_Value::Retry_Event;
    default:
      Reset();
      p_showerhandler->CleanUp();
      interface->CleanUp();
      return Return_Value::New_Event;
    }
  case Return_Value::Nothing:
    AftermathOfNoShower(blob,bloblist);
    return Return_Value::Success;
  case Return_Value::New_Event:
    p_showerhandler->CleanUp();
    interface->CleanUp();
    return Return_Value::New_Event;
  case Return_Value::Retry_Event:
    p_showerhandler->CleanUp();
    interface->CleanUp();
    return Return_Value::Retry_Event;
  case Return_Value::Error:
    msg_Error()<<"ERROR in "<<METHOD<<":"<<std::endl
	       <<"   DefineInitialConditions yields an error for "<<std::endl<<(*blob)
	       <<"   Return 'Error' and hope for the best."<<std::endl;
    blob->SetStatus(blob_status::inactive);
    p_showerhandler->CleanUp();
    interface->CleanUp();
    return Return_Value::Error;
  default :
    msg_Error()<<"ERROR in "<<METHOD<<":"<<std::endl
	       <<"   Unexpected status of DefineInitialConditions for "<<std::endl<<(*blob)
	       <<"   Return 'Error' and hope for the best."<<std::endl;
    blob->SetStatus(blob_status::inactive);
    p_showerhandler->CleanUp();
    return Return_Value::Error;    
  }
  return Return_Value::Error;    
}


void Jet_Evolution::AftermathOfNoShower(Blob * blob,Blob_List * bloblist)
{
  Blob * myblob;
  for (int i=0;i<2;i++) {
    myblob = new Blob();
    myblob->SetType(btp::IS_Shower);
    if (Sign(blob->InParticle(i)->Momentum()[3])==1-2*i) myblob->SetBeam(i);
    else myblob->SetBeam(1-i);
    myblob->SetStatus(blob_status::needs_beams);
    Particle * p = new Particle(*blob->InParticle(i));
    p->SetStatus(part_status::decayed);
    myblob->AddToInParticles(p);
    myblob->AddToOutParticles(blob->InParticle(i));
    blob->InParticle(i)->SetStatus(part_status::decayed);
    myblob->SetId();
    bloblist->insert(bloblist->begin(),myblob);
  }
  for (int i=0;i<blob->NOutP();i++) {
    myblob = new Blob();
    myblob->SetType(btp::FS_Shower);
    myblob->SetStatus(blob_status::needs_hadronization);
    Particle * p = new Particle(*blob->OutParticle(i));
    myblob->AddToInParticles(blob->OutParticle(i));
    blob->OutParticle(i)->SetStatus(part_status::decayed);
    myblob->AddToOutParticles(p);
    myblob->SetId();
    bloblist->push_back(myblob);
  }
  //std::cout<<METHOD<<": found a blob for status=0"<<std::endl<<(*blob)<<std::endl;
  blob->SetStatus(blob_status::inactive);
}

void Jet_Evolution::AftermathOfSuccessfulShower(Blob * blob,Blob_List * bloblist,
						Perturbative_Interface * interface)
{
  Blob * myblob;
  if (blob->NInP()==1 && 
      blob->Type()!=btp::Hadron_Decay) blob->InParticle(0)->SetInfo('h');
  interface->FillBlobs(bloblist);
  //std::cout<<METHOD<<": found a blob for status=0"<<std::endl<<(*blob)<<std::endl;
  blob->SetStatus(blob_status::inactive);
  if (blob->NInP()!=1 && !p_showerhandler->ISROn()) {
    for (int i=0;i<2;i++) {
      // new ISR Blob
      myblob = new Blob();
      myblob->SetType(btp::IS_Shower);
      myblob->SetStatus(blob_status::needs_harddecays |
			blob_status::needs_beams |
			blob_status::needs_hadronization);
      myblob->SetBeam( int( blob->InParticle(1-i)->Momentum()[3] 
			    > blob->InParticle(i)->Momentum()[3]) );
      Particle * p = new Particle(*blob->InParticle(i));
      p->SetStatus(part_status::decayed);
      myblob->AddToInParticles(p);
      myblob->AddToOutParticles(blob->InParticle(i));
      blob->InParticle(i)->SetStatus(part_status::decayed);
      myblob->SetId();
      bloblist->insert(bloblist->begin(),myblob);
    }
  }
  if (!(p_showerhandler->FSROn())) {
    for (int i=0;i<blob->NOutP();i++) {
      myblob = new Blob();
      myblob->SetType(btp::FS_Shower);
      myblob->SetBeam(i);
      myblob->SetStatus(blob_status::needs_harddecays |
			blob_status::needs_hadronization);
      Particle * p = new Particle(*blob->OutParticle(i));
      if (blob->OutParticle(i)->DecayBlob()) {
	Blob * dec  = blob->OutParticle(i)->DecayBlob();
	if (dec->Type()==btp::Hard_Decay) {
	  dec->RemoveInParticle(blob->OutParticle(i));
	  dec->AddToInParticles(p);
	}
      }
      myblob->AddToInParticles(blob->OutParticle(i));
      blob->OutParticle(i)->SetStatus(part_status::decayed);
      myblob->AddToOutParticles(p);
      myblob->SetId();
      bloblist->push_back(myblob);
    }
  }
  else {
    if (blob->Type()!=btp::Hadron_Decay) SetDecayBlobPointers(blob,bloblist);
  }
}

void Jet_Evolution::SetDecayBlobPointers(Blob * blob,Blob_List * bloblist) 
{ 
  if (m_decmap.empty()) return;
  Blob     * dec;
  Particle * partin, * partout;
  Particle_Blob_Map::iterator pbiter;
  for (int i=0;i<blob->NOutP();i++) {
    partin = blob->OutParticle(i);
    if (partin->Flav().IsStable()) continue;
    pbiter = m_decmap.find(partin);
    if (pbiter==m_decmap.end()) {
      msg_Error()<<"ERROR in "<<METHOD<<":"<<endl
		 <<"   Did not find particle in map of decay blobs."<<endl
		 <<"   Particle : "<<partin<<endl
		 <<"   Will abort the run."<<endl;
      abort();
    }
    dec     = pbiter->second;
    if (dec->Type()==btp::Hard_Decay) {
      //cout<<"Check for "<<partin->Flav()<<" / "<<partin->FinalMass()<<endl;
      partout = FollowUp(partin,dec);
      partout->SetDecayBlob(dec);
      dec->RemoveInParticle(partin);
      dec->AddToInParticles(partout);
      //cout<<"Match : "<<partin->Number()<<" -> "<<partout->Number()<<endl
      //	       <<(*dec)<<endl;
    }
  }
}

Particle * Jet_Evolution::FollowUp(Particle * partin,Blob * dec) 
{
  Blob * current = partin->DecayBlob();
  Particle * partout;
  if (current==0 || current==dec) return partin;
  //cout<<"Current : "<<current->Id()<<" <- "<<partin->Number()
  //	   <<" / "<<partin->Flav()<<endl;
  for (int i=0;i<current->NOutP();++i) {
    partout = current->OutParticle(i);
    if (partout->Flav()==partin->Flav() &&
	partout->FinalMass()==partin->FinalMass()) return FollowUp(partout,dec);
  }
  msg_Error()<<"ERROR in JetEvolution::FollowUp:"<<endl
	     <<"   Did not find a suitable particle to follow up decay initiator through blob list."<<endl
	     <<"   Particle = "<<partin<<endl
	     <<"   Will abort the run."<<endl;
  abort();
}

void Jet_Evolution::FillDecayBlobMap(Blob * blob,Blob_List * bloblist) 
{
  if (blob->Type()==btp::Hadron_Decay) {
    return;
  }
  //  cout<<"Which blob ? "<<endl<<(*blob)<<endl;
  for (int i=0;i<blob->NOutP();++i) {
    if (blob->OutParticle(i)->DecayBlob()) 
      m_decmap.insert(make_pair(blob->OutParticle(i),
				     blob->OutParticle(i)->DecayBlob()));
  }
  //  cout<<"length of map : "<<m_decmap.size()<<endl;
}

void Jet_Evolution::CleanUp() 
{ 
  m_decmap.clear();
  p_showerhandler->CleanUp();
}

void Jet_Evolution::Reset() 
{
  p_showerhandler->GetISRHandler()->Reset(0);
  p_showerhandler->GetISRHandler()->Reset(1);
}

bool Jet_Evolution::DefineInitialConditions(const Blob *blob,
					    const Blob_List *bloblist) 
{ 
  Reset();
  for (::Blob_List::const_iterator blit=bloblist->begin();
       blit!=bloblist->end();++blit) {
    if (((*blit)->Type()==::btp::Signal_Process ||
	 (*blit)->Type()==::btp::Hard_Collision) && *blit!=blob) {
      Update(*blit,0);
      Update(*blit,1);
    }
  }
  return true;
}

void Jet_Evolution::Update(const Blob *blob,const size_t beam) 
{ 
  for (size_t i=0;i<(size_t)blob->NInP();++i) {
    const Particle *cur=blob->ConstInParticle(i);
    if (i==beam || blob->NInP()<=1) {
      if (cur->ProductionBlob()!=NULL) Update(cur->ProductionBlob(),beam);
      else {
	p_showerhandler->GetISRHandler()->Extract(cur->Flav(),cur->Momentum()[0],beam);
	return;
      }
    }
  }
}

void Jet_Evolution::Finish(const string &) 
{
}

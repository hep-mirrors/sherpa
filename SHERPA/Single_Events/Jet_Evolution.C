#include "SHERPA/Single_Events/Jet_Evolution.H"

#include "SHERPA/PerturbativePhysics/Perturbative_Interface.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Jet_Evolution::Jet_Evolution(MEHandlersMap *_mehandlers,HDHandlersMap *_hdhandlers,
			     MI_Handler *_mihandler,Shower_Handler *_showerhandler) :
  p_showerhandler(_showerhandler)
{
  m_name      = string("Jet_Evolution:")+p_showerhandler->ShowerGenerator();
  m_type      = eph::Perturbative;

  Perturbative_Interface * interface;
  MEHandlerIter            meIter;
  for (meIter=_mehandlers->begin();meIter!=_mehandlers->end();++meIter) {
    interface = new Perturbative_Interface(meIter->second,p_showerhandler);
    if (interface!=NULL) m_interfaces.insert(make_pair(meIter->first,interface));
  }
  HDHandlersIter            hditer;
  for (hditer=_hdhandlers->begin();hditer!=_hdhandlers->end();++hditer) {
    interface = new Perturbative_Interface(hditer->second,p_showerhandler);
    if (interface!=NULL) m_interfaces.insert(make_pair("HadronDecays",interface));
    break;
  }
  if (_mihandler) {
    interface = new Perturbative_Interface(_mihandler,p_showerhandler);
    if (interface!=NULL) m_interfaces.insert(make_pair("MPIs",interface));
  }
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
	  case (int(btp::Hard_Collision)) : tag = string("MPIs"); break;
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
	  if (piIter->second->MEHandler()) weight *= piIter->second->Weight();
	  break;
	case Return_Value::New_Event  : return Return_Value::New_Event;
	case Return_Value::Nothing    : return Return_Value::Nothing;
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
  if (hit) {
    // enable shower generator independent FS QED correction to ME
    // TODO: check first, whether shower did FS QED
    bloblist->FindLast(btp::Shower)->AddStatus(blob_status::needs_extraQED);
    return Return_Value::Success;
  }
  return Return_Value::Nothing;
}

Return_Value::code Jet_Evolution::AttachShowers(Blob * blob,Blob_List * bloblist,
						Perturbative_Interface * interface) 
{
  if (!p_showerhandler->On()) {
    AftermathOfNoShower(blob,bloblist);
    return Return_Value::Nothing;
  }
  int shower(0);
  Return_Value::code stat(interface->DefineInitialConditions(blob));
  if (stat==Return_Value::New_Event) return stat;
  if (blob->Type()!=::btp::Hadron_Decay) {
    msg_Debugging()<<METHOD<<"(): Setting scale for MI {\n";
    double scale(0.0);
    Cluster_Amplitude *ampl(interface->Amplitude());
    while (ampl->Next()) ampl=ampl->Next();
    if (msg_LevelIsDebugging()) ampl->Print();
    for (size_t i(ampl->NIn());i<ampl->Legs().size();++i)
      if (ampl->Leg(i)->Flav().Strong()) 
	scale=Max(scale,ampl->Leg(i)->Mom().PPerp());
    if (scale==0.0) scale=(ampl->Leg(0)->Mom()+ampl->Leg(1)->Mom()).Mass();
    blob->AddData("MI_Scale",new Blob_Data<double>(scale));
    msg_Debugging()<<"} -> p_T = "<<scale<<"\n";
  }
  switch (stat) {
  case Return_Value::Success:
    if (blob->Type()!=::btp::Hadron_Decay) DefineInitialConditions(blob,bloblist);
    if (blob->NInP()==1) shower = interface->PerformDecayShowers();
    if (blob->NInP()==2) shower = interface->PerformShowers();
    switch (shower) {
    case 1: 
      Reset();
      AftermathOfSuccessfulShower(blob,bloblist,interface);    
      return Return_Value::Success;
    case 0:
      Reset();
      p_showerhandler->CleanUp();
      interface->CleanUp();
      return Return_Value::New_Event;
    default:
      THROW(fatal_error,"Invalid return value from shower");
    }
  case Return_Value::Nothing:
    AftermathOfNoShower(blob,bloblist);
    return Return_Value::Success;
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
    myblob->SetType(btp::Shower);
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
    myblob->SetType(btp::Shower);
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
  if (!p_showerhandler->On()) {
    if (blob->NInP()!=1)
    for (int i=0;i<2;i++) {
      // new ISR Blob
      myblob = new Blob();
      myblob->SetType(btp::Shower);
      myblob->SetStatus(blob_status::needs_beams);
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
    for (int i=0;i<blob->NOutP();i++) {
      myblob = new Blob();
      myblob->SetType(btp::Shower);
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
  msg_Debugging()<<METHOD<<"(): {\n";
  for (::Blob_List::const_iterator blit=bloblist->begin();
       blit!=bloblist->end();++blit) 
    if ((*blit)->Type()==::btp::Shower) {
      Update(*blit,0);
      Update(*blit,1);
    }
  msg_Debugging()<<"}\n";
  return true;
}

void Jet_Evolution::Update(const Blob *blob,const size_t beam) 
{ 
  size_t cbeam=0;
  for (int i=0;i<blob->NInP();++i) {
    const Particle *cur=blob->ConstInParticle(i);
    if (!cur->Flav().Strong() || cur->ProductionBlob()) continue;
    if (cbeam==beam) {
      msg_Debugging()<<"  "<<*cur<<"\n";
      p_showerhandler->GetISRHandler()->Extract(cur->Flav(),cur->Momentum()[0],beam);
      return;
    }
    ++cbeam;
  }
}

void Jet_Evolution::Finish(const string &) 
{
}

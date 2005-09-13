#include "Jet_Evolution.H"
#include "SimpleXS_Apacic_Interface.H"
#include "Amegic_Apacic_Interface.H"
#ifdef USING__Adicic    
#include "SimpleXS_Adicic_Interface.H"
#endif
#ifdef USING__CSS    
#include "SimpleXS_CSS_Interface.H"
#endif
#ifdef USING__Amisic    
#include "MI_Base.H"
#endif

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
  m_name      = std::string("Jet_Evolution:")+p_showerhandler->ShowerGenerator();
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
#ifdef USING__Adicic    
    if (meIter->second->Name()==string("SimpleXS") &&
      p_showerhandler->ShowerGenerator()==string("Adicic")) 
      interface = new SimpleXS_Adicic_Interface(meIter->second,p_showerhandler);
#endif
#ifdef USING__CSS
    if (meIter->second->Name()==string("SimpleXS") &&
      p_showerhandler->ShowerGenerator()==string("CSS")) 
      interface = new SimpleXS_CSS_Interface(meIter->second,p_showerhandler);
#endif
    if (interface!=NULL) m_interfaces.insert(make_pair(meIter->first,interface));
  }
}

Jet_Evolution::~Jet_Evolution() 
{ 
  while (m_interfaces.size()>0) {
    if (m_interfaces.begin()->second!=NULL) delete m_interfaces.begin()->second;
    m_interfaces.erase(m_interfaces.begin());
  }
}


bool Jet_Evolution::Treat(Blob_List * _bloblist, double & weight)
{
  PROFILE_LOCAL("Jet_Evolution::Treat");
  if (_bloblist->empty()) {
    msg.Error()<<"Potential error in Jet_Evolution::Treat."<<endl
	       <<"   Incoming blob list contains "<<_bloblist->size()<<" entries."<<endl
	       <<"   Continue and hope for the best."<<endl;
    return 0;
  }
  PertInterfaceIter piIter;
  std::string tag;
  bool found = 1;
  bool hit   = 0;
  Blob * blob;
  while (found) {
    found = 0;
    for (size_t i=0;i<_bloblist->size();++i) {
      blob = (*_bloblist)[i];
      if (blob->Status()==1 && blob->Type()==btp::Signal_Process) {
	FillDecayBlobMap(blob,_bloblist);
	piIter = m_interfaces.find(string("SignalMEs"));
	if (piIter==m_interfaces.end()) {
	  msg.Error()<<"Error in Jet_Evolution::Treat :"<<endl
		     <<"   No Perturbative_Interface found for type : SignalMEs."<<endl
		     <<"   Abort the run."<<endl;
	  abort();
	}	
	found   = AttachShowers(blob,_bloblist,piIter->second);
	weight *= piIter->second->Weight();
      }  
      if (blob->Status()==1 && blob->Type()==btp::Hard_Collision) {
	FillDecayBlobMap(blob,_bloblist);
	piIter = m_interfaces.find(string("MIMEs"));
	if (piIter==m_interfaces.end()) {
	  msg.Error()<<"Error in Jet_Evolution::Treat :"<<endl
		     <<"   No Perturbative_Interface found for type : MIMEs."<<endl
		     <<"   Abort the run."<<endl;
	  abort();
	}	
	found   = AttachShowers(blob,_bloblist,piIter->second);
	weight *= piIter->second->Weight();
      }  
      if ((blob->Status()==1 && blob->Type()==btp::Soft_Collision) ||
	  (blob->Status()==5 && blob->Type()==btp::Signal_Process)) {
	found = AttachShowers(blob,_bloblist,piIter->second);
      }  
      if (blob->Status()==1 && blob->Type()==btp::Hard_Decay) {
	FillDecayBlobMap(blob,_bloblist);
	piIter = m_interfaces.find(string("HardDecays"));
	if (piIter==m_interfaces.end()) {
	  msg.Error()<<"Error in Jet_Evolution::Treat :"<<endl
		     <<"   No Perturbative_Interface found for type : SignalMEs."<<endl
		     <<"   Abort the run."<<endl;
	  abort();
	}
	found   = AttachShowers(blob,_bloblist,piIter->second);
	weight *= piIter->second->Weight();
      } 
    }
    if (found) hit = 1;
    Reset();
  }
  return hit;
}

int Jet_Evolution::AttachShowers(Blob * _blob,Blob_List * _bloblist,
				 Perturbative_Interface * interface) 
{
  if ((_blob->Type()==btp::Hard_Collision && !p_showerhandler->ShowerMI()) ||
      _blob->Type()==btp::Soft_Collision || 
      (_blob->Type()==btp::Signal_Process && _blob->Status()==5)) {
    Blob * myblob;
    for (int i=0;i<2;i++) {
      myblob = new Blob();
      myblob->SetType(btp::IS_Shower);
      if (Sign(_blob->InParticle(i)->Momentum()[3])==1-2*i) myblob->SetBeam(i);
      else myblob->SetBeam(1-i);
      myblob->SetStatus(1);
      Particle * p = new Particle(*_blob->InParticle(i));
      p->SetStatus(2);
      myblob->AddToInParticles(p);
      myblob->AddToOutParticles(_blob->InParticle(i));
      _blob->InParticle(i)->SetStatus(2);
      myblob->SetId();
      _bloblist->insert(_bloblist->begin(),myblob);
    }
    for (int i=0;i<_blob->NOutP();i++) {
      myblob = new Blob();
      myblob->SetType(btp::FS_Shower);
      myblob->SetBeam(i);
      myblob->SetStatus(1);
      Particle * p = new Particle(*_blob->OutParticle(i));
      myblob->AddToInParticles(_blob->OutParticle(i));
      _blob->OutParticle(i)->SetStatus(2);
      myblob->AddToOutParticles(p);
      myblob->SetId();
      _bloblist->push_back(myblob);
    }
    _blob->SetStatus(0);
    return 1;
  }
  bool decayblob   = (_blob->NInP()==1);
  int shower, stat = interface->DefineInitialConditions(_blob);
  if (stat==3) {
    _blob->SetStatus(-1);
    p_showerhandler->CleanUp();
    return 1;
  }
  if (stat==1) {
    DefineInitialConditions(_blob,_bloblist);
    if (!decayblob) shower = interface->PerformShowers();
               else shower = interface->PerformDecayShowers();
    if (shower==1) {
      Blob * myblob;
      if (decayblob) _blob->InParticle(0)->SetInfo('h');
      interface->FillBlobs(_bloblist);
      p_showerhandler->FillBlobs(_bloblist); // BUG !!!!
      _blob->SetStatus(0);
      if ((!decayblob) && (!p_showerhandler->ISROn())) {
	for (int i=0;i<2;i++) {
	  // new ISR Blob
	  myblob = new Blob();
	  myblob->SetType(btp::IS_Shower);
	  if (Sign(_blob->InParticle(i)->Momentum()[3])==1-2*i) myblob->SetBeam(i);
	  else myblob->SetBeam(1-i);
	  myblob->SetStatus(1);
	  Particle * p = new Particle(*_blob->InParticle(i));
	  p->SetStatus(2);
	  myblob->AddToInParticles(p);
	  myblob->AddToOutParticles(_blob->InParticle(i));
	  _blob->InParticle(i)->SetStatus(2);
	  myblob->SetId();
	  _bloblist->insert(_bloblist->begin(),myblob);
	}
      }
      if (!(p_showerhandler->FSROn())) {
	for (int i=0;i<_blob->NOutP();i++) {
	  myblob = new Blob();
	  myblob->SetType(btp::FS_Shower);
	  myblob->SetBeam(i);
	  myblob->SetStatus(1);
	  Particle * p = new Particle(*_blob->OutParticle(i));
	  if (_blob->OutParticle(i)->DecayBlob()) {
	    Blob * dec  = _blob->OutParticle(i)->DecayBlob();
	    if (dec->Type()==btp::Hard_Decay) {
	      dec->RemoveInParticle(_blob->OutParticle(i));
	      dec->AddToInParticles(p);
	    }
	  }
	  myblob->AddToInParticles(_blob->OutParticle(i));
	  _blob->OutParticle(i)->SetStatus(2);
	  myblob->AddToOutParticles(p);
	  myblob->SetId();
	  _bloblist->push_back(myblob);
	}
      }
      else SetDecayBlobPointers(_blob,_bloblist);
    }
    else  if (shower==-1) {
      p_showerhandler->CleanUp();
      _blob->SetStatus(-1);
      // delete all meps blobs
      //interface->CleanBlobList(_bloblist,_blob->Type());
    }
    else {
      msg.Error()<<"Jet_Evolution::AttachShowers(..): Shower failure."<<std::endl;
      _blob->SetStatus(2);
      p_showerhandler->CleanUp();
    }
    return 1;
  }
  _blob->SetStatus(2);
  p_showerhandler->CleanUp();
  return 0;
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
      msg.Error()<<"ERROR in Jet_Evolution::SetDecayBlobPointers:"<<std::endl
		 <<"   Did not find particle in map of decay blobs."<<std::endl
		 <<"   Particle : "<<partin<<std::endl
		 <<"   Will abort the run."<<std::endl;
      abort();
    }
    dec     = pbiter->second;
    if (dec->Type()==btp::Hard_Decay) {
      //std::cout<<"Check for "<<partin->Flav()<<" / "<<partin->FinalMass()<<std::endl;
      partout = FollowUp(partin,dec);
      partout->SetDecayBlob(dec);
      dec->RemoveInParticle(partin);
      dec->AddToInParticles(partout);
      //std::cout<<"Match : "<<partin->Number()<<" -> "<<partout->Number()<<std::endl
      //	       <<(*dec)<<std::endl;
    }
  }
}

Particle * Jet_Evolution::FollowUp(Particle * partin,Blob * dec) 
{
  Blob * current = partin->DecayBlob();
  Particle * partout;
  if (current==0 || current==dec) return partin;
  //std::cout<<"Current : "<<current->Id()<<" <- "<<partin->Number()
  //	   <<" / "<<partin->Flav()<<std::endl;
  for (int i=0;i<current->NOutP();++i) {
    partout = current->OutParticle(i);
    if (partout->Flav()==partin->Flav() &&
	partout->FinalMass()==partin->FinalMass()) return FollowUp(partout,dec);
  }
  msg.Error()<<"ERROR in JetEvolution::FollowUp:"<<std::endl
	     <<"   Did not find a suitable particle to follow up decay initiator through blob list."<<std::endl
	     <<"   Particle = "<<partin<<std::endl
	     <<"   Will abort the run."<<std::endl;
  abort();
}

void Jet_Evolution::FillDecayBlobMap(Blob * blob,Blob_List * bloblist) 
{
  //  std::cout<<"Which blob ? "<<std::endl<<(*blob)<<std::endl;
  for (int i=0;i<blob->NOutP();++i) {
    if (blob->OutParticle(i)->DecayBlob()) 
      m_decmap.insert(std::make_pair(blob->OutParticle(i),
				     blob->OutParticle(i)->DecayBlob()));
  }
  //  std::cout<<"length of map : "<<m_decmap.size()<<std::endl;
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

bool Jet_Evolution::DefineInitialConditions(const ATOOLS::Blob *blob,
					    const ATOOLS::Blob_List *bloblist) 
{ 
  Reset();
  for (ATOOLS::Blob_List::const_iterator blit=bloblist->begin();
       blit!=bloblist->end();++blit) {
    if (((*blit)->Type()==ATOOLS::btp::Signal_Process ||
	 (*blit)->Type()==ATOOLS::btp::Hard_Collision) && *blit!=blob) {
      Update(*blit,0);
      Update(*blit,1);
    }
  }
  return true;
}

void Jet_Evolution::Update(const ATOOLS::Blob *blob,const size_t beam) 
{ 
  for (size_t i=0;i<(size_t)blob->NInP();++i) {
    const ATOOLS::Particle *cur=blob->ConstInParticle(i);
    if (i==beam || blob->NInP()<=1) {
      if (cur->ProductionBlob()!=NULL) Update(cur->ProductionBlob(),beam);
      else {
	p_showerhandler->GetISRHandler()->Extract(cur->Flav(),cur->Momentum()[0],beam);
	return;
      }
    }
  }
}

void Jet_Evolution::Finish(const std::string &) 
{
}

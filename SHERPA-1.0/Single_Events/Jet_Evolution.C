#include "Jet_Evolution.H"
#include "SimpleXS_Apacic_Interface.H"
#include "SimpleXS_Adicic_Interface.H"
#include "Amegic_Apacic_Interface.H"

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
  p_showerhandler(_showerhandler),
  m_initializedhardest(-1)
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
    if (meIter->second->Name()==string("SimpleXS") &&
	p_showerhandler->ShowerGenerator()==string("Adicic")) 
      interface = new SimpleXS_Adicic_Interface(meIter->second,p_showerhandler);
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
	piIter = m_interfaces.find(string("SignalMEs"));
	if (piIter==m_interfaces.end()) {
	  msg.Error()<<"Error in Jet_Evolution::Treat :"<<endl
		     <<"   No Perturbative_Interface found for type : SignalMEs."<<endl
		     <<"   Abort the run."<<endl;
	  abort();
	}	
	if (m_initializedhardest==-1) {
	  m_initializedhardest = piIter->second->DefineInitialConditions(blob);
	  piIter->second->FillBlobs(_bloblist);
	  return true;
	}
	else {
	  DefineInitialConditions(blob,_bloblist);
	}
	found   = AttachShowers(blob,_bloblist,piIter->second);
	weight *= piIter->second->Weight();
      }  
      if (blob->Status()==1 && blob->Type()==btp::Hard_Collision) {
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
      if (blob->Status()==1 && blob->Type()==btp::Hard_Decay) {
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
    ResetInterfaces();
  }
  return hit;
}

int Jet_Evolution::AttachShowers(Blob * _blob,Blob_List * _bloblist,
				  Perturbative_Interface * interface) 
{
  bool decayblob   = (_blob->NInP()==1);
  int shower, stat = m_initializedhardest;
  if (_blob->Type()!=ATOOLS::btp::Signal_Process) {
    stat = interface->DefineInitialConditions(_blob);
  }
  if (stat==3) {
    _blob->SetStatus(-1);
    p_showerhandler->CleanUp();
    return 1;
  }
  if (stat==1) {
    if (m_initializedhardest==-1) {
      interface->FillBlobs(_bloblist);
      DefineInitialConditions(_blob,_bloblist);
    }
    else m_initializedhardest=-1;
    if (!decayblob) shower = interface->PerformShowers();
               else shower = interface->PerformDecayShowers();  
    if (shower==1) {
      Blob * myblob;
      if (decayblob) _blob->InParticle(0)->SetInfo('h');
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
	  Particle * p = new Particle(_blob->InParticle(i));
	  p->SetStatus(2);
	  myblob->AddToInParticles(p);
	  myblob->AddToOutParticles(_blob->InParticle(i));
	  _blob->InParticle(i)->SetStatus(2);
	  myblob->SetId(_bloblist->size());
	  _bloblist->insert(_bloblist->begin(),myblob);
	}
      }
      if (!(p_showerhandler->FSROn())) {
	for (int i=0;i<_blob->NOutP();i++) {
	  myblob = new Blob();
	  myblob->SetType(btp::FS_Shower);
	  myblob->SetBeam(i);
	  myblob->SetStatus(1);
	  Particle * p = new Particle(_blob->OutParticle(i));
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
	  myblob->SetId(_bloblist->size());
	  _bloblist->push_back(myblob);
	}
      }
    }
    else  if (shower==3) {
      _blob->SetStatus(-1);
      p_showerhandler->CleanUp();
      // delete all meps blobs
      interface->CleanBlobList(_bloblist,_blob->Type());
    }
    else {
      _blob->SetStatus(2);
      p_showerhandler->CleanUp();
    }
    return 1;
  }
  _blob->SetStatus(2);
  p_showerhandler->CleanUp();
  return 0;
}

void Jet_Evolution::CleanUp() 
{ 
  p_showerhandler->CleanUp();
}

void Jet_Evolution::ResetInterfaces() 
{
  for (PertInterfaceMap::iterator piiter=m_interfaces.begin();
       piiter!=m_interfaces.end();++piiter) {
    piiter->second->Reset();
  }
  p_showerhandler->GetISRHandler()->Reset(0);
  p_showerhandler->GetISRHandler()->Reset(1);
}

bool Jet_Evolution::DefineInitialConditions(const ATOOLS::Blob *blob,
					    const ATOOLS::Blob_List *bloblist) 
{ 
  ResetInterfaces();
  // fetch maximum remnant energy from mi handler
  p_showerhandler->GetISRHandler()->Extract(ATOOLS::Flavour(ATOOLS::kf::none),2.,0);
  p_showerhandler->GetISRHandler()->Extract(ATOOLS::Flavour(ATOOLS::kf::none),2.,1);
  for (ATOOLS::Blob_List::const_iterator blit=bloblist->begin();
       blit!=bloblist->end();++blit) {
    if (((*blit)->Type()==ATOOLS::btp::Signal_Process ||
	 (*blit)->Type()==ATOOLS::btp::Hard_Collision) && *blit!=blob) {
      UpdateInterfaces(*blit,0);
      UpdateInterfaces(*blit,1);
    }
  }
  return true;
}

void Jet_Evolution::UpdateInterfaces(const ATOOLS::Blob *blob,const size_t beam) 
{ 
  for (size_t i=0;i<(size_t)blob->NInP();++i) {
    const ATOOLS::Particle *cur=blob->ConstInParticle(i);
    if (i==beam || blob->NInP()<=1) {
      if (cur->ProductionBlob()!=NULL) {
	UpdateInterfaces(cur->ProductionBlob(),beam);
      }
      else {
	for (PertInterfaceMap::iterator piiter=m_interfaces.begin();
	     piiter!=m_interfaces.end();++piiter) {
	  piiter->second->UpdateEnergy(cur->Momentum()[0],i);
	}
	p_showerhandler->GetISRHandler()->Extract(cur->Flav(),cur->Momentum()[0],beam);
	return;
      }
    }
  }
}

void Jet_Evolution::Finish(const std::string &) 
{
}

#include "Jet_Evolution.H"
#include "SimpleXS_Apacic_Interface.H"
#include "Amegic_Apacic_Interface.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Jet_Evolution::Jet_Evolution(MEHandlersMap *_mehandlers,Shower_Handler *_showerhandler) :
  p_showerhandler(_showerhandler)
{
  m_name      = string("Jet_Evolution : ")+p_showerhandler->ShowerGenerator();
  m_type      = string("Perturbative");

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
    m_interfaces.insert(make_pair(meIter->first,interface));
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
  bool take;
  size_t pos;
  Blob * blob;
  while (found) {
  int b=0;
    found = 0;
    for (size_t i=0;i<_bloblist->size();++i) {
      blob = (*_bloblist)[i];
      // Search for active blobs of type "Signal Process :"
      pos = blob->Type().find(std::string("Signal Process :"));
      if (blob->Status()==1 && pos!=std::string::npos) {
	piIter = m_interfaces.find(string("SignalMEs"));
	if (piIter==m_interfaces.end()) {
	  msg.Error()<<"Error in Jet_Evolution::Treat :"<<endl
		     <<"   No Perturbative_Interface found for type : SignalMEs."<<endl
		     <<"   Abort the run."<<endl;
	  abort();
	}	
	found   = AttachShowers(blob,_bloblist,piIter->second);
	weight *= piIter->second->GetWeight();
      }  

      // Search for active blobs of type "Hard Subprocess :"
      pos = blob->Type().find(std::string("Hard Subprocess :"));
      if (blob->Status()==1 && pos!=std::string::npos) {
	piIter = m_interfaces.find(string("MIMEs"));
	if (piIter==m_interfaces.end()) {
	  msg.Error()<<"Error in Jet_Evolution::Treat :"<<endl
		     <<"   No Perturbative_Interface found for type : MIMEs."<<endl
		     <<"   Abort the run."<<endl;
	  abort();
	}	
	found   = AttachShowers(blob,_bloblist,piIter->second);
	weight *= piIter->second->GetWeight();
      }  

      // Search for active blobs of type "Hard Decay :"
      pos = blob->Type().find(std::string("Hard decay :"));
      if (blob->Status()==1 && pos!=std::string::npos) {
	piIter = m_interfaces.find(string("HardDecays"));
	if (piIter==m_interfaces.end()) {
	  msg.Error()<<"Error in Jet_Evolution::Treat :"<<endl
		     <<"   No Perturbative_Interface found for type : SignalMEs."<<endl
		     <<"   Abort the run."<<endl;
	  abort();
	}
	found   = AttachShowers(blob,_bloblist,piIter->second);
	weight *= piIter->second->GetWeight();
      } 
    }
    if (found) hit = 1;
  }
  return hit;
}

bool Jet_Evolution::AttachShowers(Blob * _blob,Blob_List * _bloblist,
				  Perturbative_Interface * interface) 
{
  size_t colon= _blob->Type().find(":");
  std::string type = _blob->Type().substr(0,colon+2);
  bool decayblob   = (_blob->NInP()==1);
  int shower,stat = interface->DefineInitialConditions(_blob);
  if (stat==3) {
    _blob->SetType(type);
    _blob->SetStatus(-1);
    p_showerhandler->CleanUp();
    return true;
  }
  if (stat) {
    interface->FillBlobs(_bloblist);
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
	  myblob->SetType(string("IS Shower (none)"));
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
	  myblob->SetType(string("FS Shower (none)"));
	  myblob->SetBeam(i);
	  myblob->SetStatus(1);
	  Particle * p = new Particle(_blob->OutParticle(i));
	  if (_blob->OutParticle(i)->DecayBlob()) {
	    Blob * dec  = _blob->OutParticle(i)->DecayBlob();
	    size_t pos = dec->Type().find(std::string("Hard decay"));
	    if (pos!=std::string::npos) {
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
    else if (shower==3) {
      _blob->SetType(type);
      _blob->SetStatus(-1);
      p_showerhandler->CleanUp();
      // delete all meps blobs
      interface->CleanBlobs(_bloblist);
    }
    else {
      _blob->SetType(type);
      _blob->SetStatus(0);
      p_showerhandler->CleanUp();
    }
    return true;
  }
  _blob->SetType(type);
  _blob->SetStatus(0);
  p_showerhandler->CleanUp();
  return false;
}

void Jet_Evolution::CleanUp() 
{ 
  p_showerhandler->CleanUp();
}



#include "Jet_Evolution.H"
#include "SimpleXS_Apacic_Interface.H"
#include "Amegic_Apacic_Interface.H"


using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Jet_Evolution::Jet_Evolution(MEHandlersMap * _mehandlers,
			     Shower_Handler * _showerhandler) :
  p_showerhandler(_showerhandler)
{
  m_name      = string("Jet_Evolution : ")+p_showerhandler->ShowerGenerator();
  m_type      = string("Perturbative");

  Perturbative_Interface * interface;
  MEHandlerIter            meIter;
  for (meIter=_mehandlers->begin();meIter!=_mehandlers->end();++meIter) {
    if (meIter->second->Name()==string("Amegic") &&
	p_showerhandler->ShowerGenerator()==string("Apacic")) 
      interface = new Amegic_Apacic_Interface(meIter->second,p_showerhandler);
    if (meIter->second->Name()==string("Simple X-section") &&
	p_showerhandler->ShowerGenerator()==string("Apacic")) 
      interface = new SimpleXS_Apacic_Interface(meIter->second,p_showerhandler);
    m_interfaces.insert(make_pair(meIter->first,interface));
  }
}

Jet_Evolution::~Jet_Evolution() 
{ 
  m_interfaces.clear();
}


bool Jet_Evolution::Treat(Blob_List * _bloblist, double & weight)
{
  if (_bloblist->empty()) {
    msg.Error()<<"Potential error in Jet_Evolution::Treat."<<endl
	       <<"   Incoming blob list contains "<<_bloblist->size()<<" entries."<<endl
	       <<"   Continue and hope for the best."<<endl;
    return 0;
  }
  
  PertInterfaceIter        piIter;

  bool hit = 0,found = 1;
  int  pos;
  while (found) {
    found = 0;
    for (Blob_Iterator blit=_bloblist->begin();blit!=_bloblist->end();++blit) {
      pos = (*blit)->Type().find(string("Signal Process :"));
      if ((*blit)->Status()==1 && pos>-1) {
	piIter = m_interfaces.find(string("SignalMEs"));
	if (piIter==m_interfaces.end()) {
	  msg.Error()<<"Error in Jet_Evolution::Treat :"<<endl
		     <<"   No Perturbative_Interface found for type : SignalMEs."<<endl
		     <<"   Abort the run."<<endl;
	  abort();
	}	
	found   = AttachShowers((*blit),_bloblist,piIter->second);
	weight *= piIter->second->GetWeight();
      }  // Search for active blobs of type "Signal Process :"
      pos = (*blit)->Type().find(string("Hard decay :"));
      if ((*blit)->Status()==1 && pos>-1) {
	piIter = m_interfaces.find(string("HardDecays"));
	if (piIter==m_interfaces.end()) {
	  msg.Error()<<"Error in Jet_Evolution::Treat :"<<endl
		     <<"   No Perturbative_Interface found for type : SignalMEs."<<endl
		     <<"   Abort the run."<<endl;
	  abort();
	}
	found   = AttachShowers((*blit),_bloblist,piIter->second);
	weight *= piIter->second->GetWeight();
      } // Search for active blobs of type "Hard Decay :"
    }
    if (found) hit = 1;
  }
  return hit;
}

bool Jet_Evolution::AttachShowers(Blob * _blob,Blob_List * _bloblist,
				  Perturbative_Interface * interface) 
{
  std::string type = _blob->Type();
  bool decayblob   = (_blob->NInP()==1);
  Matrix_Element_Handler * mehandler = interface->GetMEHandler();

  int shower,stat = interface->DefineInitialConditions(_blob);
  if (stat) {
    interface->FillBlobs(_bloblist);
    if (!decayblob) {
      shower = p_showerhandler->PerformShowers(p_showerhandler->MaxJetNumber()!=mehandler->Nout(),
					       mehandler->GetISR_Handler()->X1(),
					       mehandler->GetISR_Handler()->X2());
    }
    else {
      shower = p_showerhandler->PerformDecayShowers(0);
    }      

    if (shower==1) {
      Blob * myblob;
      p_showerhandler->FillBlobs(_bloblist); // BUG !!!!
      _blob->SetStatus(0);
      if ((!decayblob) && (!p_showerhandler->ISROn())) {
	for (int i=0;i<2;i++) {
	  // new ISR Blob
	  myblob = new Blob();
	  myblob->SetType(string("IS Shower (none)"));
	  myblob->SetBeam(i);
	  myblob->SetStatus(1);
	  Particle * p = new Particle(_blob->InParticle(i));
	  p->SetProductionBlob(NULL);
	  p->SetDecayBlob(myblob);
	  p->SetStatus(2);
	  myblob->AddToInParticles(p);
	  myblob->AddToOutParticles(_blob->InParticle(i));
	  _blob->InParticle(i)->SetProductionBlob(myblob);
	  _blob->InParticle(i)->SetStatus(2);
	  myblob->SetId(_bloblist->size());
	  _bloblist->insert(_bloblist->begin(),myblob);
	}
      }
      if (!(p_showerhandler->FSROn())) {
	for (int i=0;i<_blob->NOutP();i++) {
	  myblob = new Blob();
	  myblob->SetType(string("FS Shower (none)"));
	  myblob->SetBeam(-1);
	  myblob->SetStatus(1);
	  myblob->SetId(_bloblist->size());
	  myblob->AddToInParticles(_blob->OutParticle(i));
	  Particle * p = new Particle(_blob->OutParticle(i));
	  p->SetProductionBlob(myblob);
	  p->SetDecayBlob(_blob->OutParticle(i)->DecayBlob());
	  _blob->OutParticle(i)->SetDecayBlob(myblob);
	  _blob->OutParticle(i)->SetStatus(2);
	  myblob->AddToOutParticles(p);
	  _bloblist->push_back(myblob);
	}
      }
    }
    else if (shower==3) {
      _blob->SetType(type);
      _blob->SetStatus(-1);
      p_showerhandler->CleanUp();
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



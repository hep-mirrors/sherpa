#include "Jet_Evolution.H"
#include "SimpleXS_Apacic_Interface.H"
#include "Amegic_Apacic_Interface.H"


using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Jet_Evolution::Jet_Evolution(Matrix_Element_Handler * _mehandler,
			     Shower_Handler * _showerhandler) :
  p_mehandler(_mehandler), p_showerhandler(_showerhandler), p_interface(NULL)
{
  m_name      = string("Jet_Evolution : ")+p_showerhandler->ShowerGenerator();
  m_type      = string("Perturbative");

  if (p_mehandler->SignalGenerator()==string("Amegic") &&
      p_showerhandler->ShowerGenerator()==string("Apacic"))
    p_interface = new Amegic_Apacic_Interface(p_mehandler,p_showerhandler);
  if (p_mehandler->SignalGenerator()==string("Internal") &&
      p_showerhandler->ShowerGenerator()==string("Apacic")) 
    p_interface = new SimpleXS_Apacic_Interface(p_mehandler,p_showerhandler);
}

Jet_Evolution::~Jet_Evolution() 
{ 
  if (p_interface) delete p_interface;
  // me_handler and shower_handler are deleted in initialization_handler
}


bool Jet_Evolution::Treat(Blob_List * _bloblist, double & weight)
{
  if (_bloblist->empty()) {
    msg.Error()<<"Potential error in Jet_Evolution::Treat."<<endl
	       <<"   Incoming blob list contains "<<_bloblist->size()<<" entries."<<endl
	       <<"   Continue and hope for the best."<<endl;
    return 0;
  }
  
  Blob * myblob;
  bool found = 1;
  bool hit   = 0;
  int  pos,shower;
  while (found) {
    found = 0;
    for (Blob_Iterator blit=_bloblist->begin();blit!=_bloblist->end();++blit) {
      pos = (*blit)->Type().find(string("Signal Process :"));
      if ((*blit)->Status()==1 && pos>-1) {
	myblob = (*blit);
	int stat = p_interface->DefineInitialConditions(myblob);
	weight  *= p_interface->GetWeight();
	if (stat) {
	  p_interface->FillBlobs(_bloblist);
	  shower = p_showerhandler->PerformShowers(p_showerhandler->MaxJetNumber()!=p_mehandler->Nout(),
						   p_mehandler->GetISR_Handler()->X1(),p_mehandler->GetISR_Handler()->X2());
	  if (shower==1) {
	    p_showerhandler->FillBlobs(_bloblist);
	    (*blit)->SetStatus(0);
	    if (!(p_showerhandler->ISROn())) {
	      for (int i=0;i<2;i++) {
		// new ISR Blob
		myblob = new Blob();
		myblob->SetType(string("IS Shower (none)"));
		myblob->SetBeam(i);
		myblob->SetStatus(1);
 		Particle * p = new Particle((*blit)->InParticle(i));
		p->SetProductionBlob(NULL);
		p->SetDecayBlob(myblob);
		p->SetStatus(2);
		myblob->AddToInParticles(p);
		myblob->AddToOutParticles((*blit)->InParticle(i));
		(*blit)->InParticle(i)->SetProductionBlob(myblob);
		(*blit)->InParticle(i)->SetStatus(2);
		myblob->SetId(_bloblist->size());
		_bloblist->insert(_bloblist->begin(),myblob);
	      }
	    }
	    if (!(p_showerhandler->FSROn())) {
	      for (int i=0;i<(*blit)->NOutP();i++) {
		myblob = new Blob();
		myblob->SetType(string("FS Shower (none)"));
		myblob->SetBeam(i);
		myblob->SetStatus(1);
 		Particle * p = new Particle((*blit)->OutParticle(i));
		p->SetProductionBlob(myblob);
		p->SetDecayBlob(NULL);
		myblob->AddToInParticles((*blit)->OutParticle(i));
		(*blit)->OutParticle(i)->SetDecayBlob(myblob);
		(*blit)->OutParticle(i)->SetStatus(2);
		myblob->AddToOutParticles(p);
		myblob->SetId(_bloblist->size());
		_bloblist->push_back(myblob);
	      }
	    }
	  }
	  else if (shower==3) {
	    myblob->SetType(string("Signal Process : "));
	    myblob->SetStatus(-1);
	    p_showerhandler->CleanUp();
	  }
	  else {
	    myblob->SetType(string("Signal Process : "));
	    myblob->SetStatus(0);
	    p_showerhandler->CleanUp();
	  }
	  found = hit = 1;
	}
	else {
	  myblob->SetType(string("Signal Process : "));
	  myblob->SetStatus(0);
	  p_showerhandler->CleanUp();
	}
      }
    }
  }
  return hit;
}

void Jet_Evolution::CleanUp() 
{ 
  p_showerhandler->CleanUp();
}



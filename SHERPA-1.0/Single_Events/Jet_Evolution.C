#include "Jet_Evolution.H"
#include "SimpleXS_Apacic_Interface.H"
#include "Amegic_Apacic_Interface.H"


using namespace SHERPA;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace AORGTOOLS;
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

Jet_Evolution::~Jet_Evolution() { }


bool Jet_Evolution::Treat(Blob_List * _bloblist)
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
      cout<<"Still in loop "<<(*blit)->Id()<<endl;
      pos = (*blit)->Type().find(string("Signal Process :"));
      msg.Debugging()<<"Found blob to deal with "<<(*blit)->Type()<<" "<<pos<<" "<<(*blit)->Status()
		     <<"  ("<<_bloblist->size()<<")"<<endl; 
      if ((*blit)->Status()==1 && pos>-1) {
	msg.Debugging()<<"Found blob to deal with "<<(*blit)<<" -> "<<_bloblist->size()<<endl<<(*blit)<<endl; 
	myblob = (*blit);
	int stat = p_interface->DefineInitialConditions(myblob);
	if (stat) {
	  cout<<" jetmax = "<<p_showerhandler->MaxJetNumber()<<"  nout="<<p_mehandler->Nout()<<endl;
	  shower = p_showerhandler->PerformShowers(p_showerhandler->MaxJetNumber()!=p_mehandler->Nout());
	  if (shower==1) {
	    p_showerhandler->FillBlobs(_bloblist);
	    (*blit)->SetStatus(0);
	    if (!(p_showerhandler->ISROn())) {
	      for (int i=0;i<2;i++) {
		myblob = new Blob();
		myblob->SetType(string("IS Shower (none)"));
		myblob->SetBeam(i);
		myblob->SetStatus(1);
		myblob->AddToInPartons((*blit)->InParton(i));
		myblob->AddToOutPartons((*blit)->InParton(i));
		myblob->SetId(_bloblist->size());
		blit = _bloblist->insert(_bloblist->begin(),myblob);
	      }
	    }
	    if (!(p_showerhandler->FSROn())) {
	      for (int i=0;i<(*blit)->NOutP();i++) {
		myblob = new Blob();
		myblob->SetType(string("FS Shower (none)"));
		myblob->SetBeam(i);
		myblob->SetStatus(1);
		myblob->AddToInPartons((*blit)->InParton(i));
		myblob->AddToOutPartons((*blit)->InParton(i));
		myblob->SetId(_bloblist->size());
		_bloblist->push_back(myblob);
	      }
	    }
	  }
	  else if (shower==3) {
	    msg.Tracking()<<" ASK for SAME EVENT "<<endl;
	    myblob->SetType(string("Signal Process : "));
	    myblob->SetStatus(-1);
	    p_showerhandler->CleanUp();
	  }
	  else {
	    msg.Tracking()<<" ASK for NEW EVENT (error)"<<endl;
	    myblob->SetType(string("Signal Process : "));
	    myblob->SetStatus(0);
	    p_showerhandler->CleanUp();
	  }
	  found = hit = 1;
	}
	else {
	  // probably sudakov rejection
	  msg.Tracking()<<" ASK for NEW EVENT (sudakov) "<<endl;
	  myblob->SetType(string("Signal Process : "));
	  myblob->SetStatus(0);
	  p_showerhandler->CleanUp();
	}
      }
      cout<<"Else for found blob and status/pos are good ..."<<endl;
    }
  }
  cout<<"Return "<<hit<<endl;
  return hit;
}

void Jet_Evolution::CleanUp() 
{ 
  p_showerhandler->CleanUp();
}



#include "Signal_Processes.H"


using namespace SHERPA;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace AORGTOOLS;
using namespace std;

Signal_Processes::Signal_Processes(Matrix_Element_Handler * _mehandler) :
  p_mehandler(_mehandler)
{
  m_name      = string("Signal_Processes : ")+p_mehandler->Name();
  m_type      = string("Perturbative");
}

Signal_Processes::~Signal_Processes()
{
  if (p_mehandler) { delete p_mehandler; p_mehandler = NULL; }
}


bool Signal_Processes::Treat(Blob_List * _bloblist, double & weight)
{
  if (_bloblist->size()>1) return 0;
  if (_bloblist->empty()) {
    msg.Error()<<"Potential error in Signal_Processes::Treat."<<endl
	       <<"   Incoming blob list contains "<<_bloblist->size()<<" entries."<<endl
	       <<"   Continue and hope for the best."<<endl;
    return 0;
  }
  
  Blob * myblob;
  bool found = 1;
  bool hit   = 0;
  
  while (found) {
    found = 0;
    for (Blob_Iterator blit=_bloblist->begin();blit!=_bloblist->end();++blit) {
      if (((*blit)->Type()==string("Signal Process : ")) && ((*blit)->Status()==0)) {
	myblob = (*blit);
	found  = 1;
	if (p_mehandler->GenerateOneEvent()) {
	  FillBlob(myblob);
	  weight=p_mehandler->Weight();
	  hit = 1;
	}
      }
      else if (((*blit)->Type()==string("Signal Process : ")) && ((*blit)->Status()==-1)) {
	myblob = (*blit);
	found  = 1;
	if (p_mehandler->GenerateSameEvent()) {
	  FillBlob(myblob);
	  weight=p_mehandler->Weight();
	  hit = 1;
	}
      }
    }
  }
  return hit;
}

void Signal_Processes::CleanUp() { return; }

void Signal_Processes::FillBlob(Blob * _blob)
{
  _blob->SetPosition(Vec4D(0.,0.,0.,0.));
  _blob->SetType(_blob->Type()+p_mehandler->ProcessName());
  _blob->SetStatus(1);

  Vec4D cms = Vec4D(0.,0.,0.,0.);
  for (int i=0;i<p_mehandler->Nin();i++) cms += p_mehandler->Momenta()[i];
  _blob->SetCMS(cms);
  _blob->SetBeam(-1);

  // make sure that blob is empty
  _blob->DeleteOwnedPartons();

  Parton * parton;
  for (int i=0;i<p_mehandler->Nin();i++) {
    parton = new Parton(i,p_mehandler->Flavs()[i],p_mehandler->Momenta()[i]);
    parton->SetNumber(int(parton));
    parton->SetDecayBlob(_blob);
    parton->SetStatus(2);
    parton->SetInfo('G');
    _blob->AddToInPartons(parton);
  }
  for (int i=p_mehandler->Nin();i<p_mehandler->Nin()+p_mehandler->Nout();i++) {
    parton = new Parton(i,p_mehandler->Flavs()[i],p_mehandler->Momenta()[i]);
    parton->SetProductionBlob(_blob);
    parton->SetStatus(1);
    parton->SetInfo('H');
    _blob->AddToOutPartons(parton);
  }
}

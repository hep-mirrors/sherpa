#include "Hard_Decays.H"
#include "Message.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Hard_Decays::Hard_Decays(Hard_Decay_Handler * _dechandler) :
  p_dechandler(_dechandler)
{
  m_name      = string("Hard_Decays : ")+p_dechandler->Name();
  m_type      = string("Perturbative");
}

Hard_Decays::~Hard_Decays() 
{
}

bool Hard_Decays::Treat(ATOOLS::Blob_List * _bloblist, double & weight) 
{
  if (_bloblist->empty()) {
    msg.Error()<<"Potential error in Hard_Decays::Treat."<<endl
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
      if ((*blit)->Status()==1) {
	myblob = (*blit);
	// found  = 1;
	msg.Out()<<"Found a blob for decay : "<<endl<<myblob<<endl;
      }
    }
  }
  return hit;
}

void Hard_Decays::CleanUp() { return; }

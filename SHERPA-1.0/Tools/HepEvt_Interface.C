#include "HepEvt_Interface.H"
#include "Message.H"

using namespace SHERPA;
using namespace ATOOLS;


extern "C" {
  void inhepevt_(int &,int &,int *,int *,int *,int *,double *,double *);
  void outhepevt_();
}


HepEvt_Interface::HepEvt_Interface() : m_evtnumber(0) {}


HepEvt_Interface::~HepEvt_Interface() {}


void HepEvt_Interface::Sherpa2HepEvt(Blob_List * _blobs) {
  phep     = new double[5*maxentries];
  vhep     = new double[4*maxentries];
  jmohep   = new int[2*maxentries];
  jdahep   = new int[2*maxentries];
  isthep   = new int[maxentries];
  idhep    = new int[maxentries];

  m_evtnumber++;
  int nhep = 0;
 
  ISBlobs2HepEvt(_blobs,nhep);
  HardBlob2HepEvt(_blobs,nhep);
  FSBlobs2HepEvt(_blobs,nhep);
  FragmentationBlob2HepEvt(_blobs,nhep);
  HadronDecayBlobs2HepEvt(_blobs,nhep);

  inhepevt_(m_evtnumber,nhep,isthep,idhep,jmohep,jdahep,phep,vhep);

  delete [] jmohep;
  delete [] jdahep;
  delete [] isthep;
  delete [] idhep; 
  delete [] phep;  
  delete [] vhep;  
}

void HepEvt_Interface::ISBlobs2HepEvt(Blob_List * _blobs,int & _nhep) {
  int pos;
  for (int beam=0;beam<2;beam++) {
    for (Blob_List::const_iterator bit=_blobs->begin(); bit!=_blobs->end();++bit) {
      if ((*bit)->Type()==std::string("Bunch") && (*bit)->Beam()==beam) {
	if ((*bit)->NInP()!=1) {
	  msg.Error()<<"Error in HepEvt_Interface::ISBlobs2HepEvt."<<endl
		     <<"   Bunch blob with more than one incoming particle !"<<endl
		     <<(*bit)<<endl;
	  abort();
	}
	if ((*bit)->NOutP()>1) {
	  Parton2HepEvt((*bit)->InParton(0),_nhep);
	  for (int j=0;j<(*bit)->NOutP();j++) Parton2HepEvt((*bit)->OutParton(j),_nhep);
	  EstablishRelations((*bit));
	}
      }
      if ((*bit)->Type()==std::string("Beam Remnant") && (*bit)->Beam()==beam) {
	if ((*bit)->NInP()!=1) {
	  msg.Error()<<"Error in HepEvt_Interface::ISBlobs2HepEvt."<<endl
		     <<"   Beam Remnant blob with more than one incoming particle !"<<endl
		     <<(*bit)<<endl;
	  abort();
	}
	if ((*bit)->NOutP()>1) {
	  Parton2HepEvt((*bit)->InParton(0),_nhep);
	  for (int j=0;j<(*bit)->NOutP();j++) Parton2HepEvt((*bit)->OutParton(j),_nhep);
	  EstablishRelations((*bit));
	}
      }
      pos = (*bit)->Type().find(string("IS"));
      if ((pos>-1) && ((*bit)->Beam()==beam)) {
	if ((*bit)->NInP()!=1) {
	  msg.Error()<<"Error in HepEvt_Interface::ISBlobs2HepEvt."<<endl
		     <<"   IS blob with more than one incoming particle !"<<endl
		     <<(*bit)<<endl;
	  abort();
	}
	if ((*bit)->NOutP()>1) {
	  Parton2HepEvt((*bit)->InParton(0),_nhep);
	  for (int j=0;j<(*bit)->NOutP();j++) Parton2HepEvt((*bit)->OutParton(j),_nhep);
	  EstablishRelations((*bit));
	}
	else {
	  Parton2HepEvt((*bit)->InParton(0),_nhep);
	  for (int j=0;j<(*bit)->NOutP();j++) Parton2HepEvt((*bit)->OutParton(j),_nhep);
	  EstablishRelations((*bit));
	}
      }
    }
  }
}


void HepEvt_Interface::HardBlob2HepEvt(Blob_List * _blobs,int & _nhep) {
  int pos;
  for (Blob_List::const_iterator bit=_blobs->begin(); bit!=_blobs->end();++bit) {
    pos = (*bit)->Type().find(string("Signal Process :")); 
    if (pos>-1) {
      if ((*bit)->NInP()!=2) {
	msg.Error()<<"Error in HepEvt_Interface::HardBlob2HepEvt."<<endl
		   <<"   Hard ME blob with other than 2 incoming particles !"<<endl
		   <<(*bit)<<endl;
	abort();
      }
      if ((*bit)->NOutP()>=2) {
	Parton2HepEvt((*bit)->InParton(0),_nhep);
	Parton2HepEvt((*bit)->InParton(1),_nhep);
	for (int j=0;j<(*bit)->NOutP();j++) Parton2HepEvt((*bit)->OutParton(j),_nhep);
	EstablishRelations((*bit));
      }
    }
  }
}

void HepEvt_Interface::FSBlobs2HepEvt(Blob_List * _blobs,int & _nhep) {
  Blob   * interface;
  Parton * seed, *compare;
  Flavour  flav;
  int      cols[2],number;
  for (Blob_List::const_iterator bit=_blobs->begin(); bit!=_blobs->end();++bit) {
    if ((*bit)->Type()==std::string("ME PS Interface (Sherpa, FS)")) {
      interface=(*bit);
      for (int i=0;i<interface->NInP();i++) {
	seed    = interface->InParton(i);
	flav    = seed->Flav();
	cols[0] = seed->GetFlow(1);
	cols[1] = seed->GetFlow(2);
	number  = seed->Number();
	for (Blob_List::const_iterator bit=_blobs->begin(); bit!=_blobs->end();++bit) {
	  if ((*bit)->Type()==std::string("FS Shower (APACIC++2.0)") && (*bit)->NInP()==1) {
	    compare = (*bit)->InParton(0);
	    if (compare->Flav()==flav && compare->Number()==number &&
		compare->GetFlow(1)==cols[0] && compare->GetFlow(2)==cols[1]) {
	      for (int j=0;j<(*bit)->NOutP();j++) Parton2HepEvt((*bit)->OutParton(j),_nhep);
	      EstablishRelations(seed,(*bit));
	    }
	  } 
	}
      }
    }
  }
}

void HepEvt_Interface::FragmentationBlob2HepEvt(Blob_List * _blobs,int & _nhep) {
  Blob * fragmentation;
  for (Blob_List::const_iterator bit=_blobs->begin(); bit!=_blobs->end();++bit) {
    if ((*bit)->Type()==std::string("Fragmentation (Lund : Pythia 6.163)")) {
      String2HepEvt((*bit),_nhep);;
    }
  }
}

void HepEvt_Interface::HadronDecayBlobs2HepEvt(Blob_List * _blobs,int & _nhep) {
  for (Blob_List::const_iterator bit=_blobs->begin(); bit!=_blobs->end();++bit) {
    if ((*bit)->Type()==std::string("Hadron decay")) {
      if ((*bit)->NInP()!=1) {
	msg.Error()<<"Error in HepEvt_Interface::HadronDecays2HepEvt."<<endl
		   <<"   Decay blob with other than 1 incoming particles !"<<endl
		   <<(*bit)<<endl;
	abort();
      }
      if ((*bit)->NOutP()>=2) {
	Parton2HepEvt((*bit)->InParton(0),_nhep);
	for (int j=0;j<(*bit)->NOutP();j++) Parton2HepEvt((*bit)->OutParton(j),_nhep);
	EstablishRelations((*bit));
      }      
    }
  }
}

void HepEvt_Interface::Parton2HepEvt(Parton * _part,int & _nhep)
{
  int number = m_connect.count(_part);
  if (number>0) return;
  idhep[_nhep]    = _part->Flav().HepEvt();
  jmohep[2*_nhep] = jmohep[2*_nhep+1] = 0; 
  jdahep[2*_nhep] = jdahep[2*_nhep+1] = 0; 
        
  for (short int j=1; j<4; ++j) phep[(j-1)+_nhep*5] = _part->Momentum()[j];
  phep[3+_nhep*5] = _part->Momentum()[0];
  double pabs = (_part->Momentum()).Abs2();
  if (pabs<0) phep[4+_nhep*5] = 0.;
         else phep[4+_nhep*5] = sqrt(pabs);
  if (_part->ProductionBlob()!=NULL) {
    for (short int j=1; j<4; ++j) vhep[(j-1)+_nhep*4] = _part->XProd()[j];
    vhep[3+_nhep*4]     = _part->XProd()[0];
  }
  else {
    for (short int j=1; j<4; ++j) vhep[(j-1)+_nhep*4] = 0.;
    vhep[3+_nhep*4]     = 0.;
  }
  if (_part->DecayBlob()!=NULL) isthep[_nhep] = 2;
                           else isthep[_nhep] = 1;

  m_connect.insert(std::make_pair(_part,_nhep));
  _nhep++;
}

void HepEvt_Interface::String2HepEvt(Blob * _string,int & _nhep)
{
  idhep[_nhep]    = 92;
  jmohep[2*_nhep] = jmohep[2*_nhep+1] = 0; 
  jdahep[2*_nhep] = jdahep[2*_nhep+1] = 0; 
  for (short int j=1; j<4; ++j) phep[(j-1)+_nhep*5] = _string->CMS()[j];
  phep[3+_nhep*5] = _string->CMS()[0];
  double pabs = (_string->CMS()).Abs2();
  if (pabs<0) phep[4+_nhep*5] = 0.;
         else phep[4+_nhep*5] = sqrt(pabs);
  for (short int j=1;j<4;j++) vhep[(j-1)+_nhep*4] = _string->Position()[j];
  vhep[3+_nhep*4] = _string->Position()[0];
  isthep[_nhep]   = 2;

  Parton * incoming, * outgoing;
  int number, stringnumber = _nhep;
  for (int i=0;i<_string->NInP();i++) {
    incoming = _string->InParton(i);
    number   = m_connect[incoming];
    if (i==0)                 jmohep[2*_nhep]   = number+1;
    if (i==_string->NInP()-1) jmohep[2*_nhep+1] = number+1;
    jdahep[2*number] = jdahep[2*number+1]       = _nhep+1; 
  }
  
  _nhep++;
  for (int i=0;i<_string->NOutP();i++) {
    outgoing = _string->OutParton(i);
    Parton2HepEvt(outgoing,_nhep);
    if (i==0)                  jdahep[2*stringnumber]   = _nhep;
    if (i==_string->NOutP()-1) jdahep[2*stringnumber+1] = _nhep;
    jmohep[2*m_connect[outgoing]] = jmohep[2*m_connect[outgoing]+1] = stringnumber+1;
  }
}



void HepEvt_Interface::EstablishRelations(Blob * _blob) {
  int mothers[2];
  int daughters[2];
  mothers[0]  = mothers[1]   = 0;
  for (int i=0;i<_blob->NInP();i++) mothers[i] = m_connect[_blob->InParton(i)];
  if (_blob->NOutP()>0) {
    daughters[0] = m_connect[_blob->OutParton(0)];
    if (_blob->NOutP()>1) daughters[1] = m_connect[_blob->OutParton(_blob->NOutP()-1)];
    else daughters[1] = m_connect[_blob->OutParton(0)];
  }
  else daughters[0] = daughters[1] = 0;
  
  if (_blob->NInP()>0) {
    for (int i=0;i<_blob->NInP();i++) {
      jdahep[2*mothers[i]]   = daughters[0]+1;
      jdahep[2*mothers[i]+1] = daughters[1]+1;
    }
  }
  if (_blob->NOutP()>0) {
    for (int i=0;i<_blob->NOutP();i++) {
      for (int j=0;j<_blob->NInP();j++){ 
	jmohep[2*m_connect[_blob->OutParton(i)]+j] = mothers[j]+1; 
      }
    }
  }
}

void HepEvt_Interface::EstablishRelations(Parton * _mother,Blob * _blob) {
  int mother;
  int daughters[2];
  mother = m_connect[_mother];
  if (_blob->NOutP()>0) {
    daughters[0] = m_connect[_blob->OutParton(0)];
    if (_blob->NOutP()>1) daughters[1] = m_connect[_blob->OutParton(_blob->NOutP()-1)];
    else daughters[1] = m_connect[_blob->OutParton(0)];
  }
  else daughters[0]   = daughters[1] = 0;
  
  jdahep[2*mother]    = daughters[0]+1;
  jdahep[2*mother+1]  = daughters[1]+1;
  if (_blob->NOutP()>0) {
    for (int i=0;i<_blob->NOutP();i++) {
      jmohep[2*m_connect[_blob->OutParton(i)]]   = mother+1;
      jmohep[2*m_connect[_blob->OutParton(i)]+1] = 0;
    }
  }
}

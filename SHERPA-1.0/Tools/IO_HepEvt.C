#include "IO_HepEvt.H"
#include "Message.H"

using namespace SHERPA;
using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace AMATOOLS;

/*
extern "C" {
  void f2parton_(int &, int *, double *);
  void f2hepevt_(int &, int *, int *, int *, int *, double *, double *);
  
  void pyhepc_(int &);
  void pylist_(int &);
}
*/


IO_HepEvt::IO_HepEvt()
{
  phep   = new double[5*maxentries];
  vhep   = new double[4*maxentries];
  jmohep = new int[2*maxentries];
  jdahep = new int[2*maxentries];
  isthep = new int[maxentries];
  idhep  = new int[maxentries];
}


IO_HepEvt::~IO_HepEvt() {
  if (jmohep) { delete jmohep; jmohep = 0; }
  if (jdahep) { delete jdahep; jdahep = 0; }
  if (isthep) { delete isthep; isthep = 0; }
  if (idhep)  { delete idhep;  idhep  = 0; }
  if (phep)   { delete phep;   phep   = 0; }
  if (vhep)   { delete vhep;   vhep   = 0; }
}


void IO_HepEvt::Blobs2HepEvt(Blob_List * blobs, std::string type, int & _nhep) {
  nhep = _nhep;
  if (nhep==0) mypl.clear();

  int position;
  for (Blob_List::const_iterator bit=blobs->begin(); bit!=blobs->end();++bit) {
    position = (*bit)->Type().find(type);
    if (position > -1) Blob2HepEvt((*bit),0);
  }

  for (Blob_List::const_iterator bit=blobs->begin(); bit!=blobs->end();++bit) {
    position = (*bit)->Type().find(type);
    if (position > -1) Blob2HepEvt((*bit),1);
  }

  //mypl.clear();

  // f2hepevt_(nhep, isthep, idhep, jmohep, jdahep, phep, vhep);

  int line = -1;
  _nhep = nhep;
}


void IO_HepEvt::Blob2HepEvt(Blob * blob, int mode) {
  msg.Debugging()<<"In Blob2HepEvt("<<blob->Type()<<","<<mode<<")  nhep = "<<nhep<<endl;

  int * mohep = new int[2];
  mohep[0] = mohep[1] = 0;
  int * dahep = new int[2];
  dahep[0] = dahep[1] = 0;
  if (mode==0){
    for (int i=0;i<blob->NInP();++i) Parton2HepEvt(blob->InParton(i),mohep,dahep,1-mode);
    return;
  }

  Parton * mo;
  bool     id;
  for (int i=0;i<blob->NInP();i++) {
    mo = blob->InParton(i);
    for (int j=0;j<nhep;j++) {
      if (mo->Flav().HepEvt() == idhep[j]) {
	id = 1;
	for (short int k=1;k<4;k++) {
	  if ( phep[(k-1)+j*5] != mo->Momentum()[k] ) { id = 0; break; }
	}
	if (id) {
	  if ( phep[3+j*5] != mo->Momentum()[0] ) id = 0; 
	}
	if (id) {
	  mohep[i] = j+1;
	  continue;
	}
      }
    }
  }
  if (blob->NInP() < 2) mohep[1] = mohep[0];

  dahep[0] = dahep[1] = -1;
  msg.Debugging()<<blob->Type()<<" "<<mode<<" "<<mohep[0]<<" "<<mohep[1]<<endl;
  for (int i=0;i<blob->NOutP();++i) Parton2HepEvt(blob->OutParton(i),mohep,dahep,mode);
}


void IO_HepEvt::Parton2HepEvt(Parton * parton,int mode) { 
  if (nhep >= maxentries) {
    --nhep;
    msg.Error()<<" ERROR :in void IO_HepEvt::Parton2HepEvt()"<<endl
	       <<"   nehp>=maxentries "<<endl;
    abort();
  }
        
  idhep[nhep]    = parton->Flav().HepEvt();
  isthep[nhep]   = 1;
  if ((parton->Info() == 'I') || (parton->Info() == 'B')) isthep[nhep] = 2;

  for (short int j=1; j<4; ++j) phep[(j-1)+nhep*5] = parton->Momentum()[j];
  phep[3+nhep*5] = parton->Momentum()[0];
  double pabs    = (parton->Momentum()).Abs2();
  if (pabs<0) phep[4+nhep*5] = 0.;
         else phep[4+nhep*5] = sqrt(pabs);

  if (parton->Prod() != 0) {
    for (short int j=1; j<4; ++j) {
      vhep[(j-1)+nhep*4] = parton->XProd()[j];
      vhep[3+nhep*4]     = parton->XProd()[j];
    }
  }
  else {
    for (short int j=0; j<4; ++j) vhep[j+nhep*4] = 0.;
  }

  jmohep[nhep*2]     = 0;
  jmohep[1+nhep*2]   = 0;
  jdahep[nhep*2]     = 0;
  jdahep[1+nhep*2]   = 0;

  if (mode) { mypl.push_back(parton); ++nhep; }
}



void IO_HepEvt::Parton2HepEvt(Parton * parton,int * mohep,int * dahep,int mode) { 
  if (nhep >= maxentries) {
    --nhep;
    msg.Error()<<" ERROR :in void IO_HepEvt::Parton2HepEvt()"<<endl
	       <<"   nehp>=maxentries "<<endl;
    abort();
  }

  for (Parton_List::iterator pit=mypl.begin();pit!=mypl.end();++pit) {
    if ( (parton->Flav()     == (*pit)->Flav()) &&
	 (parton->Momentum() == (*pit)->Momentum()) ) return;
  }

  Parton2HepEvt(parton,0);

  if (mohep[0] == -1) {
    mohep[0] = nhep+1;
    if (mohep[1] == 0) mohep[1] = mohep[0];
  }
  else {
    if (mohep[1] < 0) mohep[1] = nhep+1;
    else { 
      jmohep[nhep*2]         = mohep[0];
      jmohep[1+nhep*2]       = mohep[1];
      if (dahep[0]==-1) {
	dahep[0]             = nhep+1;
	jdahep[2*mohep[0]-2] = dahep[0];
	jdahep[2*mohep[1]-2] = dahep[0];
      }
      dahep[1]               = nhep+1;
      jdahep[2*mohep[0]-1]   = dahep[1];
      jdahep[2*mohep[1]-1]   = dahep[1];
    }
  }
  if (mode) { mypl.push_back(parton); ++nhep; }
}




void IO_HepEvt::Entry2HepEvt(Flavour flav,Vec4D mom,Vec4D pos,int _nhep) {
  nhep = _nhep;
  if (nhep==0) {
    idhep[nhep] = flav.HepEvt();
        
    for (short int j=1; j<4; ++j) phep[(j-1)+nhep*5] = mom[j];
    phep[3+nhep*5] = mom[0];
    double pabs = mom.Abs2();
    if (pabs<0) phep[4+nhep*5] = 0.;
           else phep[4+nhep*5] = sqrt(pabs);
    
    for (short int j=1; j<4; ++j) {
      vhep[(j-1)+nhep*4] = pos[j];
      vhep[3+nhep*4]     = pos[j];
    }

    isthep[nhep] = 1;
    jmohep[nhep*2] = 0;
    jmohep[1+nhep*2] = 0;
    
    jdahep[nhep*2] = 0;
    jdahep[1+nhep*2] = 0;
  }
  ++nhep;
}

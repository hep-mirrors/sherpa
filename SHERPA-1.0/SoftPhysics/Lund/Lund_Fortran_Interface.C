//#include "MyStrStream.H"
#include "Lund_Fortran_Interface.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "Blob.H"

using namespace SHERPA;
using namespace AORGTOOLS;
using namespace AMATOOLS;
using namespace APHYTOOLS;

extern "C" {
  void fhawface_(int&, int*, int*, int*, double*, double*);
  void apyinit_(const double&,double&,double&,double&,const int&,const int&);
  void finterf_(int&, int*, int*, int*, int*, double*, double*);
  void pylist_(int&);
};

Lund_Fortran_Interface::Lund_Fortran_Interface(double _a,double _b, double _sigma) :
  a(_a), b(_b), sigma(_sigma)
{ 
  Init();
  apyinit_(AORGTOOLS::rpa.gen.Ecms(),a,b,sigma,1,0);
  //  apyinit_(sqr(AORGTOOLS::rpa.gen.Ecms()),a,b,sigma,1,0);
  pyinit = 1;
} 

Lund_Fortran_Interface::~Lund_Fortran_Interface() {
  if (phep)   delete phep;
  if (vhep)   delete vhep;
  if (jmohep) delete jmohep;
  if (jdahep) delete jdahep;
  if (isthep) delete isthep;
  if (idhep)  delete idhep;
}

void Lund_Fortran_Interface::Init()
{
  phep   = new double[5*maxentries];
  vhep   = new double[4*maxentries];
  jmohep = new int[2*maxentries];
  jdahep = new int[2*maxentries];
  isthep = new int[maxentries];
  idhep  = new int[maxentries];
}

bool Lund_Fortran_Interface::Hadronize(APHYTOOLS::Blob * blob,
				       APHYTOOLS::Blob_List * bloblist,
				       APHYTOOLS::Parton_List * pl) {
  blob->SetType(std::string("Fragmentation (Lund : Pythia 6.163)"));
  int nhep = 0;

  if (nhep==0) {
    idhep[nhep] = Flavour(kf::photon).HepEvt();
        
    for (short int j=1; j<4; ++j) phep[(j-1)+nhep*5] = blob->CMS()[j];
    phep[3+nhep*5] = blob->CMS()[0];
    double pabs = (blob->CMS()).Abs2();
    if (pabs<0) phep[4+nhep*5] = 0.;
           else phep[4+nhep*5] = sqrt(pabs);
    
    for (short int j=0; j<4; ++j) { vhep[j+nhep*4] = 0.; }
    isthep[nhep] = 1;
    jmohep[nhep*2] = 0;
    jmohep[1+nhep*2] = 0;
    
    jdahep[nhep*2] = 0;
    jdahep[1+nhep*2] = 0;
  }

  for (int i=0;i<blob->NInP();++i) {
    AddPartonToString(blob->InParton(i),nhep);
  }

  for (int i=0;i<nhep;i++) {
    msg.Debugging()<<"  "<<i<<" : "<<std::endl
		   <<"  "<<idhep[i]<<" / ("
		   <<phep[5*i+3]<<","<<phep[5*i+0]<<","
		   <<phep[5*i+1]<<","<<phep[5*i+2]<<")"<<std::endl<<std::endl;
  }

  
  int dummy=1;
  if (rpa.gen.Events()) {
    msg.Out()<<"before hadronisation"<<std::endl;
    pylist_(dummy);
    msg.Out()<<"~~~~~~~~~~~~~~~~~~~~"<<std::endl;
  }
  finterf_(nhep, isthep, idhep, jmohep, jdahep, phep, vhep);
  if (rpa.gen.Events()) {
    msg.Out()<<"after hadronisation"<<std::endl;
    pylist_(dummy);
    msg.Out()<<std::endl<<std::endl;
  }


  FillPrimaryHadronsInBlob(blob,bloblist,pl);
}


void Lund_Fortran_Interface::AddPartonToString(Parton * parton,int & nhep)
{
  idhep[nhep] = parton->Flav().HepEvt();
        
  for (short int j=1; j<4; ++j) phep[(j-1)+nhep*5] = parton->Momentum()[j];
  phep[3+nhep*5] = parton->Momentum()[0];
  double pabs = (parton->Momentum()).Abs2();
  if (pabs<0) phep[4+nhep*5] = 0.;
         else phep[4+nhep*5] = sqrt(pabs);

  for (short int j=1; j<4; ++j) {
    vhep[(j-1)+nhep*4] = parton->XProd()[j];
    vhep[3+nhep*4]     = parton->XProd()[j];
  }

  isthep[nhep] = 1;
  jmohep[nhep*2] = 0;
  jmohep[1+nhep*2] = 0;
  
  jdahep[nhep*2] = 0;
  jdahep[1+nhep*2] = 0;

  nhep++;
}

void Lund_Fortran_Interface::FillPrimaryHadronsInBlob(APHYTOOLS::Blob * blob,
						      APHYTOOLS::Blob_List * bloblist,
						      APHYTOOLS::Parton_List * pl)
{
  pjet      = new double[4*maxentries];
  xjet      = new double[4*maxentries];
  kfjet     = new int[maxentries];
  mothers   = new int[maxentries];
  daughters = new int[2*maxentries];
    
  fhawface_(nk,kfjet,mothers,daughters,pjet,xjet);
    
  Blob    * decay;
  Parton  * parton;
  Flavour   flav;
  Vec4D     momentum,position;

  int number;
  for (int i=0; i<nk; ++i) {
    flav = Flavour(kf::code(abs(*(kfjet+i))));
    if (flav==Flavour(kf::code(92))) {
      msg.Debugging()<<"Found string !"<<i<<std::endl
		     <<" -> "<<(*(daughters+2*i))<<" / "<<(*(daughters+2*i+1))<<std::endl;
      for (int j=(*(daughters+2*i))-1;j<(*(daughters+2*i+1));j++) {
	msg.Debugging()<<"String offspring !"<<j<<std::endl
		     <<" -> "<<(*(daughters+2*j))<<" / "<<(*(daughters+2*j+1))<<std::endl
		     <<" <- "<<(*(mothers+j))<<std::endl;
	flav = Flavour(kf::code(abs(*(kfjet+j))));
	if (!flav.IsHadron()) continue;

	if ((*(kfjet+j))<0)    flav        = flav.Bar();
	for(int k=0; k<4; ++k) momentum[k] = *(pjet+j+k*2000);
	for(int k=0; k<4; ++k) position[k] = *(xjet+j+k*2000);
	parton         = new Parton(-1,flav,momentum);
	if (pl) number = pl->size();
	else number    = int(parton);
	parton->SetNumber(number);
	parton->SetStatus(1);
	parton->SetProductionBlob(blob);      
	blob->SetPosition(position);
	if (pl) pl->push_back(parton);
	blob->AddToOutPartons(parton);
	if (*(daughters+2*j)!=0 && *(daughters+2*j+1)!=0) {
	  msg.Debugging()<<"Decayed particle : "<<j<<" "<<flav<<" "<<momentum<<std::endl
			 <<"                          "<<position<<std::endl
			 <<" -> "<<(*(daughters+2*j))<<" / "<<(*(daughters+2*j+1))<<std::endl
			 <<" <- "<<(*(mothers+j))<<std::endl;
	  decay = new Blob();
	  decay->SetStatus(1);
	  decay->SetType(std::string("Hadron decay"));
	  decay->SetId(bloblist->size());
	  decay->AddToInPartons(parton);
	  parton->SetDecayBlob(decay);
	  bloblist->push_back(decay);
	  FillSecondaryHadronsInBlob(decay,bloblist,(*(daughters+2*j))-1,(*(daughters+2*j+1)),pl);
	}
	else {
	  msg.Debugging()<<"Non - Decayed particle : "<<i<<" "<<flav<<" "<<momentum<<std::endl
			 <<"                          "<<position<<std::endl
			 <<" -> "<<(*(daughters+i))<<" / "<<(*(daughters+2000+i))<<std::endl
			 <<" <- "<<(*(mothers+i))<<std::endl;
	}
      }
    }
  }
  blob->SetStatus(0);
    
  delete kfjet;
  delete pjet;
  delete mothers;
  delete daughters;
}

void Lund_Fortran_Interface::FillSecondaryHadronsInBlob(APHYTOOLS::Blob * blob,
							APHYTOOLS::Blob_List * bloblist,
							int daughter1,int daughter2,
							APHYTOOLS::Parton_List * pl) 
{
  Blob    * decay;
  Parton  * parton;
  Flavour   flav;
  Vec4D     momentum, position;
  int number;
  for (int i=daughter1; i<daughter2; ++i) {
    flav = Flavour(kf::code(abs(*(kfjet+i))));
    if ((*(kfjet+i))<0)    flav        = flav.Bar();
    for(int k=0; k<4; ++k) momentum[k] = *(pjet+i+k*2000);
    for(int k=0; k<4; ++k) position[k] = *(xjet+i+k*2000);

    parton         = new Parton(-1,flav,momentum);
    if (pl) number = pl->size();
    else number    = int(parton);
    parton->SetNumber(number);
    parton->SetStatus(1);
    parton->SetProductionBlob(blob);  
    blob->SetPosition(position);
    if (pl) pl->push_back(parton);
    blob->AddToOutPartons(parton);
    if (*(daughters+2*i)!=0 && *(daughters+2*i+1)!=0) {
      msg.Debugging()<<"Decayed particle : "<<i<<" "<<flav<<" "<<momentum<<std::endl
		     <<"                          "<<position<<std::endl
		     <<" -> "<<(*(daughters+2*i))<<" / "<<(*(daughters+2*i+1))<<std::endl
		     <<" <- "<<(*(mothers+i))<<std::endl;
      decay = new Blob();
      decay->SetStatus(1);
      decay->SetType(std::string("Hadron decay"));
      decay->SetId(bloblist->size());
      decay->AddToInPartons(parton);
      parton->SetDecayBlob(decay);
      bloblist->push_back(decay);
      FillSecondaryHadronsInBlob(decay,bloblist,(*(daughters+2*i))-1,(*(daughters+2*i+1)),pl);
    }
    else {
      msg.Debugging()<<"Non - Decayed particle : "<<i<<" "<<flav<<" "<<momentum<<std::endl
		     <<"                          "<<position<<std::endl
		     <<" -> "<<(*(daughters+2*i))<<" / "<<(*(daughters+2*i+1))<<std::endl
		     <<" <- "<<(*(mothers+i))<<std::endl;
    }
  }
}

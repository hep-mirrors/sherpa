//#include "MyStrStream.H"
#include "Lund_Fortran_Interface.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "Blob.H"

using namespace MOCAIC;
using namespace AORGTOOLS;
using namespace AMATOOLS;
using namespace APHYTOOLS;

extern "C" {
  void fhawface_(int&, int*, double*);
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

void Lund_Fortran_Interface::Hadronize(APHYTOOLS::Parton_List * pl,
				       APHYTOOLS::Blob * blob) {
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
    msg.Out()<<"before hadronisation"<<endl;
    pylist_(dummy);
    msg.Out()<<"~~~~~~~~~~~~~~~~~~~~"<<endl;
  }
  finterf_(nhep, isthep, idhep, jmohep, jdahep, phep, vhep);
  if (rpa.gen.Events()) {
    msg.Out()<<"after hadronisation"<<endl;
    pylist_(dummy);
    msg.Out()<<endl<<endl;
  }


  FillHadronsInBlob(pl,blob);
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
    vhep[3+nhep*4] = parton->XProd()[j];
  }

  isthep[nhep] = 1;
  jmohep[nhep*2] = 0;
  jmohep[1+nhep*2] = 0;
  
  jdahep[nhep*2] = 0;
  jdahep[1+nhep*2] = 0;

  nhep++;
}

void Lund_Fortran_Interface::FillHadronsInBlob(APHYTOOLS::Parton_List * pl,
					       APHYTOOLS::Blob * blob)
{
  Parton  * parton;
  Flavour   flav;
  Vec4D     momentum;
  
  double  * pjet;
  int       nk, * kfjet;
    
  kfjet = new int[maxentries];
  pjet  = new double[4*maxentries];
    
  fhawface_(nk,kfjet,pjet);
    
  for (int i=0; i<nk; ++i) {
    flav = Flavour(kf::code(abs(*(kfjet+i))));
    if ((*(kfjet+i))<0) flav = flav.Bar();
    for(int j=0; j<4; ++j) momentum[j] = *(pjet+i+j*2000);

    parton = new Parton(pl->size(),flav,momentum);
    parton->SetStatus(1);
    parton->SetProd(blob);

    pl->push_back(parton);
    blob->AddToOutPartons(parton);
  }
  blob->SetStatus(0);
    
  delete kfjet;
  delete pjet;
}

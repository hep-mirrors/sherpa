//#include "MyStrStream.H"
#include "Lund_Fortran_Interface.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "Blob.H"

using namespace SHERPA;
using namespace ATOOLS;

extern "C" {
  void fhawface_(int&, int*, int*, int*, double*, double*);
  void apyinit_(const double&,double&,double&,double&,const int&,const int&);
  void finterf_(int&, int*, int*, int*, int*, double*, double*);
  void pylist_(int&);
};

Lund_Fortran_Interface::Lund_Fortran_Interface(double _a,double _b, double _sigma) 
{ 
  Init();
  apyinit_(ATOOLS::rpa.gen.Ecms(),_a,_b,_sigma,1,0);
} 

Lund_Fortran_Interface::~Lund_Fortran_Interface() {
  if (phep)   delete [] phep;
  if (vhep)   delete [] vhep;
  if (jmohep) delete [] jmohep;
  if (jdahep) delete [] jdahep;
  if (isthep) delete [] isthep;
  if (idhep)  delete [] idhep;

  if (daughters) delete [] daughters;
  if (mothers)   delete [] mothers;
  if (kfjet)     delete [] kfjet;
  if (xjet)      delete [] xjet;
  if (pjet)      delete [] pjet;
}

void Lund_Fortran_Interface::Init()
{
  phep   = new double[5*maxentries];
  vhep   = new double[4*maxentries];
  jmohep = new int[2*maxentries];
  jdahep = new int[2*maxentries];
  isthep = new int[maxentries];
  idhep  = new int[maxentries];

  pjet      = new double[4*maxentries];
  xjet      = new double[4*maxentries];
  kfjet     = new int[maxentries];
  mothers   = new int[maxentries];
  daughters = new int[2*maxentries];
}

bool Lund_Fortran_Interface::Hadronize(ATOOLS::Blob * blob,
				       ATOOLS::Blob_List * bloblist,
				       ATOOLS::Particle_List * pl) {
  blob->SetType(btp::Fragmentation);
  blob->SetTypeSpec("Lund : Pythia 6.214");
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
    AddPartonToString(blob->InParticle(i),nhep);
  }

  int dummy=2;
  finterf_(nhep, isthep, idhep, jmohep, jdahep, phep, vhep);
  //if (rpa.gen.Debugging()) {
    msg.Out()<<"after hadronisation"<<std::endl;
    pylist_(dummy);
    msg.Out()<<std::endl<<std::endl;
    //}


  FillPrimaryHadronsInBlob(blob,bloblist,pl);
  return 1;
}


void Lund_Fortran_Interface::AddPartonToString(Particle * parton,int & nhep)
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

void Lund_Fortran_Interface::FillPrimaryHadronsInBlob(ATOOLS::Blob * blob,
						      ATOOLS::Blob_List * bloblist,
						      ATOOLS::Particle_List * pl)
{
  
  fhawface_(nk,kfjet,mothers,daughters,pjet,xjet);
    
  Blob    * decay;
  Particle  * particle;
  Flavour   flav;
  Vec4D     momentum,position;

  int number;
  for (int i=0; i<nk; ++i) {
    flav.FromHepEvt((*(kfjet+i)));
    if (flav==Flavour(kf::string) || flav==Flavour(kf::cluster)) {
      for (int j=(*(daughters+2*i))-1;j<(*(daughters+2*i+1));j++) {
	//	flav = Flavour(kf::code(abs(*(kfjet+j))));
	flav.FromHepEvt((*(kfjet+j)));
	if (!flav.IsHadron()) continue;

	// if ((*(kfjet+j))<0)    flav        = flav.Bar();
	for(int k=0; k<4; ++k) momentum[k] = *(pjet+j+k*2000);
	for(int k=0; k<4; ++k) position[k] = *(xjet+j+k*2000);
	particle         = new Particle(-1,flav,momentum);
	if (pl) number = pl->size();
	else number    = (long int)(particle);
	particle->SetNumber(number);
	particle->SetStatus(1);
	particle->SetInfo('P');
	blob->SetPosition(position);
	if (pl) pl->push_back(particle);
	blob->AddToOutParticles(particle);
	if (*(daughters+2*j)!=0 && *(daughters+2*j+1)!=0) {
	  decay = new Blob();
	  decay->SetStatus(1);
	  decay->SetType(btp::Hadron_Decay);
	  decay->SetTypeSpec("Lund Pythia 6.214");
	  decay->SetId();
	  decay->AddToInParticles(particle);
	  if (particle->Info()=='P') particle->SetInfo('p');
	  if (particle->Info()=='D') particle->SetInfo('d');
	  bloblist->push_back(decay);
	  FillSecondaryHadronsInBlob(decay,bloblist,(*(daughters+2*j))-1,(*(daughters+2*j+1)),pl);
	}
      }
    }
  }
  blob->SetStatus(0);
    
}

void Lund_Fortran_Interface::FillSecondaryHadronsInBlob(ATOOLS::Blob * blob,
							ATOOLS::Blob_List * bloblist,
							int daughter1,int daughter2,
							ATOOLS::Particle_List * pl) 
{
  Blob    * decay;
  Particle  * particle;
  Flavour   flav;
  Vec4D     momentum, position;
  int number;
  for (int i=daughter1; i<daughter2; ++i) {
    flav.FromHepEvt((*(kfjet+i)));
    //flav = Flavour(kf::code(abs(*(kfjet+i))));
    //if ((*(kfjet+i))<0)    flav        = flav.Bar();
    for(int k=0; k<4; ++k) momentum[k] = *(pjet+i+k*2000);
    for(int k=0; k<4; ++k) position[k] = *(xjet+i+k*2000);

    particle         = new Particle(-1,flav,momentum);
    if (pl) number = pl->size();
    else number    = (long int)(particle);
    particle->SetNumber(number);
    particle->SetStatus(1);
    particle->SetInfo('D');
    blob->SetPosition(position);
    if (pl) pl->push_back(particle);
    blob->AddToOutParticles(particle);
    if (*(daughters+2*i)!=0 && *(daughters+2*i+1)!=0) {
      decay = new Blob();
      decay->SetStatus(1);
      decay->SetType(btp::Hadron_Decay);
      decay->SetTypeSpec("Lund Pythia 6.214");
      decay->SetId();
      decay->AddToInParticles(particle);
      if (particle->Info()=='P') particle->SetInfo('p');
      if (particle->Info()=='D') particle->SetInfo('d');
      bloblist->push_back(decay);
      FillSecondaryHadronsInBlob(decay,bloblist,(*(daughters+2*i))-1,(*(daughters+2*i+1)),pl);
    }
  }
}

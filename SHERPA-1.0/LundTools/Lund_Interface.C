#include "Lund_Interface.H"

#include "Lund_Wrapper.H"
#include "Data_Reader.H"
#include "Particle.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "Exception.H"
#ifndef NO_EXPORT__AlphaS
#include "Running_AlphaS.H"
#endif
#include "Scaling.H"

using namespace SHERPA;

size_t Lund_Interface::s_errors=0;
size_t Lund_Interface::s_maxerrors=0;

Lund_Interface::Lund_Interface(std::string _m_path,std::string _m_file):
  m_path(_m_path),
  m_file(_m_file)
{
  double win;
  std::string beam[2], frame("CMS");
  ATOOLS::Flavour flav[2];
  for (size_t i=0;i<2;++i) flav[i]=ATOOLS::rpa.gen.Bunch(i);
  if (flav[0]==ATOOLS::kf::e && flav[1]==ATOOLS::kf::p_plus) {
    beam[0]="e-";
    beam[1]="p+";
  }
  else if (flav[0]==ATOOLS::kf::p_plus && flav[1]==ATOOLS::kf::e) {
    beam[0]="p+";
    beam[1]="e-";
  }
  else if (flav[0]==ATOOLS::kf::e && flav[1]==ATOOLS::kf::photon) {
    if (flav[0].IsAnti()) beam[0]="e+"; else beam[0]="e-";
    beam[1]="gamma";
    pysubs.msub[33]=1;    
  }
  else if (flav[0]==ATOOLS::kf::photon && flav[1]==ATOOLS::kf::e) {
    beam[0]="gamma";
    if (flav[1].IsAnti()) beam[1]="e+"; else beam[1]="e-";
    pysubs.msub[33]=1;    
  }
  else if (flav[0]==ATOOLS::kf::photon && flav[1]==ATOOLS::kf::photon) {
    for (size_t i=0;i<2;++i) beam[i]="gamma";
    pysubs.msub[57]=1;    
  }
  else if (flav[0]==ATOOLS::kf::e && flav[1]==ATOOLS::kf::e) {
    for (size_t i=0;i<2;++i) if (flav[i].IsAnti()) beam[i]="e+"; else beam[i]="e-";
    pysubs.msub[0]=1;    
    pypars.mstp[47]=1;
    pydat1.mstj[100]=5;
  }
  else {
    for (size_t i=0;i<2;++i) if (flav[i].IsAnti()) beam[i]="p-"; else beam[i]="p+";
    pysubs.msub[0]=1;    
    pypars.mstp[47]=1;
    pydat1.mstj[100]=5;
  }
  win=ATOOLS::rpa.gen.Ecms();
  s_maxerrors=ATOOLS::rpa.gen.NumberOfEvents();
  std::vector<std::vector<double> > help;
  ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader("=",";","!");
  reader->SetMatrixType(reader->MTransposed);
  reader->SetInputPath(m_path);
  reader->SetInputFile(m_file);
  reader->AddIgnore("(");
  reader->AddIgnore(")");
  reader->AddIgnore(",");
  // if (!reader->ReadFromFile(pysubs.msel,"MSEL")) pysubs.msel=1;
  pysubs.msel=0;
  reader->MatrixFromFile(help,"MSUB");
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>1) if ((int)help[i][0]>0) pysubs.msub[(int)help[i][0]-1]=(int)help[i][1];
  }
  reader->MatrixFromFile(help,"KFIN");
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>2) {
      if (((int)help[i][0]>0)&&((int)help[i][1]>-41)) {
	pysubs.kfin[(int)help[i][1]+40][(int)help[i][0]-1]=(int)help[i][2];
      }
    }
  }
  reader->MatrixFromFile(help,"CKIN");
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>1) if ((int)help[i][0]>0) pysubs.ckin[(int)help[i][0]-1]=help[i][1];
  }
  reader->MatrixFromFile(help,"MSTJ");
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>1) if ((int)help[i][0]>0) pydat1.mstj[(int)help[i][0]-1]=(int)help[i][1];
  }
  reader->MatrixFromFile(help,"MSTP");
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>1) if ((int)help[i][0]>0) pypars.mstp[(int)help[i][0]-1]=(int)help[i][1];
  }
  reader->MatrixFromFile(help,"MSTU");
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>1) if ((int)help[i][0]>0) pydat1.mstu[(int)help[i][0]-1]=(int)help[i][1];
  }
  reader->MatrixFromFile(help,"PARP");
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>1) if ((int)help[i][0]>0) pypars.parp[(int)help[i][0]-1]=help[i][1];
  }
  reader->MatrixFromFile(help,"PARJ");
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>1) if ((int)help[i][0]>0) pydat1.parj[(int)help[i][0]-1]=help[i][1];
  }
  reader->MatrixFromFile(help,"PARU");
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>1) if ((int)help[i][0]>0) pydat1.paru[(int)help[i][0]-1]=help[i][1];
  }
  // the next lines replace the apyinit_ call
  hepevt.nhep=100;
  for (int i=pydat3.mdcy[23-1][2-1];i<pydat3.mdcy[23-1][2-1]+pydat3.mdcy[23-1][3-1];++i) {
    if (abs(pydat3.kfdp[i-1][1-1])>=2) pydat3.mdme[i-1][1-1]=ATOOLS::Min(0,pydat3.mdme[i-1][1-1]);
  }
  // replacement ends here
  if (ATOOLS::msg.Level()&8) ListLundParameters();
  pyinit(frame.c_str(),beam[0].c_str(),beam[1].c_str(),win);
  delete reader;
}

Lund_Interface::~Lund_Interface()
{
}

bool Lund_Interface::ConvertParticles(std::map<int,ATOOLS::Particle*> &converted)
{
  for (int i=0;i<hepevt.nhep;++i) {
    ATOOLS::msg.Debugging()<<i<<" "<<hepevt.isthep[i]<<" "<<hepevt.idhep[i]<<" "
			   <<hepevt.jmohep[i][0]<<" "
			   <<hepevt.jdahep[i][0]<<" "<<hepevt.jdahep[i][1]<<" "
			   <<"("<<hepevt.phep[i][0]<<","<<hepevt.phep[i][1]<<","
			   <<hepevt.phep[i][2]<<","<<hepevt.phep[i][3]<<") "
			   <<"("<<hepevt.vhep[i][0]<<","<<hepevt.vhep[i][1]<<","
			   <<hepevt.vhep[i][2]<<","<<hepevt.vhep[i][3]<<")"<<std::endl;
    ATOOLS::Flavour flavour;
    flavour.FromHepEvt(hepevt.idhep[i]);
    ATOOLS::Particle *newpart = new ATOOLS::Particle(i+1,flavour,
						     ATOOLS::Vec4D(hepevt.phep[i][3],hepevt.phep[i][0],
								   hepevt.phep[i][1],hepevt.phep[i][2]));
    newpart->SetStatus(hepevt.isthep[i]);
    converted[i]=newpart;
  }
  return true;
}

bool Lund_Interface::ConstructBlobs(ATOOLS::Blob_List *blobs,std::map<int,ATOOLS::Particle*> &converted)
{
  for (int i=0;i<hepevt.nhep;++i) {
    ATOOLS::Blob *productionblob=converted[i]->ProductionBlob();
    ATOOLS::Particle *mother=converted[hepevt.jmohep[i][0]-1];
    if (productionblob==NULL) {
      bool createblob=(mother==NULL);
      if (!createblob) createblob=createblob||(mother->DecayBlob()==NULL); 
      if (createblob) {
	  productionblob = new ATOOLS::Blob();
	  productionblob->SetPosition(ATOOLS::Vec4D(hepevt.vhep[i][0],hepevt.vhep[i][1],
						    hepevt.vhep[i][2],hepevt.vhep[i][3]));
	  productionblob->AddToOutParticles(converted[i]);
	  if (mother!=NULL) productionblob->AddToInParticles(mother);
	  productionblob->SetId(blobs->size());
	  blobs->push_back(productionblob);
      }
      else {
	mother->DecayBlob()->AddToOutParticles(converted[i]);
	productionblob=converted[i]->ProductionBlob();
      }
    }
    ATOOLS::Blob *decayblob=converted[i]->DecayBlob();
    for (unsigned int j=0;j<2;++j) {
      if (decayblob==NULL) {
	ATOOLS::Particle *daughter=converted[hepevt.jdahep[i][j]-1];
	if (daughter!=NULL) {
	  if (daughter->ProductionBlob()!=NULL) {
	    daughter->ProductionBlob()->AddToInParticles(converted[i]);
	  }
	}
      }
    }
  }
  for (int i=0;i<hepevt.nhep;++i) {
    if (converted[i]->DecayBlob()==NULL) {
      ATOOLS::Particle *daughter=converted[hepevt.jdahep[i][0]-1];
      if (daughter!=NULL) daughter->ProductionBlob()->AddToInParticles(converted[i]);
    }
  }
  return true;
}

bool Lund_Interface::SetTypes(ATOOLS::Blob_List *blobs)
{
  for (ATOOLS::Blob_Iterator bit=blobs->begin();bit!=blobs->end();++bit) {
    bool final=true;
    for (int i=0;i<(*bit)->NOutP();++i) {
      if (((*bit)->OutParticle(i)->DecayBlob()!=NULL)||
	  ((*bit)->OutParticle(i)->Status()!=1)) final=false;
    }
    if (final) (*bit)->SetStatus(1);
    else (*bit)->SetStatus(0);
  }
  return true;
}

bool Lund_Interface::DeleteObsolete(std::map<int,ATOOLS::Particle*> &converted)
{
  for (int i=0;i<hepevt.nhep;++i) {
    if ((converted[i]->ProductionBlob()==NULL)&&
	(converted[i]->DecayBlob()==NULL)) {
      delete converted[i];
    }
  }
  return true;
}

bool Lund_Interface::ConvertEvent(ATOOLS::Blob_List *bloblist)
{
  if (ATOOLS::msg.Level()>2) pylist(1);
  bool success=true;
  std::map<int,ATOOLS::Particle*> converted;
  bloblist->clear();
  success=success&&ConvertParticles(converted);
  success=success&&ConstructBlobs(bloblist,converted);
  success=success&&SetTypes(bloblist);
  success=success&&DeleteObsolete(converted);
  return success;
}

bool Lund_Interface::Hadronize(ATOOLS::Blob *blob,ATOOLS::Blob_List *bloblist,
			       ATOOLS::Particle_List *pl) 
{
  int nhep = 0;
  blob->SetType(ATOOLS::btp::Fragmentation);
  blob->SetTypeSpec("Pythia_v6.214");
  if (nhep==0) {
    hepevt.idhep[nhep]=ATOOLS::Flavour(ATOOLS::kf::photon).HepEvt();
    for (short int j=1;j<4;++j) hepevt.phep[nhep][j-1]=blob->CMS()[j];
    hepevt.phep[nhep][3]=blob->CMS()[0];
    double pabs=(blob->CMS()).Abs2();
    if (pabs<0) hepevt.phep[nhep][4]=0.0;
    else hepevt.phep[nhep][4]=sqrt(pabs);
    for (short int j=0;j<4;++j) hepevt.vhep[nhep][j]=0.0;
    hepevt.isthep[nhep]=1;
    hepevt.jmohep[nhep][0]=0;
    hepevt.jmohep[nhep][1]=0;
    hepevt.jdahep[nhep][0]=0;
    hepevt.jdahep[nhep][1]=0;
  }
  for (int i=0;i<blob->NInP();++i) {
    AddPartonToString(blob->InParticle(i),nhep);
  }
  int dummy=2;
  // the next lines replace the finterf_ call
  hepevt.nevhep=0;
  hepevt.nhep=nhep;
  pyhepc(2);
  pydat1.mstu[70-1]=1;
  pydat1.mstu[71-1]=hepevt.nhep;
  pyexec();
  pyhepc(1);
  pydat1.mstu[70-1]=2;
  pydat1.mstu[72-1]=hepevt.nhep;
  // replacement ends here
  if (ATOOLS::msg.LevelIsDebugging()) {
    ATOOLS::msg.Tracking()<<"Lund_Interface::Hadronize(..): Passed hadronisation (pythia 6.214)."<<std::endl;
    pylist(dummy);
    ATOOLS::msg.Out()<<std::endl<<std::endl;
  }
  FillPrimaryHadronsInBlob(blob,bloblist,pl);
  return true;
}

void Lund_Interface::AddPartonToString(ATOOLS::Particle *parton,int &nhep)
{
  hepevt.idhep[nhep]=parton->Flav().HepEvt();
  for (short int j=1; j<4; ++j) hepevt.phep[nhep][j-1]=parton->Momentum()[j];
  hepevt.phep[nhep][3]=parton->Momentum()[0];
  double pabs=(parton->Momentum()).Abs2();
  if (pabs<0) hepevt.phep[nhep][4]=0.0;
  else hepevt.phep[nhep][4]=sqrt(pabs);
  for (short int j=1;j<4;++j) {
    hepevt.vhep[nhep][j-1]=parton->XProd()[j];
    hepevt.vhep[nhep][3]=parton->XProd()[j];
  }
  hepevt.isthep[nhep]=1;
  hepevt.jmohep[nhep][0]=0;
  hepevt.jmohep[nhep][1]=0;
  hepevt.jdahep[nhep][0]=0;
  hepevt.jdahep[nhep][1]=0;
  nhep++;
}

void Lund_Interface::FillPrimaryHadronsInBlob(ATOOLS::Blob *blob,ATOOLS::Blob_List *bloblist,
					      ATOOLS::Particle_List *pl)
{
  ATOOLS::Blob *decay;
  ATOOLS::Particle *particle;
  ATOOLS::Flavour flav;
  ATOOLS::Vec4D momentum, position;
  int number;
  for (int i=0;i<hepevt.nhep;++i) {
    if ((hepevt.isthep[i]!=2)&&(hepevt.isthep[i]!=1)&&(hepevt.isthep[i]!=149)) continue;
    flav.FromHepEvt(hepevt.idhep[i]);
    if (flav==ATOOLS::Flavour(ATOOLS::kf::string) || 
	flav==ATOOLS::Flavour(ATOOLS::kf::cluster)) {
      for (int j=hepevt.jdahep[i][0]-1;j<hepevt.jdahep[i][1];j++) {
	// flav=Flavour(ATOOLS::kf::code(abs(hepevt.idhep[j])));
	// if (hepevt.idhep[j]<0) flav=flav.Bar();
	flav.FromHepEvt(hepevt.idhep[j]);
	if (!flav.IsHadron()) continue;
	momentum=ATOOLS::Vec4D(hepevt.phep[j][3],hepevt.phep[j][0],
			       hepevt.phep[j][1],hepevt.phep[j][2]);
	position=ATOOLS::Vec4D(hepevt.vhep[j][3],hepevt.vhep[j][0],
			       hepevt.vhep[j][1],hepevt.vhep[j][2]);
	particle = new ATOOLS::Particle(-1,flav,momentum);
	if (pl) number=pl->size();
	else number=(long int)(particle);
	particle->SetNumber(number);
	particle->SetStatus(1);
	particle->SetInfo('P');
	blob->SetPosition(position);
	if (pl) pl->push_back(particle);
	blob->AddToOutParticles(particle);
	if (hepevt.jdahep[j][0]!=0 && hepevt.jdahep[j][1]!=0) {
	  decay = new ATOOLS::Blob();
	  decay->SetStatus(1);
	  decay->SetType(ATOOLS::btp::Hadron_Decay);
	  decay->SetTypeSpec("Pythia_v6.214");
	  decay->SetId(bloblist->size());
	  decay->AddToInParticles(particle);
	  if (particle->Info()=='P') particle->SetInfo('p');
	  if (particle->Info()=='D') particle->SetInfo('d');
	  bloblist->push_back(decay);
	  FillSecondaryHadronsInBlob(decay,bloblist,hepevt.jdahep[j][0]-1,hepevt.jdahep[j][1],pl);
	}
      }
    }
  }
  blob->SetStatus(0);
}

void Lund_Interface::FillSecondaryHadronsInBlob(ATOOLS::Blob *blob,ATOOLS::Blob_List *bloblist,
						int daughter1,int daughter2,ATOOLS::Particle_List *pl) 
{
  ATOOLS::Blob *decay;
  ATOOLS::Particle *particle;
  ATOOLS::Flavour flav;
  ATOOLS::Vec4D momentum, position;
  int number;
  for (int i=daughter1;i<daughter2;++i) {
    //flav=Flavour(ATOOLS::kf::code(abs(hepevt.idhep[i])));
    //if ((*(kfjet+i))<0) flav=flav.Bar();
    flav.FromHepEvt(hepevt.idhep[i]);
    momentum=ATOOLS::Vec4D(hepevt.phep[i][3],hepevt.phep[i][0],
			   hepevt.phep[i][1],hepevt.phep[i][2]);
    position=ATOOLS::Vec4D(hepevt.vhep[i][3],hepevt.vhep[i][0],
			   hepevt.vhep[i][1],hepevt.vhep[i][2]);
    particle = new ATOOLS::Particle(-1,flav,momentum);
    if (pl) number=pl->size();
    else number=(long int)(particle);
    particle->SetNumber(number);
    particle->SetStatus(1);
    particle->SetInfo('D');
    blob->SetPosition(position);
    if (pl) pl->push_back(particle);
    blob->AddToOutParticles(particle);
    if (hepevt.jdahep[i][0]!=0 && hepevt.jdahep[i][1]!=0) {
      decay = new ATOOLS::Blob();
      decay->SetStatus(1);
      decay->SetType(ATOOLS::btp::Hadron_Decay);
      decay->SetTypeSpec("Pythia_v6.214");
      decay->SetId(bloblist->size());
      decay->AddToInParticles(particle);
      if (particle->Info()=='P') particle->SetInfo('p');
      if (particle->Info()=='D') particle->SetInfo('d');
      bloblist->push_back(decay);
      FillSecondaryHadronsInBlob(decay,bloblist,hepevt.jdahep[i][0]-1,hepevt.jdahep[i][1],pl);
    }
  }
}

void Lund_Interface::Error(const int error)
{
  ++s_errors;
  if (s_errors>s_maxerrors) {
    throw(ATOOLS::Exception(ATOOLS::ex::critical_error,std::string("Pythia calls PYERRM(")+
			    ATOOLS::ToString(error)+std::string(")"),"Lund_Interface","Error"));
  }
  else {
    ATOOLS::msg.Error()<<"Lund_Interface::Error("<<error<<") "<<ATOOLS::om::red
		       <<"Pythia calls PYERRM("<<error<<")."<<ATOOLS::om::reset<<std::endl;
    if (ATOOLS::msg.LevelIsError()) pylist(2);
  }
}

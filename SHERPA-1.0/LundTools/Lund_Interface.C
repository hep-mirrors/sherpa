#include "Lund_Interface.H"

#include "Lund_Wrapper.H"
#include "Data_Reader.H"
#include "Particle.H"
#include "Run_Parameter.H"
#include "Message.H"

using namespace SHERPA;

inline void MakeFortranString(char *output,std::string input,unsigned int length)
{
  for (unsigned int i=0;i<length;++i) output[i]=(char)32;
  for (size_t j=0;j<input.length();++j) output[j]=(char)input[j];
}

Lund_Interface::Lund_Interface(std::string _m_path,std::string _m_file):
  m_path(_m_path),
  m_file(_m_file)
{
  SetType("Perturbative");
  SetName("The Lund Monte Carlo");
  double win;
  std::string frame, beam, target;
  std::vector<std::vector<double> > help;
  ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader("=",";","!");
  reader->SetMatrixType(reader->MTransposed);
  reader->SetInputPath(m_path);
  reader->SetInputFile(m_file);
  reader->AddIgnore("(");
  reader->AddIgnore(")");
  reader->AddIgnore(",");
  if (!reader->ReadFromFile(frame,"FRAME")) frame=std::string("CMS");
  if (!reader->ReadFromFile(beam,"BEAM")) beam=std::string("P+");
  if (!reader->ReadFromFile(target,"TARGET")) target=std::string("PBAR-");
  if (!reader->ReadFromFile(win,"WIN")) win=1800.0;
  ATOOLS::rpa.gen.SetEcms(win);
  if (!reader->ReadFromFile(pysubs.msel,"MSEL")) pysubs.msel=1;
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
  pyinit(frame.c_str(),beam.c_str(),target.c_str(),win);
  if (ATOOLS::msg.Level()>2) ListLundParameters();
  delete reader;
}

Lund_Interface::~Lund_Interface()
{
}

void Lund_Interface::CleanUp()
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

bool Lund_Interface::ConvertEvent(ATOOLS::Blob_List *blobs)
{
  if (ATOOLS::msg.Level()>2) pylist(1);
  bool success=true;
  std::map<int,ATOOLS::Particle*> converted;
  blobs->clear();
  success=success&&ConvertParticles(converted);
  success=success&&ConstructBlobs(blobs,converted);
  success=success&&SetTypes(blobs);
  success=success&&DeleteObsolete(converted);
  return success;
}

bool Lund_Interface::Treat(ATOOLS::Blob_List *blobs,double &weight)
{
  pyevnt();
  pyhepc(1);
  weight*=pypars.pari[10];
  bool success=!ConvertEvent(blobs);
  return success;
} 

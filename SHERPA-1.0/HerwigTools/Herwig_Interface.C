#include "Herwig_Interface.H"

#include "Herwig_Wrapper.H"
#include "ISR_Handler.H"
#include "Data_Reader.H"
#include "Particle.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "Exception.H"
#include "Running_AlphaS.H"
#include "Scaling.H"

using namespace SHERPA;

bool Herwig_Interface::s_exportas=false;
bool Herwig_Interface::s_exportpdf=false;

size_t Herwig_Interface::s_errors=0;
size_t Herwig_Interface::s_maxerrors=0;

ATOOLS::Blob_List *Herwig_Interface::s_bloblist=NULL; 
PDF::ISR_Handler *Herwig_Interface::s_isrhandler=NULL; 

inline void MakeFortranString(char *output,std::string input,unsigned int length)
{
  for (unsigned int i=0;i<length;++i) output[i]=(char)32;
  for (size_t j=0;j<input.length();++j) output[j]=(char)input[j];
}


Herwig_Interface::Herwig_Interface(std::string _m_path,std::string _m_file,bool sherpa):
  m_path(_m_path),m_file(_m_file),
  p_hepevt(NULL), 
  m_compress(true),m_writeout(false),
  p_phep(new double[5*4000]),
  p_vhep(new double[4*4000]),
  p_jmohep(new int[2*4000]),
  p_jdahep(new int[2*4000])
{
  std::string beam1, beam2, pdfgroup;
  int pdfset;
  std::vector<std::vector<double> > help;
  ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader("=",";","!");
  reader->SetMatrixType(reader->MTransposed);
  reader->SetInputPath(m_path);
  reader->SetInputFile(m_file);
  reader->AddIgnore("(");
  reader->AddIgnore(")");
  reader->AddIgnore(",");
  if (!reader->ReadFromFile(beam1,"BEAM_1"))         beam1 = std::string("P");
  if (!reader->ReadFromFile(beam2,"BEAM_2"))         beam2 = std::string("PBAR");
  if (!reader->ReadFromFile(hwproc.pbeam1,"PBEAM1")) hwproc.pbeam1 = 980.0;
  if (!reader->ReadFromFile(hwproc.pbeam2,"PBEAM2")) hwproc.pbeam2 = 980.0;
  if (!reader->ReadFromFile(hwproc.iproc,"IPROC"))   hwproc.iproc  = 1471;
  if (!reader->ReadFromFile(pdfgroup,"AUTPDF"))      pdfgroup=std::string("");
  if (!reader->ReadFromFile(pdfset,"MODPDF"))        pdfset        = 1;
  MakeFortranString(hwbmch.part1,beam1,8);
  MakeFortranString(hwbmch.part2,beam2,8);
  ATOOLS::rpa.gen.SetEcms(hwproc.pbeam1+hwproc.pbeam2);
  if (pdfgroup!=std::string("")) {
    MakeFortranString(hwprch.autpdf[0],pdfgroup,20);
    MakeFortranString(hwprch.autpdf[1],pdfgroup,20);
    hwpram.modpdf[0]=pdfset;
    hwpram.modpdf[1]=pdfset;
  }
  else {
    hwpram.modpdf[0] = -1;
    hwpram.modpdf[1] = -1;
  }
  
  hwproc.maxev=ATOOLS::rpa.gen.NumberOfEvents();
#ifdef EXPORT__AlphaS
  int orderas;
  double asmz, asdef, mz;  
  reader->SetInputFile("Model.dat");
  if (!reader->ReadFromFile(orderas,"ORDER_ALPHAS"))  orderas = 0;
  if (!reader->ReadFromFile(asmz,"ALPHAS(MZ)"))       asmz    = 0.1188;
  if (!reader->ReadFromFile(asdef,"ALPHAS(default)")) asdef   = asmz;
  reader->SetInputFile("Particle.dat");
  if (!reader->ReadFromFile(mz,"24    ")) mz=91.188;
  MODEL::as = new MODEL::Running_AlphaS(asmz,mz*mz,orderas);
  MODEL::as->SetDefault(asdef);
#endif
  hwigin();
  hwuinc();
  hwaini();
  std::string outputname;
  if (reader->ReadFromFile(m_outfilename,"OUTPUT_FILE")) {
    if (!reader->ReadFromFile(m_evtsperfile,"EVENTS_PER_FILE")) m_evtsperfile=1000;
    NextFile(true);
  }
  int helper;
  if (!reader->ReadFromFile(helper,"COMPRESS")) helper=1;
  m_compress=(bool)helper;
  if (!sherpa) {
    p_hepevt = new ATOOLS::HepEvt_Interface(ATOOLS::gtp::Herwig);
  }
  delete reader;
}

Herwig_Interface::~Herwig_Interface()
{
  NextFile(false);
  if (p_hepevt) { 
    p_hepevt->SetNhep(0);
    p_hepevt->SetIsthep(NULL);
    p_hepevt->SetIdhep(NULL);
    p_hepevt->SetJmohep(NULL);
    p_hepevt->SetJdahep(NULL);
    p_hepevt->SetPhep(NULL);
    p_hepevt->SetVhep(NULL);
    delete p_hepevt; p_hepevt = NULL; 
  }
  if (p_jmohep) { delete p_jmohep; p_jmohep = NULL; }
  if (p_jdahep) { delete p_jdahep; p_jdahep = NULL; }
  if (p_phep)   { delete p_phep;   p_phep   = NULL; }
  if (p_vhep)   { delete p_vhep;   p_vhep   = NULL; }
}

void Herwig_Interface::NextFile(const bool newfile) 
{
  if (!m_writeout) return; 
  std::string oldfile;
  bool oldfileexists=false;
  std::ofstream *outfile=p_hepevt->GetOutStream();
  if (outfile!=NULL) {
    oldfileexists=true;
    oldfile=m_outfilename+ATOOLS::ToString(m_curfile)+std::string(".evts");
    if (newfile) 
      (*outfile)<<(m_outfilename+ATOOLS::ToString(++m_curfile)+std::string(".evts"))<<std::endl;
    if (m_compress) {
      system((std::string("gzip ")+oldfile+std::string(".gz ")+oldfile).c_str());
      system((std::string("rm ")+oldfile).c_str());
    }
  }
  if (!newfile) {
    if (p_hepevt) { 
      p_hepevt->SetNhep(0);
      p_hepevt->SetIsthep(NULL);
      p_hepevt->SetIdhep(NULL);
      p_hepevt->SetJmohep(NULL);
      p_hepevt->SetJdahep(NULL);
      p_hepevt->SetPhep(NULL);
      p_hepevt->SetVhep(NULL);
      delete p_hepevt; 
      p_hepevt = NULL; 
    }
    return;
  }
  std::string file = std::string(m_outfilename+ATOOLS::ToString(m_curfile)+std::string(".evts"));
  p_hepevt->ChangeOutStream(file,m_evtsperfile);
}

bool Herwig_Interface::OneEvent(ATOOLS::Blob_List * const blobs,double &weight)
{
  hwuine();
  hwepro();
  hwbgen();
  hwdhob();
  hwcfor();
  hwcdec();
  hwdhad();
  hwdhvy();
  hwmevt();
  hwufne();
  weight=1.; 
  for (int i=0;i<hepevt.nhep;i++) {
    for (int j=0;j<2;j++) {
      p_jmohep[2*i+j] = hepevt.jmohep[i][j]; 
      p_jdahep[2*i+j] = hepevt.jdahep[i][j];
    } 
    for (int j=0;j<5;j++) p_phep[5*i+j] = hepevt.phep[i][j];
    for (int j=0;j<4;j++) p_vhep[4*i+j] = hepevt.vhep[i][j];
  }
  p_hepevt->SetNhep(hepevt.nhep);
  p_hepevt->SetIsthep(hepevt.isthep);
  p_hepevt->SetIdhep(hepevt.idhep);
  p_hepevt->SetJmohep(p_jmohep);
  p_hepevt->SetJdahep(p_jdahep);
  p_hepevt->SetPhep(p_phep);
  p_hepevt->SetVhep(p_vhep);
  //if (ATOOLS::msg.LevelIsDebugging()) pylist(3);
  if (p_hepevt->HepEvt2Sherpa(blobs)) return true; 
  return false;
} 

bool Herwig_Interface::ConvertEvent(ATOOLS::Blob_List *bloblist)
{
  //if (ATOOLS::msg.LevelIsDebugging()) pylist(1);
  bool success=true;
  std::map<int,ATOOLS::Particle*> converted;
  bloblist->clear();
  success=success&&ConvertParticles(converted);
  success=success&&ConstructBlobs(bloblist,converted);
  success=success&&SetTypes(bloblist);
  success=success&&DeleteObsolete(converted);
  return success;
}


bool Herwig_Interface::ConvertParticles(std::map<int,ATOOLS::Particle*> &converted)
{
  for (int i=0;i<hepevt.nhep;++i) {
    if (ATOOLS::msg.LevelIsTracking()) {
      ATOOLS::msg.Debugging()<<i<<" "<<hepevt.isthep[i]<<" "<<hepevt.idhep[i]<<" "
			     <<hepevt.jmohep[i][0]<<" "
			     <<hepevt.jdahep[i][0]<<" "<<hepevt.jdahep[i][1]<<" "
			     <<"("<<hepevt.phep[i][0]<<","<<hepevt.phep[i][1]<<","
			     <<hepevt.phep[i][2]<<","<<hepevt.phep[i][3]<<") "
			     <<"("<<hepevt.vhep[i][0]<<","<<hepevt.vhep[i][1]<<","
			     <<hepevt.vhep[i][2]<<","<<hepevt.vhep[i][3]<<")"<<std::endl;
    }
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

bool Herwig_Interface::ConstructBlobs(ATOOLS::Blob_List *blobs,std::map<int,ATOOLS::Particle*> &converted)
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
	  productionblob->SetId();
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

bool Herwig_Interface::SetTypes(ATOOLS::Blob_List *blobs)
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

bool Herwig_Interface::DeleteObsolete(std::map<int,ATOOLS::Particle*> &converted)
{
  for (int i=0;i<hepevt.nhep;++i) {
    if ((converted[i]->ProductionBlob()==NULL)&&
	(converted[i]->DecayBlob()==NULL)) {
      delete converted[i];
    }
  }
  return true;
}

void Herwig_Interface::Error(const int error)
{
  ++s_errors;
  if (s_errors>s_maxerrors) {
    throw(ATOOLS::Exception(ATOOLS::ex::critical_error,std::string("Herwig calls PYERRM(")+
			    ATOOLS::ToString(error)+std::string(")"),"Herwig_Interface","Error"));
  }
  else {
    ATOOLS::msg.Error()<<"Herwig_Interface::Error("<<error<<") "<<ATOOLS::om::red
		       <<"Herwig calls PYERRM("<<error<<") in event "
		       <<ATOOLS::rpa.gen.NumberOfDicedEvents()<<"."
		       <<ATOOLS::om::reset<<std::endl;
    if (ATOOLS::msg.LevelIsDebugging()) {
      ATOOLS::msg.Tracking()<<*s_bloblist<<std::endl;
      //pylist(2);
    }
  }
}

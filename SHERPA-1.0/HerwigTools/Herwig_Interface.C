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

Herwig_Interface::Herwig_Interface(std::string _m_path,std::string _m_file,bool sherpa):
  m_path(_m_path),m_file(_m_file),
  p_hepevt(NULL), 
  m_compress(true),m_writeout(false),
  p_phep(new double[5*4000]),
  p_vhep(new double[4*4000]),
  p_jmohep(new int[2*4000]),
  p_jdahep(new int[2*4000])
{
  ReadInTheParameters();
  if (!sherpa) {
    p_hepevt = new ATOOLS::HepEvt_Interface(ATOOLS::gtp::Herwig);
  }
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

void Herwig_Interface::ReadInTheParameters()
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

  // HWBMCH common block
  if (!reader->ReadFromFile(beam1,"BEAM_1"))           beam1 = std::string("P");
  if (!reader->ReadFromFile(beam2,"BEAM_2"))           beam2 = std::string("PBAR");
  MakeFortranString(hwbmch.part1,beam1,8);
  MakeFortranString(hwbmch.part2,beam2,8);

  // HWPROC common block
  if (!reader->ReadFromFile(hwproc.pbeam1,"PBEAM1"))   hwproc.pbeam1 = 980.0;
  if (!reader->ReadFromFile(hwproc.pbeam2,"PBEAM2"))   hwproc.pbeam2 = 980.0;
  if (!reader->ReadFromFile(hwproc.iproc,"IPROC"))     hwproc.iproc  = 1471;
  hwproc.maxev=ATOOLS::rpa.gen.NumberOfEvents();
  ATOOLS::rpa.gen.SetEcms(hwproc.pbeam1+hwproc.pbeam2);

  // HWPRAM common block
  if (!reader->ReadFromFile(hwpram.qcdlam,"QCDLAM"))   hwpram.qcdlam   = 0.18;
  if (!reader->ReadFromFile(hwpram.vqcut,"VQCUT"))     hwpram.vqcut    = 0.48;
  if (!reader->ReadFromFile(hwpram.vgcut,"VGCUT"))     hwpram.vgcut    = 0.10;
  if (!reader->ReadFromFile(hwpram.vpcut,"VPCUT"))     hwpram.vpcut    = 0.40;
  if (!reader->ReadFromFile(hwpram.clmax,"CLMAX"))     hwpram.clmax    = 3.35;
  if (!reader->ReadFromFile(hwpram.clpow,"CLPOW"))     hwpram.clpow    = 2.;
  if (!reader->ReadFromFile(hwpram.qdiqk,"QDIQK"))     hwpram.qdiqk    = 0.;
  if (!reader->ReadFromFile(hwpram.pdiqk,"PDIQK"))     hwpram.pdiqk    = 5.;
  if (!reader->ReadFromFile(hwpram.qspac,"QSPAC"))     hwpram.qspac    = 2.5;
  if (!reader->ReadFromFile(hwpram.ptrms,"PTRMS"))     hwpram.ptrms    = 0.;
  if (!reader->ReadFromFile(hwpram.psplt[0],"PSPLT1")) hwpram.psplt[0] = 1.;
  if (!reader->ReadFromFile(hwpram.psplt[1],"PSPLT2")) hwpram.psplt[1] = hwpram.psplt[0];

  if (!reader->ReadFromFile(pdfgroup,"PDFGROUP_MC"))   pdfgroup=std::string("");
  if (!reader->ReadFromFile(pdfset,"PDFSET_MC"))       pdfset        = 1;
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

  // HWHARD common block 
  if (!reader->ReadFromFile(hwhard.ibrn[0],"HWSEED1")) hwhard.ibrn[0]  = 1246579;
  if (!reader->ReadFromFile(hwhard.ibrn[1],"HWSEED2")) hwhard.ibrn[1]  = 8447766;
  if (!reader->ReadFromFile(hwhard.ptmin,"PTMIN"))     hwhard.ptmin    = 20.;

  // HWPROP common block
  if (!reader->ReadFromFile(hwprop.rmass[2],"M_DOWN"))  hwprop.rmass[2]   = 0.32;
  if (!reader->ReadFromFile(hwprop.rmass[3],"M_UP"))    hwprop.rmass[3]   = 0.32;
  if (!reader->ReadFromFile(hwprop.rmass[4],"M_STRN"))  hwprop.rmass[4]   = 0.50;
  if (!reader->ReadFromFile(hwprop.rmass[5],"M_CHRM"))  hwprop.rmass[5]   = 1.55;
  if (!reader->ReadFromFile(hwprop.rmass[6],"M_BOTT"))  hwprop.rmass[6]   = 4.95;
  if (!reader->ReadFromFile(hwprop.rmass[7],"M_TOP"))   hwprop.rmass[7]   = 174.3;
  if (!reader->ReadFromFile(hwprop.rmass[14],"M_GLUE")) hwprop.rmass[14]  = 0.75;
  if (!reader->ReadFromFile(hwprop.rmass[199],"M_W"))   hwprop.rmass[199] = 80.42;
  if (!reader->ReadFromFile(hwprop.rmass[201],"M_Z"))   hwprop.rmass[201] = 91.188;
  if (!reader->ReadFromFile(hwprop.rmass[202],"M_H"))   hwprop.rmass[202] = 115.;

  // HWEVNT common block
  if (!reader->ReadFromFile(hwevnt.maxpr,"MAXPR"))      hwevnt.maxpr    = 0; 

  
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

  std::string outputname;
  if (reader->ReadFromFile(m_outfilename,"OUTPUT_FILE")) {
    if (!reader->ReadFromFile(m_evtsperfile,"EVENTS_PER_FILE")) m_evtsperfile=1000;
    NextFile(true);
  }
  int helper;
  if (!reader->ReadFromFile(helper,"COMPRESS")) helper=1;
  m_compress=(bool)helper;
  delete reader;
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

bool Herwig_Interface::Initialize() 
{ 
  hwcint(); 
  return true;
}

bool Herwig_Interface::OneEvent(ATOOLS::Blob_List * const blobs,double &weight)
{
  for (int trials=0;trials<5;++trials) {
    hwcgse();
    weight = hwevnt.evwgt; 
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
    if (p_hepevt->HepEvt2Sherpa(blobs)) {
      for (ATOOLS::Blob_List::const_iterator bit=blobs->begin(); bit!=blobs->end();++bit) {
	if ((*bit)->Type()==ATOOLS::btp::Signal_Process) { 
	  (*bit)->AddData("ME_Weight",new ATOOLS::Blob_Data<double>(weight));
	  return true;
	} 
      }
    }
    else {
      if (!blobs->empty()) {
	for (ATOOLS::Blob_List::const_iterator blit=blobs->begin();blit!=blobs->end();++blit) delete (*blit);
	blobs->clear();
      }
      ATOOLS::msg.Error()<<"Error in Herwig_Interface::OneEvent."<<std::endl
			 <<"   Could not translate HEPEVT common block into blobs."<<std::endl
			 <<"   Trials so far : "<<trials;
      if (trials<4) ATOOLS::msg.Error()<<", will continue."<<std::endl;
      else ATOOLS::msg.Error()<<", will return false, unknown result."<<std::endl;
    }
  }
  return false;
} 

void Herwig_Interface::Error(const int error)
{
  ++s_errors;
  if (s_errors>s_maxerrors) {
    throw(ATOOLS::Exception(ATOOLS::ex::critical_error,std::string("Herwig calls HWWARN(")+
			    ATOOLS::ToString(error)+std::string(")"),"Herwig_Interface","Error"));
  }
  else {
    ATOOLS::msg.Error()<<"Herwig_Interface::Error("<<error<<") "<<ATOOLS::om::red
		       <<"Herwig calls HWWARN("<<error<<") in event "
		       <<ATOOLS::rpa.gen.NumberOfDicedEvents()<<"."
		       <<ATOOLS::om::reset<<std::endl;
    if (ATOOLS::msg.LevelIsDebugging()) {
      ATOOLS::msg.Tracking()<<*s_bloblist<<std::endl;
    }
  }
}

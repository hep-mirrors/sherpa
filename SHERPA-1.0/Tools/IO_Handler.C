#include "IO_Handler.H"
#include "Message.H"
#include "Run_Parameter.H"
#include "Data_Read.H"
#include "Exception.H"
#include "MyStrStream.H"
#include <stdio.h>

using namespace SHERPA;
using namespace ATOOLS;

extern "C" {
  void outhepevt_();
}

const iotype::code SHERPA::operator|(const iotype::code code1,const iotype::code code2)
{
  return (iotype::code)((int)code1|(int)code2);
}
  
const iotype::code SHERPA::operator&(const iotype::code code1,const iotype::code code2)
{
  return (iotype::code)((int)code1&(int)code2);
}  

IO_Handler::IO_Handler(const std::string _mode):
  m_on(true), m_io(1), m_outtype(iotype::Unknown), m_intype(iotype::Unknown), 
#ifdef CLHEP_SUPPORT
  p_hepmc(NULL), 
#endif
  p_hepevt(NULL), 
  p_instream(NULL) 
{
  if (_mode==std::string("Sherpa")) {
    m_outtype = iotype::Sherpa;
  }
  else if (_mode==std::string("HepEvt")) {
    m_outtype = iotype::HepEvt;
    p_hepevt  = new HepEvt_Interface(true,10);
  }    
#ifdef CLHEP_SUPPORT
  else if (_mode==std::string("HepMC")) {
    m_outtype = iotype::HepMC; 
    p_hepmc   = new HepMC_Interface();
  }
#endif
  else {
    msg.Events()<<"Potential Error in IO_Handler::IO_Handler("<<m_outtype<<")"<<std::endl
		<<"   No output format specified. Continue run."<<std::endl;
    msg.LogFile()<<"Potential Error in IO_Handler::IO_Handler("<<m_outtype<<")"<<std::endl
		 <<"   No output format specified. Continue run."<<std::endl;
    m_on = false;
    m_io = 0;
  }
}

IO_Handler::IO_Handler(const std::vector<std::string> & outfiles,
		       const std::vector<std::string> & infiles,
		       const std::string _path, const int _filesize):
  m_on(true), m_outtype(iotype::Unknown), m_intype(iotype::Unknown), 
#ifdef CLHEP_SUPPORT
  p_hepmc(NULL), 
#endif
  p_hepevt(NULL), 
  p_instream(NULL),
  m_path(_path), 
  m_filesize(_filesize), 
  m_evtnumber(0), 
  m_evtcount(0)
{
  for (size_t i=0;i<infiles.size();++i) {
    m_intype=m_intype|(iotype::code)(pow(2.,int(i))*(infiles[i]!=std::string("")));
  }
  for (size_t i=0;i<outfiles.size();++i) {
    m_outtype=m_outtype|(iotype::code)(pow(2.,int(i))*(outfiles[i]!=std::string("")));
  }

  if (m_outtype!=iotype::Unknown && m_intype!=iotype::Unknown) {
    m_on=false;
    return; 
  }

  m_io=2*(m_outtype>0)+4*(m_intype>0);
  switch (m_outtype) {
  case iotype::Sherpa:
    m_file     = outfiles[0];
    m_filename = m_path+std::string("/")+m_file+std::string(".evts");
    m_outstream.open(m_filename.c_str());
    if (!m_outstream.good()) { 
      msg.Error()<<"ERROR in IO_Handler."<<std::endl
		 <<"   Could not open event file "<<m_filename<<"."<<std::endl
		 <<"   Will abort the run."<<std::endl;
      abort();
    }
    break;
#ifdef CLHEP_SUPPORT
  case iotype::HepMC: 
    p_hepmc = new HepMC_Interface();
    break;
#endif
  case iotype::HepEvt:
    p_hepevt = new HepEvt_Interface(true,1,m_path,outfiles[2],m_filesize);
    break;
  default :
    msg.LogFile()<<"Potential Error in IO_Handler::IO_Handler("<<m_outtype<<")"<<std::endl
		 <<"   No output format specified. Continue run."<<std::endl;
    break;
  }

  std::string gentype;
  switch (m_intype) {
  case iotype::Sherpa:
    m_filename = m_path+std::string("/")+infiles[0]+std::string(".evts"); 
    p_instream = new std::ifstream(m_filename.c_str()); 
    if (!p_instream->good()) {
      msg.Error()<<"ERROR in IO_Handler."<<std::endl
		 <<"   Event file "<<m_filename<<" not found."<<std::endl
		 <<"   Will abort the run."<<std::endl;
      abort();
    }
    (*p_instream)>>gentype>>m_filesize;
    if (gentype!=std::string("Sherpa")) {
      msg.Error()<<"ERROR in HepEvt_Interface."<<std::endl
		 <<"   Generator type:"<<gentype<<" not SHERPA."<<std::endl
		 <<"   Cannot guarantee i/o operations. Will abort the run."<<std::endl;
      abort();
    }
    break;
#ifdef CLHEP_SUPPORT
  case iotype::HepMC: 
    if (p_hepmc==NULL) p_hepmc  = new HepMC_Interface();
    break;
#endif
  case iotype::HepEvt:
    // Obacht !!!!!
    if (p_hepevt==NULL) p_hepevt = new HepEvt_Interface(false,1,m_path,infiles[2]);
    else abort();  // Schlamassel !!!!
    break;
  default :
    break;
  }
}

IO_Handler::~IO_Handler() 
{
#ifdef CLHEP_SUPPORT
  if (p_hepmc)  { delete p_hepmc;  p_hepmc  = NULL; }
#endif
  if (p_hepevt!=NULL) { delete p_hepevt; p_hepevt=NULL; }
  if (p_instream) {
    p_instream->close();
    delete p_instream; p_instream=NULL;
  }
}

void IO_Handler::AddOutputMode(const iotype::code c1) 
{
  if (m_on==false) return;
  msg.Error()<<"Error in IO_Handler::AddOutputMode("<<c1<<")"<<std::endl
	     <<"   Method not yet implemented. Continue run."<<std::endl;
}

void IO_Handler::AddInputMode(const iotype::code c1) 
{
  if (m_on==false) return;
  msg.Error()<<"Error in IO_Handler::AddInputMode("<<c1<<")"<<std::endl
	     <<"   Method not yet implemented. Continue run."<<std::endl;
}

void IO_Handler::PrintEvent() {
  if (m_on==false) return;
  switch (m_outtype) {
  case iotype::Sherpa:
    break;
#ifdef CLHEP_SUPPORT
  case iotype::HepMC: 
    p_hepmc->PrintHepMCEvent();
    break;
#endif 
  case iotype::HepEvt: 
    break;
  default:
    msg.Error()<<"Error in IO_Handler::OutputToFormat."<<std::endl
	       <<"   Unknown Output format : "<<m_outtype<<std::endl
	       <<"   No output, continue run ..."<<std::endl;
    break;
  }
}

bool IO_Handler::OutputToFormat(ATOOLS::Blob_List *const blobs,const double weight) 
{
  //std::cout<<"Output : "<<m_io<<" : "<<m_outtype<<std::endl;
  if (m_on==false)          return false;
  if (!(m_io&1)&&!(m_io&2)) return false;
  if (m_io&1 && msg.LevelIsEvents()) {
    if (!blobs->empty()) {
      msg.Out()<<"  -------------------------------------------------  "<<std::endl;
      for (Blob_Iterator blit=blobs->begin();blit!=blobs->end();++blit) 
	msg.Out()<<*(*blit)<<std::endl;
      msg.Out()<<"  -------------------------------------------------  "<<std::endl;
    }
    else {
      msg.Out()<<"  ******** Empty event ********  "<<std::endl;
      return false;
    }
  }

  for (int i=1;i<(int)iotype::size;i*=2) {
    if (m_outtype&i) {
      switch (m_outtype) {
      case iotype::Sherpa:
	if (m_io&2) SherpaOutput(blobs,weight); 
	break;
#ifdef CLHEP_SUPPORT
      case iotype::HepMC: 
	p_hepmc->Sherpa2HepMC(blobs);
	if (m_io&1 && msg.LevelIsDebugging()) {
	  p_hepmc->PrintHepMCEvent();
	}
	break;
#endif 
      case iotype::HepEvt: 
	p_hepevt->Sherpa2HepEvt(blobs); 
	if (m_io&1 && msg.LevelIsDebugging()) {
	  p_hepevt->PrintHepEvtEvent(p_hepevt->Nhep());
	}
	break;
      default:
	msg.Error()<<"Error in IO_Handler::OutputToFormat."<<std::endl
		   <<"   Unknown Output format : "<<m_outtype<<std::endl
		   <<"   No output, continue run ..."<<std::endl;
	break;
      }
    }
  }
  return false;
}

bool IO_Handler::InputFromFormat(ATOOLS::Blob_List *const blobs) 
{
  if (m_on==false) return false;
  if (!(m_io&4)) return false;
  switch (m_intype) {
  case iotype::Sherpa: return SherpaInput(blobs); 
#ifdef CLHEP_SUPPORT
  case iotype::HepMC: 
    throw(ATOOLS::Exception(ATOOLS::ex::not_implemented,
			    "Reading input from HepMC is not yet possible.",
			    "IO_Handler","InputFromFormat"));
#endif
  case iotype::HepEvt: return p_hepevt->HepEvt2Sherpa(blobs); 
  default:
    msg.Error()<<"Error in IO_Handler::InputFromFormat."<<std::endl
	       <<"   Unknown Input format : "<<m_intype<<std::endl
	       <<"   No input, continue run ... ."<<std::endl;
    break;
  }
  return false;
}

/*------------------------------------------------------------------
  Sherpa-specific I/O methods : Output is ASCII-format
  ------------------------------------------------------------------*/

void IO_Handler::SherpaOutput(ATOOLS::Blob_List *const blobs,const double weight) 
{ 
  if (m_on==false) return;
  ATOOLS::Particle_Int_Map           P2I;
  ATOOLS::Particle_Int_Map::iterator P2Iiter;
  
  for (Blob_Iterator blit=blobs->begin();blit!=blobs->end();++blit) {
    for (int i=0;i<(*blit)->NInP();i++) {
      if (P2I.find((*blit)->InParticle(i))==P2I.end()) 
	P2I.insert(std::make_pair((*blit)->InParticle(i),P2I.size()+1));
    }
    for (int i=0;i<(*blit)->NOutP();i++) {
      if (P2I.find((*blit)->OutParticle(i))==P2I.end()) 
	P2I.insert(std::make_pair((*blit)->OutParticle(i),P2I.size()+1));
    }
  }
  m_evtnumber++;
  m_evtcount++;
  m_outstream<<m_evtnumber<<" "<<P2I.size()<<" "<<blobs->size()<<" "<<weight<<std::endl;
  Particle * part;
  int kfc;
  for (P2Iiter=P2I.begin();P2Iiter!=P2I.end();P2Iiter++) {
    part = P2Iiter->first;
    kfc  = part->Flav().Kfcode(); if (part->Flav().IsAnti()) kfc=-kfc;
    m_outstream<<P2Iiter->second<<" "<<part->Status()<<" "<<part->Info()<<" "<<kfc<<" "
	       <<" "<<part->Momentum()[0]<<" "<<part->Momentum()[1]
	       <<" "<<part->Momentum()[2]<<" "<<part->Momentum()[3]<<" \n";
  }
  for (Blob_Iterator blit=blobs->begin();blit!=blobs->end();++blit) {
    m_outstream<<(*blit)->Id()<<" "<<(*blit)->Status()<<" "<<(int)(*blit)->Type()<<" "<<(*blit)->TypeSpec()
	       <<" "<<(*blit)->NInP()<<" "<<(*blit)->NOutP()<<" \n"
	       <<" "<<(*blit)->Position()[0]<<" "<<(*blit)->Position()[1]
	       <<" "<<(*blit)->Position()[2]<<" "<<(*blit)->Position()[3]<<" \n";
    for (int i=0;i<(*blit)->NInP();i++)  m_outstream<<P2I.find((*blit)->InParticle(i))->second<<" ";
    for (int i=0;i<(*blit)->NOutP();i++) m_outstream<<P2I.find((*blit)->OutParticle(i))->second<<" ";
    m_outstream<<" \n";
  }
  if (m_evtcount%m_filesize==0) {
    m_evtcount = 0;
    m_filename = m_file+ToString(int(m_evtnumber/m_filesize));
    m_outstream<<m_filename<<" \n";
    m_outstream.close();
    m_filename = m_path+std::string("/")+m_filename+std::string(".evts");
    m_outstream.open(m_filename.c_str());
    if (!m_outstream.good()) { 
      msg.Error()<<"ERROR in IO_Handler."<<std::endl
		 <<"   Could not open event file "<<m_filename<<"."<<std::endl
		 <<"   Will abort the run."<<std::endl;
      abort();
    }
  }
}

bool IO_Handler::SherpaInput(ATOOLS::Blob_List *const blobs) 
{ 
  if (m_on==false) return false;
  blobs->clear();

  m_evtcount++;
  int panumber, blnumber, weight;
  (*p_instream)>>m_evtnumber>>panumber>>blnumber>>weight;

  ATOOLS::Int_Particle_Map           I2P;
  ATOOLS::Int_Particle_Map::iterator I2Piter;

  int        paid, status, kfc;
  double     mom[4];
  Vec4D      momentum;
  char       info;
  Flavour    flav;
  Particle * part;
  
  for (int i=0;i<panumber;i++) {
    (*p_instream)>>paid>>status>>info>>kfc>>mom[0]>>mom[1]>>mom[2]>>mom[3];
    flav     = Flavour(kf::code(abs(kfc))); if (kfc<0) flav=flav.Bar();
    momentum = Vec4D(mom[0],mom[1],mom[2],mom[3]);
    part     = new Particle(paid,flav,momentum);
    part->SetStatus(status);
    part->SetInfo(info);
    I2P.insert(std::make_pair(paid,part));
  }

  int         type, ninp, noutp, blid;
  std::string typespec;
  double      pos[4];
  Vec4D       position;
  Blob      * blob;

  for (int i=0;i<blnumber;i++) {
    (*p_instream)>>blid>>status>>type>>typespec>>ninp>>noutp;
    (*p_instream)>>pos[0]>>pos[1]>>pos[2]>>pos[3];

    position = Vec4D(pos[0],pos[1],pos[2],pos[3]);
    blob     = new Blob(position,blid);
    blob->SetStatus(status);
    blob->SetType((btp::code)type);
    blob->SetTypeSpec(typespec);
    for (int i=0;i<ninp;i++) {
      (*p_instream)>>paid;
      I2Piter = I2P.find(paid);
      if (I2Piter==I2P.end()) {
	msg.Error()<<"Error in IO_Handler::SherpaInput."<<std::endl
		   <<"   Particle with number "<<paid<<" not found in event."<<std::endl
		   <<"   Will return false, continue & hope for the best."<<std::endl;
	blobs->clear();
	return false;
      }
      blob->AddToInParticles(I2Piter->second);
    }
    for (int i=0;i<noutp;i++) {
      (*p_instream)>>paid;
      I2Piter = I2P.find(paid);
      if (I2Piter==I2P.end()) {
	msg.Error()<<"Error in IO_Handler::SherpaInput."<<std::endl
		   <<"   Particle with number "<<paid<<" not found in event."<<std::endl
		   <<"   Will return false, continue & hope for the best."<<std::endl;
	blobs->clear();
	return false;
      }
      blob->AddToOutParticles(I2Piter->second);
    }

    blobs->push_back(blob);
  }

  if (m_evtcount%m_filesize==0) {
    std::string file, filename;
    (*p_instream)>>file;
    p_instream->close();
    filename =  m_path+std::string("/")+file+std::string(".evts"); 
    delete p_instream;
    p_instream = new std::ifstream(filename.c_str()); 
    if (!p_instream->good()) {
      msg.Error()<<"ERROR in HepEvt_Interface."<<std::endl
		 <<"   Event file "<<filename<<" not found."<<std::endl
		 <<"   Will abort the run."<<std::endl;
      abort();
    }
    std::string gentype;
    (*p_instream)>>gentype>>m_filesize;
    if (gentype!=std::string("Sherpa")) {
      msg.Error()<<"ERROR in HepEvt_Interface."<<std::endl
		 <<"   Generator type:"<<gentype<<" not SHERPA."<<std::endl
		 <<"   Cannot guarantee i/o operations. Will abort the run."<<std::endl;
      abort();
    }
    m_evtcount=0;
  }
  return true;
}

#include "Output_Handler.H"
#include "Message.H"
#include "Run_Parameter.H"
#include "Data_Read.H"

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

Output_Handler::Output_Handler():
  m_io(0), m_outtype(iotype::Unknown), m_intype(iotype::Unknown), 
#ifdef _USE_HEPMC_
  p_hepmc(NULL), p_event(NULL),
#endif
  p_hepevt(NULL), p_instream(NULL) {}

Output_Handler::Output_Handler(const std::vector<std::string> & outfiles,
			       const std::vector<std::string> & infiles):
  m_outtype(iotype::Unknown), m_intype(iotype::Unknown), 
#ifdef _USE_HEPMC_
  p_hepmc(NULL), p_event(NULL),
#endif
  p_hepevt(NULL), p_instream(NULL)
{
  for (size_t i=0;i<infiles.size();++i) {
    m_intype=m_intype|(iotype::code)(pow(2,i)*(infiles[i]!=std::string("")));
  }
  for (size_t i=0;i<outfiles.size();++i) {
    m_outtype=m_outtype|(iotype::code)(pow(2,i)*(outfiles[i]!=std::string("")));
  }

  m_io=(m_outtype>0)+2*(m_intype>0);
  switch (m_outtype) {
  case iotype::Sherpa:
    m_path     = std::string(".");
    m_filename = m_path+std::string("/")+outfiles[0]; 
    m_outstream.open(m_filename.c_str());
    if (!m_outstream.good()) { 
      msg.Error()<<"ERROR in Output_Handler."<<std::endl
		 <<"   Could not open event file "<<m_filename<<"."<<std::endl
		 <<"   Will abort the run."<<std::endl;
      abort();
    }
    break;
#ifdef _USE_HEPMC_
  case iotype::HepMC: 
    p_hepmc  = new HepMC_Interface();
    break;
#endif
  case iotype::HepEvt:
    p_hepevt = new HepEvt_Interface(true,1,std::string("."),outfiles[2]);
    break;
  default :
    msg.Error()<<"Potential Error in Output_Handler::Output_Handler("<<m_outtype<<")"<<std::endl
	       <<"   No output format specified. Continue run."<<std::endl;
    break;
  }
  switch (m_intype) {
  case iotype::Sherpa:
    std::cout<<"Try to init Sherpa-read in."<<std::endl;
    m_path     = std::string(".");
    m_filename = m_path+std::string("/")+infiles[0]; 
    p_instream = new std::ifstream(m_filename.c_str()); 
    if (!p_instream->good()) {
      msg.Error()<<"ERROR in Output_Handler."<<std::endl
		 <<"   Event file "<<m_filename<<" not found."<<std::endl
		 <<"   Will abort the run."<<std::endl;
      abort();
    }
    break;
#ifdef _USE_HEPMC_
  case iotype::HepMC: 
    if (p_hepmc==NULL) p_hepmc  = new HepMC_Interface();
    break;
#endif
  case iotype::HepEvt:
    // Obacht !!!!!
    if (p_hepevt==NULL) p_hepevt = new HepEvt_Interface(false,1,std::string("."),infiles[2]);
    else abort();  // Schlamassel !!!!
    break;
  default :
    break;
  }
}

Output_Handler::~Output_Handler() 
{
#ifdef _USE_HEPMC_
  if (p_hepmc!=NULL)  { delete p_hepmc;  p_hepmc=NULL;  }
  if (p_event!=NULL)  { delete p_event;  p_event=NULL;  }
#endif
  if (p_hepevt!=NULL) { delete p_hepevt; p_hepevt=NULL; }
  if (p_instream) {
    p_instream->close();
    delete p_instream; p_instream=NULL;
  }
}

void Output_Handler::AddOutputMode(const iotype::code c1) 
{
  msg.Error()<<"Error in Output_Handler::AddOutputMode("<<c1<<")"<<std::endl
	     <<"   Method not yet implemented. Continue run."<<std::endl;
}

void Output_Handler::AddInputMode(const iotype::code c1) 
{
  msg.Error()<<"Error in Output_Handler::AddInputMode("<<c1<<")"<<std::endl
	     <<"   Method not yet implemented. Continue run."<<std::endl;
}


bool Output_Handler::OutputToFormat(ATOOLS::Blob_List *const blobs,const double weight) 
{
  if (!(m_io&1)) return 0;
  for (int i=1;i<(int)iotype::size;i*=2) {
    if (m_outtype&i) {
      switch (m_outtype) {
      case iotype::Sherpa:
	if (!blobs->empty()) { SherpaOutput(blobs,weight); return 1; }
	else { 
	  msg.Error()<<"Error in Output_Handler::OutputToFormat."<<std::endl
		     <<"   empty bloblist."<<std::endl
		     <<"   No output, continue run ..."<<std::endl;
	  break;
	}
#ifdef _USE_HEPMC_
      case iotype::HepMC: 
	if (p_event) { delete p_event; p_event = NULL; }
	p_event = new HepMC::GenEvent();
	p_hepmc->Sherpa2HepMC(_blobs,p_event);
	p_event->print();
	return 1;
#endif
      case iotype::HepEvt: 
	p_hepevt->Sherpa2HepEvt(blobs); return 1;
      default:
	msg.Error()<<"Error in Output_Handler::OutputToFormat."<<std::endl
		   <<"   Unknown Output format : "<<m_outtype<<std::endl
		   <<"   No output, continue run ..."<<std::endl;
	break;
      }
    }
  }
  return 0;
}

bool Output_Handler::InputFromFormat(ATOOLS::Blob_List *const blobs) 
{
  if (!(m_io&2)) return 0;
  switch (m_intype) {
  case iotype::Sherpa: return SherpaInput(blobs); 
#ifdef _USE_HEPMC_
  case iotype::HepMC: 
    //return p_hepmc->Ouput(blobs); 
#endif
  case iotype::HepEvt: return p_hepevt->HepEvt2Sherpa(blobs); 
  default:
    msg.Error()<<"Error in Output_Handler::InputFromFormat."<<std::endl
	       <<"   Unknown Input format : "<<m_intype<<std::endl
	       <<"   No input, continue run ... ."<<std::endl;
    break;
  }
  return 0;
}

/*------------------------------------------------------------------
  Sherpa-specific I/O methods : Output is ASCII-format
  ------------------------------------------------------------------*/

void Output_Handler::SherpaOutput(ATOOLS::Blob_List *const blobs,const double weight) 
{ 
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
  //  for (P2Iiter=P2I.begin();P2Iiter!=P2I.end();P2Iiter++) {
  //    std::cout<<P2Iiter->first->Flav()<<" -> "<<P2Iiter->second<<std::endl;
  //  }
  m_outstream<<P2I.size()<<" "<<blobs->size()<<" "<<weight<<std::endl;
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

  m_outstream.close();
  p_instream = new std::ifstream(m_filename.c_str()); 
  SherpaInput(blobs);
}

bool Output_Handler::SherpaInput(ATOOLS::Blob_List *const blobs) 
{ 
  blobs->clear();

  int panumber, blnumber, weight;
  (*p_instream)>>panumber>>blnumber>>weight;

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
	msg.Error()<<"Error in Output_Handler::SherpaInput."<<std::endl
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
	msg.Error()<<"Error in Output_Handler::SherpaInput."<<std::endl
		   <<"   Particle with number "<<paid<<" not found in event."<<std::endl
		   <<"   Will return false, continue & hope for the best."<<std::endl;
	blobs->clear();
	return false;
      }
      blob->AddToOutParticles(I2Piter->second);
    }

    blobs->push_back(blob);
  }
  std::cout<<(*blobs)<<std::endl;
}

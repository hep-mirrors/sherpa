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
  p_hepevt(NULL) {}

Output_Handler::Output_Handler(const std::vector<std::string> & outfiles,
			       const std::vector<std::string> & infiles):
  m_outtype(iotype::Unknown), m_intype(iotype::Unknown), 
#ifdef _USE_HEPMC_
  p_hepmc(NULL), p_event(NULL),
#endif
  p_hepevt(NULL)
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
    // new sherpa output
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
    // new sherpa output
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

bool Output_Handler::OutputToFormat(ATOOLS::Blob_List *const blobs) 
{
  if (!(m_io&1)) return 0;
  for (int i=1;i<(int)iotype::size;i*=2) {
    if (m_outtype&i) {
      switch (m_outtype) {
      case iotype::Sherpa:
	if (!blobs->empty()) {
	  for (Blob_Iterator blit=blobs->begin();blit!=blobs->end();++blit) {
	    msg.Events()<<(*blit)<<std::endl;
	  }
	}
#ifdef _USE_HEPMC_
      case iotype::HepMC: 
	if (p_event) { delete p_event; p_event = NULL; }
	p_event = new HepMC::GenEvent();
	p_hepmc->Sherpa2HepMC(_blobs,p_event);
	p_event->print();
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
  case iotype::Sherpa:
    //return Sherpa_Output(blobs); 
#ifdef _USE_HEPMC_
  case iotype::HepMC: 
    //return p_hepmc->Ouput(blobs); 
#endif
  case iotype::HepEvt:
    //if (ATOOLS::msg.Level()>=1) outhepevt_();
    return p_hepevt->HepEvt2Sherpa(blobs); 
  default:
    msg.Error()<<"Error in Output_Handler::InputFromFormat."<<std::endl
	       <<"   Unknown Input format : "<<m_intype<<std::endl
	       <<"   No input, continue run ... ."<<std::endl;
    break;
  }
  return 0;
}

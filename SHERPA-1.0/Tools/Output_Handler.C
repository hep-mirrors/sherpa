#include "Output_Handler.H"
#include "Message.H"

using namespace SHERPA;
using namespace ATOOLS;

extern "C" {
  void outhepevt_();
}

Output_Handler::Output_Handler(int type) :
#ifdef _USE_HEPMC_
  p_hepmc(NULL), 
  p_event(NULL),  
#endif
  p_hepevt(NULL), 
  m_active(0),
  m_type(type)
{
  switch (m_type) {
  case 1: 
#ifdef _USE_HEPMC_
    p_hepmc  = new HepMC_Interface();
    m_active = 1;
#endif
    return;
  case 2:
    p_hepevt = new HepEvt_Interface();
    m_active = 1;
    return;
  default :
    msg.Error()<<"Potential Error in Output_Handler::Output_Handler("<<type<<")"<<std::endl
	       <<"   No output format specified. Continue run."<<std::endl;
  }
}

Output_Handler::~Output_Handler() {
#ifdef _USE_HEPMC_
  if (p_hepmc)  { delete p_hepmc;  p_hepmc  = NULL; }
#endif
  if (p_hepevt) { delete p_hepevt; p_hepevt = NULL; }
}

void Output_Handler::OutputToFormat(Blob_List * _blobs)
{
  if (!m_active) return;
  switch (m_type) {
  case 1: 
#ifdef _USE_HEPMC_
    if (p_event) { delete p_event; p_event = NULL; }
    p_event = new HepMC::GenEvent();
    p_hepmc->Sherpa2HepMC(_blobs,p_event);
    p_event->print();
#endif
    return;
  case 2:
    p_hepevt->Sherpa2HepEvt(_blobs);
    if (ATOOLS::msg.Level()>=1) outhepevt_();
    return;
  default:
    msg.Error()<<"Potential Error in Output_Handler::OutputToFormat("<<m_type<<")"<<std::endl
	       <<"   No output format specified. Continue run."<<std::endl;
  }
}


#include "Output_Handler.H"
#include "Message.H"

using namespace SHERPA;
using namespace APHYTOOLS;
using namespace AORGTOOLS;

Output_Handler::Output_Handler(int type) :
  m_type(type), p_hepmc(NULL), p_event(NULL), m_active(0)
{
  switch (m_type) {
  case 1: 
    p_hepmc  = new HepMC_Interface();
    m_active = 1;
    return;
  default :
    msg.Error()<<"Potential Error in Output_Handler::Output_Handler("<<type<<")"<<endl
	       <<"   No output format specified. Continue run."<<endl;
  }
}

Output_Handler::~Output_Handler() {
  if (p_hepmc) { delete p_hepmc; p_hepmc = NULL; }
}

void Output_Handler::OutputToFormat(Blob_List * _blobs)
{
  if (!m_active) return;
  switch (m_type) {
  case 1: 
    if (p_event) { delete p_event; p_event = NULL; }
    p_event = new HepMC::GenEvent();
    p_hepmc->Sherpa2HepMC(_blobs,p_event);
    p_event->print();
    return;
  default:
    msg.Error()<<"Potential Error in Output_Handler::OutputToFormat("<<m_type<<")"<<endl
	       <<"   No output format specified. Continue run."<<endl;
  }
}

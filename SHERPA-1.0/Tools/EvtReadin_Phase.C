#include "EvtReadin_Phase.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;



EvtReadin_Phase::EvtReadin_Phase(IO_Handler * _iohandler) :
  p_iohandler(_iohandler), p_evtreader(NULL), m_path(".")
{
  m_type = eph::Read_In;
}    

EvtReadin_Phase::EvtReadin_Phase(Event_Reader * _evtreader) :
  p_iohandler(NULL), p_evtreader(_evtreader), m_path(_evtreader->GetPath())
{
  m_type = eph::Read_In;
}    

bool EvtReadin_Phase::Treat(ATOOLS::Blob_List * _bloblist, double & weight) 
{
  if (!m_read) {
    if (p_iohandler) {
      p_iohandler->InputFromFormat(_bloblist);
      //weight = p_iohandler->Weight();
    }
    if (p_evtreader) {
      p_evtreader->FillBlobs(_bloblist);
      weight = p_evtreader->Weight();
    }
    m_read = true;
  }
  return (!m_read);
}


void EvtReadin_Phase::CleanUp()                   { m_read = false; }
void EvtReadin_Phase::Finish(const std::string &) {}


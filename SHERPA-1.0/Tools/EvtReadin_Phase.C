#include "EvtReadin_Phase.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;



EvtReadin_Phase::EvtReadin_Phase(Event_Reader * _evtreader) :
  p_evtreader(_evtreader), m_path(_evtreader->GetPath())
{
  m_type = std::string("Read-in");
}    

bool EvtReadin_Phase::Treat(ATOOLS::Blob_List * _bloblist, double & weight) 
{
  if (!m_read) {
    //std::cout<<"In EvtReadin_Phase::Treat"<<std::endl;
    p_evtreader->FillBlobs(_bloblist);
    //what about the weight ???
    weight = p_evtreader->Weight();
    m_read = true;
  }
  return (!m_read);
}


void EvtReadin_Phase::CleanUp()                   { m_read = false; }
void EvtReadin_Phase::Finish(const std::string &) {}


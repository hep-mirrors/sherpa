#include "SHERPA/Single_Events/EvtReadin_Phase.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;



EvtReadin_Phase::EvtReadin_Phase(Event_Reader_Base * evtreader) :
  p_evtreader(evtreader), m_path(evtreader->GetPath())
{
  m_type = eph::Read_In;
}    

Return_Value::code EvtReadin_Phase::Treat(ATOOLS::Blob_List * bloblist, 
					  double & weight) 
{
  if (!m_read) {
    if (p_evtreader) {
      p_evtreader->FillBlobs(bloblist);
      weight = p_evtreader->Weight();
    }
    m_read = true;
    return Return_Value::Success;
  }
  return Return_Value::Nothing;
}


void EvtReadin_Phase::CleanUp(const size_t & mode) { m_read = false; }
void EvtReadin_Phase::Finish(const std::string &)  {}


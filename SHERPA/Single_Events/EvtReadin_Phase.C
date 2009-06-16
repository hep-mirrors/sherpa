#include "SHERPA/Single_Events/EvtReadin_Phase.H"
#include "SHERPA/Tools/Input_Output_Handler.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;



EvtReadin_Phase::EvtReadin_Phase(Input_Output_Handler * iohandler) :
  p_iohandler(iohandler), p_evtreader(NULL), m_path("."), m_read(false)
{
  m_type = eph::Read_In;
}    

EvtReadin_Phase::EvtReadin_Phase(Event_Reader * evtreader) :
  p_iohandler(NULL), p_evtreader(evtreader), m_path(evtreader->GetPath())
{
  m_type = eph::Read_In;
}    

Return_Value::code EvtReadin_Phase::Treat(ATOOLS::Blob_List * bloblist, 
					  double & weight) 
{
  if (!m_read) {
    if (p_iohandler) p_iohandler->InputFromFormat(bloblist);
    if (p_evtreader) {
      p_evtreader->FillBlobs(bloblist);
      weight = p_evtreader->Weight();
    }
    m_read = true;
    return Return_Value::Success;
  }
  return Return_Value::Nothing;
}


void EvtReadin_Phase::CleanUp()                   { m_read = false; }
void EvtReadin_Phase::Finish(const std::string &) {}


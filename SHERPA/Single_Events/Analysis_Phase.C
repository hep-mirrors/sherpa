#include "SHERPA/Single_Events/Analysis_Phase.H"

#include "SHERPA/Tools/Analysis_Interface.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

#include <limits>

using namespace SHERPA;
using namespace ATOOLS;

Analysis_Phase::Analysis_Phase(Analysis_Interface *const analysis):
  Event_Phase_Handler(std::string("Analysis")), 
  p_analysis(analysis), m_wit(std::numeric_limits<size_t>::max()), m_init(0)
{
  m_type=eph::Analysis;
  Data_Reader read(" ",";","!","=");
  double wit;
  if (read.ReadFromFile(wit,"ANALYSIS_WRITEOUT_INTERVAL")) {
    if (wit<1.0) {
      if (wit*rpa.gen.NumberOfEvents()>1.0)
        m_wit=(size_t)(wit*rpa.gen.NumberOfEvents());
    }
    else m_wit=(size_t)(wit);
    msg_Info()<<METHOD<<"(): Set writeout interval "<<m_wit<<" events.\n";
  }
}

Return_Value::code Analysis_Phase::Treat(Blob_List *bloblist,double &weight) 
{
  if (!m_init) m_init=p_analysis->Init();
  if (!bloblist->empty()) p_analysis->Run(bloblist);
  if (rpa.gen.NumberOfDicedEvents()%m_wit==0 &&
      rpa.gen.NumberOfDicedEvents()<rpa.gen.NumberOfEvents()) 
    p_analysis->WriteOut();
  return Return_Value::Nothing;
}

void Analysis_Phase::CleanUp() 
{
  p_analysis->CleanUp();
}

void Analysis_Phase::Finish(const std::string &path)
{
  p_analysis->Finish();
}

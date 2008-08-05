#include "Analysis_Phase.H"

#include "Analysis_Handler.H"
#include "Run_Parameter.H"
#include "Data_Reader.H"

#include <limits>

using namespace SHERPA;
using namespace ATOOLS;

Analysis_Phase::Analysis_Phase(ANALYSIS::Analysis_Handler *const analysis):
  Event_Phase_Handler(std::string("Analysis")), 
  p_analysis(analysis), m_wit(std::numeric_limits<size_t>::max())
{
  m_type=eph::Analysis;
  Data_Reader read(" ",";","!","=");
  double wit;
  if (read.ReadFromFile(wit,"ANALYSIS_WRITEOUT_INTERVAL")) {
    if (wit<1.0) m_wit=(size_t)(wit*rpa.gen.NumberOfEvents());
    else m_wit=(size_t)(wit);
    msg_Info()<<METHOD<<"(): Set writeout interval "<<m_wit<<" events.\n";
  }
}

Return_Value::code Analysis_Phase::Treat(Blob_List *bloblist,double &weight) 
{
  if (!bloblist->empty()) p_analysis->DoAnalysis(bloblist,weight);
  if (rpa.gen.NumberOfDicedEvents()%m_wit==0 &&
      rpa.gen.NumberOfDicedEvents()<rpa.gen.NumberOfEvents()) 
    p_analysis->WriteOut();
  return Return_Value::Nothing;
}

void Analysis_Phase::CleanUp() 
{
  p_analysis->Clear();
}

void Analysis_Phase::Finish(const std::string &path)
{
  p_analysis->Finish(path);
}

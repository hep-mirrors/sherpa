#include "Analysis_Phase.H"
#include "Analysis_Handler.H"

#ifdef PROFILE__all
#define PROFILE__Analysis_Phase
#endif
#ifdef PROFILE__Analysis_Phase
#include "prof.hh"
#else 
#define PROFILE_HERE 
#define PROFILE_LOCAL(LOCALNAME) 
#endif

using namespace SHERPA;
using namespace ATOOLS;

Analysis_Phase::Analysis_Phase(ANALYSIS::Analysis_Handler *const analysis):
  Event_Phase_Handler(std::string("Analysis")), 
  p_analysis(analysis) 
{
  m_type=eph::Analysis;
}

Return_Value::code Analysis_Phase::Treat(Blob_List *bloblist,double &weight) 
{
  PROFILE_HERE;
  if (!bloblist->empty()) p_analysis->DoAnalysis(bloblist,weight);
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

#include "Analysis_Phase.H"
#include "MyStrStream.H"
#include "Run_Parameter.H"
#include "Message.H"
#include <vector>

#include "Jet_Observables.H"
#include "One_Particle_Observables.H"
#include "Two_Particle_Observables.H"
#include "Four_Particle_Observables.H"

#include "Primitive_Detector.H"
#include "Final_Selector.H"

#ifdef PROFILE__Analysis_Phase
#include "prof.hh"
#else 
#define PROFILE_HERE {}
#define PROFILE_LOCAL(LOCALNAME) {}
#endif
using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Analysis_Phase::Analysis_Phase() :
  Event_Phase_Handler(std::string("Analysis")), p_analysis(NULL) 
{ 
  m_type = eph::Analysis;
}

Analysis_Phase::Analysis_Phase(Sample_Analysis * ana,const std::string  & iter) :
  Event_Phase_Handler(std::string("Analysis:")+ana->Phase()), 
  p_analysis(ana) 
{
  if (iter!=std::string("")) m_name += std::string("_")+iter;
  m_type = eph::Analysis;
}

bool Analysis_Phase::Treat(ATOOLS::Blob_List * blist, double & weight) 
{
  PROFILE_LOCAL("Analysis_Phase::Treat");
  if (!blist->empty()) p_analysis->DoAnalysis(blist,weight);
  return 0;
}

void  Analysis_Phase::CleanUp() 
{
  p_analysis->Clear();
}

void Analysis_Phase::Finish(const std::string & dirname)
{
  p_analysis->Finish();
}

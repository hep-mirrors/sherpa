#include "Sherpa.H"
#include "Analysis_Phase.H"
#include "EvtReadin_Phase.H"
#include "Signal_Processes.H"
#include "Hard_Decays.H"
#include "Multiple_Interactions.H"
#include "Jet_Evolution.H"
#include "Hadronization.H"
#include "Message.H"

#ifdef PROFILE__Sherpa
#include "prof.hh"
#else 
#define PROFILE_HERE {}
#define PROFILE_LOCAL(LOCALNAME) {}
#endif

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

namespace SHERPA {
  Sherpa generator;
}

extern "C" {
  void apainit_() {
    SHERPA::generator.InitializeTheRun(0,NULL);
    SHERPA::generator.InitializeTheEventHandler();
  }
  void aparun_() {
    SHERPA::generator.GenerateOneEvent();
  }
  void apaend_() {
    SHERPA::generator.SummarizeRun();
  }
}

Sherpa::Sherpa() :
  p_inithandler(NULL), p_eventhandler(NULL), p_output(NULL), p_analysis(NULL)  
{
  PROFILE_HERE;
  m_errors = 0;
  m_trials = 100;
}

Sherpa::~Sherpa() 
{
  PROFILE_HERE;
  if (p_analysis)     { delete p_analysis;     p_analysis     = NULL; }
  if (p_eventhandler) { delete p_eventhandler; p_eventhandler = NULL; }
  if (p_inithandler)  { delete p_inithandler;  p_inithandler  = NULL; }
}

bool Sherpa::InitializeTheRun(int argc,char * argv[]) 
{ 
  PROFILE_HERE;
  m_path = std::string("./");
  p_inithandler  = new Initialization_Handler(argc, argv);
  int mode = p_inithandler->Mode();  
  if (mode==14) {
    return PerformScan();
  }
  else {
    if (p_inithandler->InitializeTheFramework()) {
      return p_inithandler->CalculateTheHardProcesses();
    }
  }
  msg.Error()<<"Error in Sherpa::InitializeRun("<<m_path<<")"<<endl
	     <<"   Did not manage to initialize the framework."<<endl
	     <<"   Try to run nevertheless ... ."<<endl;
  return 0;
}


bool Sherpa::InitializeTheEventHandler() 
{
  p_output       = p_inithandler->GetOutputHandler();
  int mode       = p_inithandler->Mode();
  p_eventhandler = new Event_Handler();
  if (mode==9999) p_eventhandler->AddEventPhase(new EvtReadin_Phase(p_output));
  else {
      cout<<" ============================================ "<<std::endl;
      p_eventhandler->AddEventPhase(new Signal_Processes(p_inithandler->GetMatrixElementHandler(std::string("SignalMEs")),
							 p_inithandler->GetHardDecayHandler()));
      //p_eventhandler->AddEventPhase(new Hard_Decays(p_inithandler->GetHardDecayHandler()));
      p_eventhandler->AddEventPhase(new Multiple_Interactions(p_inithandler->GetMIHandler()));
      p_eventhandler->AddEventPhase(new Jet_Evolution(p_inithandler->GetMatrixElementHandlers(),
						      p_inithandler->GetShowerHandler()));
      p_eventhandler->AddEventPhase(new Hadronization(p_inithandler->GetBeamRemnantHandler(),
						      p_inithandler->GetFragmentationHandler()));
  }
  AnalysesMap * analyses = p_inithandler->GetSampleAnalyses();
  for (AnalysesIter ana=analyses->begin();ana!=analyses->end();ana++) {
    p_eventhandler->AddEventPhase(new Analysis_Phase(ana->second,ana->first));
  }

  return 1;
}


bool Sherpa::GenerateOneEvent() 
{
  PROFILE_HERE;
  for (int i=0;i<m_trials;i++) {
    if (p_eventhandler->GenerateEvent(p_inithandler->Mode())) {
      if (p_output->OutputOn()) {
	p_output->OutputToFormat(p_eventhandler->GetBlobs());
      }
      return 1;
    }
    m_errors++;
  }
  return 0;
}



bool Sherpa::SummarizeRun() { p_eventhandler->Finish(); }

void Sherpa::DrawLogo() { }

bool Sherpa::PerformScan() 
{
  int np = p_inithandler->NumberOfSteps();
  for (int i=0;i<=np;++i) {
    if (p_inithandler->InitializeTheFramework(i)) {
      if (!p_inithandler->CalculateTheHardProcesses()) {
	msg.Error()<<"ERROR in Sherpa::PerfomScan("<<i<<")"<<endl;
      }
    }
    else {
      msg.Error()<<"ERROR in Sherpa::InitializeRun("<<m_path<<")"<<endl;
    }
  }

  return 1;
}

#include "Sherpa.H"
#include "Signal_Processes.H"
#include "Jet_Evolution.H"
#include "Hadronization.H"

using namespace SHERPA;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace std;

Sherpa::Sherpa()
{
  m_errors = 0;
  m_trials = 100;
  p_inithandler  = NULL;
  p_eventhandler = NULL;
  p_analysis     = NULL;
}

Sherpa::~Sherpa() 
{
  if (p_analysis)     { delete p_analysis;     p_analysis     = NULL; }
  if (p_eventhandler) { delete p_eventhandler; p_eventhandler = NULL; }
  if (p_inithandler)  { delete p_inithandler;  p_inithandler  = NULL; }
}

bool Sherpa::InitializeTheRun(std::string _path) { 
  cout<<"Initialization path = "<<_path<<endl;
  m_path = _path;
  p_inithandler  = new Initialization_Handler(m_path);
  if (p_inithandler->InitializeTheFramework()) {
    return p_inithandler->CalculateTheHardProcesses();
  }
  msg.Error()<<"Error in Sherpa::InitializeRun("<<_path<<")"<<endl
	     <<"   Did not manage to initialize the framework."<<endl
	     <<"   Try to run nevertheless ... ."<<endl;
  return 0;
}


bool Sherpa::InitializeTheEventHandler() {
  p_eventhandler = new Event_Handler();
  p_eventhandler->AddEventPhase(new Signal_Processes(p_inithandler->GetMatrixElementHandler()));
  p_eventhandler->AddEventPhase(new Jet_Evolution(p_inithandler->GetMatrixElementHandler(),
						  p_inithandler->GetShowerHandler()));
  //  p_eventhandler->AddEventPhase(new Hadronization(p_inithandler->GetBeamRemnantHandler(),
  //						  p_inithandler->GetFragmentationHandler()));
  return 1;
}


bool Sherpa::GenerateOneEvent() {
  msg.Out()<<"Starting event generation now. "<<std::endl;
  for (int i=0;i<m_trials;i++) {
    if (p_eventhandler->GenerateEvent()) return 1;
    m_errors++;
  }
  return 0;
}



bool Sherpa::SummarizeRun() {
  // p_analysis->Finish();
}

void Sherpa::DrawLogo() { }

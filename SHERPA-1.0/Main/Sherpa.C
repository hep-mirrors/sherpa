#include "Sherpa.H"
#include "Signal_Processes.H"
#include "Jet_Evolution.H"
#include "Hadronization.H"
#include "Message.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Sherpa::Sherpa() :
  p_inithandler(NULL), p_eventhandler(NULL), p_analysis(NULL)  
{
  m_errors = 0;
  m_trials = 100;
}

Sherpa::~Sherpa() 
{
  if (p_analysis)     { delete p_analysis;     p_analysis     = NULL; }
  if (p_eventhandler) { delete p_eventhandler; p_eventhandler = NULL; }
  if (p_inithandler)  { delete p_inithandler;  p_inithandler  = NULL; }
}

bool Sherpa::InitializeTheRun(std::string _path) { 
  m_path = _path;
  p_inithandler  = new Initialization_Handler(m_path);
  if (p_inithandler->InitializeTheFramework()) {
    p_output     = new Output_Handler(1);
    return p_inithandler->CalculateTheHardProcesses();
  }
  msg.Error()<<"Error in Sherpa::InitializeRun("<<_path<<")"<<endl
	     <<"   Did not manage to initialize the framework."<<endl
	     <<"   Try to run nevertheless ... ."<<endl;
  return 0;
}


bool Sherpa::InitializeTheEventHandler() {
  p_analysis        = new Sample_Analysis();
  p_eventhandler    = new Event_Handler();
  p_eventhandler->AddEventPhase(new Signal_Processes(p_inithandler->GetMatrixElementHandler()));
  p_eventhandler->AddEventPhase(new Analysis_Phase(p_analysis,1));
  p_eventhandler->AddEventPhase(new Jet_Evolution(p_inithandler->GetMatrixElementHandler(),
						  p_inithandler->GetShowerHandler()));
  p_eventhandler->AddEventPhase(new Analysis_Phase(p_analysis,2));
  p_eventhandler->AddEventPhase(new Hadronization(p_inithandler->GetBeamRemnantHandler(),
  						  p_inithandler->GetFragmentationHandler()));
  return 1;
}


bool Sherpa::GenerateOneEvent() {
  for (int i=0;i<m_trials;i++) {
    if (p_eventhandler->GenerateEvent()) {
      if (p_output->Active()) p_output->OutputToFormat(p_eventhandler->GetBlobs());
      return 1;
    }
    m_errors++;
  }
  return 0;
}



bool Sherpa::SummarizeRun() {
  if (p_analysis) p_analysis->Finish();
  return 1;
}

void Sherpa::DrawLogo() { }

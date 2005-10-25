#include "Sherpa.H"
#include "Initialization_Handler.H"
#include "Event_Handler.H"
#include "Analysis_Phase.H"
#include "Analysis_Handler.H"
#include "Input_Output_Handler.H"
#include "EvtReadin_Phase.H"
#include "Signal_Processes.H"
#include "Hard_Decays.H"
#include "Multiple_Interactions.H"
#include "Jet_Evolution.H"
#include "Hadronization.H"
#include "Hadron_Decays.H"
#include "MC_Interface.H"
#include "Message.H"
#include "Scaling.H"

#ifdef PROFILE__all
#define PROFILE__Sherpa
#endif
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
  p_inithandler(NULL), p_eventhandler(NULL), p_iohandler(NULL)
{
  PROFILE_HERE;
  m_errors = 0;
  m_trials = 100;
}

Sherpa::~Sherpa() 
{
  PROFILE_HERE;
  if (p_eventhandler) { delete p_eventhandler; p_eventhandler = NULL; }
  if (p_inithandler)  { delete p_inithandler;  p_inithandler  = NULL; }
}

bool Sherpa::InitializeTheRun(int argc,char * argv[]) 
{ 
  PROFILE_HERE;
  m_path = std::string("./");
  p_inithandler  = new Initialization_Handler(argc, argv);
  DrawLogo();

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
  p_iohandler     = p_inithandler->GetIOHandler();
  int mode        = p_inithandler->Mode();
  p_eventhandler  = new Event_Handler();
  
  std::string sme = std::string("SignalMEs");
  switch (mode) {
  case 9000:
    p_eventhandler->AddEventPhase(new MC_Interface(p_inithandler->GetPythiaInterface())); 
    break;
#ifdef USING__MCatNLO
  case 9001:
    p_eventhandler->AddEventPhase(new MC_Interface(p_inithandler->GetHerwigInterface())); 
    break;
  case 9002:
    p_eventhandler->AddEventPhase(new MC_Interface(p_inithandler->GetMCatNLOInterface())); 
    break;
#endif
  case 9999: 
    p_eventhandler->AddEventPhase(new EvtReadin_Phase(p_inithandler->GetEventReader())); 
    break;
  default:
    p_eventhandler->AddEventPhase(new Signal_Processes(p_inithandler->GetMatrixElementHandler(sme),
													   p_inithandler->GetHardDecayHandler()));
    p_eventhandler->AddEventPhase(new Hard_Decays(p_inithandler->GetHardDecayHandler()));
    p_eventhandler->AddEventPhase(new Jet_Evolution(p_inithandler->GetMatrixElementHandlers(),
													p_inithandler->GetShowerHandler()));
    p_eventhandler->AddEventPhase(new Multiple_Interactions(p_inithandler->GetMIHandler()));
    p_eventhandler->AddEventPhase(new Hadronization(p_inithandler->GetBeamRemnantHandler(),
													p_inithandler->GetFragmentationHandler()));
    p_eventhandler->AddEventPhase(new Hadron_Decays(p_inithandler->GetHadronDecayHandlers()));
    break;
  }
  ANALYSIS::Analysis_Handler * ana = p_inithandler->GetSampleAnalysis();
  if (ana) p_eventhandler->AddEventPhase(new Analysis_Phase(ana));
  p_eventhandler->PrintGenericEventStructure();
  return 1;
}


bool Sherpa::GenerateOneEvent() 
{
  PROFILE_HERE;
  ATOOLS::rpa.gen.SetNumberOfDicedEvents(ATOOLS::rpa.gen.NumberOfDicedEvents()+1);
  for (int i=0;i<m_trials;i++) {
    if (p_eventhandler->GenerateEvent(p_inithandler->Mode())) {
      if (p_iohandler->OutputOn()) {
		p_iohandler->OutputToFormat(p_eventhandler->GetBlobs());
      }
      return 1;
    }
    m_errors++;
  }
  return 0;
}



bool Sherpa::SummarizeRun() 
{ 
  p_eventhandler->Finish(); 
  return true; 
}

void Sherpa::DrawLogo() 
{ 
  msg_Info()<<"-----------------------------------------------------------------------------"<<std::endl;
  msg.Out()<<"-----------    Event generation run with SHERPA started .......   -----------"<<std::endl;
  msg_Info()<<"-----------------------------------------------------------------------------"<<std::endl
	    <<"................................................ |       +                   "<<std::endl
	    <<"................................................ ||  |       +  +            "<<std::endl
	    <<"...................................        ....  | |         /   +           "<<std::endl
	    <<"................. ................   _,_ |  ....  ||         +|  +  +        "<<std::endl
	    <<"...............................  __.'  ,\\|  ...  ||    /    +|          +    "<<std::endl
	    <<".............................. (´  \\  ´  \\   ...  | |  |   + + \\         +   "<<std::endl
	    <<".............................  (    \\   -/  .... ||       +    |          +  "<<std::endl
	    <<"........ ...................  <S   /()))))~~~~~~~~##     +     /\\    +       "<<std::endl
	    <<"............................ (!H   (~~)))))~~~~~~#/     +  +    |  +         "<<std::endl
	    <<"................ ........... (!E   (~~~)))))     /|/    +         +          "<<std::endl
	    <<"............................ (!R   (~~~)))))   |||   + +            +        "<<std::endl
	    <<"..... ...................... (!P    (~~~~)))   /|  + +          +            "<<std::endl
	    <<"............................ (!A>    (~~~~~~~~~##        + +        +        "<<std::endl
	    <<"............................. ~~(!    '~~~~~~~ \\       +     + +      +      "<<std::endl
	    <<"...............................  `~~~QQQQQDb //   |         + + +        +   "<<std::endl
	    <<"........................ ..........   IDDDDP||     \\  + + + + +             +"<<std::endl
	    <<"....................................  IDDDI||       \\                      + "<<std::endl
	    <<".................................... IHD HD||         \\ + +  + + + + +      +"<<std::endl
	    <<"...................................  IHD ##|            :-) + +\\          +  "<<std::endl
	    <<"......... ............... ......... IHI ## /      /   +  + + + +\\       +    "<<std::endl
	    <<"................................... IHI/ /       /      + + + +        +     "<<std::endl
	    <<"................................... ## | | /    / + +      + + /      +      "<<std::endl
	    <<"....................... /TT\\ .....  ##/ ///  / + + + + + + +/       +        "<<std::endl
	    <<"......................./TTT/T\\ ... /TT\\/\\\\\\ / + + + + + + +/   \\         +   "<<std::endl
	    <<"version 1.0.5 ......../TTT/TTTT\\...|TT/T\\\\\\/   +    ++  + /                  "<<std::endl
	    <<"-----------------------------------------------------------------------------"<<std::endl
	    <<std::endl
	    <<"          SHERPA version 1.0.5.                                              "<<std::endl
	    <<"                                                                             "<<std::endl
	    <<"          AUTHORS: Tanju Gleisberg, Stefan Hoeche, Frank Krauss,             "<<std::endl
	    <<"               Andreas Schaelicke, Steffen Schumann, Jan Winter              "<<std::endl
	    <<"                                                                             "<<std::endl
	    <<"                                                                             "<<std::endl
	    <<"          This program uses a lot of genuine and original research           "<<std::endl
	    <<"          work by other people. Users are encouraged to refer to             "<<std::endl
	    <<"          the various original publications.                                 "<<std::endl
	    <<"                                                                             "<<std::endl
	    <<"          Users are kindly asked to refer to the documentation               "<<std::endl
	    <<"          published under JHEP 0402 (2004) 056.                              "<<std::endl
	    <<"                                                                             "<<std::endl
	    <<"          Please visit also our homepage                                     "<<std::endl
	    <<"          http://www.physik.tu-dresden.de/~krauss/hep/index.html             "<<std::endl
	    <<"          for news, bugreports, updates and new releases.                    "<<std::endl
	    <<"                                                                             "<<std::endl
	    <<"-----------------------------------------------------------------------------"<<std::endl
	    <<std::endl;
}

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

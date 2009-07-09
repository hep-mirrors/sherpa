#include "SHERPA/Main/Sherpa.H"
#include "SHERPA/Initialization/Initialization_Handler.H"
#include "SHERPA/Single_Events/Event_Handler.H"
#include "SHERPA/Single_Events/Analysis_Phase.H"
#include "SHERPA/Tools/Input_Output_Handler.H"
#include "SHERPA/Single_Events/EvtReadin_Phase.H"
#include "SHERPA/Single_Events/Signal_Processes.H"
#include "SHERPA/Single_Events/Hard_Decays.H"
#include "SHERPA/Single_Events/Multiple_Interactions.H"
#include "SHERPA/Single_Events/Jet_Evolution.H"
#include "SHERPA/Single_Events/Beam_Remnants.H"
#include "SHERPA/Single_Events/Hadronization.H"
#include "SHERPA/Single_Events/Hadron_Decays.H"
#ifdef USING__PYTHIA
#include "SHERPA/Single_Events/MC_Interface.H"
#endif
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Data_Writer.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include <cstring>

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Sherpa::Sherpa() :
  p_inithandler(NULL), p_eventhandler(NULL), p_iohandler(NULL)
{
  ATOOLS::s_loader = new Library_Loader();
  m_errors = 0;
  m_trials = 100;
  m_debuginterval = 0;
  m_debugstep = 0;
  exh->AddTerminatorObject(this);
}

Sherpa::~Sherpa() 
{
  if (p_eventhandler) { delete p_eventhandler; p_eventhandler = NULL; }
  if (p_inithandler)  { delete p_inithandler;  p_inithandler  = NULL; }
  exh->RemoveTerminatorObject(this);
  delete ATOOLS::s_loader;
}

bool Sherpa::InitializeTheRun(int argc,char * argv[]) 
{ 
  m_path = std::string("./");
  int oldc(argc);
  char **oldargs(NULL);
  std::string statuspath;
  for (int i(1);i<argc;++i) {
    std::string cur(argv[i]);
    size_t pos(cur.find("STATUS_PATH"));
    if (pos==0 && cur.length()>11 && cur[11]=='=') {
      statuspath=cur.substr(12);
      if (statuspath=="") continue;
      if (statuspath[statuspath.length()-1]!='/') statuspath+=std::string("/");
      Data_Reader reader;
      reader.SetInputFile(statuspath+"cmd");
      String_Matrix args;
      reader.MatrixFromFile(args);
      oldc=argc;
      oldargs=argv;
      argc+=args.size();
      argv = new char*[argc];
      argv[0] = new char[strlen(oldargs[0])+1];
      strcpy(argv[0],oldargs[0]);
      for (int j(0);j<(int)args.size();++j) {
	std::string cur(args[j].front());
	for (size_t k(1);k<args[j].size();++k) cur+=args[j][k];
	argv[j+1] = new char[cur.length()+1];
	strcpy(argv[j+1],cur.c_str());
      }
      for (int j(1);j<oldc;++j) {
	argv[args.size()+j] = new char[strlen(oldargs[j])+1];
	strcpy(argv[args.size()+j],oldargs[j]);
      }
      break;
    }
  }

  p_inithandler  = new Initialization_Handler(argc, argv);
  DrawLogo();

  if (p_inithandler->InitializeTheFramework()) {
    if (!p_inithandler->CalculateTheHardProcesses()) return false;
    if (statuspath!="") {
      bool res(exh->ReadInStatus(statuspath));
      if (oldargs) {
        for (int i(0);i<argc;++i) delete [] argv[i];
        delete [] argv;
      }
      return res;
    }

    Data_Reader read(" ",";","!","=");
    long int debuginterval(0);
    if (read.ReadFromFile(debuginterval,"DEBUG_INTERVAL")) {
      m_debuginterval=debuginterval;
      msg_Info()<<"Setting debug interval to "<<m_debuginterval<<std::endl;
    }
    long int debugstep(0);
    if (read.ReadFromFile(debugstep,"DEBUG_STEP")) {
      m_debugstep=debugstep;
      ran.ReadInStatus(("random."+ToString(m_debugstep)+".dat").c_str());
    }
    return true;
  }
  msg_Error()<<"Error in Sherpa::InitializeRun("<<m_path<<")"<<endl
	     <<"   Did not manage to initialize the framework."<<endl
	     <<"   Try to run nevertheless ... ."<<endl;
  
  return 0;
}


bool Sherpa::InitializeTheEventHandler() 
{
  p_iohandler     = p_inithandler->GetIOHandler();
  int mode        = p_inithandler->Mode();
  p_eventhandler  = new Event_Handler();
  p_iohandler->SetEventHandler(p_eventhandler);
  Analysis_Interface *ana(p_inithandler->GetAnalysis());
  if (ana) p_inithandler->GetAnalysis()->SetEventHandler(p_eventhandler);
  
  std::string sme = std::string("SignalMEs");
  switch (mode) {
#ifdef USING__PYTHIA
  case 9000:
    p_eventhandler->AddEventPhase(new MC_Interface(p_inithandler->GetPythiaInterface())); 
    break;
#endif
  case 9999: 
    p_eventhandler->AddEventPhase(new EvtReadin_Phase(p_inithandler->GetEventReader())); 
    break;
  default:
    p_eventhandler->AddEventPhase(new Signal_Processes(p_inithandler->GetMatrixElementHandler(sme)));
    p_eventhandler->AddEventPhase(new Hard_Decays(p_inithandler->GetHardDecayHandler()));
    p_eventhandler->AddEventPhase(new Jet_Evolution(p_inithandler->GetMatrixElementHandlers(),
						    p_inithandler->GetHadronDecayHandlers(),
						    p_inithandler->GetMIHandler(),
						    p_inithandler->GetShowerHandler()));
    p_eventhandler->AddEventPhase(new Multiple_Interactions(p_inithandler->GetMIHandler()));
    p_eventhandler->AddEventPhase(new Beam_Remnants(p_inithandler->GetBeamRemnantHandler()));
    p_eventhandler->AddEventPhase(new Hadronization(p_inithandler->GetFragmentationHandler()));
    p_eventhandler->AddEventPhase(new Hadron_Decays(p_inithandler->GetHadronDecayHandlers(),
                                                    p_inithandler->GetSoftPhotonHandler()));

    break;
  }
  if (ana) p_eventhandler->AddEventPhase(new Analysis_Phase(ana));
  p_eventhandler->PrintGenericEventStructure();
  return 1;
}


bool Sherpa::GenerateOneEvent(bool reset) 
{
  ATOOLS::rpa.gen.SetNumberOfDicedEvents(ATOOLS::rpa.gen.NumberOfDicedEvents()+1);
  for (int i=0;i<m_trials;i++) {
    if(m_debuginterval>0 && rpa.gen.NumberOfDicedEvents()%m_debuginterval==0){
      std::string fname=ToString(rpa.gen.NumberOfDicedEvents())+".dat";
      ran.WriteOutStatus(("random."+fname).c_str());
    }
    if (m_debugstep>0) {
      ran.ReadInStatus(("random."+ToString(m_debugstep)+".dat").c_str());
    }
    if (reset) p_eventhandler->Reset();
    if (p_eventhandler->GenerateEvent(p_inithandler->Mode())) {
      if(m_debuginterval>0 && rpa.gen.NumberOfDicedEvents()%m_debuginterval==0){
        std::string fname=ToString(rpa.gen.NumberOfDicedEvents())+".dat";
        std::ofstream eventout(("refevent."+fname).c_str());
        eventout<<*p_eventhandler->GetBlobs()<<std::endl;
        eventout.close();
      }
      if (m_debugstep>0) {
        std::ofstream event(("event."+ToString(m_debugstep)+".dat").c_str());
        event<<*p_eventhandler->GetBlobs()<<std::endl;
        event.close();
        THROW(normal_exit,"Debug event written.");
      }
      p_iohandler->OutputToFormat(p_eventhandler->GetBlobs());
      p_iohandler->PrintEvent(p_eventhandler->GetBlobs());
      return 1;
    }
    m_errors++;
  }
  return 0;
}

void Sherpa::PrepareTerminate()
{
  SummarizeRun();
  exh->RemoveTerminatorObject(this);
}

void Sherpa::WriteCitationInfo()
{
  std::string refname("Sherpa_References.tex");
  Data_Writer writer;
  writer.SetComment("%%");
  writer.SetOutputPath(rpa.gen.Variable("SHERPA_RUN_PATH")+"/");
  writer.SetOutputFile(refname);
  writer.SetVectorType(vtc::vertical);
  writer.WriteComment("Citation summary file generated by Sherpa");
  writer.WriteComment("PID "+ToString(getpid())+" on "+rpa.gen.Timer().TimeString(0));
  writer.WriteToFile(std::string("\n\\documentclass{article}\n\n\\begin{document}\n"));
  writer.VectorToFile(rpa.gen.Citations());
  writer.WriteToFile(std::string("\n\\end{document}"));
  msg_Out()<<METHOD<<"(): {\n";
  msg_Out()<<"  "<<om::bold<<"Please cite the publications listed in '"
	   <<om::red<<refname<<om::reset<<om::bold<<"'."<<om::reset
	   <<"\n  Extract the bibtex list by running 'bin/SHERPA-MC/get_bibtex "
	   <<refname<<"'\n  or email the file to 'slaclib2@slac.stanford.edu'"
	   <<", subject 'generate'.\n";
  msg_Out()<<"}"<<std::endl;
}

bool Sherpa::SummarizeRun() 
{ 
  p_eventhandler->Finish(); 
  WriteCitationInfo();
  return true; 
}

void Sherpa::DrawLogo() 
{ 
  msg_Info()<<"-----------------------------------------------------------------------------"<<std::endl;
  if (msg->Level()>0) msg_Out()<<"-----------    Event generation run with SHERPA started .......   -----------"<<std::endl;
  msg_Info()<<"-----------------------------------------------------------------------------"<<std::endl
	    <<"................................................ |       +                   "<<std::endl
	    <<"................................................ ||  |       +  +            "<<std::endl
	    <<"...................................        ....  | |         /   +           "<<std::endl
	    <<"................. ................   _,_ |  ....  ||         +|  +  +        "<<std::endl
	    <<"...............................  __.'  ,\\|  ...  ||    /    +|          +    "<<std::endl
	    <<".............................. (  \\    \\   ...  | |  |   + + \\         +   "<<std::endl
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
	    <<"version "<<SHERPA_VERSION<<"."<<SHERPA_SUBVERSION<<" ......../TTT/TTTT\\...|TT/T\\\\\\/   +    ++  + /                  "<<std::endl
	    <<"-----------------------------------------------------------------------------"<<std::endl
	    <<std::endl
	    <<"     SHERPA version "<<SHERPA_VERSION<<"."<<SHERPA_SUBVERSION
	    <<"                                              "<<std::endl
	    <<"                                                                             "<<std::endl
	    <<"     Authors:        Tanju Gleisberg, Stefan Hoeche, Frank Krauss,           "<<std::endl
	    <<"                     Marek Schoenherr, Steffen Schumann, Frank Siegert,      "<<std::endl
            <<"                     Jan Winter."<<std::endl
	    <<"     Former Authors: Timo Fischer, Ralf Kuhn, Thomas Laubrich,               "<<std::endl
	    <<"                     Andreas Schaelicke                                      "<<std::endl
	    <<"                                                                             "<<std::endl
	    <<"     This program uses a lot of genuine and original research work           "<<std::endl
	    <<"     by other people. Users are encouraged to refer to                       "<<std::endl
	    <<"     the various original publications.                                      "<<std::endl
	    <<"                                                                             "<<std::endl
	    <<"     Users are kindly asked to refer to the documentation                    "<<std::endl
	    <<"     published under JHEP 02(2009)007                                        "<<std::endl
	    <<"                                                                             "<<std::endl
	    <<"     Please visit also our homepage                                          "<<std::endl
	    <<"                                                                             "<<std::endl
	    <<"       http://www.sherpa-mc.de                                               "<<std::endl
	    <<"                                                                             "<<std::endl
	    <<"     for news, bugreports, updates and new releases.                         "<<std::endl
	    <<"                                                                             "<<std::endl
	    <<"-----------------------------------------------------------------------------"<<std::endl
	    <<std::endl;
  rpa.gen.AddCitation
    (0,"The complete Sherpa package is published under \\cite{Gleisberg:2008ta}.");
}

#include "SHERPA/Main/Sherpa.H"
#include "SHERPA/Initialization/Initialization_Handler.H"
#include "SHERPA/Single_Events/Event_Handler.H"
#include "SHERPA/Single_Events/Analysis_Phase.H"
#include "SHERPA/Single_Events/Output_Phase.H"
#include "SHERPA/Single_Events/EvtReadin_Phase.H"
#include "SHERPA/Single_Events/Signal_Processes.H"
#include "SHERPA/Single_Events/Hard_Decays.H"
#include "SHERPA/Single_Events/Minimum_Bias.H"
#include "SHERPA/Single_Events/Multiple_Interactions.H"
#include "SHERPA/Single_Events/Jet_Evolution.H"
#include "SHERPA/Single_Events/Signal_Process_FS_QED_Correction.H"
#include "SHERPA/Single_Events/Beam_Remnants.H"
#include "SHERPA/Single_Events/Hadronization.H"
#include "SHERPA/Single_Events/Hadron_Decays.H"
#include "SHERPA/PerturbativePhysics/Hard_Decay_Handler.H"
#include "SHERPA/Tools/HepMC2_Interface.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include <cstring>

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Sherpa::Sherpa() :
  p_inithandler(NULL), p_eventhandler(NULL), p_hepmc2(NULL)
{
  ATOOLS::mpi = new My_MPI();
  ATOOLS::exh = new Exception_Handler();
  ATOOLS::msg = new Message();
  ATOOLS::ran = new Random(1234);
  ATOOLS::rpa = new Run_Parameter();
  ATOOLS::s_loader = new Library_Loader();
  m_trials = 100;
  m_debuginterval = 0;
  m_debugstep = -1;
  m_displayinterval = 100;
  m_evt_starttime = -1.0;
  exh->AddTerminatorObject(this);
}

Sherpa::~Sherpa() 
{
  if (msg_LevelIsInfo()) Return_Value::PrintStatistics(msg->Out());
  rpa->gen.WriteCitationInfo();
  if (p_eventhandler) { delete p_eventhandler; p_eventhandler = NULL; }
  if (p_inithandler)  { delete p_inithandler;  p_inithandler  = NULL; }
  exh->RemoveTerminatorObject(this);
  delete ATOOLS::s_loader;
  delete ATOOLS::rpa;
  delete ATOOLS::ran;
#ifdef USING__MPI
  int dummy=0;
  MPI::COMM_WORLD.Bcast(&dummy,1,MPI::INT,0);
#endif  
  delete ATOOLS::msg;
  delete ATOOLS::exh;
  delete ATOOLS::mpi;
  for (KF_Table::const_iterator kfit(s_kftable.begin());kfit!=s_kftable.end();++kfit)
    delete kfit->second;
  ATOOLS::s_kftable.clear();
}

bool Sherpa::InitializeTheRun(int argc,char * argv[]) 
{ 
  m_path = std::string("");
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

  mpi->SetUpSendRecv(p_inithandler->DataReader());

  DrawLogo(p_inithandler->DataReader()->GetValue("PRINT_VERSION_INFO",0));

  if (p_inithandler->InitializeTheFramework()) {
    Data_Reader read(" ",";","!","=");
    int initonly=read.GetValue<int>("INIT_ONLY",0);
    if (initonly==1) THROW(normal_exit,"Initialization complete.");
    if (initonly==2) return true;
    if (!p_inithandler->CalculateTheHardProcesses()) return false;
    bool res(true);
    if (statuspath!="") {
      res=exh->ReadInStatus(statuspath);
      if (oldargs) {
        for (int i(0);i<argc;++i) delete [] argv[i];
        delete [] argv;
      }
    }

    long int debuginterval(0);
    if (read.ReadFromFile(debuginterval,"DEBUG_INTERVAL")) {
      m_debuginterval=debuginterval;
      msg_Info()<<"Setting debug interval to "<<m_debuginterval<<std::endl;
    }
    long int debugstep(-1);
    if (read.ReadFromFile(debugstep,"DEBUG_STEP")) {
      m_debugstep=debugstep;
      ran->ReadInStatus(("random."+ToString(m_debugstep)+".dat").c_str());
    }

    m_displayinterval=read.GetValue<int>("EVENT_DISPLAY_INTERVAL",100);
    
    return res;
  }
  msg_Error()<<"Error in Sherpa::InitializeRun("<<m_path<<")"<<endl
	     <<"   Did not manage to initialize the framework."<<endl
	     <<"   Try to run nevertheless ... ."<<endl;
  
  return 0;
}


bool Sherpa::InitializeTheEventHandler() 
{
  eventtype::code mode = p_inithandler->Mode();
  p_eventhandler  = new Event_Handler();
  Output_Vector *outs(p_inithandler->GetOutputs());
  Analysis_Vector *anas(p_inithandler->GetAnalyses());
  for (Analysis_Vector::iterator it=anas->begin(); it!=anas->end(); ++it) {
    (*it)->SetEventHandler(p_eventhandler);
  }
  
  if (mode==eventtype::EventReader) {
    p_eventhandler->AddEventPhase(new EvtReadin_Phase(p_inithandler->GetEventReader())); 
    p_eventhandler->AddEventPhase(new Beam_Remnants(p_inithandler->GetBeamRemnantHandler()));
  }
  else {
    p_eventhandler->AddEventPhase(new Signal_Processes(p_inithandler->GetMatrixElementHandler(),
                                                       p_inithandler->GetVariations()));
    p_eventhandler->AddEventPhase(new Hard_Decays(p_inithandler->GetHardDecayHandler()));
    p_eventhandler->AddEventPhase(new Jet_Evolution(p_inithandler->GetMatrixElementHandler(),
                                                    p_inithandler->GetHardDecayHandler(),
						    p_inithandler->GetHDHandler(),
						    p_inithandler->GetMIHandler(),
						    p_inithandler->GetSoftCollisionHandler(),
						    p_inithandler->GetShowerHandlers()));
    p_eventhandler->AddEventPhase(new Signal_Process_FS_QED_Correction(p_inithandler->GetMatrixElementHandler(),
                                                                       p_inithandler->GetSoftPhotonHandler()));
    p_eventhandler->AddEventPhase(new Multiple_Interactions(p_inithandler->GetMIHandler()));
    p_eventhandler->AddEventPhase(new Minimum_Bias(p_inithandler->GetSoftCollisionHandler()));
    p_eventhandler->AddEventPhase(new Beam_Remnants(p_inithandler->GetBeamRemnantHandler()));
    p_eventhandler->AddEventPhase(new Hadronization(p_inithandler->GetFragmentationHandler()));
    p_eventhandler->AddEventPhase(new Hadron_Decays(p_inithandler->GetHDHandler()));

  }
  if (!anas->empty()) p_eventhandler->AddEventPhase(new Analysis_Phase(anas));
  if (!outs->empty()) p_eventhandler->AddEventPhase(new Output_Phase(outs,p_eventhandler));
  p_eventhandler->PrintGenericEventStructure();

  msg->SetLevel(p_inithandler->DataReader()->GetValue<int>("EVT_OUTPUT",msg->Level()));

  return 1;
}


bool Sherpa::GenerateOneEvent(bool reset) 
{
    if(m_debuginterval>0 && rpa->gen.NumberOfGeneratedEvents()%m_debuginterval==0){
      std::string fname=ToString(rpa->gen.NumberOfGeneratedEvents())+".dat";
      ran->WriteOutStatus(("random."+fname).c_str());
    }
    if (m_debugstep>=0) {
      ran->ReadInStatus(("random."+ToString(m_debugstep)+".dat").c_str());
    }

    if (m_evt_starttime<0.0) m_evt_starttime=rpa->gen.Timer().RealTime();
    
    if (reset) p_eventhandler->Reset();
    if (p_eventhandler->GenerateEvent(p_inithandler->Mode())) {
      if(m_debuginterval>0 && rpa->gen.NumberOfGeneratedEvents()%m_debuginterval==0){
        std::string fname=ToString(rpa->gen.NumberOfGeneratedEvents())+".dat";
        std::ofstream eventout(("refevent."+fname).c_str());
        eventout<<*p_eventhandler->GetBlobs()<<std::endl;
        eventout.close();
      }
      if (m_debugstep>=0) {
        std::ofstream event(("event."+ToString(m_debugstep)+".dat").c_str());
        event<<*p_eventhandler->GetBlobs()<<std::endl;
        event.close();
        THROW(normal_exit,"Debug event written.");
      }
      rpa->gen.SetNumberOfGeneratedEvents(rpa->gen.NumberOfGeneratedEvents()+1);
      Blob_List *blobs(p_eventhandler->GetBlobs());
      if (msg_LevelIsEvents()) {
	if (!blobs->empty()) {
	  msg_Out()<<"  -------------------------------------------------  "<<std::endl;
	  for (Blob_List::iterator blit=blobs->begin();blit!=blobs->end();++blit) 
	    msg_Out()<<*(*blit)<<std::endl;
	  msg_Out()<<"  -------------------------------------------------  "<<std::endl;
	}
	else msg_Out()<<"  ******** Empty event ********  "<<std::endl;
      }

      for (Blob_List::const_iterator bit=blobs->begin(); bit!=blobs->end();++bit) {
          double currQ = (*bit)->CheckChargeConservation();
          if (fabs(currQ)>1e-12) {
              msg_Error() << "Charge conservation failed, aborting: " << currQ << "\n";
              msg_Error() << (**bit) << "\n";
              abort();
          }
      }

      int i=rpa->gen.NumberOfGeneratedEvents();
      int nevt=rpa->gen.NumberOfEvents();
      msg_Events()<<"Sherpa : Passed "<<i<<" events."<<std::endl;
      int exp;
      for (exp=5; i/int(pow(10,exp))==0; --exp) {}
      if (((rpa->gen.BatchMode()&4 && i%m_displayinterval==0) ||
           (!(rpa->gen.BatchMode()&4) && i%int(pow(10,exp))==0)) &&
          i<rpa->gen.NumberOfEvents()) {
        double diff=rpa->gen.Timer().RealTime()-m_evt_starttime;
        msg_Info()<<"  Event "<<i<<" ( "
                  <<FormatTime(size_t(diff))<<" elapsed / "
                  <<FormatTime(size_t((nevt-i)/(double)i*diff))
                  <<" left ) -> ETA: "<<rpa->gen.Timer().
          StrFTime("%a %b %d %H:%M",time_t((nevt-i)/(double)i*diff))<<"  ";
        double xs(GetEventHandler()->TotalXSMPI());
        double err(GetEventHandler()->TotalErrMPI());
        if (!(rpa->gen.BatchMode()&2)) msg_Info()<<"\n  ";
        msg_Info()<<"XS = "<<xs<<" pb +- ( "<<err<<" pb = "
                  <<((int(err/xs*10000))/100.0)<<" % )  ";
        if (!(rpa->gen.BatchMode()&2))
          msg_Info()<<mm(1,mm::up);
        if (rpa->gen.BatchMode()&2) { msg_Info()<<std::endl; }
        else { msg_Info()<<bm::cr<<std::flush; }
      }
      
      return 1;
    }
    return 0;
}

void Sherpa::FillHepMCEvent(HepMC::GenEvent& event)
{
#ifdef USING__HEPMC2
  if (p_hepmc2==NULL) p_hepmc2 = new SHERPA::HepMC2_Interface();
  ATOOLS::Blob_List* blobs=GetEventHandler()->GetBlobs();
  p_hepmc2->Sherpa2HepMC(blobs, event, blobs->Weight());
#else
  THROW(fatal_error, "HepMC not linked.");
#endif
}

double Sherpa::TotalXS()
{
  return p_eventhandler->TotalXS();
}

double Sherpa::TotalErr()
{
  return p_eventhandler->TotalErr();
}

std::string Sherpa::PDFInfo()
{
  std::string pdf="Unknown";
  PDF::ISR_Handler* isr=GetInitHandler()->GetISRHandler(PDF::isr::hard_process);
  if (isr) {
    if (isr->PDF(0)) {
      pdf=isr->PDF(0)->Type();
      if (isr->PDF(1) && isr->PDF(1)->Type()!=pdf) {
        pdf="Unknown";
      }
    }
  }
  return pdf;
}

void Sherpa::PrepareTerminate()
{
  SummarizeRun();
  exh->RemoveTerminatorObject(this);
}

bool Sherpa::SummarizeRun() 
{
  if (p_eventhandler) {
    msg_Info()<<"  Event "<<rpa->gen.NumberOfGeneratedEvents()<<" ( "
              <<size_t(rpa->gen.Timer().RealTime()-m_evt_starttime)
              <<" s total ) = "
              << rpa->gen.NumberOfGeneratedEvents()*3600*24/
                 ((size_t) rpa->gen.Timer().RealTime()-m_evt_starttime)
              <<" evts/day                    "<<std::endl;
  }
  p_eventhandler->Finish(); 
  return true; 
}

long int Sherpa::NumberOfEvents() const
{
  return rpa->gen.NumberOfEvents();
}

const Blob_List &Sherpa::GetBlobList() const
{
  return *p_eventhandler->GetBlobs();
}

double Sherpa::GetMEWeight(const Cluster_Amplitude &ampl,const int mode) const
{
  return p_inithandler->GetMatrixElementHandler()->
    GetWeight(ampl,PHASIC::nlo_type::lo,mode);
}

void Sherpa::DrawLogo(const int mode) 
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
	    <<"....................../TTT/TTTT\\...|TT/T\\\\\\/   +    ++  + /              "<<std::endl
	    <<"-----------------------------------------------------------------------------"<<std::endl
	    <<std::endl
	    <<"     SHERPA version "<<SHERPA_VERSION<<"."<<SHERPA_SUBVERSION<<" ("<<SHERPA_NAME<<")"<<std::endl
	    <<"                                                                             "<<std::endl
	    <<"     Authors:        Stefan Hoeche, Frank Krauss, Silvan Kuttimalai,         "<<std::endl
	    <<"                     Marek Schoenherr, Holger Schulz, Steffen Schumann,      "<<std::endl
	    <<"                     Frank Siegert, Korinna Zapp."<<std::endl
	    <<"     Former Authors: Timo Fischer, Tanju Gleisberg, Hendrik Hoeth,           "<<std::endl
	    <<"                     Ralf Kuhn, Thomas Laubrich, Andreas Schaelicke,         "<<std::endl
	    <<"                     Jan Winter                                              "<<std::endl
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
	    <<"       http://sherpa.hepforge.org                                            "<<std::endl
	    <<"                                                                             "<<std::endl
	    <<"     for news, bugreports, updates and new releases.                         "<<std::endl
	    <<"                                                                             "<<std::endl
	    <<"-----------------------------------------------------------------------------"<<std::endl
	    <<std::endl;
  rpa->gen.PrintSVNVersion(msg->Info(),mode);
  rpa->gen.AddCitation
    (0,"The complete Sherpa package is published under \\cite{Gleisberg:2008ta}.");
}

#include "Sherpa.H"
#include "Message.H"
#include "prof.hh"
#include "Random.H"
#include "Exception.H"
#include "Run_Parameter.H"

#ifdef PROFILE__all
#include "prof.hh"
#endif
#ifdef ROOT_SUPPORT
#include "My_Root.H"
#endif

using namespace SHERPA;

int main(int argc,char* argv[]) 
{  
#ifdef ROOT_SUPPORT
  MYROOT::myroot = new MYROOT::My_Root(argc,argv);
  ATOOLS::Exception_Handler::AddTerminatorObject(MYROOT::myroot);
#endif
  std::set_terminate(ATOOLS::Exception_Handler::Terminate);
  std::set_unexpected(ATOOLS::Exception_Handler::Terminate);
  signal(SIGSEGV,ATOOLS::Exception_Handler::SignalHandler);
  signal(SIGINT,ATOOLS::Exception_Handler::SignalHandler);
  signal(SIGBUS,ATOOLS::Exception_Handler::SignalHandler);
  signal(SIGFPE,ATOOLS::Exception_Handler::SignalHandler);
  signal(SIGABRT,ATOOLS::Exception_Handler::SignalHandler);
  signal(SIGTERM,ATOOLS::Exception_Handler::SignalHandler);
  signal(SIGXCPU,ATOOLS::Exception_Handler::SignalHandler);
  try {
#ifdef PROFILE__all    
    set_prof();
#endif
    Sherpa Generator;
    Generator.InitializeTheRun(argc,argv);
    int nevt=ATOOLS::rpa.gen.NumberOfEvents();
    if (nevt>0) {
      ATOOLS::msg.Out()<<"=========================================================================="<<std::endl
		       <<"Sherpa will start event generation now : "
		       <<nevt<<" events"<<std::endl
		       <<"=========================================================================="<<std::endl;
      Generator.InitializeTheEventHandler();
      for (int i=1;i<=nevt;i++) {
	if (i%100==0) {
	  msg_Info()<<"  Event "<<i<<" ( "
		    <<ATOOLS::rpa.gen.Timer().UserTime()
		    <<" s )  "<<ATOOLS::bm::cr<<std::flush; 
	}
	if (Generator.GenerateOneEvent()) msg_Events()<<"Sherpa : Passed "<<i<<" events."<<std::endl;
      }
      msg_Info()<<std::endl;      
      Generator.SummarizeRun();
    }
    ATOOLS::msg.Out()<<"=========================================================================="<<std::endl
		     <<"Sherpa finished its simulation run with "
		     <<Generator.NumberOfErrors()<<" errors."<<std::endl
		     <<"=========================================================================="<<std::endl;
#ifdef ROOT_SUPPORT
    delete MYROOT::myroot;
#endif
#ifdef PROFILE__all    
    std::ofstream *output = new std::ofstream("profile.out",std::ios::out);
    print_profile(*output);
    delete output;
#endif
    return 0;
  }
  catch (ATOOLS::Exception exception) {
    exception.UpdateLogFile();
    ATOOLS::msg.Error()<<exception<<std::endl;
    std::terminate();
  }
  catch (std::exception exception) {
    std::cout<<"Sherpa: throws std::exception "<<exception.what()<<" ..."<<std::endl;
    std::terminate();
  }
}




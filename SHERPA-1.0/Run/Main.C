#include "Sherpa.H"
#include "Message.H"
#include "prof.hh"
#include "Random.H"
#include "Exception.H"
#include "Run_Parameter.H"
#ifdef USING__CLHEP
#include "Input_Output_Handler.H"
#endif

#ifdef PROFILE__all
#include "prof.hh"
#endif
#ifdef USING__ROOT
#include "My_Root.H"
#endif

#ifdef TRACE_malloc
#include <mcheck.h>
#endif

using namespace SHERPA;

int main(int argc,char* argv[]) 
{
#ifdef TRACE_malloc
  setenv("MALLOC_TRACE","malloc_trace.log",1);
  mtrace();  
#endif
#ifdef USING__ROOT
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
      double starttime=ATOOLS::rpa.gen.Timer().UserTime();
      for (int i=1;i<=ATOOLS::rpa.gen.NumberOfEvents();i++) {
	//if (i==259) ATOOLS::msg.SetLevel(15); 
	if (i%100==0) {
	  double diff=ATOOLS::rpa.gen.Timer().UserTime()-starttime;
	  msg_Info()<<"  Event "<<i<<" ( "
		    <<int(diff)<<" s elapsed / "
		    <<int((nevt-i)/(double)i*diff)
		    <<" s left / "<<int(nevt/(double)i*diff)
		    <<" s total )   "<<ATOOLS::bm::cr<<std::flush;
	}
	if (Generator.GenerateOneEvent()) msg_Events()<<"Sherpa : Passed "<<i<<" events."<<std::endl;
#ifdef USING__CLHEP
// 	Generator.GetIOHandler()->GetHepMCInterface()->PrintHepMCEvent();
#endif
      }
      msg_Info()<<std::endl;      
      Generator.SummarizeRun();
    }
    ATOOLS::msg.Out()<<"=========================================================================="<<std::endl
		     <<"Sherpa finished its simulation run with "
		     <<Generator.NumberOfErrors()<<" errors."<<std::endl
		     <<"=========================================================================="<<std::endl;
#ifdef USING__ROOT
    delete MYROOT::myroot;
#endif
#ifdef PROFILE__all    
    std::ofstream *output = new std::ofstream("profile.out",std::ios::out);
    print_profile(*output);
    delete output;
#endif
#ifdef TRACE_malloc
    muntrace();  
#endif
    return 0;
  }
  catch (ATOOLS::Exception exception) {
    exception.UpdateLogFile();
    ATOOLS::msg.Error()<<exception<<std::endl;
#ifdef TRACE_malloc
    muntrace();  
#endif
    std::terminate();
  }
  catch (std::exception exception) {
    std::cout<<"Sherpa: throws std::exception "<<exception.what()<<" ..."<<std::endl;
#ifdef TRACE_malloc
    muntrace();  
#endif
    std::terminate();
  }
}




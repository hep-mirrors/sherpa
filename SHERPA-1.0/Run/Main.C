#include "Sherpa.H"
#include "Message.H"
#include "prof.hh"
#include "Random.H"
#include "Exception.H"

using namespace SHERPA;

int main(int argc,char* argv[]) 
{  
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
    Sherpa Generator;
    Generator.InitializeTheRun(argc,argv);
    int nevt=ATOOLS::rpa.gen.NumberOfEvents();
    if (nevt>0) {
      ATOOLS::msg.Out()<<"=========================================================================="<<std::endl
		       <<"Sherpa will start event generation now : "
		       <<ATOOLS::rpa.gen.NumberOfEvents()<<" events"<<std::endl
		       <<"=========================================================================="<<std::endl;
      Generator.InitializeTheEventHandler();
      for (int i=1;i<=nevt;i++) {
	if (i%500==0) {
	  ATOOLS::msg.Out()<<" Event "<<i<<std::endl;      
	}
	if (Generator.GenerateOneEvent()) ATOOLS::msg.Events()<<"Sherpa : Passed "<<i<<" events."<<std::endl;
      }
      Generator.SummarizeRun();
    }
    ATOOLS::msg.Out()<<"=========================================================================="<<std::endl
		     <<"Sherpa finished its simulation run with "
		     <<Generator.NumberOfErrors()<<" errors."<<std::endl
		     <<"=========================================================================="<<std::endl;
    return 0;
  }
  catch (ATOOLS::Exception exception) {
    exception.UpdateLogFile();
    std::cout<<exception<<std::endl;
    std::terminate();
  }
  catch (std::exception exception) {
    std::cout<<"Sherpa: throws std::exception "<<exception.what()<<" ..."<<std::endl;
    std::terminate();
  }
}




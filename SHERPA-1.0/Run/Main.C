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
    ATOOLS::msg.Out()<<" Process initialization started "<<std::endl;
    Generator.InitializeTheRun(argc,argv);
    if (ATOOLS::rpa.gen.NumberOfEvents()>0) {
      ATOOLS::msg.Out()<<"generate "<<ATOOLS::rpa.gen.NumberOfEvents()<<" events"<<std::endl;
      Generator.InitializeTheEventHandler();
      int nevt=ATOOLS::rpa.gen.NumberOfEvents();
      if (nevt>0) ATOOLS::msg.Out()<<"Starting event generation now. "<<std::endl;
      for (int i=1;i<=nevt;i++) {
	if (i%500==0) {
	  ATOOLS::msg.Out()<<" Event "<<i<<std::endl;      
	}
	if (Generator.GenerateOneEvent()) ATOOLS::msg.Events()<<"Sherpa : Passed "<<i<<" events."<<std::endl;
      }
      Generator.SummarizeRun();
      ATOOLS::msg.Events()<<"Sherpa did "<<nevt<<" with "<<Generator.NumberOfErrors()<<" errors."<<std::endl;
    }
    ATOOLS::msg.Out()<<" Simulation finished "<<std::endl;
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




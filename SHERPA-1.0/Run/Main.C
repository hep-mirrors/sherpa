#include "Sherpa.H"
#include "Message.H"
#include "prof.hh"
#include "Random.H"
#include "Exception.H"

using namespace std;

extern "C" {
  void apainit_();
  void aparun_();
}

using namespace SHERPA;

int main(int argc,char* argv[]) 
{  
  std::set_terminate(ATOOLS::Exception::Terminate);
  std::set_unexpected(ATOOLS::Exception::Terminate);
  signal(SIGSEGV,ATOOLS::Exception::SEGVHandler);
  signal(SIGINT,ATOOLS::Exception::INTHandler);
  signal(SIGFPE,ATOOLS::Exception::FPEHandler);
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
	  ATOOLS::msg.Out()<<" Event "<<i<<endl;      
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
    if (exception.ApproveTerminate()) terminate();
  }
  catch (std::exception exception) {
    std::cout<<"Sherpa: throws std::exception "<<exception.what()<<" ..."<<std::endl;
    terminate();
  }
}




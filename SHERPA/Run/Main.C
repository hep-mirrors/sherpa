#include "Sherpa.H"
#include "Message.H"
#include "prof.hh"
#include "Random.H"
#include "Exception.H"
#include "Run_Parameter.H"
#include "CXXFLAGS.H"
#include "CXXFLAGS_PACKAGES.H"

#ifdef USING__ROOT
#include "My_Root.H"
#endif

using namespace SHERPA;

#ifdef F77_MAIN
extern "C" int F77_MAIN(int argc,char* argv[]) 
#else
int main(int argc,char* argv[]) 
#endif
{
  ATOOLS::exh->Init();
#ifdef USING__ROOT
  MYROOT::myroot = new MYROOT::My_Root(argc,argv);
  ATOOLS::exh->AddTerminatorObject(MYROOT::myroot);
#endif
  std::set_terminate(ATOOLS::Terminate);
  std::set_unexpected(ATOOLS::Terminate);
  signal(SIGSEGV,ATOOLS::SignalHandler);
  signal(SIGINT,ATOOLS::SignalHandler);
  signal(SIGPIPE,ATOOLS::SignalHandler);
  signal(SIGBUS,ATOOLS::SignalHandler);
  signal(SIGFPE,ATOOLS::SignalHandler);
  signal(SIGABRT,ATOOLS::SignalHandler);
  signal(SIGTERM,ATOOLS::SignalHandler);
  signal(SIGXCPU,ATOOLS::SignalHandler);
  try {
    Sherpa Generator;
    Generator.InitializeTheRun(argc,argv);
    int nevt=ATOOLS::rpa.gen.NumberOfEvents();
    if (nevt>0) {
      msg_Events()<<"=========================================================================="<<std::endl
		       <<"Sherpa will start event generation now : "
		       <<nevt<<" events"<<std::endl
		       <<"=========================================================================="<<std::endl;
      Generator.InitializeTheEventHandler();
      double starttime=ATOOLS::rpa.gen.Timer().UserTime();
      for (int i=1;i<=ATOOLS::rpa.gen.NumberOfEvents();i++) {
	if (i%100==0) {
	  double diff=ATOOLS::rpa.gen.Timer().UserTime()-starttime;
	  msg_Info()<<"  Event "<<i<<" ( "
		    <<int(diff)<<" s elapsed / "
		    <<int((nevt-i)/(double)i*diff)
		    <<" s left / "<<int(nevt/(double)i*diff)
		    <<" s total )   "<<ATOOLS::bm::cr<<std::flush;
	}
	if (Generator.GenerateOneEvent()) msg_Events()<<"Sherpa : Passed "<<i<<" events."<<std::endl;
      }
      msg_Info()<<std::endl;      
      Generator.SummarizeRun();
    }
    msg_Events()<<"=========================================================================="<<std::endl
		     <<"Sherpa finished its simulation run with "
		     <<Generator.NumberOfErrors()<<" errors."<<std::endl
		     <<"=========================================================================="<<std::endl;
#ifdef USING__ROOT
    delete MYROOT::myroot;
#endif
    return 0;
  }
  catch (ATOOLS::Exception exception) {
    exception.UpdateLogFile();
    msg_Error()<<exception<<std::endl;
    std::terminate();
  }
  catch (std::exception exception) {
    std::cout<<"Sherpa: throws std::exception "<<exception.what()<<" ..."<<std::endl;
    std::terminate();
  }
}




#include "SHERPA/Main/Sherpa.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"

#ifdef USING__ROOT
#include "ATOOLS/Org/My_Root.H"
#endif

using namespace SHERPA;

#ifdef FC_DUMMY_MAIN
extern "C" int FC_DUMMY_MAIN() { return 1; }
#endif

int main(int argc,char* argv[])
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
    Sherpa *Generator(new Sherpa());
    Generator->InitializeTheRun(argc,argv);
    int nevt=ATOOLS::rpa.gen.NumberOfEvents();
    if (nevt>0) {
      msg_Events()<<"=========================================================================="<<std::endl
		       <<"Sherpa will start event generation now : "
		       <<nevt<<" events"<<std::endl
		       <<"=========================================================================="<<std::endl;
      Generator->InitializeTheEventHandler();
      double starttime=ATOOLS::rpa.gen.Timer().RealTime();
      for (int i=1;i<=ATOOLS::rpa.gen.NumberOfEvents();i++) {
	if (i%100==0) {
	  double diff=ATOOLS::rpa.gen.Timer().RealTime()-starttime;
	  msg_Info()<<"  Event "<<i<<" ( "
		    <<int(diff)<<" s elapsed / "
		    <<int((nevt-i)/(double)i*diff)
		    <<" s left ) -> ETA: "<<ATOOLS::rpa.gen.Timer().
	    StrFTime("%a %b %d %H:%M",time_t((nevt-i)/(double)i*diff))<<"  ";
	  if (ATOOLS::rpa.gen.BatchMode()&2) { msg_Info()<<std::endl; }
	  else { msg_Info()<<ATOOLS::bm::cr<<std::flush; }
	}
	if (Generator->GenerateOneEvent()) msg_Events()<<"Sherpa : Passed "<<i<<" events."<<std::endl;
      }
      msg_Info()<<std::endl;      
      Generator->SummarizeRun();
//       std::string sigmas(ATOOLS::rpa.gen.Variable("TOTAL_CROSS_SECTION"));
//       if (sigmas=="") msg_Error()<<"Error: xs info not available"<<std::endl;
//       double sigma(ATOOLS::ToType<double>(sigmas));
    }
    msg_Events()<<"=========================================================================="<<std::endl
		     <<"Sherpa finished its simulation run with "
		     <<Generator->NumberOfErrors()<<" errors."<<std::endl
		     <<"=========================================================================="<<std::endl;
    delete Generator;
#ifdef USING__ROOT
    delete MYROOT::myroot;
#endif
    delete ATOOLS::exh;
    return 0;
  }
  catch (ATOOLS::Exception exception) {
    exception.UpdateLogFile();
    msg_Error()<<exception<<std::endl;
    std::terminate();
  }
}




#include "SHERPA/Main/Sherpa.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/MyTiming.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Org/Data_Reader.H"

using namespace SHERPA;

#ifdef FC_DUMMY_MAIN
extern "C" int FC_DUMMY_MAIN() { return 1; }
#endif

int main(int argc,char* argv[])
{
  try {
    Sherpa *Generator(new Sherpa(argc,argv));
    int nevt=ATOOLS::rpa->gen.NumberOfEvents();
    if (nevt>0) {
      msg_Events()<<"=========================================================================="<<std::endl
		       <<"Sherpa will start event generation now : "
		       <<nevt<<" events"<<std::endl
		       <<"=========================================================================="<<std::endl;
      Generator->InitializeTheEventHandler();
      ATOOLS::Data_Reader read(" ",";","!","=");
      ATOOLS::msg->SetLevel(read.GetValue<int>("EVT_OUTPUT",ATOOLS::msg->Level()));
      double starttime=ATOOLS::rpa->gen.Timer().RealTime();
      for (int i=1;i<=ATOOLS::rpa->gen.NumberOfEvents();i++) {
	if (i%100==0 && i<ATOOLS::rpa->gen.NumberOfEvents()) {
	  double diff=ATOOLS::rpa->gen.Timer().RealTime()-starttime;
	  msg_Info()<<"  Event "<<i<<" ( "
		    <<ATOOLS::FormatTime(size_t(diff))<<" elapsed / "
		    <<ATOOLS::FormatTime(size_t((nevt-i)/(double)i*diff))
		    <<" left ) -> ETA: "<<ATOOLS::rpa->gen.Timer().
	    StrFTime("%a %b %d %H:%M",time_t((nevt-i)/(double)i*diff))<<"  ";
	  if (ATOOLS::rpa->gen.BatchMode()&2) { msg_Info()<<std::endl; }
	  else { msg_Info()<<ATOOLS::bm::cr<<std::flush; }
	}
	if (Generator->GenerateOneEvent()) msg_Events()<<"Sherpa : Passed "<<i<<" events."<<std::endl;
      }
      msg_Info()<<"  Event "<<ATOOLS::rpa->gen.NumberOfEvents()<<" ( "
		<<size_t(ATOOLS::rpa->gen.Timer().RealTime()-starttime)
		<<" s total )                                         "<<std::endl;      
      Generator->SummarizeRun();
    }
    msg_Events()<<"=========================================================================="<<std::endl
		     <<"Sherpa finished its simulation run with "
		     <<Generator->NumberOfErrors()<<" errors."<<std::endl
		     <<"=========================================================================="<<std::endl;
    delete Generator;
    return 0;
  }
  catch (ATOOLS::Exception exception) {
    msg_Error()<<exception<<std::endl;
    std::terminate();
  }
}




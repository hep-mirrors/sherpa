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
#include "SHERPA/Single_Events/Event_Handler.H"
#ifdef USING__MPI
#include "mpi.h"
#endif

using namespace SHERPA;
using namespace ATOOLS;

#ifdef FC_DUMMY_MAIN
extern "C" int FC_DUMMY_MAIN() { return 1; }
#endif

int main(int argc,char* argv[])
{
#ifdef USING__MPI
  MPI::Init(argc,argv);
#endif
  Sherpa *Generator(new Sherpa());
  try {
    Generator->InitializeTheRun(argc,argv);
    int nevt=rpa->gen.NumberOfEvents();
    if (nevt>0) {
      Generator->InitializeTheEventHandler();
      Data_Reader read(" ",";","!","=");
      msg->SetLevel(read.GetValue<int>("EVT_OUTPUT",msg->Level()));
      double starttime=rpa->gen.Timer().RealTime();
      for (size_t i=1;i<=rpa->gen.NumberOfEvents();) {
        if (Generator->GenerateOneEvent()) {
          msg_Events()<<"Sherpa : Passed "<<i<<" events."<<std::endl;
          int exp;
          for (exp=5; i/int(pow(10,exp))==0; --exp) {}
          if (((rpa->gen.BatchMode()&4 && i%100==0) ||
               (!(rpa->gen.BatchMode()&4) && i%int(pow(10,exp))==0)) &&
              i<rpa->gen.NumberOfEvents()) {
            double diff=rpa->gen.Timer().RealTime()-starttime;
            msg_Info()<<"  Event "<<i<<" ( "
                      <<FormatTime(size_t(diff))<<" elapsed / "
                      <<FormatTime(size_t((nevt-i)/(double)i*diff))
                      <<" left ) -> ETA: "<<rpa->gen.Timer().
              StrFTime("%a %b %d %H:%M",time_t((nevt-i)/(double)i*diff))<<"  ";
            double xs(Generator->GetEventHandler()->TotalXSMPI());
            double err(Generator->GetEventHandler()->TotalErrMPI());
            if (!(rpa->gen.BatchMode()&2)) msg_Info()<<"\n  ";
            msg_Info()<<"XS = "<<xs<<" pb +- ( "<<err<<" pb = "
                      <<((int(err/xs*10000))/100.0)<<" % )  ";
            if (!(rpa->gen.BatchMode()&2))
              msg_Info()<<mm(1,mm::up);
            if (rpa->gen.BatchMode()&2) { msg_Info()<<std::endl; }
            else { msg_Info()<<bm::cr<<std::flush; }
          }
          ++i;
        }
      }
      msg_Info()<<"  Event "<<rpa->gen.NumberOfEvents()<<" ( "
		<<size_t(rpa->gen.Timer().RealTime()-starttime)
		<<" s total )                                         "<<std::endl;      
      Generator->SummarizeRun();
    }
  }
  catch (Exception exception) {
    msg_Error()<<exception<<std::endl;
    std::terminate();
  }
  delete Generator;
#ifdef USING__MPI
  MPI::Finalize();
#endif
  return 0;
}




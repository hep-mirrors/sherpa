#include "SHERPA/Main/Sherpa.H"
#include "ATOOLS/Org/Terminator_Objects.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Org/My_MPI.H"

using namespace SHERPA;
using namespace ATOOLS;

#ifdef FC_DUMMY_MAIN
extern "C" int FC_DUMMY_MAIN() { return 1; }
#endif

int main(int argc,char* argv[])
{

#ifdef USING__MPI
  MPI_Init(&argc, &argv);
#endif

  try {
    Sherpa* Generator = new Sherpa(argc, argv);
    Generator->InitializeTheRun();
    int nevt=rpa->gen.NumberOfEvents();
    if (nevt>0) {
      Generator->InitializeTheEventHandler();
      for (size_t i=1;i<=rpa->gen.NumberOfEvents();) {
        if (Generator->GenerateOneEvent()) ++i;
      }
      Generator->SummarizeRun();
    }
    delete Generator;
  }
  catch (const normal_exit& exception) {
    msg_Error() << exception << std::endl;
    exh->Terminate(0);
  }
  catch (const Exception& exception) {
    msg_Error() << exception << std::endl;
    exh->Terminate(1);
  }
  catch (const std::exception& exception) {
    msg_Error() << exception.what() << std::endl;
    exh->Terminate(1);
  }

#ifdef USING__MPI
  MPI_Finalize();
#endif

  return 0;
}

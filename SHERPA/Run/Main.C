#include "SHERPA/Main/Sherpa.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
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
#ifdef USING__Threading
  int provided=MPI::Init_thread(argc,argv,MPI_THREAD_SERIALIZED);
  if (provided<MPI_THREAD_SERIALIZED) {
    printf("MPI library does not provide required thread support\n");
    MPI::Finalize();
    return 1;
  }
#else
  MPI::Init(argc,argv);
#endif
#endif
  Sherpa *Generator(new Sherpa());
  try {
    Generator->InitializeTheRun(argc,argv);
    int nevt=rpa->gen.NumberOfEvents();
    if (nevt>0) {
      Generator->InitializeTheEventHandler();
      for (size_t i=1;i<=rpa->gen.NumberOfEvents();) {
        if (Generator->GenerateOneEvent()) ++i;
      }
      Generator->SummarizeRun();
    }
  }
  catch (Exception exception) {
    std::terminate();
  }
  delete Generator;
#ifdef USING__MPI
  MPI::Finalize();
#endif
  return 0;
}

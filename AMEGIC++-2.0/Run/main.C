#include "Amegic.H"
#include "MyTiming.H"
#include "Environment.H"


#include "Message.H"

//#define _USE_MPI_
#ifdef _USE_MPI_
#include <mpi++.h>
#endif

using namespace AMEGIC;

int main(int argc,char* argv[]) 
{    
  
#ifdef _USE_MPI_
  MPI::Init(argc, argv);
#endif
  std::string path("./");
  Environment environment(path,std::string("Run.dat"));
  environment.InitializeTheEnvironment();

  Amegic generator(path,environment.GetMEFile(),environment.GetModel());

  /*    
	if (generator.InitializeDecays()) {
	generator.CalculateBranchingWidths();
	}
  */

  if (generator.InitializeProcesses(environment.GetBeamSpectraHandler(),
				    environment.GetISRHandler())) { 
    generator.CalculateTotalXSec();
    generator.SingleEvents();
  }
  return 1;
  
#ifdef _USE_MPI_
  MPI::Finalize();
#endif

}














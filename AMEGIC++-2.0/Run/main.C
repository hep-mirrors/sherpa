#include "Amegic.H"
#include "MyTiming.H"



#include "Run_Parameter.H"
#include "Running_AlphaS.H"
#include "Running_AlphaQED.H"
#include "Flow.H"
#include "Message.H"

//#define _USE_MPI_
#ifdef _USE_MPI_
#include <mpi++.h>
#endif

using namespace AMEGIC;
using namespace APHYTOOLS;
using namespace AORGTOOLS;


int main(int argc,char* argv[]) 
{    
  
#ifdef _USE_MPI_
  MPI::Init(argc, argv);
#endif
  
  int runmode = AMEGIC::AMPLITUDE_MODE; 

  std::string path("Testrun");
  if (argc==2) {
    if (std::string(argv[1])==std::string("-XS"))         runmode = AMEGIC::XS_MODE;
    else if (std::string(argv[1])==std::string("-P")) runmode = AMEGIC::PATCH_MODE;
    else path = std::string(argv[1]);
  }

  ParticleInit(path);
  rpa.Init(path);
  

  if (!as)   as   = new Running_AlphaS;
  if (!aqed) aqed = new Running_AlphaQED;
  
  Amegic generator(path,0);
  
  generator.DecCalc();

  if (generator.InitializeProcesses(runmode)) { 
    generator.CalculateTotalXSec();
    generator.SingleEvents();
  }

  msg.Tracking()<<"AMEGIC is deleting Running_AlphaS ..."<<std::endl;
  if (as)   delete as; 
  if (aqed) delete aqed;
  msg.Tracking()<<"Deleted couplings in the wrapper "<<std::endl;
  return 1;

#ifdef _USE_MPI_
  MPI::Finalize();
#endif


}














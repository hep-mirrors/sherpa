#include "Amegic.H"
#include "MyTiming.H"



#include "Run_Parameter.H"
#include "Running_AlphaS.H"
#include "Running_AlphaQED.H"
#include "Flow.H"
#include "Message.H"
//#include <mpi++.h>


using namespace AMEGIC;
using namespace APHYTOOLS;
using namespace AORGTOOLS;



int main(int argc,char* argv[]) 
{    
  std::string path("Testrun");
  if (argc==2) path = std::string(argv[1]);

  APHYTOOLS::particle_init(path);
  AORGTOOLS::rpa.Init(path);
  

  if (!as)   as   = new Running_AlphaS;
  if (!aqed) aqed = new Running_AlphaQED;
  
  Amegic Test(path,0);

  if (Test.InitializeProcesses()) { 
    Test.CalculateTotalXSec();
    Test.SingleEvents();
  }

  msg.Tracking()<<"AMEGIC is deleting Running_AlphaS ..."<<endl;
  if (as)   delete as; 
  if (aqed) delete aqed;
  msg.Tracking()<<"Deleted couplings in the wrapper "<<std::endl;
  return 1;
}














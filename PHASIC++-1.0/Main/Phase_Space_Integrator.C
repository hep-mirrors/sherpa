#include "Phase_Space_Integrator.H"
#include "Phase_Space_Handler.H"
#include "Run_Parameter.H"
#include "Message.H"

#include "Random.H"

//#define _USE_MPI_
#ifdef _USE_MPI_
#include <mpi++.h>
#endif

using namespace PHASIC;
using namespace AORGTOOLS;
using namespace AMATOOLS;
using namespace std;

Phase_Space_Integrator::Phase_Space_Integrator() {
  nmax      = 100000000;
  // nmax      = 100;
}


double Phase_Space_Integrator::Calculate(Phase_Space_Handler * psh,double maxerror) 
{
  msg.Tracking()<<"Starting the calculation. Lean back and enjoy ... ."<<endl; 
  if (maxerror >= 1.) nmax = 1;

  result = max = error = 0.;

  int numberofchannels = 1;

  msg.Tracking()<<"Integrators : "<<psh->BeamIntegrator()<<" / "
		<<psh->ISRIntegrator()<<" / "<<psh->FSRIntegrator()<<endl;


  if ((psh->BeamIntegrator())) {
    (psh->BeamIntegrator())->Reset();
    numberofchannels *= psh->NumberOfBeamIntegrators();
    msg.Tracking()<<"   found "<<psh->NumberOfBeamIntegrators()<<" beam integrators."<<endl;
  }
  if ((psh->ISRIntegrator())) {
    (psh->ISRIntegrator())->Reset();
    numberofchannels *= psh->NumberOfISRIntegrators();
    msg.Tracking()<<"   found "<<psh->NumberOfISRIntegrators()<<" isr integrators."<<endl;
  }

  (psh->FSRIntegrator())->Reset();
  numberofchannels *= psh->NumberOfFSRIntegrators();
  msg.Tracking()<<"   found "<<psh->NumberOfFSRIntegrators()<<" fsr integrators."<<endl;

  //iter      = Max(1+int(numberofchannels/maxerror),20000);
  //nopt      = Max(1+int(numberofchannels/maxerror),10);

  iter      = Max(1+10*int(numberofchannels),20000);
  nopt      = 10; 

  maxopt    = iter*nopt;

  bool flag = 0;
  long int  n;
  int       endopt = 0;
  double    value;

  /*
  int nran= Ran.WriteOutStatus("Random_C.dat");
  cout<< " Written: "<< nran <<endl;

  Ran.ReadInStatus("Random_A.dat",0);
  */

  cout.precision(10);

  for (n=1;n<=nmax;n++) {
    value = psh->Differential();

    if ((psh->BeamIntegrator())) (psh->BeamIntegrator())->AddPoint(value);    
    if ((psh->ISRIntegrator()))  (psh->ISRIntegrator())->AddPoint(value);    
    (psh->FSRIntegrator())->AddPoint(value);    
    psh->AddPoint(value);

    if (value>max) max = value;
    if ((!(n%iter)) || (n==maxopt)) {
      msg.Tracking()<<" n="<<n<<"  iter="<<iter<<"  maxopt="<<maxopt<<endl;
      if ((n<=maxopt) && (endopt<2)) {
	if ((psh->BeamIntegrator())) (psh->BeamIntegrator())->Optimize(maxerror);
	if ((psh->ISRIntegrator()))  (psh->ISRIntegrator())->Optimize(maxerror);
	(psh->FSRIntegrator())->Optimize(maxerror);
      }
      if ((n==maxopt) && (endopt<2)) {
	if ((psh->BeamIntegrator())) (psh->BeamIntegrator())->EndOptimize(maxerror);
	if ((psh->ISRIntegrator()))  (psh->ISRIntegrator())->EndOptimize(maxerror);
	(psh->FSRIntegrator())->EndOptimize(maxerror);
	iter   *= Min(int(1./maxerror),(psh->FSRIntegrator())->Number());
	maxopt += 4*iter;
	endopt++;
      }
      if (!((psh->FSRIntegrator())->Result()>0) && !((psh->FSRIntegrator())->Result()<0)) {
	msg.Error()<<"FS - Channel result is a NaN. Knockout!!!!"<<endl;
	break;
      }
      if ( AMATOOLS::IsZero((psh->FSRIntegrator())->Result()) ) break;
      error = (psh->FSRIntegrator())->Variance()/(psh->FSRIntegrator())->Result() * 
	(psh->FSRIntegrator())->N();
      msg.Tracking()<<(psh->FSRIntegrator())->Result()/(psh->FSRIntegrator())->N() * 
	rpa.Picobarn()<<" pb"
		    <<" +- ("<<(psh->FSRIntegrator())->Variance()*rpa.Picobarn()
		    <<" pb = "<<error*100<<"% )."<<endl;
      if (error<maxerror) break;
    }
  }
  result = (psh->FSRIntegrator())->Result() / (psh->FSRIntegrator())->N() * rpa.Picobarn();
  return result;
}



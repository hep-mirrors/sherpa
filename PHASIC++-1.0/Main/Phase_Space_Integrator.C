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
using namespace ATOOLS;
using namespace std;

long int Phase_Space_Integrator::nmax=10000000;                  

double Phase_Space_Integrator::Calculate(Phase_Space_Handler * psh,double maxerror) 
{
  msg.Tracking()<<"Starting the calculation. Lean back and enjoy ... ."<<endl; 
  if (maxerror >= 1.) nmax = 1;

  result = max = error = 0.;


#ifdef _USE_MPI_
  int rank = MPI::COMM_WORLD.Get_rank();
  int size = MPI::COMM_WORLD.Get_size();
  
  if (psh->ISRIntegrator() || psh->BeamIntegrator()) {
    cerr<<" Parallel ISR Integration not supported yet. Sorry."<<endl;
    if (size>1) abort();
  }
#else
  //MPI off
  
  int size = 1;
  int rank = 0;
#endif

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

  iter      = Max(1+10*int(numberofchannels),20000);
  nopt      = 10; 

  maxopt    = iter*nopt;

  long int  n;
  int       endopt = 1;
  double    value;

  ran.SetSeed(-100*rank);
  
  // ------ total sums for MPI ---
  long int alln;
  int iterall;
  double allsum  = 0.;
  double allsum2 = 0.;
  // ------- local sums ------
  struct {
    int    id;
    int    n;
    double sum;
    double sum2;
    double max;
  } local, message;
  
#ifdef _USE_MPI_
  
  local.id   = rank;
  local.n    = 0;
  local.sum  = 0.;
  local.sum2 = 0.;
  local.max  = 0.;

  MPI::Datatype Mes_Type;
  MPI::Datatype type[5]     = {MPI::INT,MPI::INT,MPI::DOUBLE,MPI::DOUBLE,MPI::DOUBLE};
  int           blocklen[5] = {1,1,1,1,1};
  MPI::Aint     disp[5];

  disp[0] = MPI::Get_address(&message.id);
  disp[1] = MPI::Get_address(&message.n);
  disp[2] = MPI::Get_address(&message.sum);
  disp[3] = MPI::Get_address(&message.sum2);
  disp[4] = MPI::Get_address(&message.max);
    
  for (short int i=4;i>=0;i--) disp[i] -= disp[0];

  Mes_Type = MPI::Datatype::Create_struct(5,blocklen,disp,type);
  Mes_Type.Commit();
  // --- end define Message ---

  if (rank==0) {
    for (short int i=1;i<size;i++) MPI::COMM_WORLD.Send(&iterall, 1, MPI::INT, i, 10);   
    msg.Out()<<"Node "<<rank<<"starts calculation with "<<iterall<<" steps."<<endl;
  }
  else {
    MPI::COMM_WORLD.Recv(&iterall, 1, MPI::INT, 0, MPI::ANY_TAG);
    msg.Out()<<"Node "<<rank<<"starts calculation with "<<iterall<<" steps."<<endl;
  }
#endif

  int over = 0;
  int saveiter = 0;

  for (n=1;n<=nmax;n++) {
    value = psh->Differential();

    if ((psh->BeamIntegrator())) (psh->BeamIntegrator())->AddPoint(value);    
    if ((psh->ISRIntegrator()))  (psh->ISRIntegrator())->AddPoint(value);    
    (psh->FSRIntegrator())->AddPoint(value);    
    psh->AddPoint(value);

    local.sum+=value;
    local.sum2+=sqr(value);
    local.n++; 


    if (value>max) message.max = max = value;
#ifdef _USE_MPI_
    if (rank>0) {  // slaves
              
      if (!(n%iterall)) {
	// check if last message was recieved!
	// create message
	message=local;
	// send
	MPI::COMM_WORLD.Isend(&message, 1, Mes_Type, 0, 1);  // to master
	local.sum=0.;
	local.sum2=0.;
	local.n=0;

	double local_result = (psh->FSRIntegrator())->Result() / (psh->FSRIntegrator())->N() * rpa.Picobarn();
	double local_error = (psh->FSRIntegrator())->Variance()/(psh->FSRIntegrator())->Result() * 
	  (psh->FSRIntegrator())->N();

	msg.Out()<<n<<" Node "<<rank<<" : Result : "<<local_result<<" pb +- "
		 <<local_error*100<<" % "<<"steps "<<n<<endl;
      }
      // change iteration steps or end integrations
      if (MPI::COMM_WORLD.Iprobe(0, 10)) {	
	int olditerall=iterall;
	MPI::COMM_WORLD.Recv(&iterall, 1, MPI::INT, 0, 10);
	if (iterall<=n/olditerall)   iterall=n/olditerall+1;
	msg.Out()<<"Slave "<<rank<<" changed iteration steps to "<<iterall<<endl; 
	if (iterall==0) break;  
      }
      // check for optimizations
      if (MPI::COMM_WORLD.Iprobe(0, 5)) {
	//Tag 5 means optimization
	int opt;
	MPI::COMM_WORLD.Recv(&opt, 1, MPI::INT, 0, 5);
	msg.Out()<<"Slave "<<rank<<" received opt-tag: "<<opt<<endl;
	if (opt==1) psh->FSRIntegrator()->MPIOptimize(maxerror);
	if (opt==2) psh->FSRIntegrator()->EndOptimize(maxerror);
      }
    }
    else { // master

      if (!(n%iterall)) {
	saveiter+= local.n;
	alln    += local.n;
	allsum  += local.sum;
	allsum2 += local.sum2;
	local.sum=0.;
	local.sum2=0.;
	local.n=0;

  	double local_result = (psh->FSRIntegrator())->Result() / (psh->FSRIntegrator())->N() * rpa.Picobarn();
	double local_error = (psh->FSRIntegrator())->Variance()/(psh->FSRIntegrator())->Result() * 
	(psh->FSRIntegrator())->N();
	msg.Out()<<n<<" Node 0 : Result : "<<local_result<<" +- "<<local_error*100.<<" % "<<endl;
	
	// total 
	double total_result=allsum/alln*rpa.Picobarn();
	double total_error=sqrt((allsum2/alln-(sqr(allsum)-allsum2)/alln/(alln-1))/alln)/allsum*alln;
	
	// abort conditions
	if (total_error<maxerror || alln>=nmax) over = 1;	
      }
      // receive data from slaves:

      while (MPI::COMM_WORLD.Iprobe(MPI::ANY_SOURCE, MPI::ANY_TAG)) {
	MPI::COMM_WORLD.Recv(&message, 1, Mes_Type, MPI::ANY_SOURCE, MPI::ANY_TAG);	
	msg.Out()<<" received Data from Knot "<<message.id<<endl;

	saveiter+= message.n;
	alln    += message.n;
	allsum  += message.sum;
	allsum2 += message.sum2;
	if (message.max>max) max=message.max;

	int olditerall = iterall;
	
	if (saveiter>=iter) { // check if optimization steps neccessary
	  saveiter = 0;
	  if (alln<=maxopt) {
	    msg.Out()<<"Start optimization"<<endl;
	    //Send Optimize Message with tag 5
	    int opt = 1;
	    for (short int i=1;i<size;i++) MPI::COMM_WORLD.Isend(&opt, 1, MPI::INT, i, 5);   
	    psh->FSRIntegrator()->MPIOptimize(maxerror);
	  }
	  if (alln>=maxopt && alln<maxopt+iter) {
	    msg.Out()<<"End optimization"<<endl;
	    int opt = 2;
	    for (short int i=1;i<size;i++) MPI::COMM_WORLD.Isend(&opt, 1, MPI::INT, i, 5);   
	    psh->FSRIntegrator()->EndOptimize(maxerror);
	    iterall = 20000;
	    iter    = 20000;
	  }
	  if (alln+iterall*size>nmax) {
	    iterall=(nmax-alln)/size+1;
	  }
	}
       
	double total_result=allsum/alln*rpa.Picobarn();
	double total_error=sqrt((allsum2/alln-(sqr(allsum)-allsum2)/alln/(alln-1))/alln)/allsum*alln;
	msg.Out()<<alln<<" Master : Result : "<<total_result
		 <<" pb +- "<<total_error*100.<<" % "<<endl;
	
	if (total_error<maxerror || alln>=nmax) over = 1;	
 
	if (iterall!=olditerall) {
	  for (short int i=1;i<size;i++) MPI::COMM_WORLD.Isend(&iterall, 1, MPI::INT, i, 10);   
	}
      }     
      if (over) break;
    }
#endif

    if ((!(n%iter)) || (n==maxopt)) {
#ifndef _USE_MPI_ // non MPI mode
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
	iter   *= 2;
	maxopt += 4*iter;
	endopt++;
      }

      if (!((psh->FSRIntegrator())->Result()>0.) && !((psh->FSRIntegrator())->Result()<0.)
	  && !((psh->FSRIntegrator())->Result()==0.)) {
	msg.Error()<<"FS - Channel result is a NaN. Knockout!!!!"<<endl;
	break;
      }
      if ( ATOOLS::IsZero((psh->FSRIntegrator())->Result()) ) break;
      error = (psh->FSRIntegrator())->Variance()/(psh->FSRIntegrator())->Result() * 
	(psh->FSRIntegrator())->N();
      msg.Tracking()<<(psh->FSRIntegrator())->Result()/(psh->FSRIntegrator())->N() * rpa.Picobarn()<<" pb"
		    <<" +- ("<<(psh->FSRIntegrator())->Variance()*rpa.Picobarn()
		    <<" pb = "<<error*100<<"% )."<<endl;
      if (error<maxerror) break;
#endif

    }
  }

  result = (psh->FSRIntegrator())->Result() / (psh->FSRIntegrator())->N() * rpa.Picobarn();
  return result;
  
}

double Phase_Space_Integrator::CalculateDecay(Phase_Space_Handler* psh,double mass, double maxerror) 
{ 
  msg.Tracking()<<"Starting the calculation for a decay. Lean back and enjoy ... ."<<endl; 
  
  iter      = 20000;
  nopt      = 10; 

  maxopt    = iter*nopt;

  // double flux = 1./(2.*mass);
  long int n;
  double value;
  double max = 0.;
  double error;
  
  (psh->FSRIntegrator())->Reset();

  double oldvalue = 0.;

  for (n=1;n<=nmax;n++) {
    do { value = psh->Differential(); }
    while (value > 1./ATOOLS::Accu());
    (psh->FSRIntegrator())->AddPoint(value);
    
    //new SS
    psh->AddPoint(value);
    
    if (value>max) max = value;

    if (value!=0. && value==oldvalue) break;
    oldvalue = value;
    
    //Nan Check
    if (!((psh->FSRIntegrator())->Result()>0.) && !((psh->FSRIntegrator())->Result()<0.)) {
      msg.Out()<<"NaN knockout!!!!"<<endl;
      break;
    }
    
    if (!(n%iter)) {
      if (n<=maxopt) psh->FSRIntegrator()->Optimize(maxerror);
      if (n==maxopt) {
	psh->FSRIntegrator()->EndOptimize(maxerror);
	iter = 50000;
      }
      //Nan Check
      if (!((psh->FSRIntegrator())->Result()>0.) && !((psh->FSRIntegrator())->Result()<0.)) {
	msg.Out()<<"NaN knockout!!!!"<<endl;
	break;
      }
      if (ATOOLS::IsZero((psh->FSRIntegrator())->Result())) break;
      
      msg.Out()<<n<<". Result : "<<(psh->FSRIntegrator())->Result()/(psh->FSRIntegrator())->N()<<" GeV";
      error = (psh->FSRIntegrator())->Variance() / (psh->FSRIntegrator())->Result()*(psh->FSRIntegrator())->N();
      msg.Out()<<" +- "<<error*100<<"%"<<endl;
      if (error<maxerror) break;
    }
  }
  return (psh->FSRIntegrator())->Result() / (psh->FSRIntegrator())->N();
}

long int Phase_Space_Integrator::MaxPoints()                  
{ return nmax; };

void     Phase_Space_Integrator::SetMaxPoints(long int _nmax) 
{ nmax=_nmax;  };


#include "Phase_Space_Integrator.H"
#include "Phase_Space_Handler.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "Single_Channel.H"
#include "PI_Interface.H"

#include "Random.H"

//#define _USE_MPI_
#ifdef _USE_MPI_
#include <mpi++.h>
#endif

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

long int Phase_Space_Integrator::nmax=10000000;                  

double Phase_Space_Integrator::Calculate(Phase_Space_Handler * psh,double maxerror, int fin_opt) 
{
  p_psh=psh;
  msg_Info()<<"Starting the calculation. Lean back and enjoy ... ."<<endl; 
  if (maxerror >= 1.) nmax = 1;

  result = max = error = 0.;


#ifdef _USE_MPI_
  int rank = MPI::COMM_WORLD.Get_rank();
  int size = MPI::COMM_WORLD.Get_size();
  
  if (psh->ISRIntegrator() || psh->BeamIntegrator()) {
    cerr<<" Parallel ISR Integration not supported yet. Sorry."<<endl;
    if (size>1) abort();
  }
  int size = 1;
#else
  //MPI off
#endif

  int numberofchannels = 1;

  msg_Tracking()<<"Integrators : "<<psh->BeamIntegrator()<<" / "
		<<psh->ISRIntegrator()<<" / "<<psh->FSRIntegrator()
		<<" / "<<psh->KMRZIntegrator()<<" / "<<psh->KMRKPIntegrator()<<endl;
  
   if ((psh->BeamIntegrator())) {
     (psh->BeamIntegrator())->Reset();
     numberofchannels = psh->NumberOfBeamIntegrators();
     msg.Tracking()<<"   Found "<<psh->NumberOfBeamIntegrators()<<" Beam Integrators."<<endl;
   }
   if ((psh->ISRIntegrator())) {
     (psh->ISRIntegrator())->Reset();
     numberofchannels += psh->NumberOfISRIntegrators();
     msg.Tracking()<<"   Found "<<psh->NumberOfISRIntegrators()<<" ISR Integrators."<<endl;
   }

  if ((psh->KMRZIntegrator())) {
    (psh->KMRZIntegrator())->Reset();
    numberofchannels += psh->NumberOfKMRZIntegrators();
    msg.Tracking()<<"   Found "<<psh->NumberOfKMRZIntegrators()<<" KMR z Integrators."<<endl;
  }
  if ((psh->KMRKPIntegrator())) {
    (psh->KMRKPIntegrator())->Reset();
    numberofchannels += psh->NumberOfKMRKPIntegrators();
    msg.Tracking()<<"   Found "<<psh->NumberOfKMRKPIntegrators()<<" KMR k_\\perp Integrators."<<endl;
  }



  (psh->FSRIntegrator())->Reset();
  numberofchannels += psh->NumberOfFSRIntegrators();
  msg.Tracking()<<"   Found "<<psh->NumberOfFSRIntegrators()<<" FSR integrators."<<endl;

  iter = iter0 = Max(20*int(numberofchannels),5000);
  iter1      = Max(100*int(numberofchannels),10000);
  if (iter1>50000) iter1=Max(iter0,50000);
  int hlp = (iter1-1)/iter0+1;
  iter1   = hlp*iter0;
  nopt      = 25; 

  maxopt    = (5/hlp+21)*iter1;
  int ncontrib = psh->FSRIntegrator()->ValidN();
  if (ncontrib/iter0>=5) iter=iter1;

  long int  n;
  int       endopt = 1;
  double    value;
  int nlo=0;
  if (ncontrib>maxopt) endopt=2;
  
#ifdef _USE_MPI_
  // ------ total sums for MPI ---
  ran.SetSeed(-100*rank);
  long int alln;
  int iterall;
  double allsum  = 0.;
  double allsum2 = 0.;
#endif
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

#ifdef _USE_MPI_
  int over = 0;
  int saveiter = 0;
#endif

  for (n=ATOOLS::Max(psh->Process()->Points(),(long int)0)+1;n<=nmax;n++) {
    if (!rpa.gen.CheckTime()) {
      ATOOLS::msg.Error()<<ATOOLS::om::bold
			 <<"\nPhase_Space_Integrator::Calculate(): "
			 <<ATOOLS::om::reset<<ATOOLS::om::red
			 <<"Timeout. Interrupt integration."
			 <<ATOOLS::om::reset<<std::endl;
      kill(getpid(),SIGINT);
    }

    value = psh->Differential();

    if ((psh->BeamIntegrator())) (psh->BeamIntegrator())->AddPoint(value);    
    if ((psh->ISRIntegrator()))  (psh->ISRIntegrator())->AddPoint(value);    
    if ((psh->KMRZIntegrator()))  (psh->KMRZIntegrator())->AddPoint(value);    
    if ((psh->KMRKPIntegrator()))  (psh->KMRKPIntegrator())->AddPoint(value);    
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
    ncontrib = psh->FSRIntegrator()->ValidN();
    if ( ncontrib!=nlo && ((ncontrib%iter)==0 || ncontrib==maxopt)) {
      nlo=ncontrib;
#ifndef _USE_MPI_ // non MPI mode
      msg_Tracking()<<" n="<<ncontrib<<"  iter="<<iter<<"  maxopt="<<maxopt<<endl;
      if ((ncontrib<=maxopt) && (endopt<2)) {
	if ((psh->BeamIntegrator())) (psh->BeamIntegrator())->Optimize(maxerror);
	if ((psh->ISRIntegrator()))  (psh->ISRIntegrator())->Optimize(maxerror);
	if ((psh->KMRZIntegrator()))  (psh->KMRZIntegrator())->Optimize(maxerror);
	if ((psh->KMRKPIntegrator()))  (psh->KMRKPIntegrator())->Optimize(maxerror);
	(psh->FSRIntegrator())->Optimize(maxerror);
	(psh->Process())->ResetMax(2);
	if (ncontrib%iter1==0) (psh->Process())->OptimizeResult();
      }
      else {
	(psh->Process())->ResetMax(0);
      }
      if ((ncontrib==maxopt) && (endopt<2)) {
	if ((psh->BeamIntegrator())) (psh->BeamIntegrator())->EndOptimize(maxerror);
	if ((psh->ISRIntegrator()))  (psh->ISRIntegrator())->EndOptimize(maxerror);
	if ((psh->KMRZIntegrator()))  (psh->KMRZIntegrator())->EndOptimize(maxerror);
	if ((psh->KMRKPIntegrator()))  (psh->KMRKPIntegrator())->EndOptimize(maxerror);
	(psh->FSRIntegrator())->EndOptimize(maxerror);
	iter   *= 2;
	maxopt += 4*iter;
	endopt++;
	(psh->Process())->ResetMax(1);
	(psh->Process())->InitWeightHistogram();
      }

      if (!((psh->Process())->TotalResult()>0.) && !((psh->Process())->TotalResult()<0.)
	  && !((psh->Process())->TotalResult()==0.)) {
	msg.Error()<<"FS - Channel result is a NaN. Knockout!!!!"<<endl;
	break;
      }
      error=(psh->Process())->TotalVar()/(psh->Process())->TotalResult();
      msg_Info()<<om::blue
		<<(psh->Process())->TotalResult()*rpa.Picobarn()
		<<" pb"<<om::reset<<" +- ( "<<om::red
		<<(psh->Process())->TotalVar()*rpa.Picobarn()
		<<" pb = "<<error*100<<" %"<<om::reset<<" )."<<endl;
      if (ncontrib/iter0==5) iter=iter1;
      bool allowbreak = true;
      if (fin_opt==1 && (endopt<2||ncontrib<maxopt)) allowbreak = false;
      if (p_psh->PI()!=NULL && ncontrib/iter<10) allowbreak = false;
      if (error<maxerror && allowbreak) break;
      if (ncontrib/iter0==5 && p_psh->UsePI()>0 && p_psh->PI()==NULL) {
	p_psh->CreatePI();
	if (p_psh->PI()==NULL) THROW(fatal_error,"Cannot initialize PI.");
	p_psh->PI()->Initialize();
      }
#endif

    }
  }

  result = (psh->Process())->TotalResult() * rpa.Picobarn();
  return result;
  
}

double Phase_Space_Integrator::CalculateDecay(Phase_Space_Handler* psh,double mass, double maxerror) 
{ 
  msg_Info()<<"Starting the calculation for a decay. Lean back and enjoy ... ."<<endl; 
  
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
    if (!((psh->Process())->TotalResult()>0.) && !((psh->Process())->TotalResult()<0.)) {
      msg.Out()<<"NaN knockout!!!!"<<endl;
      break;
    }
    
    if (!(n%iter)) {
      if (n<=maxopt) {
	psh->FSRIntegrator()->Optimize(maxerror);
	(psh->Process())->OptimizeResult();
      }
      if (n==maxopt) {
	psh->FSRIntegrator()->EndOptimize(maxerror);
	iter = 50000;
      }
      //Nan Check
      if (!((psh->Process())->TotalResult()>0.) && !((psh->Process())->TotalResult()<0.)) {
	msg.Out()<<"NaN knockout!!!!"<<endl;
	break;
      }
      if ((psh->Process())->TotalResult()==0.) break;
      
      error = (psh->Process())->TotalVar() / (psh->Process())->TotalResult();
      msg_Info()<<n<<". Result : "<<(psh->Process())->TotalResult()
		<<" GeV"<<" +- "<<error*100<<"%"<<endl;
      if (error<maxerror) break;
    }
  }
  return (psh->Process())->TotalResult();
}

long int Phase_Space_Integrator::MaxPoints()                  
{ return nmax; };

void     Phase_Space_Integrator::SetMaxPoints(long int _nmax) 
{ nmax=_nmax;  };


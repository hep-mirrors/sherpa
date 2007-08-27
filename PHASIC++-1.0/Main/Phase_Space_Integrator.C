#include "Phase_Space_Integrator.H"
#include "Phase_Space_Handler.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "Single_Channel.H"
#include "PI_Interface.H"
#include "Algebra_Interpreter.H"
#include "MyStrStream.H"
#include "Helicity_Integrator.H"
#include "Color_Integrator.H"
#include "Data_Reader.H"

#include "Random.H"
#include <unistd.h>

//#define _USE_MPI_
#ifdef _USE_MPI_
#include <mpi++.h>
#endif

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

struct TDouble: public Term {
  double m_value;
};// end of struct Double

struct TVec4D: public Term {
  Vec4D m_value;
  TVec4D(const Vec4D &value): m_value(value) {}
};// end of struct Vec4D

long unsigned int Phase_Space_Integrator::nmax=std::numeric_limits<long unsigned int>::max();
Phase_Space_Integrator::Phase_Space_Integrator():
  p_interpreter(NULL) 
{
  local.id=local.n=0;
  local.sum=local.sum2=local.max=0.0;
  addtime=0.0;
}

Phase_Space_Integrator::~Phase_Space_Integrator()
{
  if (p_interpreter!=NULL) delete p_interpreter;
}

double Phase_Space_Integrator::Calculate(Phase_Space_Handler *_psh,double _maxerror, int _fin_opt) 
{
  local.sum=local.sum2=0.0;
  local.n=0; 
  maxerror=_maxerror;
  fin_opt=_fin_opt;
  psh=_psh;
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
     msg_Tracking()<<"   Found "<<psh->NumberOfBeamIntegrators()<<" Beam Integrators."<<endl;
   }
   if ((psh->ISRIntegrator())) {
     (psh->ISRIntegrator())->Reset();
     numberofchannels += psh->NumberOfISRIntegrators();
     msg_Tracking()<<"   Found "<<psh->NumberOfISRIntegrators()<<" ISR Integrators."<<endl;
   }

  if ((psh->KMRZIntegrator())) {
    (psh->KMRZIntegrator())->Reset();
    numberofchannels += psh->NumberOfKMRZIntegrators();
    msg_Tracking()<<"   Found "<<psh->NumberOfKMRZIntegrators()<<" KMR z Integrators."<<endl;
  }
  if ((psh->KMRKPIntegrator())) {
    (psh->KMRKPIntegrator())->Reset();
    numberofchannels += psh->NumberOfKMRKPIntegrators();
    msg_Tracking()<<"   Found "<<psh->NumberOfKMRKPIntegrators()<<" KMR k_\\perp Integrators."<<endl;
  }



  (psh->FSRIntegrator())->Reset();
  numberofchannels += psh->NumberOfFSRIntegrators();
  msg_Tracking()<<"   Found "<<psh->NumberOfFSRIntegrators()<<" FSR integrators."<<endl;
  iter = iter0 = Max(20*int(numberofchannels),5000);
  iter1      = Max(100*int(numberofchannels),10000);
  if (iter1>50000) iter1=Max(iter0,50000);
  int hlp = (iter1-1)/iter0+1;
  iter1   = hlp*iter0;
  nopt      = 25;

  maxopt    = (5/hlp+21)*iter1;
  ncontrib = psh->FSRIntegrator()->ValidN();
  if (ncontrib/iter0>=5) iter=iter1;

  endopt = 1;
  nlo=0;
  if (ncontrib>maxopt) endopt=2;

  addtime = 0.0;
  lotime = starttime = ATOOLS::rpa.gen.Timer().UserTime();
  if (psh->Stats().size()>0)
    addtime=psh->Stats().back()[6];
  totalopt  = maxopt+8.*iter1;

  nstep = ncstep = 0;
  
#ifdef _USE_MPI_
  // ------ total sums for MPI ---
  ran.SetSeed(-100*rank);
  allsum  = 0.;
  allsum2 = 0.;
#endif
  
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
    msg_Out()<<"Node "<<rank<<"starts calculation with "<<iterall<<" steps."<<endl;
  }
  else {
    MPI::COMM_WORLD.Recv(&iterall, 1, MPI::INT, 0, MPI::ANY_TAG);
    msg_Out()<<"Node "<<rank<<"starts calculation with "<<iterall<<" steps."<<endl;
  }
#endif

#ifdef _USE_MPI_
  int over = 0;
  int saveiter = 0;
#endif

  for (n=ATOOLS::Max(psh->Process()->Points(),(long int)0)+1;n<=nmax;n++) {
    if (!rpa.gen.CheckTime()) {
      msg_Error()<<ATOOLS::om::bold
			 <<"\nPhase_Space_Integrator::Calculate(): "
			 <<ATOOLS::om::reset<<ATOOLS::om::red
			 <<"Timeout. Interrupt integration."
			 <<ATOOLS::om::reset<<std::endl;
      kill(getpid(),SIGINT);
    }

    value = psh->Differential();
    if (AddPoint(value)) break;
    
  }
  
  msg_Info()<<om::blue
	    <<(psh->Process())->TotalResult()*rpa.Picobarn()
	    <<" pb"<<om::reset<<" +- ( "<<om::red
	    <<(psh->Process())->TotalVar()*rpa.Picobarn()
	    <<" pb = "<<error*100<<" %"<<om::reset<<" ) "
	    <<ncontrib<<" ( "<<(ncontrib*1000/n)/10.0<<" % )"<<endl;
  result = (psh->Process())->TotalResult() * rpa.Picobarn();
  return result;
  
}

bool Phase_Space_Integrator::AddPoint(const double value)
{
  nstep++;
  if (value>0.) ncstep++;
  Integrable_Base *proc=psh->Process();
  std::string func=proc->EnhanceFunction();
  double enhance=1.0;
  if (func!="1") {
    if (p_interpreter==NULL) {
      p_interpreter = new Algebra_Interpreter();
      p_interpreter->SetTagReplacer(this);
      for (size_t i=0;i<proc->NIn()+proc->NOut();++i) 
	p_interpreter->AddTag("p["+ToString(i)+"]",ToString(psh->Point()[i]));
      p_interpreter->Interprete(func);
    }
    enhance=((TDouble*)p_interpreter->Calculate())->m_value;
  }
  
  if (psh->ColorIntegrator()==NULL || psh->ColorIntegrator()->ValidPoint()) {
    if ((psh->BeamIntegrator())) (psh->BeamIntegrator())->AddPoint(value*enhance);    
    if ((psh->ISRIntegrator()))  (psh->ISRIntegrator())->AddPoint(value*enhance);    
    if ((psh->KMRZIntegrator()))  (psh->KMRZIntegrator())->AddPoint(value*enhance);    
    if ((psh->KMRKPIntegrator()))  (psh->KMRKPIntegrator())->AddPoint(value*enhance);    
    if (psh->HelicityIntegrator()) psh->HelicityIntegrator()->AddPoint(value*enhance);
    (psh->FSRIntegrator())->AddPoint(value*enhance);    
  }
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

	msg_Out()<<n<<" Node "<<rank<<" : Result : "<<local_result<<" pb +- "
		 <<local_error*100<<" % "<<"steps "<<n<<endl;
      }
      // change iteration steps or end integrations
      if (MPI::COMM_WORLD.Iprobe(0, 10)) {	
	int olditerall=iterall;
	MPI::COMM_WORLD.Recv(&iterall, 1, MPI::INT, 0, 10);
	if (iterall<=n/olditerall)   iterall=n/olditerall+1;
	msg_Out()<<"Slave "<<rank<<" changed iteration steps to "<<iterall<<endl; 
	if (iterall==0) return true;  
      }
      // check for optimizations
      if (MPI::COMM_WORLD.Iprobe(0, 5)) {
	//Tag 5 means optimization
	int opt;
	MPI::COMM_WORLD.Recv(&opt, 1, MPI::INT, 0, 5);
	msg_Out()<<"Slave "<<rank<<" received opt-tag: "<<opt<<endl;
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
	msg_Out()<<n<<" Node 0 : Result : "<<local_result<<" +- "<<local_error*100.<<" % "<<endl;
	
	// total 
	double total_result=allsum/alln*rpa.Picobarn();
	double total_error=sqrt((allsum2/alln-(sqr(allsum)-allsum2)/alln/(alln-1))/alln)/allsum*alln;
	
	// abort conditions
	if (total_error<maxerror || alln>=nmax) over = 1;	
      }
      // receive data from slaves:

      while (MPI::COMM_WORLD.Iprobe(MPI::ANY_SOURCE, MPI::ANY_TAG)) {
	MPI::COMM_WORLD.Recv(&message, 1, Mes_Type, MPI::ANY_SOURCE, MPI::ANY_TAG);	
	msg_Out()<<" received Data from Knot "<<message.id<<endl;

	saveiter+= message.n;
	alln    += message.n;
	allsum  += message.sum;
	allsum2 += message.sum2;
	if (message.max>max) max=message.max;

	int olditerall = iterall;
	
	if (saveiter>=iter) { // check if optimization steps neccessary
	  saveiter = 0;
	  if (alln<=maxopt) {
	    msg_Out()<<"Start optimization"<<endl;
	    //Send Optimize Message with tag 5
	    int opt = 1;
	    for (short int i=1;i<size;i++) MPI::COMM_WORLD.Isend(&opt, 1, MPI::INT, i, 5);   
	    psh->FSRIntegrator()->MPIOptimize(maxerror);
	  }
	  if (alln>=maxopt && alln<maxopt+iter) {
	    msg_Out()<<"End optimization"<<endl;
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
	msg_Out()<<alln<<" Master : Result : "<<total_result
		 <<" pb +- "<<total_error*100.<<" % "<<endl;
	
	if (total_error<maxerror || alln>=nmax) over = 1;	
 
	if (iterall!=olditerall) {
	  for (short int i=1;i<size;i++) MPI::COMM_WORLD.Isend(&iterall, 1, MPI::INT, i, 10);   
	}
      }     
      if (over) return true;
    }
#endif
    ncontrib = psh->FSRIntegrator()->ValidN();
    if ( ncontrib!=nlo && ncontrib>0 && ((ncontrib%iter)==0 || ncontrib==maxopt)) {
      nlo=ncontrib;
      bool optimized=false;
#ifndef _USE_MPI_ // non MPI mode
      bool fotime = false;
      msg_Tracking()<<" n="<<ncontrib<<"  iter="<<iter<<"  maxopt="<<maxopt<<endl;
      if ((ncontrib<=maxopt) && (endopt<2)) {
	if ((psh->BeamIntegrator())) (psh->BeamIntegrator())->Optimize(maxerror);
	if ((psh->ISRIntegrator()))  (psh->ISRIntegrator())->Optimize(maxerror);
	if ((psh->KMRZIntegrator()))  (psh->KMRZIntegrator())->Optimize(maxerror);
	if ((psh->KMRKPIntegrator()))  (psh->KMRKPIntegrator())->Optimize(maxerror);
	if ((psh->HelicityIntegrator()))  (psh->HelicityIntegrator())->Optimize();
	(psh->FSRIntegrator())->Optimize(maxerror);
	(psh->Process())->ResetMax(2);
	if (ncontrib%iter1==0) {
	  (psh->Process())->OptimizeResult();
	  if ((psh->Process())->SPoints()==0) lotime = ATOOLS::rpa.gen.Timer().UserTime();
	}
	fotime = true;
	if ((psh->FSRIntegrator())->OptimizationFinished()) { 
	  if (!(psh->ISRIntegrator()) || ncontrib/iter1>=8) { 
	    maxopt=ncontrib;
	  }
	}
	optimized=true;
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
	psh->UpdateIntegrators();
	iter   = iter0;
	maxopt += 4*iter;
	endopt++;
	(psh->Process())->ResetMax(1);
	(psh->Process())->InitWeightHistogram();
  lotime = ATOOLS::rpa.gen.Timer().UserTime();
	return false;
      }

      if (!((psh->Process())->TotalResult()>0.) && !((psh->Process())->TotalResult()<0.)
	  && !((psh->Process())->TotalResult()==0.)) {
	msg_Error()<<"FS - Channel result is a NaN. Knockout!!!!"<<endl;
	return true;
      }
      
      double time = ATOOLS::rpa.gen.Timer().UserTime();
      double timeest=0.;
      timeest = totalopt/double(ncontrib)*(time-starttime);
      if (!fotime) {
	if (fin_opt==1)
	  timeest = ATOOLS::Max(timeest,(psh->Process())->RemainTimeFactor(maxerror)*
				(time-lotime)+lotime-starttime);
	else timeest = (psh->Process())->RemainTimeFactor(maxerror)*
	       (time-lotime)+lotime-starttime;
      }
      error=(psh->Process())->TotalVar()/(psh->Process())->TotalResult();
      msg_Info()<<om::blue
		<<(psh->Process())->TotalResult()*rpa.Picobarn()
		<<" pb"<<om::reset<<" +- ( "<<om::red
		<<(psh->Process())->TotalVar()*rpa.Picobarn()
		<<" pb = "<<error*100<<" %"<<om::reset<<" ) "
		<<ncontrib<<" ( "<<(ncstep*1000/nstep)/10.0
		<<" % )"<<endl;
      if (optimized) nstep = ncstep = 0;
      if (fotime) {
	msg_Info()<<"full optimization: ";
      }
      else msg_Info()<<"integration time: ";
      msg_Info()<<" ( "<<int(time-starttime)<<" s elapsed / "
		    <<int(timeest)-int((time-starttime))
		    <<" s left / "<<int(timeest)
		    <<" s total )   "<<endl;

      std::vector<double> stats(6);
      stats[0]=psh->Process()->TotalResult()*rpa.Picobarn();
      stats[1]=psh->Process()->TotalVar()*rpa.Picobarn();
      stats[2]=error;
      stats[3]=ncontrib;
      stats[4]=ncontrib/(double)n;
      stats[5]=time-starttime+addtime;
      psh->AddStats(stats);

      if (ncontrib/iter0==5) iter=iter1;
      bool allowbreak = true;
      if (fin_opt==1 && (endopt<2||ncontrib<maxopt)) allowbreak = false;
      if (error<maxerror && allowbreak) return true;
#endif

    }
    return false;
}

double Phase_Space_Integrator::CalculateDecay(Phase_Space_Handler* psh,double mass, double maxerror) 
{ 
  msg_Info()<<"Starting the calculation for a decay. Lean back and enjoy ... ."<<endl; 
  
  iter      = 20000;
  nopt      = 10; 

  maxopt    = iter*nopt;

  // double flux = 1./(2.*mass);
  long unsigned int n;
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
      msg_Out()<<"NaN knockout!!!!"<<endl;
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
	msg_Out()<<"NaN knockout!!!!"<<endl;
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
{ return nmax; }

void     Phase_Space_Integrator::SetMaxPoints(long int _nmax) 
{ nmax=_nmax;  }

std::string Phase_Space_Integrator::ReplaceTags(std::string &expr) const
{
  return p_interpreter->ReplaceTags(expr);
}

Term *Phase_Space_Integrator::ReplaceTags(Term *term) const
{
  int i=ATOOLS::ToType<int>(term->m_tag.substr(2,term->m_tag.length()-3));
  ((TVec4D*)term)->m_value=psh->Point()[i];
  return term;
}

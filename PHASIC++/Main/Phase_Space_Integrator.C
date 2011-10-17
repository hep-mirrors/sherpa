#include "PHASIC++/Main/Phase_Space_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "PHASIC++/Channels/Single_Channel.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "PHASIC++/Main/Process_Integrator.H"

#include "ATOOLS/Math/Random.H"
#include <unistd.h>
#ifdef USING__MPI
#include "mpi.h"
#endif

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

long unsigned int Phase_Space_Integrator::nmax=
  std::numeric_limits<long unsigned int>::max();
Phase_Space_Integrator::Phase_Space_Integrator()
{
  Data_Reader read(" ",";","!","=");
  read.AddComment("#");
  read.AddWordSeparator("\t");
  read.SetInputPath(rpa->GetPath());
  read.SetInputFile(rpa->gen.Variable("INTEGRATION_DATA_FILE"));
  if (!read.ReadFromFile(nmax,"PSI_NMAX")) 
    nmax=std::numeric_limits<long unsigned int>::max();
  else msg_Info()<<METHOD<<"(): Set n_{max} = "<<nmax<<".\n";
  if (!read.ReadFromFile(itmin,"PSI_ITMIN")) itmin=5000;
  else msg_Info()<<METHOD<<"(): Set n_{iter} = "<<itmin<<".\n";
  if (!read.ReadFromFile(wadjust,"PSI_ADJUST_POINTS")) wadjust=1;
  addtime=0.0;
}

Phase_Space_Integrator::~Phase_Space_Integrator()
{
}

void Phase_Space_Integrator::MPISync()
{
#ifdef USING__MPI
  double nrtime=ATOOLS::rpa->gen.Timer().RealTime();
  exh->MPISync();
  psh->MPISync();
  int size=MPI::COMM_WORLD.Get_size(), nact=1;
  if (size>1) {
    int rank=MPI::COMM_WORLD.Get_rank();
    double values[4];
    if (rank==0) {
      double trtime=(nrtime-lrtime)/mn;
      std::vector<double> times(size,nrtime-lrtime);
      for (int tag=1;tag<size;++tag) {
	if (!exh->MPIStat(tag)) continue;
	MPI::COMM_WORLD.Recv(&values,4,MPI::DOUBLE,MPI::ANY_SOURCE,tag);
	mn+=values[0];
	mnstep+=values[1];
	mncstep+=values[2];
	trtime+=times[tag]=values[3]/values[0];
	++nact;
      }
      int sum=0, max=0, min=std::numeric_limits<int>::max();
      for (int tag=1;tag<size;++tag) {
	if (!exh->MPIStat(tag)) continue;
	if (!wadjust) sum+=times[tag]=(int)iter/size;
	else sum+=times[tag]=(int)Max(10.0,iter*times[tag]/trtime);
	min=Min(min,(int)times[tag]);
	max=Max(max,(int)times[tag]);
      }
      optiter=iter-sum;
      min=Min(min,optiter);
      max=Max(max,optiter);
      if (wadjust) {
	msg_Info()<<"MPI point range: "
		  <<min<<" .. "<<max<<" ("<<nact<<" nodes)"<<std::endl;
	if (msg_LevelIsTracking() || wadjust>1) {
	msg_Info()<<"New weights {\n  master: "<<optiter<<"\n";
	for (int tag=1;tag<size;++tag)
	  msg_Info()<<"  node "<<tag<<": "<<times[tag]<<"\n";
	msg_Info()<<"}"<<std::endl;
	}
      }
      values[0]=mn;
      values[1]=mnstep;
      values[2]=mncstep;
      for (int tag=1;tag<size;++tag) {
	if (!exh->MPIStat(tag)) continue;
	values[3]=times[tag];
	MPI::COMM_WORLD.Send(&values,4,MPI::DOUBLE,tag,size+tag);
      }
    }
    else {
      values[0]=mn;
      values[1]=mnstep;
      values[2]=mncstep;
      values[3]=nrtime-lrtime;
      MPI::COMM_WORLD.Send(&values,4,MPI::DOUBLE,0,rank);
      MPI::COMM_WORLD.Recv(&values,4,MPI::DOUBLE,0,size+rank);
      mn=values[0];
      mnstep=values[1];
      mncstep=values[2];
      optiter=values[3];
    }
  }
  n+=mn;
  nstep+=mnstep;
  ncstep+=mncstep;
  mn=mnstep=mncstep=0;
  ncontrib=psh->FSRIntegrator()->ValidN();
  nlo=0;
#else
  nlo=psh->FSRIntegrator()->ValidN();
#endif
  lrtime=ATOOLS::rpa->gen.Timer().RealTime();
}

double Phase_Space_Integrator::Calculate(Phase_Space_Handler *_psh,double _maxerror, int _fin_opt) 
{
  mn=mnstep=mncstep=0;
  maxerror=_maxerror;
  fin_opt=_fin_opt;
  psh=_psh;
  msg_Info()<<"Starting the calculation. Lean back and enjoy ... ."<<endl; 
  if (maxerror >= 1.) nmax = 1;

  int numberofchannels = 1;

  msg_Tracking()<<"Integrators : "<<psh->BeamIntegrator()<<" / "
		<<psh->ISRIntegrator()<<" / "<<psh->FSRIntegrator()<<endl;
  
   if ((psh->BeamIntegrator())) {
     (psh->BeamIntegrator())->Reset();
     numberofchannels = psh->BeamIntegrator()->NChannels();
     msg_Tracking()<<"   Found "<<psh->BeamIntegrator()->NChannels()<<" Beam Integrators."<<endl;
   }
   if ((psh->ISRIntegrator())) {
     (psh->ISRIntegrator())->Reset();
     numberofchannels += psh->ISRIntegrator()->NChannels();
     msg_Tracking()<<"   Found "<<psh->ISRIntegrator()->NChannels()<<" ISR Integrators."<<endl;
   }

  (psh->FSRIntegrator())->Reset();
  numberofchannels += psh->FSRIntegrator()->NChannels();
  msg_Tracking()<<"   Found "<<psh->FSRIntegrator()->NChannels()<<" FSR integrators."<<endl;
  iter = iter0 = Max(itmin,Max((int)psh->Process()->ItMin(),Max(20*int(numberofchannels),5000)));
  iter1      = Max(2*itmin,Max(2*(int)psh->Process()->ItMin(),Max(100*int(numberofchannels),10000)));
  if (iter1>50000) iter1=Max(iter0,50000);
  int hlp = (iter1-1)/iter0+1;
  iter1   = hlp*iter0;
  nopt      = 25;

  maxopt    = (5/hlp+21)*iter1;
  ncontrib = psh->FSRIntegrator()->ValidN();
  if (ncontrib/iter0>=5) iter=iter1;

  endopt = 1;
#ifdef USING__MPI
  nlo=0;
#else
  nlo=psh->FSRIntegrator()->ValidN();
#endif
  if (ncontrib>maxopt) endopt=2;

  addtime = 0.0;
#ifdef USING__Threading
  rlotime = rstarttime = ATOOLS::rpa->gen.Timer().RealTime();
#endif
  lotime = starttime = ATOOLS::rpa->gen.Timer().UserTime();
  if (psh->Stats().size()>0)
    addtime=psh->Stats().back()[6];
  totalopt  = maxopt+8.*iter1;

  nstep = ncstep = 0;

  lrtime = ATOOLS::rpa->gen.Timer().RealTime();
  optiter=iter;
#ifdef USING__MPI
  int size = MPI::COMM_WORLD.Get_size();
  optiter /= size;
  if (MPI::COMM_WORLD.Get_rank()==0) optiter+=iter-(iter/size)*size;
#endif
  
  for (n=psh->Process()->Points();n<=nmax;) {
    if (!rpa->gen.CheckTime()) {
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
  
  return psh->Process()->TotalResult() * rpa->Picobarn();
  
}

bool Phase_Space_Integrator::AddPoint(const double value)
{
  if (IsBad(value)) {
    msg_Error()<<METHOD<<"(): value = "<<value<<". Skip."<<endl;
    return false;
  }
      
#ifdef USING__MPI
  ++mn;
  mnstep++;
  if (value!=0.) mncstep++;
#else
  ++n;
  nstep++;
  if (value!=0.) ncstep++;
#endif
  
    psh->AddPoint(value);

#ifdef USING__MPI
    ncontrib = psh->FSRIntegrator()->ValidMN();
#else
    ncontrib = psh->FSRIntegrator()->ValidN();
#endif
    if ( ncontrib!=nlo && ncontrib>0 && ((ncontrib%optiter)==0 || ncontrib==maxopt)) {
      MPISync();
      bool optimized=false;
      bool fotime = false;
      msg_Tracking()<<" n="<<ncontrib<<"  iter="<<iter<<"  maxopt="<<maxopt<<endl;
      if ((ncontrib<=maxopt) && (endopt<2)) {
	psh->Optimize();
	if (ncontrib%iter1==0) {
	  (psh->Process())->OptimizeResult();
	  if ((psh->Process())->SPoints()==0) {
#ifdef USING__Threading
	    rlotime = ATOOLS::rpa->gen.Timer().RealTime();
#endif
	    lotime = ATOOLS::rpa->gen.Timer().UserTime();
	  }
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
      if ((ncontrib>=maxopt) && (endopt<2)) {
	psh->EndOptimize();
	int oiter=iter;
	if (psh->UpdateIntegrators()) iter=iter0;
	else iter*=2;
	optiter*=iter/oiter;
	maxopt += 4*iter;
	endopt++;
	(psh->Process())->ResetMax(1);
	(psh->Process())->InitWeightHistogram();
#ifdef USING__Threading
	rlotime = ATOOLS::rpa->gen.Timer().RealTime();
#endif
	lotime = ATOOLS::rpa->gen.Timer().UserTime();
	return false;
      }

#ifdef USING__Threading
      double rtime = ATOOLS::rpa->gen.Timer().RealTime();
      double rtimeest=0.;
      rtimeest = totalopt/double(ncontrib)*(rtime-rstarttime);
#endif
      double time = ATOOLS::rpa->gen.Timer().UserTime();
      double timeest=0.;
      timeest = totalopt/double(ncontrib)*(time-starttime);
      if (!fotime) {
	if (fin_opt==1) {
	  timeest = ATOOLS::Max(timeest,(psh->Process())->RemainTimeFactor(maxerror)*
				(time-lotime)+lotime-starttime);
#ifdef USING__Threading
	  rtimeest = ATOOLS::Max(rtimeest,(psh->Process())->RemainTimeFactor(maxerror)*
				(rtime-rlotime)+rlotime-rstarttime);
#endif
	}
	else {
	  timeest = (psh->Process())->RemainTimeFactor(maxerror)*
	    (time-lotime)+lotime-starttime;
#ifdef USING__Threading
	  rtimeest = (psh->Process())->RemainTimeFactor(maxerror)*
	    (rtime-rlotime)+rlotime-rstarttime;
#endif
	}
      }
      double error=dabs(psh->Process()->TotalVar()/psh->Process()->TotalResult());
      msg_Info()<<om::blue
		<<(psh->Process())->TotalResult()*rpa->Picobarn()
		<<" pb"<<om::reset<<" +- ( "<<om::red
		<<(psh->Process())->TotalVar()*rpa->Picobarn()
		<<" pb = "<<error*100<<" %"<<om::reset<<" ) "
		<<ncontrib<<" ( "<<n<<" -> "<<(ncstep*1000/nstep)/10.0
		<<" % )"<<endl;
      if (optimized) nstep = ncstep = 0;
      if (fotime) {
	msg_Info()<<"full optimization: ";
      }
      else msg_Info()<<"integration time: ";
#ifdef USING__Threading
      msg_Info()<<" ( "<<FormatTime(size_t(rtime-rstarttime+0.5))<<"("
		<<FormatTime(size_t(time-starttime+0.5))<<") elapsed / "
		<<FormatTime(size_t(rtimeest+0.5)
			     -size_t((rtime-rstarttime+0.5)))<<"("
		<<FormatTime(size_t(timeest+0.5)
			     -size_t((time-starttime+0.5)))
		<<") left )   "<<endl;
#else
      msg_Info()<<" ( "<<FormatTime(size_t(time-starttime))<<" elapsed / " 
		<<FormatTime(size_t(timeest)-size_t((time-starttime))) 
		<<" left )   "<<endl; 
#endif
      std::vector<double> stats(6);
      stats[0]=psh->Process()->TotalResult()*rpa->Picobarn();
      stats[1]=psh->Process()->TotalVar()*rpa->Picobarn();
      stats[2]=error;
      stats[3]=ncontrib;
      stats[4]=ncontrib/(double)n;
      stats[5]=time-starttime+addtime;
      psh->AddStats(stats);
      psh->Process()->StoreResults(1);
      if (ncontrib/iter0==5) iter=iter1;
      bool allowbreak = true;
      if (fin_opt==1 && (endopt<2||ncontrib<maxopt)) allowbreak = false;
      if (dabs(error)<maxerror && allowbreak) return true;

    }
    return false;
}

double Phase_Space_Integrator::CalculateDecay(Phase_Space_Handler* psh,
                                              double maxerror) 
{ 
  mn=mnstep=mncstep=0;
  msg_Info()<<"Starting the calculation for a decay. Lean back and enjoy ... ."<<endl; 
  
  optiter = iter = 20000;
  nopt      = 10; 

  maxopt    = iter*nopt;

  long unsigned int n;
  double value;
  double max = 0.;
  double error;
  
  (psh->FSRIntegrator())->Reset();

  double oldvalue = 0.;

  for (n=1;n<=nmax;n++) {
    do { value = psh->Differential(); }
    while (dabs(value) > 1./ATOOLS::Accu());
    
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
	psh->Optimize();
	(psh->Process())->OptimizeResult();
      }
      if (n==maxopt) {
	psh->EndOptimize();
	optiter = iter = 50000;
      }
      //Nan Check
      if (!((psh->Process())->TotalResult()>0.) && !((psh->Process())->TotalResult()<0.)) {
	msg_Out()<<"NaN knockout!!!!"<<endl;
	break;
      }
      if ((psh->Process())->TotalResult()==0.) break;
      
      error = (psh->Process())->TotalVar() / (psh->Process())->TotalResult();

      msg_Info()<<om::blue
                <<(psh->Process())->TotalResult()
                <<" GeV"<<om::reset<<" +- ( "<<om::red
                <<(psh->Process())->TotalVar()
                <<" GeV = "<<error*100<<" %"<<om::reset<<" ) "<<endl;
      if (error<maxerror) break;
    }
  }
  return (psh->Process())->TotalResult()*rpa->Picobarn();
}

long int Phase_Space_Integrator::MaxPoints()                  
{ return nmax; }

void     Phase_Space_Integrator::SetMaxPoints(long int _nmax) 
{ nmax=_nmax;  }


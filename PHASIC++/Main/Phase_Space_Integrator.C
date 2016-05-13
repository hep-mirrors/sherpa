#include "PHASIC++/Main/Phase_Space_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "PHASIC++/Channels/Single_Channel.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/My_MPI.H"
#include "PHASIC++/Main/Process_Integrator.H"

#include "ATOOLS/Org/RUsage.H"
#include <unistd.h>

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

long unsigned int Phase_Space_Integrator::nmax=
  std::numeric_limits<long unsigned int>::max();
Phase_Space_Integrator::Phase_Space_Integrator(Phase_Space_Handler *_psh):
  psh(_psh)
{
  Data_Reader read(" ",";","!","=");
  read.AddComment("#");
  read.AddWordSeparator("\t");
  read.SetInputPath(rpa->GetPath());
  read.SetInputFile(rpa->gen.Variable("INTEGRATION_DATA_FILE"));
  if (!read.ReadFromFile(nmax,"PSI_NMAX")) 
    nmax=std::numeric_limits<long unsigned int>::max();
  else msg_Info()<<METHOD<<"(): Set n_{max} = "<<nmax<<".\n";
  read.SetAllowUnits(true);
  if (!read.ReadFromFile(itmin,"PSI_ITMIN"))
    itmin=psh->Process()->Process()->Info().m_itmin;
  else msg_Info()<<METHOD<<"(): Set n_{it,min} = "<<itmin<<".\n";
  if (!read.ReadFromFile(itmax,"PSI_ITMAX")) itmax=100*itmin;
  else msg_Info()<<METHOD<<"(): Set n_{it,max} = "<<itmax<<".\n";
  if (!read.ReadFromFile(nopt,"PSI_NOPT")) nopt=25;
  else msg_Info()<<METHOD<<"(): Set n_{opt} = "<<nopt<<".\n";
  if (!read.ReadFromFile(maxopt,"PSI_MAXOPT")) maxopt=5;
  else msg_Info()<<METHOD<<"(): Set n_{maxopt} = "<<maxopt<<".\n";
  if (!read.ReadFromFile(ndecopt,"PSI_NDECOPT")) ndecopt=10;
  else msg_Info()<<METHOD<<"(): Set n_{opt,dec} = "<<ndecopt<<".\n";
  addtime=0.0;
  lastrss=0;
#ifdef USING__MPI
  int size=MPI::COMM_WORLD.Get_size();
  if (size) {
    int helpi;
    if (read.ReadFromFile(helpi,"PSI_ITMIN_BY_NODE")) {
      itmin=helpi*size;
      msg_Info()<<METHOD<<"(): Set n_{it,min} = "<<itmin<<".\n";
    }
    if (read.ReadFromFile(helpi,"PSI_ITMAX_BY_NODE")) {
      itmax*=helpi*size;
      msg_Info()<<METHOD<<"(): Set n_{it,max} = "<<itmax<<".\n";
    }
  }
#endif
}

Phase_Space_Integrator::~Phase_Space_Integrator()
{
}

void Phase_Space_Integrator::MPISync()
{
#ifdef USING__MPI
  psh->MPISync();
  int size=MPI::COMM_WORLD.Get_size();
  if (size>1) {
    double values[3];
    values[0]=mn;
    values[1]=mnstep;
    values[2]=mncstep;
    mpi->MPIComm()->Allreduce(MPI_IN_PLACE,values,3,MPI::DOUBLE,MPI::SUM);
    mn=values[0];
    mnstep=values[1];
    mncstep=values[2];
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

double Phase_Space_Integrator::Calculate(double _maxerror, double _maxabserror, int _fin_opt) 
{
  mn=mnstep=mncstep=0;
  maxerror=_maxerror;
  maxabserror=_maxabserror;
  fin_opt=_fin_opt;
  msg_Info()<<"Starting the calculation at "
	    <<rpa->gen.Timer().StrFTime("%H:%M:%S")<<". Lean back and enjoy ... ."<<endl; 
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
  iter = iter0 = Min(itmax,Max(itmin,Max((int)psh->Process()->ItMin(),20*int(numberofchannels))));
  iter1      = Min(2*itmax,Max(2*itmin,Max(2*(int)psh->Process()->ItMin(),100*int(numberofchannels))));
  int hlp = (iter1-1)/iter0+1;
  iter1   = hlp*iter0;

  ncontrib = psh->FSRIntegrator()->ValidN();
  if (ncontrib/iter0>=6) iter=iter1;

#ifdef USING__MPI
  nlo=0;
#else
  nlo=psh->FSRIntegrator()->ValidN();
#endif

  addtime = 0.0;
#if (defined USING__Threading)
  rlotime = rstarttime = ATOOLS::rpa->gen.Timer().RealTime();
#endif
  lotime = starttime = ATOOLS::rpa->gen.Timer().UserTime();
  if (psh->Stats().size()>0)
    addtime=psh->Stats().back()[6];

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
    if ( ncontrib!=nlo && ncontrib>0 && ((ncontrib%optiter)==0)) {
      MPISync();
      bool optimized=false;
      bool fotime = false;
      msg_Tracking()<<" n="<<ncontrib<<"  iter="<<iter<<endl;
      if (psh->Stats().size()<nopt) {
	psh->Optimize();
	if (ncontrib%iter1==0) {
	  (psh->Process())->OptimizeResult();
	  if ((psh->Process())->SPoints()==0) {
#if (defined USING__Threading)
	    rlotime = ATOOLS::rpa->gen.Timer().RealTime();
#endif
	    lotime = ATOOLS::rpa->gen.Timer().UserTime();
	  }
	}
	fotime = true;
	optimized=true;
      }
      else if (psh->Stats().size()==nopt) {
	(psh->Process())->ResetMax(0);
	psh->EndOptimize();
	int oiter=iter;
	if (psh->UpdateIntegrators()) iter=iter0;
	else iter*=2;
	optiter*=iter/(double)oiter;
	(psh->Process())->ResetMax(1);
	(psh->Process())->InitWeightHistogram();
	(psh->Process())->EndOptimize();
#if (defined USING__Threading)
	rlotime = ATOOLS::rpa->gen.Timer().RealTime();
#endif
	lotime = ATOOLS::rpa->gen.Timer().UserTime();
      }

#if (defined USING__Threading)
      double rtime = ATOOLS::rpa->gen.Timer().RealTime();
      double rtimeest=0.;
      rtimeest = (5*iter0+(nopt-5)*iter1+2*maxopt*iter1)/double(ncontrib)*(rtime-rstarttime);
#endif
      double time = ATOOLS::rpa->gen.Timer().UserTime();
      double timeest=0.;
      timeest = (5*iter0+(nopt-5)*iter1+2*maxopt*iter1)/double(ncontrib)*(time-starttime);
      if (!fotime) {
	if (fin_opt==1) {
	  timeest = ATOOLS::Max(timeest,(psh->Process())->RemainTimeFactor(maxerror)*
				(time-lotime)+lotime-starttime);
#if (defined USING__Threading)
	  rtimeest = ATOOLS::Max(rtimeest,(psh->Process())->RemainTimeFactor(maxerror)*
				(rtime-rlotime)+rlotime-rstarttime);
#endif
	}
	else {
	  timeest = (psh->Process())->RemainTimeFactor(maxerror)*
	    (time-lotime)+lotime-starttime;
#if (defined USING__Threading)
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
#if (defined USING__Threading)
      msg_Info()<<" ( "<<FormatTime(size_t(rtime-rstarttime+0.5))<<" ("
		<<FormatTime(size_t(time-starttime+0.5))<<") elapsed / "
		<<FormatTime(size_t(rtimeest+0.5)
			     -size_t((rtime-rstarttime+0.5)))<<" ("
		<<FormatTime(size_t(timeest+0.5)
			     -size_t((time-starttime+0.5)))
		<<") left ) ["<<rpa->gen.Timer().StrFTime("%H:%M:%S")<<"]   "<<endl;
#else
      msg_Info()<<" ( "<<FormatTime(size_t(time-starttime))<<" elapsed / " 
		<<FormatTime(size_t(timeest)-size_t((time-starttime))) 
		<<" left ) ["<<rpa->gen.Timer().StrFTime("%H:%M:%S")<<"]   "<<endl; 
#endif
      size_t currentrss=GetCurrentRSS();
      if (lastrss==0) lastrss=currentrss;
      else if (currentrss>lastrss+ToType<int>
	  (rpa->gen.Variable("MEMLEAK_WARNING_THRESHOLD"))) {
	msg_Error()<<METHOD<<"() {\n"<<om::bold<<"  Memory usage increased by "
		   <<(currentrss-lastrss)/(1<<20)<<" MB,"
		   <<" now "<<currentrss/(1<<20)<<" MB.\n"
		   <<om::red<<"  This might indicate a memory leak!\n"
		   <<"  Please monitor this process closely.\n"<<om::reset<<"}"<<std::endl;
	lastrss=currentrss;
      }
      std::vector<double> stats(6);
      stats[0]=psh->Process()->TotalResult()*rpa->Picobarn();
      stats[1]=psh->Process()->TotalVar()*rpa->Picobarn();
      stats[2]=error;
      stats[3]=ncontrib;
      stats[4]=ncontrib/(double)n;
      stats[5]=time-starttime+addtime;
      psh->AddStats(stats);
      psh->Process()->StoreResults(1);
      if (ncontrib/iter0==6) {
	optiter=iter=iter1;
#ifdef USING__MPI
	int size = MPI::COMM_WORLD.Get_size();
	optiter /= size;
	if (MPI::COMM_WORLD.Get_rank()==0) optiter+=iter-(iter/size)*size;
#endif
      }
      bool wannabreak = dabs(error)<maxerror ||
        dabs(psh->Process()->TotalVar()*rpa->Picobarn())<maxabserror;
      if (fin_opt==0 && nopt>psh->Stats().size() && wannabreak) nopt=psh->Stats().size();
      if (wannabreak && psh->Stats().size()>=nopt+maxopt) return true;
    }
    return false;
}

double Phase_Space_Integrator::CalculateDecay(double maxerror) 
{ 
  mn=mnstep=mncstep=0;
  msg_Info()<<"Starting the calculation for a decay. Lean back and enjoy ... ."<<endl; 
  
  optiter = iter = 20000;

  long unsigned int n;
  double value;
  double max = 0.;
  double error;
  
  (psh->FSRIntegrator())->Reset();

  for (n=1;n<=nmax;n++) {
    value = psh->Differential();
    psh->AddPoint(value);
    
    if (value>max) max = value;
    
    if (!(n%iter)) {
      MPISync();
      if (psh->Stats().size()<=ndecopt) {
	psh->Optimize();
	(psh->Process())->OptimizeResult();
      }
      if (psh->Stats().size()==ndecopt) {
	psh->EndOptimize();
	optiter = iter = 50000;
      }
      if ((psh->Process())->TotalResult()==0.) break;
      
      error = (psh->Process())->TotalVar() / (psh->Process())->TotalResult();

      msg_Info()<<om::blue
                <<(psh->Process())->TotalResult()
                <<" GeV"<<om::reset<<" +- ( "<<om::red
                <<(psh->Process())->TotalVar()
                <<" GeV = "<<error*100<<" %"<<om::reset<<" ) "<<n<<endl;
      if (error<maxerror) break;
    }
  }
  return (psh->Process())->TotalResult()*rpa->Picobarn();
}

long int Phase_Space_Integrator::MaxPoints()                  
{ return nmax; }

void     Phase_Space_Integrator::SetMaxPoints(long int _nmax) 
{ nmax=_nmax;  }


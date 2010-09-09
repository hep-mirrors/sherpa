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
  read.SetInputPath(rpa.GetPath());
  read.SetInputFile(rpa.gen.Variable("INTEGRATION_DATA_FILE"));
  if (!read.ReadFromFile(nmax,"PSI_NMAX")) 
    nmax=std::numeric_limits<long unsigned int>::max();
  else msg_Info()<<METHOD<<"(): Set n_{max} = "<<nmax<<".\n";
  local.id=local.n=0;
  local.sum=local.sum2=local.max=0.0;
  addtime=0.0;
}

Phase_Space_Integrator::~Phase_Space_Integrator()
{
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

  int numberofchannels = 1;

  msg_Tracking()<<"Integrators : "<<psh->BeamIntegrator()<<" / "
		<<psh->ISRIntegrator()<<" / "<<psh->FSRIntegrator()<<endl;
  
   if ((psh->BeamIntegrator())) {
     (psh->BeamIntegrator())->Reset();
     numberofchannels = psh->BeamIntegrator()->Number();
     msg_Tracking()<<"   Found "<<psh->BeamIntegrator()->Number()<<" Beam Integrators."<<endl;
   }
   if ((psh->ISRIntegrator())) {
     (psh->ISRIntegrator())->Reset();
     numberofchannels += psh->ISRIntegrator()->Number();
     msg_Tracking()<<"   Found "<<psh->ISRIntegrator()->Number()<<" ISR Integrators."<<endl;
   }

  (psh->FSRIntegrator())->Reset();
  numberofchannels += psh->FSRIntegrator()->Number();
  msg_Tracking()<<"   Found "<<psh->FSRIntegrator()->Number()<<" FSR integrators."<<endl;
  iter = iter0 = Max((int)psh->Process()->ItMin(),Max(20*int(numberofchannels),5000));
  iter1      = Max(2*(int)psh->Process()->ItMin(),Max(100*int(numberofchannels),10000));
  if (iter1>50000) iter1=Max(iter0,50000);
  int hlp = (iter1-1)/iter0+1;
  iter1   = hlp*iter0;
  nopt      = 25;

  maxopt    = (5/hlp+21)*iter1;
  ncontrib = psh->FSRIntegrator()->ValidN();
  if (ncontrib/iter0>=5) iter=iter1;

  endopt = 1;
  nlo=psh->FSRIntegrator()->ValidN();
  if (ncontrib>maxopt) endopt=2;

  addtime = 0.0;
#ifdef USING__Threading
  rlotime = rstarttime = ATOOLS::rpa.gen.Timer().RealTime();
#endif
  lotime = starttime = ATOOLS::rpa.gen.Timer().UserTime();
  if (psh->Stats().size()>0)
    addtime=psh->Stats().back()[6];
  totalopt  = maxopt+8.*iter1;

  nstep = ncstep = 0;
  
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
  
  result = (psh->Process())->TotalResult() * rpa.Picobarn();
  return result;
  
}

bool Phase_Space_Integrator::AddPoint(const double value)
{
  if (IsBad(value)) {
    msg_Error()<<METHOD<<"(): value = "<<value<<". Skip."<<endl;
    return false;
  }
      
  nstep++;
  if (value!=0.) ncstep++;
  double enhance=psh->EnhanceFactor();
  
  if (value!=0.0) {
    if ((psh->BeamIntegrator())) (psh->BeamIntegrator())->AddPoint(value*enhance);    
    if ((psh->ISRIntegrator()))  (psh->ISRIntegrator())->AddPoint(value*enhance);    
    (psh->FSRIntegrator())->AddPoint(value*enhance);    
  }
    psh->AddPoint(value);

    local.sum+=value;
    local.sum2+=sqr(value);
    local.n++; 


    if (value>max) message.max = max = value;
    ncontrib = psh->FSRIntegrator()->ValidN();
    if ( ncontrib!=nlo && ncontrib>0 && ((ncontrib%iter)==0 || ncontrib==maxopt)) {
      nlo=ncontrib;
      bool optimized=false;
      bool fotime = false;
      msg_Tracking()<<" n="<<ncontrib<<"  iter="<<iter<<"  maxopt="<<maxopt<<endl;
      if ((ncontrib<=maxopt) && (endopt<2)) {
	if ((psh->BeamIntegrator())) (psh->BeamIntegrator())->Optimize(maxerror);
	if ((psh->ISRIntegrator()))  (psh->ISRIntegrator())->Optimize(maxerror);
	(psh->FSRIntegrator())->Optimize(maxerror);
	psh->Process()->Optimize();
	(psh->Process())->ResetMax(2);
	if (ncontrib%iter1==0) {
	  (psh->Process())->OptimizeResult();
	  if ((psh->Process())->SPoints()==0) {
#ifdef USING__Threading
	    rlotime = ATOOLS::rpa.gen.Timer().RealTime();
#endif
	    lotime = ATOOLS::rpa.gen.Timer().UserTime();
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
      if ((ncontrib==maxopt) && (endopt<2)) {
	if ((psh->BeamIntegrator())) (psh->BeamIntegrator())->EndOptimize(maxerror);
	if ((psh->ISRIntegrator()))  (psh->ISRIntegrator())->EndOptimize(maxerror);
	(psh->FSRIntegrator())->EndOptimize(maxerror);
	psh->Process()->EndOptimize();
	if (psh->UpdateIntegrators()) iter=iter0;
	else iter*=2;
	maxopt += 4*iter;
	endopt++;
	(psh->Process())->ResetMax(1);
	(psh->Process())->InitWeightHistogram();
#ifdef USING__Threading
	rlotime = ATOOLS::rpa.gen.Timer().RealTime();
#endif
	lotime = ATOOLS::rpa.gen.Timer().UserTime();
	return false;
      }

#ifdef USING__Threading
      double rtime = ATOOLS::rpa.gen.Timer().RealTime();
      double rtimeest=0.;
      rtimeest = totalopt/double(ncontrib)*(rtime-rstarttime);
#endif
      double time = ATOOLS::rpa.gen.Timer().UserTime();
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
      error=(psh->Process())->TotalVar()/(psh->Process())->TotalResult();
      msg_Info()<<om::blue
		<<(psh->Process())->TotalResult()*rpa.Picobarn()
		<<" pb"<<om::reset<<" +- ( "<<om::red
		<<(psh->Process())->TotalVar()*rpa.Picobarn()
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
      stats[0]=psh->Process()->TotalResult()*rpa.Picobarn();
      stats[1]=psh->Process()->TotalVar()*rpa.Picobarn();
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
  msg_Info()<<"Starting the calculation for a decay. Lean back and enjoy ... ."<<endl; 
  
  iter      = 20000;
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
	psh->Process()->Optimize();
	(psh->Process())->OptimizeResult();
      }
      if (n==maxopt) {
	psh->FSRIntegrator()->EndOptimize(maxerror);
	psh->Process()->EndOptimize();
	iter = 50000;
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
  return (psh->Process())->TotalResult()*rpa.Picobarn();
}

long int Phase_Space_Integrator::MaxPoints()                  
{ return nmax; }

void     Phase_Space_Integrator::SetMaxPoints(long int _nmax) 
{ nmax=_nmax;  }


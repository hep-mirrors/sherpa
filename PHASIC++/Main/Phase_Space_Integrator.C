#include "PHASIC++/Main/Phase_Space_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "PHASIC++/Channels/Single_Channel.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Default_Reader.H"
#include "ATOOLS/Org/My_MPI.H"
#include "PHASIC++/Main/Process_Integrator.H"

#include "ATOOLS/Org/RUsage.H"
#include <unistd.h>

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

long unsigned int Phase_Space_Integrator::m_nmax(std::numeric_limits<long unsigned int>::max());
long unsigned int Phase_Space_Integrator::m_nrawmax(std::numeric_limits<long unsigned int>::max());

Phase_Space_Integrator::Phase_Space_Integrator(Phase_Space_Handler *_psh):
  m_iter(5000), m_itmin(5000), m_itmax(500000), m_itminbynode(2),
  m_nmin(0), m_nrawmin(0),
  m_n(0), m_nstep(0), m_ncstep(0), m_mn(0), m_mnstep(0), m_mncstep(0),
  m_ncontrib(0), m_maxopt(0), m_stopopt(1000), m_nlo(0), m_fin_opt(1),
  m_starttime(0.), m_lotime(0.), m_addtime(0.), m_lrtime(0.),
  m_maxerror(0.), m_maxabserror(0.),
  m_lastrss(0), p_psh(_psh)
{
  Default_Reader reader;
  reader.SetInputPath(rpa->GetPath());
  reader.SetInputFile(rpa->gen.Variable("INTEGRATION_DATA_FILE"));
  reader.SetAllowUnits(true);

  // total number of points
  m_nmax = reader.Get("PSI_NMAX", std::numeric_limits<long unsigned int>::max(), "n_{max}", METHOD);
  m_nmax = reader.Get("PSI_NMIN", 0, "n_{min}", METHOD);
  m_nrawmax = reader.Get("PSI_NRAWMAX", std::numeric_limits<long unsigned int>::max(), "n_{max,raw}", METHOD);
  m_nrawmax = reader.Get("PSI_NRAWMIN", 0, "n_{min}", METHOD);

  // number of points per iteration
  m_itmin = reader.Get("PSI_ITMIN", p_psh->Process()->Process()->Info().m_itmin, "n_{it,min,raw}", METHOD);
  m_itmax = reader.Get("PSI_ITMAX", 100 * m_itmin, "n_{it,max}", METHOD);

  // number of optimisation steps
  m_nopt = reader.Get("PSI_NOPT", 25, "n_{opt}", METHOD);
  m_maxopt = reader.Get("PSI_MAXOPT", 5, "n_{maxopt}", METHOD);
  m_stopopt = reader.Get("PSI_STOPOPT", 1000, "n_{stopopt}", METHOD);
  m_ndecopt = reader.Get("PSI_NDECOPT", 10, "n_{opt,dec}", METHOD);

  // time steps
  m_timestep = reader.Get("PSI_TIMESTEP_OFFSET", 0.0, "\\Delta t offset", METHOD);
  m_timeslope = reader.Get("PSI_TIMESTEP_SLOPE", 0.0, "\\Delta t slope", METHOD);

#ifdef USING__MPI
  int size=MPI::COMM_WORLD.Get_size();
  m_itminbynode=Max(1,Max(1000,(int)m_itmin)/size);
  if (size) {
    int helpi;
    if (reader.Read(helpi,"PSI_ITMIN_BY_NODE", 0)) {
      m_itmin=(m_itminbynode=helpi)*size;
      msg_Info()<<METHOD<<"(): Set n_{it,min} = "<<m_itmin<<".\n";
    }
    if (reader.Read(helpi,"PSI_ITMAX_BY_NODE", 0)) {
      m_itmax*=helpi*size;
      msg_Info()<<METHOD<<"(): Set n_{it,max} = "<<m_itmax<<".\n";
    }
    if (reader.Read(helpi,"PSI_IT_BY_NODE", 0)) {
      m_itmin=m_itmax=(m_itminbynode=helpi)*size;
      msg_Info()<<METHOD<<"(): Set n_{it} = "<<m_itmin<<".\n";
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
  p_psh->MPISync();
  int size=MPI::COMM_WORLD.Get_size();
  if (size>1) {
    double values[3];
    values[0]=m_mn;
    values[1]=m_mnstep;
    values[2]=m_mncstep;
    mpi->MPIComm()->Allreduce(MPI_IN_PLACE,values,3,MPI::DOUBLE,MPI::SUM);
    m_mn=values[0];
    m_mnstep=values[1];
    m_mncstep=values[2];
  }
  m_n+=m_mn;
  m_nstep+=m_mnstep;
  m_ncstep+=m_mncstep;
  m_mn=m_mnstep=m_mncstep=0;
  m_ncontrib=p_psh->FSRIntegrator()->ValidN();
  m_nlo=0;
#else
  m_nlo=p_psh->FSRIntegrator()->ValidN();
#endif
  m_lrtime=ATOOLS::rpa->gen.Timer().RealTime();
}

double Phase_Space_Integrator::Calculate(double _maxerror, double _maxabserror,
                                         int _fin_opt)
{
  m_mn=m_mnstep=m_mncstep=0;
  m_maxerror=_maxerror;
  m_maxabserror=_maxabserror;
  m_fin_opt=_fin_opt;
  msg_Info()<<"Starting the calculation at "
            <<rpa->gen.Timer().StrFTime("%H:%M:%S")
            <<". Lean back and enjoy ... ."<<endl;
  if (m_maxerror >= 1.) { m_nrawmin=0; m_nmin=0; m_nrawmax=1; m_nmax=1; }

  long unsigned int numberofchannels = 1;

  msg_Tracking()<<"Integrators : "<<p_psh->BeamIntegrator()<<" / "
                <<p_psh->ISRIntegrator()<<" / "<<p_psh->FSRIntegrator()<<endl;

   if ((p_psh->BeamIntegrator())) {
     (p_psh->BeamIntegrator())->Reset();
     numberofchannels = p_psh->BeamIntegrator()->NChannels();
     msg_Tracking()<<"   Found "<<p_psh->BeamIntegrator()->NChannels()
                   <<" Beam Integrators."<<endl;
   }
   if ((p_psh->ISRIntegrator())) {
     (p_psh->ISRIntegrator())->Reset();
     numberofchannels += p_psh->ISRIntegrator()->NChannels();
     msg_Tracking()<<"   Found "<<p_psh->ISRIntegrator()->NChannels()
                   <<" ISR Integrators."<<endl;
   }

  p_psh->FSRIntegrator()->Reset();
  numberofchannels += p_psh->FSRIntegrator()->NChannels();
  msg_Tracking()<<"   Found "<<p_psh->FSRIntegrator()->NChannels()
                <<" FSR integrators."<<endl;
  long unsigned int procitmin(p_psh->Process()->ItMin());
  m_iter = Min(m_itmax,Max(m_itmin,Max(p_psh->Process()->ItMin(),20*numberofchannels)));

  m_ncontrib = p_psh->FSRIntegrator()->ValidN();

#ifdef USING__MPI
  m_nlo=0;
#else
  m_nlo=p_psh->FSRIntegrator()->ValidN();
#endif

  m_addtime = 0.0;
  m_stepstart = m_lotime = m_starttime = ATOOLS::rpa->gen.Timer().RealTime();
  if (p_psh->Stats().size()>0)
    m_addtime=p_psh->Stats().back()[6];

  m_nstep = m_ncstep = 0;

  m_lrtime = ATOOLS::rpa->gen.Timer().RealTime();
  m_optiter=m_iter;
#ifdef USING__MPI
  int size = MPI::COMM_WORLD.Get_size();
  m_optiter /= size;
  if (MPI::COMM_WORLD.Get_rank()==0) m_optiter+=m_iter-(m_iter/size)*size;
#endif

  while (m_n<m_nrawmax && m_ncontrib<m_nmax) {
    if (!rpa->gen.CheckTime()) {
      msg_Error()<<ATOOLS::om::bold
			 <<"\nPhase_Space_Integrator::Calculate(): "
			 <<ATOOLS::om::reset<<ATOOLS::om::red
			 <<"Timeout. Interrupt integration."
			 <<ATOOLS::om::reset<<std::endl;
      kill(getpid(),SIGINT);
    }

    if (AddPoint(p_psh->Differential())) break;
  }

  return p_psh->Process()->TotalResult() * rpa->Picobarn();

}

bool Phase_Space_Integrator::AddPoint(const double value)
{
  if (IsBad(value)) {
    msg_Error()<<METHOD<<"(): value = "<<value<<". Skip."<<endl;
    return false;
  }

#ifdef USING__MPI
  ++m_mn;
  m_mnstep++;
  if (value!=0.) m_mncstep++;
#else
  ++m_n;
  m_nstep++;
  if (value!=0.) m_ncstep++;
#endif

  p_psh->AddPoint(value);

#ifdef USING__MPI
  m_ncontrib = p_psh->FSRIntegrator()->ValidMN();
#else
  m_ncontrib = p_psh->FSRIntegrator()->ValidN();
#endif
  double deltat(0.);
  double targettime(m_timestep+dabs(m_timeslope)*(p_psh->Process()->NOut()-2));
  if (m_timeslope<0.0) targettime*=p_psh->Process()->Process()->Size();
  if (m_timestep>0.0) deltat = ATOOLS::rpa->gen.Timer().RealTime()-m_stepstart;
  if ((m_timestep==0.0 && m_ncontrib!=m_nlo && m_ncontrib>0 &&
       ((m_ncontrib%m_optiter)==0)) ||
      (m_timestep>0.0 && deltat>=targettime)) {
    MPISync();
    bool optimized=false;
    bool fotime = false;
    msg_Tracking()<<" n="<<m_ncontrib<<"  iter="<<m_iter<<endl;
    if (p_psh->Stats().size()<m_nopt) {
      p_psh->Optimize();
      p_psh->Process()->OptimizeResult();
      if ((p_psh->Process())->SPoints()==0)
        m_lotime = ATOOLS::rpa->gen.Timer().RealTime();
      fotime    = true;
      optimized = true;
    }
    else if (p_psh->Stats().size()==m_nopt) {
      p_psh->Process()->ResetMax(0);
      p_psh->EndOptimize();
      p_psh->Process()->ResetMax(1);
      p_psh->Process()->InitWeightHistogram();
      p_psh->Process()->EndOptimize();
      m_lotime = ATOOLS::rpa->gen.Timer().RealTime();
    }

    double time = ATOOLS::rpa->gen.Timer().RealTime();
    double timeest=0.;
    timeest = (m_nopt*m_iter+m_maxopt*m_iter)/double(m_ncontrib)*(time-m_starttime);
    if (!fotime) {
      if (m_fin_opt==1) {
        timeest = ATOOLS::Max(timeest,
                              p_psh->Process()->RemainTimeFactor(m_maxerror)*
                              (time-m_lotime)+m_lotime-m_starttime);
      }
      else {
        timeest = p_psh->Process()->RemainTimeFactor(m_maxerror)*
                  (time-m_lotime)+m_lotime-m_starttime;
      }
    }
    double error=dabs(p_psh->Process()->TotalVar()/
                      p_psh->Process()->TotalResult());
    if (m_maxabserror>0.0) {
      msg_Info()<<om::blue
                <<p_psh->Process()->TotalResult()*rpa->Picobarn()
                <<" pb"<<om::reset<<" +- ( "<<om::red
                <<p_psh->Process()->TotalVar()*rpa->Picobarn()
                <<" pb <-> "<<m_maxabserror<<" pb"<<om::reset<<" ) "
                <<m_ncontrib<<" ( "<<m_n<<" -> "<<(m_ncstep*1000/m_nstep)/10.0
                <<" % )"<<endl;
    }
    else {
      msg_Info()<<om::blue
                <<p_psh->Process()->TotalResult()*rpa->Picobarn()
                <<" pb"<<om::reset<<" +- ( "<<om::red
                <<p_psh->Process()->TotalVar()*rpa->Picobarn()
                <<" pb = "<<error*100<<" %"<<om::reset<<" ) "
                <<m_ncontrib<<" ( "<<m_n<<" -> "<<(m_ncstep*1000/m_nstep)/10.0
                <<" % )"<<endl;
    }
    if (optimized) m_nstep = m_ncstep = 0;
    if (fotime) { msg_Info()<<"full optimization: "; }
    else        { msg_Info()<<"integration time:  "; }
    msg_Info()<<" ( "<<FormatTime(size_t(time-m_starttime))<<" elapsed / "
              <<FormatTime(size_t(timeest)-size_t((time-m_starttime)))
              <<" left ) ["<<rpa->gen.Timer().StrFTime("%H:%M:%S")<<"]"<<endl;
    size_t currentrss=GetCurrentRSS();
    if (m_lastrss==0) m_lastrss=currentrss;
    else if (currentrss>m_lastrss+ToType<int>
        (rpa->gen.Variable("MEMLEAK_WARNING_THRESHOLD"))) {
      msg_Error()<<METHOD<<"() {\n"<<om::bold<<"  Memory usage increased by "
                 <<(currentrss-m_lastrss)/(1<<20)<<" MB,"
                 <<" now "<<currentrss/(1<<20)<<" MB.\n"
                 <<om::red<<"  This might indicate a memory leak!\n"
                 <<"  Please monitor this process closely.\n"<<om::reset
                 <<"}"<<std::endl;
      m_lastrss=currentrss;
    }
    std::vector<double> stats(6);
    stats[0]=p_psh->Process()->TotalResult()*rpa->Picobarn();
    stats[1]=p_psh->Process()->TotalVar()*rpa->Picobarn();
    stats[2]=error;
    stats[3]=m_ncontrib;
    stats[4]=m_ncontrib/(double)m_n;
    stats[5]=time-m_starttime+m_addtime;
    p_psh->AddStats(stats);
    p_psh->Process()->StoreResults(1);
    m_stepstart=ATOOLS::rpa->gen.Timer().RealTime();
    if (m_n>=m_nrawmin && m_ncontrib>=m_nmin) {
      double var(p_psh->Process()->TotalVar());
      bool wannabreak = dabs(error)<m_maxerror ||
                        (var!=0. && dabs(var*rpa->Picobarn())<m_maxabserror);
      if (m_fin_opt==0 && m_nopt>p_psh->Stats().size() && wannabreak)
        m_nopt=p_psh->Stats().size();
      if (wannabreak && p_psh->Stats().size()>=m_nopt+m_maxopt) return true;
      if (p_psh->Stats().size()>=m_nopt+m_stopopt) return true;
    }
  }
  return false;
}

double Phase_Space_Integrator::CalculateDecay(double maxerror)
{
  m_mn=m_mnstep=m_mncstep=0;
  msg_Info()<<"Starting the calculation for a decay. Lean back and enjoy ... ."
            <<endl;

  m_optiter = m_iter = 20000;

  p_psh->FSRIntegrator()->Reset();

  for (long unsigned int n=1;n<=m_nrawmax;n++) {
    double value = p_psh->Differential();
    p_psh->AddPoint(value);

    if (!(n%m_iter)) {
      MPISync();
      if (p_psh->Stats().size()<=m_ndecopt) {
        p_psh->Optimize();
        p_psh->Process()->OptimizeResult();
      }
      if (p_psh->Stats().size()==m_ndecopt) {
        p_psh->EndOptimize();
        m_optiter = m_iter = 50000;
      }
      if (p_psh->Process()->TotalResult()==0.) break;

      double error = p_psh->Process()->TotalVar()/
                     p_psh->Process()->TotalResult();

      msg_Info()<<om::blue
                <<p_psh->Process()->TotalResult()
                <<" GeV"<<om::reset<<" +- ( "<<om::red
                <<p_psh->Process()->TotalVar()
                <<" GeV = "<<error*100<<" %"<<om::reset<<" ) "<<n<<endl;
      if (error<maxerror) break;
    }
  }
  return p_psh->Process()->TotalResult()*rpa->Picobarn();
}


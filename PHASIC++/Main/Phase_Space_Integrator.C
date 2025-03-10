#include "PHASIC++/Main/Phase_Space_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "PHASIC++/Channels/Single_Channel.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/RUsage.H"
#include "PHASIC++/Main/Process_Integrator.H"

#include <signal.h>
#include <unistd.h>

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

long unsigned int Phase_Space_Integrator::m_nrawmax(std::numeric_limits<long unsigned int>::max());

Phase_Space_Integrator::Phase_Space_Integrator(Phase_Space_Handler *_psh):
  m_iter(1000), m_itmin(1000), m_itmax(1000000),
  m_n(0), m_nstep(0), m_ncstep(0), m_mn(0), m_mnstep(0), m_mncstep(0),
  m_ncontrib(0), m_maxopt(0), m_stopopt(1000), m_nlo(0), m_fin_opt(true),
  m_starttime(0.), m_lotime(0.), m_addtime(0.), m_lrtime(0.),
  m_maxerror(0.), m_maxabserror(0.), m_lastrss(0), p_psh(_psh)
{
  RegisterDefaults();
  Scoped_Settings s{ Settings::GetMainSettings()["PSI"] };
  // total number of points
  m_nrawmax = s["NRAWMAX"].Get<long unsigned int>();
  // number of optimisation steps
  m_npower = s["NPOWER"].Get<double>();
  m_nopt = s["NOPT"].GetScalarWithOtherDefault
    <long unsigned int>(m_npower?10:25);
  m_maxopt = s["MAXOPT"].GetScalarWithOtherDefault
    <long unsigned int>(m_npower?1:5);
  m_stopopt = s["STOPOPT"].Get<long unsigned int>();
  // number of points per iteration
  const size_t procitmin = p_psh->Process()->Process()->Info().m_itmin;
  m_itmin = s["ITMIN"].GetScalarWithOtherDefault
    <long unsigned int>((m_npower?5:5)*procitmin);
  const size_t procitmax = p_psh->Process()->Process()->Info().m_itmax;
  m_itmax = s["ITMAX"].GetScalarWithOtherDefault
    <long unsigned int>((m_npower?5:5)*procitmax);
  // time steps
  m_timestep = s["TIMESTEP_OFFSET"].Get<double>();
  m_timeslope = s["TIMESTEP_SLOPE"].Get<double>();
#ifdef USING__MPI
  int size=mpi->Size();
  long unsigned int itminbynode=Max(1,(int)m_itmin/size),
                    itmaxbynode=Max(1,(int)m_itmax/size);
  if (size) {
    int helpi;
    if (s["ITMIN_BY_NODE"].IsSetExplicitly())
      itminbynode = s["ITMIN_BY_NODE"].Get<long unsigned int>();
    m_itmin = itminbynode * size;
    if (s["ITMAX_BY_NODE"].IsSetExplicitly())
      itmaxbynode = s["ITMAX_BY_NODE"].Get<long unsigned int>();
    m_itmax = itmaxbynode * size;
  }
#endif
  m_nexpected = m_itmin;
  for (size_t i(1);i<m_nopt;++i) m_nexpected+=m_itmin*pow(2.,i*m_npower);
  m_nexpected+=m_maxopt*m_itmin*pow(2.,m_nopt*m_npower);
  msg_Info()<<"Integration parameters: n_{min} = "<<m_itmin
            <<", n_{max} = "<<m_itmax
	    <<", N_{opt} = "<<m_nopt<<", N_{max} = "<<m_maxopt;
  if (m_npower) msg_Info()<<", exponent = "<<m_npower;
  msg_Info()<<std::endl;
}

Phase_Space_Integrator::~Phase_Space_Integrator()
{
}

void Phase_Space_Integrator::RegisterDefaults() const
{
  Scoped_Settings s{ Settings::GetMainSettings()["PSI"] };
  s["NRAWMAX"].SetDefault(std::numeric_limits<long unsigned int>::max());  // n_{max,raw}
  s["NPOWER"].SetDefault(.5);
  s["STOPOPT"].SetDefault(0);  // n_{stopopt}
  s["TIMESTEP_OFFSET"].SetDefault(0.0);  // \Delta t offset
  s["TIMESTEP_SLOPE"].SetDefault(0.0);  // \Delta t slope
  s["ITMIN_BY_NODE"].SetDefault(0);
  s["IT_BY_NODE"].SetDefault(0);
}

void Phase_Space_Integrator::MPISync()
{
#ifdef USING__MPI
  p_psh->MPISync();
  int size=mpi->Size();
  if (size>1) {
    double values[3];
    values[0]=m_mn;
    values[1]=m_mnstep;
    values[2]=m_mncstep;
    mpi->Allreduce(values,3,MPI_DOUBLE,MPI_SUM);
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
                                         bool _fin_opt)
{
  if (p_psh->Stats().size()>=m_nopt+m_maxopt+m_stopopt) return true;
  m_mn=m_mnstep=m_mncstep=0;
  m_maxerror=_maxerror;
  m_maxabserror=_maxabserror;
  m_fin_opt=_fin_opt;
  msg_Info()<<"Starting the calculation at "
            <<rpa->gen.Timer().StrFTime("%H:%M:%S")
            <<". Lean back and enjoy ... ."<<endl;

  msg_Tracking()<<"Integrators : "<<p_psh->BeamIntegrator()<<" / "
                <<p_psh->ISRIntegrator()<<" / "<<p_psh->FSRIntegrator()<<endl;

   if ((p_psh->BeamIntegrator())) {
     (p_psh->BeamIntegrator())->Reset();
     msg_Tracking()<<"   Found "<<p_psh->BeamIntegrator()->NChannels()
                   <<" Beam Integrators."<<endl;
   }
   if (p_psh->ISRIntegrator()) {
     (p_psh->ISRIntegrator())->Reset();
     msg_Tracking()<<"   Found "<<p_psh->ISRIntegrator()->NChannels()
                   <<" ISR Integrators."<<endl;
   }

  p_psh->FSRIntegrator()->Reset();
  msg_Tracking()<<"   Found "<<p_psh->FSRIntegrator()->NChannels()
                <<" FSR integrators."<<endl;

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
  m_iter = Min(m_itmin,m_itmax);
  if (p_psh->Stats().size()) {
    m_iter=p_psh->Stats().back()[4];
    if (p_psh->Stats().size()>1)
      m_iter-=p_psh->Stats()[p_psh->Stats().size()-2][4];
    m_iter*=pow(2.,m_npower);
  }
#ifdef USING__MPI
  int size = mpi->Size();
  m_iter /= size;
#endif

  while (m_n<m_nrawmax) {
    if (!rpa->gen.CheckTime()) {
      msg_Error()<<ATOOLS::om::bold
			 <<"\nPhase_Space_Integrator::Calculate(): "
			 <<ATOOLS::om::reset<<ATOOLS::om::red
			 <<"Timeout. Interrupt integration."
			 <<ATOOLS::om::reset<<std::endl;
      kill(getpid(),SIGINT);
    }

    if (AddPoint(double(p_psh->Differential(Variations_Mode::nominal_only)))) {
      break;
    }
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
       (((m_ncontrib-m_nlo)%m_iter)==0)) ||
      (m_timestep>0.0 && deltat>=targettime)) {
    MPISync();
    bool optimized=false;
    bool fotime = false;
    msg_Tracking()<<" n="<<m_ncontrib<<"  iter="<<m_iter<<endl;
    if (p_psh->Stats().size()<m_nopt-1) {
      p_psh->Optimize();
      p_psh->Process()->OptimizeResult();
      if ((p_psh->Process())->SPoints()==0)
        m_lotime = ATOOLS::rpa->gen.Timer().RealTime();
      fotime    = true;
      optimized = true;
      m_iter*=pow(2.,m_npower);
    }
    else if (p_psh->Stats().size()==m_nopt-1) {
      p_psh->Process()->ResetMax(0);
      p_psh->EndOptimize();
      p_psh->Process()->ResetMax(1);
      p_psh->Process()->InitWeightHistogram();
      p_psh->Process()->EndOptimize();
      m_lotime = ATOOLS::rpa->gen.Timer().RealTime();
      fotime    = true;
    }
    double time = ATOOLS::rpa->gen.Timer().RealTime();
    double timeest=0.;
    timeest = m_nexpected/double(m_ncontrib)*(time-m_starttime);
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
    double var(p_psh->Process()->TotalVar());
    bool wannabreak = dabs(error)<m_maxerror ||
      (var!=0. && dabs(var*rpa->Picobarn())<m_maxabserror);
    if (!m_fin_opt && wannabreak && m_nopt>p_psh->Stats().size())
      m_nopt=p_psh->Stats().size();
    if (wannabreak && p_psh->Stats().size()>=m_nopt+m_maxopt) return true;
    if (p_psh->Stats().size()>=m_nopt+m_maxopt+m_stopopt) return true;
  }
  return false;
}

double Phase_Space_Integrator::CalculateDecay(double maxerror)
{
  m_mn=m_mnstep=m_mncstep=0;
  msg_Info()<<"Starting the calculation for a decay. Lean back and enjoy ... ."
            <<endl;

  m_iter = 20000;

  p_psh->FSRIntegrator()->Reset();

  for (long unsigned int n=1;n<=m_nrawmax;n++) {
    double value = double(p_psh->Differential(Variations_Mode::nominal_only));
    p_psh->AddPoint(value);

    if (!(n%m_iter)) {
      MPISync();
      if (p_psh->Stats().size()<=m_nopt) {
        p_psh->Optimize();
        p_psh->Process()->OptimizeResult();
      }
      if (p_psh->Stats().size()==m_nopt) {
        p_psh->EndOptimize();
        m_iter = 50000;
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


#ifndef PHASIC_Main_Phase_Space_Integrator_H
#define PHASIC_Main_Phase_Space_Integrator_H

#include <stddef.h>

namespace PHASIC {class Phase_Space_Handler;}

namespace PHASIC {

  class Phase_Space_Integrator {
  private:
    static long unsigned int m_nrawmax;
    long unsigned int m_nopt, m_iter, m_itmin, m_nexpected;

    long unsigned int m_n, m_nstep, m_ncstep;
    long unsigned int m_mn, m_mnstep, m_mncstep;
    long unsigned int m_ncontrib, m_maxopt, m_nlo;

    bool m_fin_opt;

    double m_starttime, m_addtime, m_lrtime;
    double m_maxerror, m_maxabserror, m_npower;

    size_t m_lastrss;

    double m_stepstart, m_timestep, m_timeslope;

    Phase_Space_Handler *p_psh;

    void RegisterDefaults() const;

  public:

    Phase_Space_Integrator(Phase_Space_Handler*);

    ~Phase_Space_Integrator();

    double Calculate(double,double,bool);
    double CalculateDecay(double);

    void MPISync();

    bool AddPoint(const double value);

  };
  /*!
    To perform the phase space integral. Its main task is to
    call the integration channel repeatedly and to sample the process
    over the resulting points. In doing so, the Phase_Space_Integrator
    optimizes the Multi_Channel (if this type of integration is
    employed). It should be noted that the Phase_Space_Integrator
    calls the Channels via the Phase_Space_Handler, since it
    contains the links to both channels and processes.
    In that respect, the Integrator is a mere tool, initialized and
    steered by the handler.
  */
  /*! 
    Global results of the integration:
    Value of the results, its actual error and the maximal value 
    given by the combination of matrixelement and phase space point
    so far.
  */
  /*!
    Parameters for steering the integration
  */
  /*!
    Constructor initializes the steering parameters and the
    channel(s) for the integration.
  */
  /*!
    Calculate the total cross section of a 2->N process up to
    an error given externally by the Phase_Space_Handler.
  */
}
#endif

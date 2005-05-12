#ifndef Phase_Space_Integrator_h
#define Phase_Space_Integrator_h

namespace PHASIC {class Phase_Space_Handler;}
namespace ATOOLS {class Algebra_Interpreter;}

namespace PHASIC {

  class Phase_Space_Integrator {
  private:

    double          result, max, error;
    int             iter,iter0,iter1,nopt;
    long int        maxopt;
    static long int nmax;

    double starttime, lotime, totalopt;
    double value, maxerror;
    int ncontrib, endopt, nlo, fin_opt;
    long int  n;

#ifdef _USE_MPI_
    long int alln;
    int iterall;
    double allsum;
    double allsum2;
#endif
    struct {
      int    id;
      int    n;
      double sum;
      double sum2;
      double max;
    } local, message;

    Phase_Space_Handler *psh;

    ATOOLS::Algebra_Interpreter *p_interpreter;
 
  public:

    Phase_Space_Integrator();
    
    ~Phase_Space_Integrator();

    double Calculate(Phase_Space_Handler*,double,int); 
    double CalculateDecay(Phase_Space_Handler*,double,double); 

    //! Simple access methods.
    inline double GetMax()            { return max;  };
    inline void   SetMax(double _max) { max = _max;  };

    static long int MaxPoints();
    static void     SetMaxPoints(long int _nmax);

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
